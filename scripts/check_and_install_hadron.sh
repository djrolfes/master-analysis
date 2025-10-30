#!/bin/bash
# Fix: Get the project root directory (parent of scripts/)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Set R_LIBS_USER to a local directory within the project to keep it self-contained
export R_LIBS_USER="${PROJECT_ROOT}/.R/library"
mkdir -p "$R_LIBS_USER"

HADRON_DIR="${PROJECT_ROOT}/extern/hadron"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Checking hadron package installation...${NC}"
echo -e "${YELLOW}Project root: ${PROJECT_ROOT}${NC}"
echo -e "${YELLOW}R library path set to: ${R_LIBS_USER}${NC}"
echo -e "${YELLOW}Hadron directory: ${HADRON_DIR}${NC}"
echo -e "${YELLOW}Note: This script also installs minpack.lm and other dependencies for analysis scripts${NC}"

# Function to check and install required R packages
check_and_install_r_deps() {
    # Define required packages in an array
    local required_pkgs=("$@")
    # Convert bash array to a comma-separated string for R
    local pkgs_str=$(printf "'%s'," "${required_pkgs[@]}")
    pkgs_str="c(${pkgs_str%,})"

    # Capture R script output to check for our custom error message
    local output
    output=$(R --slave -e "
    # Ensure the local library path is used
    .libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths()))
    
    installed_pkgs <- installed.packages(lib.loc=.libPaths())[,'Package']
    missing_pkgs <- setdiff(${pkgs_str}, installed_pkgs)
    
    if (length(missing_pkgs) > 0) {
        cat('Installing missing R packages:', paste(missing_pkgs, collapse=', '), '\n')
        
        # Attempt to install
        install.packages(
            missing_pkgs, 
            dependencies=TRUE, 
            repos='https://cran.r-project.org',
            lib=Sys.getenv('R_LIBS_USER')
        )
        
        # After installation, check again if they are present
        installed_pkgs_after <- installed.packages(lib.loc=.libPaths())[,'Package']
        still_missing <- setdiff(missing_pkgs, installed_pkgs_after)
        
        if (length(still_missing) > 0) {
            cat('ERROR: The following R packages failed to install:', paste(still_missing, collapse=', '), '\n')
            cat('This is often due to missing system-level development libraries.\n')
            cat('On Debian/Ubuntu, try: sudo apt-get install libgsl-dev libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev\n')
            cat('On RHEL/CentOS, try: sudo dnf install gsl-devel libcurl-devel openssl-devel libxml2-devel fontconfig-devel freetype-devel harfbuzz-devel fribidi-devel libpng-devel libtiff-devel libjpeg-devel libwebp-devel\n')
            cat('Please contact your cluster administrator if you cannot install these.\n')
            q(status = 1) # Exit R with an error code
        }
    } else {
        cat('All required R packages are already installed.\n')
    }
    " 2>&1) # Capture both stdout and stderr

    # Print the output from the R script
    printf "%s\n" "$output"

    # Explicitly check for the error message and exit if found
    if printf "%s" "$output" | grep -q "ERROR: The following R packages failed to install:"; then
        printf "\n${RED}Aborting script because required R packages failed to install.${NC}\n"
        exit 1
    fi
}

# Function to check if hadron is installed
check_hadron_installed() {
    R --slave -e "
    .libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths()))
    if ('hadron' %in% installed.packages()[,'Package']) {
        cat('INSTALLED')
        cat('|')
        cat(as.character(packageVersion('hadron')))
        cat('|')
        cat(find.package('hadron'))
    } else {
        cat('NOT_INSTALLED')
    }" 2>/dev/null
}

# Function to get local version from DESCRIPTION
get_local_version() {
    if [ -f "${HADRON_DIR}/DESCRIPTION" ]; then
        grep "^Version:" "${HADRON_DIR}/DESCRIPTION" | cut -d' ' -f2
    else
        echo "unknown"
    fi
}

# Function to show R library paths
show_r_library_info() {
    echo -e "${YELLOW}R Library Information:${NC}"
    R --slave -e "
    .libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths()))
    cat('Library paths (writable are preferred for installation):\n')
    for(i in seq_along(.libPaths())) {
        cat(sprintf('  %d. %s', i, .libPaths()[i]))
        if(file.access(.libPaths()[i], 2) == 0) {
            cat(' (writable)')
        } else {
            cat(' (read-only)')
        }
        cat('\n')
    }
    "
}

# Show R library information
show_r_library_info

# --- Check for Essential System Dependencies ---
printf "\n${YELLOW}Checking for essential system dependencies (like GSL)...${NC}\n"
if ! command -v gsl-config &> /dev/null; then
    printf "${RED}ERROR: 'gsl-config' command not found.${NC}\n"
    printf "The GNU Scientific Library (GSL) is required by the 'hadron' package.\n"
    printf "On Debian/Ubuntu, try: sudo apt-get install libgsl-dev\n"
    printf "On RHEL/CentOS, try: sudo dnf install gsl-devel\n"
    exit 1
else
    printf "${GREEN}GSL (gsl-config) found.${NC}\n"
fi

# --- Check and Install Core R Dependencies ---
printf "\n${YELLOW}Checking for core R dependencies...${NC}\n"
# Note: This step requires system dependencies like libcurl-devel to be installed.
# Added minpack.lm, ggplot2, yaml, and errors for analysis scripts
check_and_install_r_deps "devtools" "roxygen2" "Rcpp" "abind" "boot" "dplyr" "R6" "stringr" "zoo" "tikzDevice" "ggplot2" "yaml" "errors" #"minpack.lm"

# Check if hadron directory exists
if [ ! -d "${HADRON_DIR}" ]; then
    echo -e "${RED}Error: ${HADRON_DIR} directory not found!${NC}"
    echo -e "${RED}Current working directory: $(pwd)${NC}"
    exit 1
fi

# Check current installation status
echo -e "\n${YELLOW}Checking for local 'hadron' package...${NC}"
INSTALL_STATUS=$(check_hadron_installed)

if [[ $INSTALL_STATUS == "NOT_INSTALLED" ]]; then
    echo -e "${YELLOW}hadron package not found. Installing locally...${NC}"
    NEED_INSTALL=true
else
    IFS='|' read -r status installed_version install_path <<< "$INSTALL_STATUS"
    local_version=$(get_local_version)
    
    echo -e "${GREEN}hadron package found:${NC}"
    echo "  Installed version: $installed_version"
    echo "  Local version: $local_version"
    echo "  Install path: $install_path"
    
    # Check if we need to update
    if [ "$local_version" != "$installed_version" ] && [ "$local_version" != "unknown" ]; then
        echo -e "${YELLOW}Version mismatch detected. Reinstalling...${NC}"
        NEED_INSTALL=true
    else
        echo -e "${GREEN}hadron package is up to date. Nothing to do.${NC}"
        NEED_INSTALL=false
    fi
fi

# Install if needed
if [ "$NEED_INSTALL" = true ]; then
    echo -e "${YELLOW}Installing hadron package from: ${HADRON_DIR}${NC}"
    
    # Remove existing installation if present (to ensure clean install)
    echo -e "${YELLOW}Removing any existing hadron installation...${NC}"
    R --slave -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); remove.packages('hadron', lib=Sys.getenv('R_LIBS_USER'))" 2>/dev/null || true
    
    # Save current directory
    ORIGINAL_DIR=$(pwd)
    
    cd "${HADRON_DIR}"
    
    # Check if install script exists and is executable
    if [ ! -f "install" ]; then
        echo -e "${RED}Error: install script not found in ${HADRON_DIR}${NC}"
        exit 1
    fi
    
    if [ ! -x "install" ]; then
        echo -e "${YELLOW}Making install script executable...${NC}"
        chmod +x install
    fi
    
    # Modify DESCRIPTION to disable lazy loading
    # This is the only way to prevent the "preparing package for lazy loading" step
    echo -e "${YELLOW}Modifying DESCRIPTION to disable lazy loading for cluster compatibility...${NC}"
    
    # Backup original DESCRIPTION
    cp DESCRIPTION DESCRIPTION.bak
    
    # Add or modify LazyData and LazyLoad settings
    if grep -q "^LazyData:" DESCRIPTION; then
        sed -i 's/^LazyData:.*/LazyData: no/' DESCRIPTION
    else
        echo "LazyData: no" >> DESCRIPTION
    fi
    
    if grep -q "^LazyLoad:" DESCRIPTION; then
        sed -i 's/^LazyLoad:.*/LazyLoad: no/' DESCRIPTION
    else
        echo "LazyLoad: no" >> DESCRIPTION
    fi
    
    echo -e "${YELLOW}Installing hadron package with lazy loading disabled...${NC}"
    echo -e "${YELLOW}Note: This prevents CPU instruction incompatibility issues on the cluster${NC}"
    
    # Set environment variables to disable JIT and compilation
    export R_ENABLE_JIT=0
    export R_COMPILE_PKGS=0
    export R_KEEP_PKG_SOURCE=no
    
    # Use R CMD INSTALL directly
    R --vanilla CMD INSTALL \
        --no-staged-install \
        --no-byte-compile \
        --no-data \
        --no-help \
        --no-html \
        --no-demo \
        --no-test-load \
        --library="${R_LIBS_USER}" \
        .
    
    # Store the result
    INSTALL_RESULT=$?
    
    # Restore original DESCRIPTION
    mv DESCRIPTION.bak DESCRIPTION
    
    # Check if installation failed
    if [ $INSTALL_RESULT -ne 0 ]; then
        echo -e "${RED}Installation command failed with exit code ${INSTALL_RESULT}${NC}"
        cd "${ORIGINAL_DIR}"
        exit 1
    fi
    
    # Return to original directory
    cd "${ORIGINAL_DIR}"
    
    # Verify installation
    VERIFY_STATUS=$(check_hadron_installed)
    
    if [[ $VERIFY_STATUS == "NOT_INSTALLED" ]]; then
        echo -e "${RED}Error: Installation failed!${NC}"
        exit 1
    else
        IFS='|' read -r status new_version install_path <<< "$VERIFY_STATUS"
        echo -e "${GREEN}Success! hadron package installed:${NC}"
        echo "  Version: $new_version"
        echo "  Path: $install_path"
    fi
fi

echo -e "${GREEN}hadron package is ready to use!${NC}"

# Optional: Test loading the package
echo -e "${YELLOW}Testing package loading...${NC}"
if R --slave -e ".libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths())); library(hadron); cat('Package loaded successfully')" 2>/dev/null; then
    echo -e "${GREEN}✓ Package loads correctly${NC}"
else
    echo -e "${RED}✗ Warning: Package installed but failed to load${NC}"
fi