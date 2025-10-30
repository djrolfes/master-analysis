#!/bin/bash
# Script to verify the cluster fix deployment

echo "=== Cluster Fix Verification Script ==="
echo ""

# Check if we're in the right directory
if [ ! -f "analysis/analysis_wilson_temp.R" ]; then
    echo "ERROR: Run this script from the master-analysis root directory"
    exit 1
fi

echo "✓ Found analysis/analysis_wilson_temp.R"

# Check if the fix is present
if grep -q "error_scale <- median(safe_errors)" analysis/analysis_wilson_temp.R; then
    echo "✓ Enhanced numerical safeguards present in analysis_wilson_temp.R"
else
    echo "✗ WARNING: Enhanced safeguards NOT found - old version may be deployed"
    echo "  Expected to find: error_scale <- median(safe_errors)"
fi

# Check if weight range logging is present
if grep -q "weight range:" analysis/analysis_wilson_temp.R; then
    echo "✓ Enhanced logging (weight range) present"
else
    echo "✗ WARNING: Enhanced logging NOT found"
fi

# Check if hadron install script includes new dependencies
if grep -q "minpack.lm" scripts/check_and_install_hadron.sh; then
    echo "✓ Updated dependency list includes minpack.lm, ggplot2, yaml"
else
    echo "✗ WARNING: Dependency installation script not updated"
fi

echo ""
echo "=== R Package Check ==="

# Check installed R packages
R --slave -e "
required_pkgs <- c('hadron', 'dplyr', 'ggplot2', 'yaml')
installed <- installed.packages()[,'Package']
missing <- setdiff(required_pkgs, installed)

if (length(missing) > 0) {
    cat('✗ MISSING packages:', paste(missing, collapse=', '), '\n')
    cat('  Run: ./scripts/check_and_install_hadron.sh\n')
} else {
    cat('✓ All required packages installed\n')
    for (pkg in required_pkgs) {
        cat('  -', pkg, 'version', as.character(packageVersion(pkg)), '\n')
    }
}

# Check optional packages
optional <- c('minpack.lm')
missing_opt <- setdiff(optional, installed)
if (length(missing_opt) > 0) {
    cat('⚠ OPTIONAL packages missing:', paste(missing_opt, collapse=', '), '\n')
    cat('  (Not critical but recommended)\n')
}
" 2>/dev/null

echo ""
echo "=== R Version ==="
R --version | head -1

echo ""
echo "=== Next Steps ==="
echo "1. If any checks failed (✗), update the files and re-run this script"
echo "2. Run: ./scripts/check_and_install_hadron.sh"
echo "3. Test with: cd analysis && Rscript analysis_wilson_temp.R <data_dir> 0"
echo "4. Check log: tail -50 <data_dir>/wilson_temp_analysis.log | grep 'weight range'"
echo ""
