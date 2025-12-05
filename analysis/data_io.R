library(yaml)

# Function to find the YAML file in the directory
find_yaml_file <- function(directory) {
  yaml_files <- list.files(directory, pattern = "\\.yaml$", full.names = TRUE)
  if (length(yaml_files) == 0) {
    stop(paste("Error: No .yaml file found in directory:", directory))
  }
  if (length(yaml_files) > 1) {
    warning(paste("Multiple .yaml files found in directory:", directory, ". Using the first one:", basename(yaml_files[1])))
  }
  yaml_files[1]
}

# Function to read the YAML configuration
read_yaml_config <- function(directory) {
  yaml_path <- find_yaml_file(directory)
  if (is.na(yaml_path) || !file.exists(yaml_path)) {
    stop(paste("Error: YAML file not found or invalid path:", yaml_path, "in directory:", directory))
  }
  yaml.load_file(yaml_path)
}

# Function to read data files
read_data_file <- function(filepath) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))

  # Read first line, remove leading # if present, use as header
  con <- file(filepath, "r")
  header_line <- readLines(con, n = 1)
  close(con)

  header_line <- gsub("^#\\s*", "", header_line)
  col_names <- strsplit(header_line, ",")[[1]]
  col_names <- trimws(col_names)

  # Make column names syntactically valid for R
  col_names <- make.names(col_names)

  # Read rest of file, skipping first line, using comma as separator and stripping whitespace
  data <- read.table(filepath, header = FALSE, sep = ",", skip = 1, comment.char = "#", stringsAsFactors = FALSE, strip.white = TRUE)

  # Ensure the number of columns match
  if (ncol(data) != length(col_names)) {
    warning(paste("Column name/data mismatch in", filepath, ". Expected", length(col_names), "got", ncol(data), ". Truncating to fit."))
    min_cols <- min(ncol(data), length(col_names))
    data <- data[, 1:min_cols, drop = FALSE]
    col_names <- col_names[1:min_cols]
  }

  colnames(data) <- col_names
  data
}

read_wilsonflow_details <- function(directory) {
  config <- read_yaml_config(directory)

  # Check if WilsonFlowParams exists
  if (is.null(config$GaugeObservableParams$WilsonFlowParams)) {
    stop("Error: WilsonFlowParams not found in YAML configuration.")
  }

  # Check if wilson_flow_filename exists
  filename <- config$GaugeObservableParams$WilsonFlowParams$wilson_flow_filename
  if (is.null(filename) || filename == "") {
    stop("Error: wilson_flow_filename not found in WilsonFlowParams.")
  }

  # Construct the filepath
  filepath <- file.path(directory, filename)

  # Check if the file exists
  if (!file.exists(filepath)) {
    stop(paste("Error: File not found:", filepath))
  }

  # Read the data file
  read_data_file(filepath)
}

# Function to read simulation log file
read_simulation_log <- function(directory, skip_steps = 0) {
  config <- read_yaml_config(directory)
  filename <- config$SimulationLoggingParams$log_filename
  read_data_file(file.path(directory, filename))
}

# Function to read plaquette data
read_data_plaquette_filename <- function(directory) {
  config <- read_yaml_config(directory)
  filename <- config$GaugeObservableParams$plaquette_filename
  read_data_file(file.path(directory, filename))
}

# Function to read action density data
read_data_action_density_filename <- function(directory) {
  config <- read_yaml_config(directory)
  filename <- config$GaugeObservableParams$action_density_filename
  read_data_file(file.path(directory, filename))
}

# Function to read W_temp data
read_data_W_temp_filename <- function(directory) {
  config <- read_yaml_config(directory)
  filename <- config$GaugeObservableParams$W_temp_filename
  read_data_file(file.path(directory, filename))
}

# Function to read Wilson flow data (e.g., topological_charge_cumulative.txt)
read_wilsonflow_data <- function(directory, filename) {
  filepath <- file.path(directory, filename)
  read_data_file(filepath)
}

# Function to read action densities data (e.g., action_densities_*.txt)
read_action_densities <- function(directory) {
  files <- list.files(directory, pattern = "action_densities_.*\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No action densities files found in the directory.")
  result <- lapply(files, read_data_file)
  names(result) <- basename(files)
  result
}

# Function to read sp_max.txt file
read_sp_max <- function(directory) {
  filepath <- file.path(directory, "sp_max.txt")
  if (!file.exists(filepath)) stop("sp_max.txt file not found in the directory.")
  read_data_file(filepath)
}

# Function to read sp_avg.txt file
read_sp_avg <- function(directory) {
  filepath <- file.path(directory, "sp_avg.txt")
  if (!file.exists(filepath)) stop("sp_avg.txt file not found in the directory.")
  read_data_file(filepath)
}

# Function to read topological_charge_cumulative.txt file
read_topological_charge_cumulative <- function(directory) {
  filepath <- file.path(directory, "topological_charge_cumulative.txt")
  if (!file.exists(filepath)) stop("topological_charge_cumulative.txt file not found in the directory.")
  read_data_file(filepath)
}

# Function to read and parse the ptbcsimulation_log.txt file
read_ptbc_simulation_log <- function(directory) {
  config <- read_yaml_config(directory)

  # Check if PTBCSimulationLoggingParams exists
  if (is.null(config$PTBCSimulationLoggingParams)) {
    stop("Error: PTBCSimulationLoggingParams not found in YAML configuration.")
  }

  # Check if log_filename exists
  filename <- config$PTBCSimulationLoggingParams$log_filename
  if (is.null(filename) || filename == "") {
    stop("Error: log_filename not found in PTBCSimulationLoggingParams.")
  }

  filepath <- file.path(directory, filename)
  if (!file.exists(filepath)) {
    stop(paste("Error: File not found:", filepath))
  }

  lines <- readLines(filepath)
  lines <- lines[!grepl("^#", lines)] # Remove comments

  if (length(lines) == 0) {
    return(data.frame())
  }

  # Helper to parse bracketed strings
  parse_vector <- function(s) {
    as.numeric(strsplit(gsub("\\[|\\]", "", s), "\\s+")[[1]])
  }

  data_list <- lapply(lines, function(line) {
    parts <- strsplit(line, ", ")[[1]]

    # Parse the 4th field - could be "ascending" (old) or "starting_defect_value" (new)
    field4 <- parts[4]
    defects <- parse_vector(parts[5])

    # Determine if this is old format (ascending: 0 or 1) or new format (starting_defect_value: actual defect value)
    # Old format: field4 is "0" or "1" (boolean)
    # New format: field4 is a defect value (e.g., "0.0", "0.137", "1.0")
    field4_numeric <- as.numeric(field4)

    # Check if field4 is a boolean (0 or 1) or a defect value
    # If it's 0 or 1, it's the old "ascending" format
    # Otherwise, it's the new "starting_defect_value" format
    if (field4_numeric %in% c(0, 1)) {
      # Old format: ascending
      ascending <- as.logical(field4_numeric)
      starting_defect_value <- NA
    } else {
      # New format: starting_defect_value
      starting_defect_value <- field4_numeric
      # Determine ascending based on whether starting value is min or max
      ascending <- (starting_defect_value == min(defects))
    }

    list(
      step = as.integer(parts[1]),
      accepts = list(parse_vector(parts[2])),
      delta_H_swap = list(parse_vector(parts[3])),
      ascending = ascending,
      starting_defect_value = starting_defect_value,
      defects = list(defects),
      prev_defects = list(parse_vector(parts[6]))
    )
  })

  # Combine list of lists into a data frame
  df <- dplyr::bind_rows(data_list)

  return(df)
}
