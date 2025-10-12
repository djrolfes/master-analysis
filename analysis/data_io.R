read_data_file <- function(filepath) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  read.table(filepath, header=FALSE, comment.char="#")
}

# Function to read simulation log file
read_simulation_log <- function(filepath, skip_steps = 0) {
  filename <- config$SimulationLoggingParams$log_filename
  read_data_file(file.path(base_dir, filename))
}

read_data_plaquette_filename <- function(config, base_dir) {
  filename <- config$GaugeObservableParams$plaquette_filename
  read_data_file(file.path(base_dir, filename))
}

read_data_action_density_filename <- function(config, base_dir) {
  filename <- config$GaugeObservableParams$action_density_filename
  read_data_file(file.path(base_dir, filename))
}

# For files not present in example_data:
read_data_W_temp_filename <- function(config, base_dir) {
  # Example data should be provided; amend as needed for actual format
  filename <- config$GaugeObservableParams$W_temp_filename
  read_data_file(file.path(base_dir, filename))
}

library(yaml)

read_yaml_config <- function(yaml_path) {
  if (!file.exists(yaml_path)) stop(paste("YAML file not found:", yaml_path))
  yaml.load_file(yaml_path)
}
