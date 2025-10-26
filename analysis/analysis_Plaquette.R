library(ggplot2)
library(dplyr)
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory



analyze_plaquette <- function(directory, skip_steps = 0) {
  # Read the plaquette data
  plaquette_data <- read_data_plaquette_filename(directory)
  print(head(plaquette_data))
  pdf(file.path(directory, "bootstrap_plaquette.pdf"))
  bootstrap_results <- hadron::bootstrap.analysis(plaquette_data$plaquette, skip = skip_steps, pl = TRUE)
  dev.off()
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analysis_plaquette.R <directory> <skip_steps>")
}
directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_steps <- args[2]
analyze_plaquette(directory, skip_steps = as.integer(skip_steps))
