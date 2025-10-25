library(ggplot2)
library(dplyr)
library(hadron)
source("data_io.R")  # relies on analysis/ being the working dir when called

# Simple per-run logger (compatible with other analysis_* scripts)
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

analyze_topological_charge <- function(directory, skip_initial = 0) {
  # Configure logfile for this run
  assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
  write_log(paste0("analyze_topological_charge: start directory=", directory, " skip_initial=", skip_initial))

  # Read YAML config to find filename
  cfg <- tryCatch({ read_yaml_config(directory) }, error = function(e) {
    write_log(paste0("ERROR reading YAML config: ", conditionMessage(e)))
    stop(e)
  })

  topo_filename <- cfg$GaugeObservableParams$topological_charge_filename
  if (is.null(topo_filename) || topo_filename == "") {
    write_log("ERROR: topological_charge_filename not found in YAML under GaugeObservableParams")
    stop("topological_charge_filename not found in YAML")
  }

  topo_path <- file.path(directory, topo_filename)
  if (!file.exists(topo_path)) {
    write_log(paste0("ERROR: topological charge file not found: ", topo_path))
    stop(paste0("Topological charge file not found: ", topo_path))
  }

  # Read the data file using the shared data_io reader
  df <- tryCatch({ read_data_file(topo_path) }, error = function(e) {
    write_log(paste0("ERROR reading topological charge file: ", conditionMessage(e)))
    stop(e)
  })

  write_log(paste0("analyze_topological_charge: read data with ", nrow(df), " rows and columns: ", paste(names(df), collapse = ", ")))

  # Detect column names: step and topological charge
  col_names <- tolower(names(df))
  step_col <- names(df)[which(grepl("step", col_names) | grepl("hmc_step", col_names))][1]
  charge_col_candidates <- names(df)[which(grepl("topo", col_names) | grepl("charge", col_names) | grepl("topological", col_names))]
  charge_col <- if(length(charge_col_candidates) >= 1) charge_col_candidates[1] else names(df)[2]

  if (is.na(step_col) || length(step_col) == 0) {
    write_log("WARNING: could not find 'step' column; using first column as step")
    step_col <- names(df)[1]
  }
  if (is.na(charge_col) || length(charge_col) == 0) {
    write_log("ERROR: could not find topological charge column")
    stop("Could not find topological charge column in file")
  }

  # Prepare a tidy data.frame with numeric values
  data <- df %>%
    transmute(
      step = as.integer(as.numeric(.data[[step_col]])),
      topo = as.numeric(.data[[charge_col]])
    ) %>%
    filter(!is.na(step) & !is.na(topo))

  write_log(paste0("analyze_topological_charge: prepared data rows = ", nrow(data)))

  # 1) Plot topological charge vs HMC step
  timeseries_plot <- ggplot(data, aes(x = step, y = topo)) +
    geom_line(color = "steelblue") +
    geom_point(size = 0.8, alpha = 0.7) +
    labs(title = "Topological Charge vs HMC Step", x = "HMC Step", y = "Topological Charge") +
    theme_minimal()

  out_ts_pdf <- file.path(directory, "topological_charge_timeseries.pdf")
  write_log(paste0("analyze_topological_charge: saving timeseries plot to ", out_ts_pdf))
  tryCatch({
    ggsave(out_ts_pdf, plot = timeseries_plot, width = 8, height = 4)
    write_log("analyze_topological_charge: timeseries plot saved")
  }, error = function(e) {
    write_log(paste0("ERROR saving timeseries plot: ", conditionMessage(e)))
  })

  # --- Bootstrap analysis of the topological charge (mean ± error) ---
  out_boot_pdf <- file.path(directory, "topological_charge_bootstrap.pdf")
  out_boot_txt <- file.path(directory, "topological_charge_bootstrap.txt")
  write_log(paste0("analyze_topological_charge: computing bootstrap.analysis and saving to ", out_boot_pdf, " and ", out_boot_txt))
  tryCatch({
    # Open PDF device so bootstrap.analysis can plot
    pdf(out_boot_pdf, width = 8, height = 6)
    # capture printed output from bootstrap.analysis into a character vector
    boot_output <- capture.output(hadron::bootstrap.analysis(data$topo, skip = skip_initial, pl = TRUE))
    dev.off()

    # Write textual output to .txt file
    writeLines(boot_output, con = out_boot_txt)
    write_log("analyze_topological_charge: bootstrap plot and text output saved")
  }, error = function(e) {
    write_log(paste0("analyze_topological_charge: ERROR running bootstrap.analysis: ", conditionMessage(e)))
  })

  # 2) Compute autocorrelation time using hadron::uwerr
  # Skip initial configs if requested
  if (skip_initial < 0) skip_initial <- 0
  if (skip_initial >= nrow(data)) {
    write_log(paste0("WARNING: skip_initial >= number of rows (", nrow(data), "); nothing to analyze for autocorr"))
    return(list(timeseries = timeseries_plot, autocorr = NULL))
  }

  ac_data <- data$topo[(skip_initial + 1):nrow(data)]
  ac_data <- as.numeric(ac_data)
  ac_data <- ac_data[!is.na(ac_data)]

  write_log(paste0("analyze_topological_charge: computing uwerr on ", length(ac_data), " samples after skipping ", skip_initial))
  if (length(ac_data) < 2) {
    write_log("Not enough data points to compute autocorrelation")
    return(list(timeseries = timeseries_plot, autocorr = NULL))
  }

  # uwerrprimary handles primary (1D) observables
  uw <- tryCatch({ hadron::uwerrprimary(ac_data, pl = TRUE) }, error = function(e) {
    write_log(paste0("ERROR computing uwerrprimary: ", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(uw)) {
    write_log("analyze_topological_charge: uwerr failed")
    return(list(timeseries = timeseries_plot, autocorr = NULL))
  }

  # Extract tauint and its error if available
  tauint <- NULL
  dtauint <- NULL
  if (!is.null(uw$tauint)) tauint <- uw$tauint
  if (!is.null(uw$dtauint)) dtauint <- uw$dtauint

  write_log(paste0("analyze_topological_charge: uwerr results - tauint=", tauint, ", dtauint=", dtauint))

  # Save the uwerr plot to PDF (plot method for uwerr)
  out_ac_pdf <- file.path(directory, "topological_charge_autocorr.pdf")
  write_log(paste0("analyze_topological_charge: saving autocorr plot to ", out_ac_pdf))
  tryCatch({
    pdf(out_ac_pdf, width = 8, height = 6)
    plot(uw)
    dev.off()
    write_log("analyze_topological_charge: autocorr plot saved")
  }, error = function(e) {
    write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
  })

  # Save summary to a small text file
  summary_file <- file.path(directory, "topological_charge_autocorr_summary.txt")
  write_log(paste0("analyze_topological_charge: writing summary to ", summary_file))
  tryCatch({
    cat(sprintf("n_points_after_skip: %d\n", length(ac_data)), file = summary_file)
    cat(sprintf("tauint: %s\n", as.character(tauint)), file = summary_file, append = TRUE)
    cat(sprintf("dtauint: %s\n", as.character(dtauint)), file = summary_file, append = TRUE)
    write_log("analyze_topological_charge: summary written")
  }, error = function(e) {
    write_log(paste0("ERROR writing summary file: ", conditionMessage(e)))
  })

  return(list(timeseries = timeseries_plot, autocorr = uw, tauint = tauint, dtauint = dtauint))
}

# CLI entrypoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analysis_topological_charge.R <directory> [skip_initial]")
}

directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_initial <- if (length(args) >= 2) as.integer(args[2]) else 0
analyze_topological_charge(directory, skip_initial = skip_initial)
