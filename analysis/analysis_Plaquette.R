library(ggplot2)
library(dplyr)
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory

# Simple per-run logger
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

analyze_plaquette <- function(directory, skip_steps = 0) {
  write_log(paste0("analyze_plaquette: start directory=", directory, " skip_steps=", skip_steps))

  # Read the plaquette data
  plaquette_data <- read_data_plaquette_filename(directory)
  write_log(paste0("analyze_plaquette: read ", nrow(plaquette_data), " rows"))
  print(head(plaquette_data))

  # Skip initial steps
  if (skip_steps > 0 && nrow(plaquette_data) > skip_steps) {
    analysis_data <- plaquette_data$plaquette[(skip_steps + 1):nrow(plaquette_data)]
  } else {
    analysis_data <- plaquette_data$plaquette
  }

  write_log(paste0("analyze_plaquette: analyzing ", length(analysis_data), " points after skip"))

  # 1. Autocorrelation analysis with uwerr (best practices)
  write_log("analyze_plaquette: computing uwerr for autocorrelation analysis")
  uw <- tryCatch(
    {
      hadron::uwerrprimary(analysis_data, pl = FALSE)
    },
    error = function(e) {
      write_log(paste0("ERROR in uwerrprimary: ", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(uw)) {
    tauint <- uw$tauint
    dtauint <- uw$dtauint
    mean_plaq <- uw$value
    error_plaq <- uw$dvalue
    n_eff <- length(analysis_data) / (2 * tauint)

    write_log(sprintf("analyze_plaquette: <P> = %.6f ± %.6f", mean_plaq, error_plaq))
    write_log(sprintf("analyze_plaquette: τ_int = %.2f ± %.2f", tauint, dtauint))
    write_log(sprintf("analyze_plaquette: n_eff = %.1f", n_eff))

    # Determine bootstrap parameters based on τ_int
    boot_R <- if (tauint > 10) 500 else 400
    boot_l <- max(1, floor(tauint / 2))

    write_log(sprintf("analyze_plaquette: using boot.R=%d, boot.l=%d based on τ_int", boot_R, boot_l))
  } else {
    # Fallback if uwerr fails
    write_log("analyze_plaquette: uwerr failed, using default bootstrap parameters")
    boot_R <- 400
    boot_l <- 2
    mean_plaq <- mean(analysis_data)
    error_plaq <- sd(analysis_data) / sqrt(length(analysis_data))
    tauint <- NA
    dtauint <- NA
    n_eff <- NA
  }

  # 2. Bootstrap analysis with plotting
  write_log("analyze_plaquette: running bootstrap.analysis")
  # Note: bootstrap.analysis with pl=TRUE creates its own plot output
  # We need to redirect it to a PDF file
  pdf_file <- file.path(directory, "bootstrap_plaquette.pdf")
  pdf(pdf_file, width = 10, height = 8)
  bootstrap_results <- hadron::bootstrap.analysis(
    plaquette_data$plaquette,
    skip = skip_steps,
    boot.R = boot_R,
    boot.l = boot_l,
    pl = TRUE
  )
  dev.off()
  write_log(paste0("analyze_plaquette: bootstrap plot saved to ", pdf_file))

  # 3. Save autocorrelation plot if uwerr succeeded
  if (!is.null(uw)) {
    write_log("analyze_plaquette: saving uwerr autocorrelation plot")
    pdf(file.path(directory, "plaquette_autocorr.pdf"), width = 8, height = 6)
    plot(uw)
    dev.off()
    write_log("analyze_plaquette: autocorrelation plot saved")
  }

  # 4. Save numerical results to file
  summary_file <- file.path(directory, "plaquette_summary.txt")
  write_log(paste0("analyze_plaquette: saving summary to ", summary_file))

  cat("Plaquette Analysis Results\n", file = summary_file)
  cat("===========================\n\n", file = summary_file, append = TRUE)
  cat(sprintf("Number of measurements: %d\n", nrow(plaquette_data)), file = summary_file, append = TRUE)
  cat(sprintf("Thermalization steps skipped: %d\n", skip_steps), file = summary_file, append = TRUE)
  cat(sprintf("Measurements analyzed: %d\n\n", length(analysis_data)), file = summary_file, append = TRUE)

  cat(sprintf("Mean plaquette: %.6f ± %.6f\n\n", mean_plaq, error_plaq), file = summary_file, append = TRUE)

  if (!is.na(tauint)) {
    cat("Autocorrelation Analysis\n", file = summary_file, append = TRUE)
    cat("------------------------\n", file = summary_file, append = TRUE)
    cat(sprintf("Integrated autocorrelation time: %.2f ± %.2f\n", tauint, dtauint), file = summary_file, append = TRUE)
    cat(sprintf("Effective number of independent measurements: %.1f\n", n_eff), file = summary_file, append = TRUE)
    cat(sprintf("Error enhancement factor: sqrt(2*τ_int) = %.2f\n\n", sqrt(2 * tauint)), file = summary_file, append = TRUE)
  }

  cat("Bootstrap Parameters\n", file = summary_file, append = TRUE)
  cat("--------------------\n", file = summary_file, append = TRUE)
  cat(sprintf("Bootstrap samples: %d\n", boot_R), file = summary_file, append = TRUE)
  cat(sprintf("Block length: %d\n", boot_l), file = summary_file, append = TRUE)

  write_log("analyze_plaquette: complete")

  return(list(
    mean = mean_plaq,
    error = error_plaq,
    tauint = tauint,
    n_eff = n_eff,
    bootstrap_results = bootstrap_results
  ))
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analysis_plaquette.R <directory> <skip_steps>")
}
directory <- args[1]
logs_dir <- file.path(directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
assign("WF_LOG_FILE", file.path(logs_dir, "analysis_Plaquette.log"), envir = .GlobalEnv)
skip_steps <- args[2]
analyze_plaquette(directory, skip_steps = as.integer(skip_steps))
write_log("=== Analysis completed successfully ===")
