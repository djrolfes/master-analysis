library(ggplot2)
library(dplyr)
library(hadron)
source("data_io.R") # relies on analysis/ being the working dir when called

# simple per-run logger (compatible with other analysis_* scripts)
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

analyze_wilsonflow_improv <- function(directory, filename = "action_densities_clover_cumulative.txt", n_boot = 400) {
  logs_dir <- file.path(directory, "logs")
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
  assign("WF_LOG_FILE", file.path(logs_dir, "analysis_wilsonflow_improv.log"), envir = .GlobalEnv)
  write_log(paste0("analyze_wilsonflow_improv: start dir=", directory, " file=", filename, " n_boot=", n_boot))

  path <- file.path(directory, filename)
  if (!file.exists(path)) {
    write_log(paste0("analyze_wilsonflow_improv: file not found: ", path))
    stop("Input file not found: ", path)
  }

  # Read CSV produced by wilsonflow_improv_testing (comma separated)
  df <- tryCatch(
    {
      read.csv(path, header = TRUE, stringsAsFactors = FALSE)
    },
    error = function(e) {
      write_log(paste0("analyze_wilsonflow_improv: ERROR reading file: ", conditionMessage(e)))
      stop(e)
    }
  )

  write_log(paste0("analyze_wilsonflow_improv: read ", nrow(df), " rows and ", ncol(df), " cols"))

  # Identify columns
  col_names <- names(df)
  # Expect hmc_step and a column like W_normal (or similar) plus factor columns (e.g. "1","2","4"...)
  normal_col_candidates <- grep("W[_ ]?normal", col_names, ignore.case = TRUE, value = TRUE)
  if (length(normal_col_candidates) >= 1) {
    normal_col <- normal_col_candidates[1]
  } else if (length(col_names) >= 2) {
    # fallback to second column if no clear "normal" column found
    normal_col <- col_names[2]
    write_log(paste0("analyze_wilsonflow_improv: WARNING - didn't find explicit W_normal column; using ", normal_col))
  } else {
    write_log("analyze_wilsonflow_improv: ERROR - unexpected file format")
    stop("Unexpected file format; cannot find normal Wilson flow column")
  }

  # Factor columns are all except hmc_step and the normal column
  factor_cols <- setdiff(col_names, c("hmc_step", normal_col))
  if (length(factor_cols) == 0) {
    write_log("analyze_wilsonflow_improv: ERROR - no improved flow columns found")
    stop("No improved flow columns found")
  }

  write_log(paste0("analyze_wilsonflow_improv: normal_col=", normal_col, " factor_cols=", paste(factor_cols, collapse = ",")))

  # Compute absolute difference per factor and obtain bootstrap mean±error
  results <- lapply(factor_cols, function(fc) {
    # convert column to numeric safely
    v_norm <- as.numeric(df[[normal_col]])
    v_imp <- as.numeric(df[[fc]])
    absdiff <- abs(v_imp - v_norm)
    # Remove NA rows
    absdiff <- absdiff[!is.na(absdiff)]
    if (length(absdiff) == 0) {
      write_log(paste0("analyze_wilsonflow_improv: WARNING - no valid values for factor ", fc))
      br <- c(mean = NA, error = NA)
    } else {
      br <- hadron::bootstrap.meanerror(absdiff, n_boot)
      # bootstrap.meanerror may return a named vector or scalar; ensure we extract mean/error
      if (!is.null(names(br)) && all(c("mean", "error") %in% names(br))) {
        br <- c(mean = as.numeric(br["mean"]), error = as.numeric(br["error"]))
      } else if (length(br) == 2) {
        br <- as.numeric(br)
        names(br) <- c("mean", "error")
      } else {
        # fallback: compute mean and sd/sqrt(N) as naive error
        br <- c(mean = mean(absdiff), error = sd(absdiff) / sqrt(length(absdiff)))
      }
    }
    # factor value parsed from header name if numeric
    fac_val <- suppressWarnings(as.numeric(fc))
    if (is.na(fac_val)) fac_val <- fc
    data.frame(factor = fac_val, mean = br["mean"], error = br["error"], stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    arrange(as.numeric(factor))

  write_log(paste0("analyze_wilsonflow_improv: aggregated ", nrow(results), " factor entries"))

  # Plot mean ± error vs factor
  p <- ggplot(results, aes(x = as.numeric(factor), y = as.numeric(mean))) +
    geom_point() +
    geom_errorbar(aes(ymin = as.numeric(mean) - as.numeric(error), ymax = as.numeric(mean) + as.numeric(error)), width = 0.05) +
    labs(
      title = "Absolute difference (improved - normal) of action density",
      x = "Improvement factor (multiplied to \u03b5) [dimensionless]",
      y = "Mean |E_improved - E_normal| \u00b1 error [lattice units]"
    ) +
    theme_minimal()

  out_pdf <- file.path(directory, "wilsonflow_improv_absdiff_bootstrap.pdf")
  write_log(paste0("analyze_wilsonflow_improv: saving plot to ", out_pdf))
  tryCatch(
    {
      ggsave(out_pdf, plot = p, width = 7, height = 5)
      write_log(paste0("analyze_wilsonflow_improv: saved ", out_pdf))
    },
    error = function(e) {
      write_log(paste0("analyze_wilsonflow_improv: ERROR saving plot: ", conditionMessage(e)))
    }
  )

  return(list(data = results, plot = p))
}

# --- CLI entrypoint ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analysis_wilsonflow_improv.R <directory> [n_boot]")
}
directory <- args[1]
logs_dir <- file.path(directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
assign("WF_LOG_FILE", file.path(logs_dir, "analysis_wilsonflow_improv.log"), envir = .GlobalEnv)
n_boot <- if (length(args) >= 2) as.integer(args[2]) else 200L
analyze_wilsonflow_improv(directory, n_boot = n_boot)
write_log("=== Analysis completed successfully ===")
