library(ggplot2)
library(dplyr)
library(hadron)
source("data_io.R") # relies on analysis/ being the working dir when called

# Simple per-run logger (compatible with other analysis_* scripts)
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

# Safe wrapper for computeacf that handles BLAS/CPU incompatibility
# Note: On some clusters, this function will cause R to abort due to BLAS incompatibility
# The error cannot be caught with tryCatch as it occurs at the C library level
safe_computeacf <- function(data, W.max, label = "data") {
  # Check if we should skip computeacf entirely (set this env var to disable)
  if (Sys.getenv("SKIP_COMPUTEACF", "0") == "1") {
    write_log(paste0("SKIPPED computeacf for ", label, ": SKIP_COMPUTEACF=1"))
    return(NULL)
  }

  result <- tryCatch(
    {
      hadron::computeacf(data, W.max)
    },
    error = function(e) {
      # Check if it's the BLAS/illegal instruction error
      err_msg <- conditionMessage(e)
      if (grepl("illegal", err_msg, ignore.case = TRUE) ||
        grepl("operand", err_msg, ignore.case = TRUE)) {
        write_log(paste0("SKIPPED computeacf for ", label, ": BLAS/CPU incompatibility detected (", err_msg, ")"))
        write_log(paste0("  This is a known issue on some clusters. Using uwerrprimary and bootstrap.analysis instead."))
      } else {
        write_log(paste0("ERROR computing computeacf for ", label, ": ", err_msg))
      }
      return(NULL)
    }
  )
  return(result)
}

analyze_topological_charge <- function(directory, skip_initial = 0, s = 2.5, s_squared = 5.0, r = 1,
                                       xlim_min = NULL, xlim_max = 120000) {
  # Apply parameters to all analyses
  s_q <- s
  s_q_rounded <- s
  s_q_squared <- s_squared
  s_q_squared_rounded <- s_squared
  r_q <- r
  r_q_rounded <- r
  r_q_squared <- r
  r_q_squared_rounded <- r

  # Configure logfile for this run
  write_log(paste0("analyze_topological_charge: start directory=", directory, " skip_initial=", skip_initial, " s=", s, " s_squared=", s_squared, " r=", r))
  output_dir <- file.path(directory, "topological_charge_scan")
  logs_dir <- file.path(output_dir, "logs")
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  param_suffix <- sprintf("_s%.2f_s2_%.2f_r%d", s, s_squared, r)
  assign("WF_LOG_FILE", file.path(logs_dir, paste0("analysis_topological_charge", param_suffix, ".log")), envir = .GlobalEnv)
  write_log(paste0("Using custom parameters: outputs in ", output_dir))


  # Read YAML config to find filename
  cfg <- tryCatch(
    {
      read_yaml_config(directory)
    },
    error = function(e) {
      write_log(paste0("ERROR reading YAML config: ", conditionMessage(e)))
      stop(e)
    }
  )

  topo_filename <- cfg$GaugeObservableParams$topological_charge_filename
  if (is.null(topo_filename) || topo_filename == "") {
    write_log("ERROR: topological_charge_filename not found in YAML under GaugeObservableParams")
    stop("topological_charge_filename not found in YAML")
  }

  # Extract measurement_interval for HMC normalization
  measurement_interval <- cfg$GaugeObservableParams$measurement_interval
  if (is.null(measurement_interval)) {
    write_log("WARNING: measurement_interval not found in YAML; using default value of 1")
    measurement_interval <- 1
  }

  # Extract number of replicas for PTBC normalization (if available)
  n_replicas <- NA
  if (!is.null(cfg$PTBCParams) && !is.null(cfg$PTBCParams$defect_values)) {
    n_replicas <- length(cfg$PTBCParams$defect_values)
    write_log(paste0("analyze_topological_charge: found ", n_replicas, " replicas in PTBCParams"))
  } else {
    write_log("analyze_topological_charge: PTBCParams not found; PTBC normalization will not be computed")
  }

  topo_path <- file.path(directory, topo_filename)
  if (!file.exists(topo_path)) {
    write_log(paste0("ERROR: topological charge file not found: ", topo_path))
    stop(paste0("Topological charge file not found: ", topo_path))
  }

  # Read the data file using the shared data_io reader
  df <- tryCatch(
    {
      read_data_file(topo_path)
    },
    error = function(e) {
      write_log(paste0("ERROR reading topological charge file: ", conditionMessage(e)))
      stop(e)
    }
  )

  write_log(paste0("analyze_topological_charge: read data with ", nrow(df), " rows and columns: ", paste(names(df), collapse = ", ")))

  # Detect column names: step and topological charge
  col_names <- tolower(names(df))
  step_col <- names(df)[which(grepl("step", col_names) | grepl("hmc_step", col_names))][1]
  charge_col_candidates <- names(df)[which(grepl("topo", col_names) | grepl("charge", col_names) | grepl("topological", col_names))]
  charge_col <- if (length(charge_col_candidates) >= 1) charge_col_candidates[1] else names(df)[2]

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

  # 1) Plot topological charge vs HMC step (base R with hadron style)
  out_ts_pdf <- file.path(output_dir, paste0("topological_charge_timeseries", param_suffix, ".pdf"))
  write_log(paste0("analyze_topological_charge: saving timeseries plot to ", out_ts_pdf))
  timeseries_plot <- NULL
  tryCatch(
    {
      pdf(out_ts_pdf, width = 8, height = 4)
      par(bty = "o") # Complete frame around plot
      plot(data$step, data$topo,
        type = "l", col = "steelblue", lwd = 1.5,
        xlab = "HMC Step",
        ylab = "Topological Charge Q",
        # main = "Topological Charge vs HMC Step"
      )
      points(data$step, data$topo, pch = 19, cex = 0.6, col = rgb(0, 0, 0, alpha = 0.7))
      abline(v = skip_initial, lty = 2, col = "red", lwd = 0.8)
      y_range <- max(data$topo, na.rm = TRUE) - min(data$topo, na.rm = TRUE)
      text(skip_initial, max(data$topo, na.rm = TRUE) - 0.05 * y_range,
        labels = sprintf("skip=%d", skip_initial),
        pos = 4, col = "red", cex = 0.9
      )
      timeseries_plot <- recordPlot()
      dev.off()
      write_log("analyze_topological_charge: timeseries saved")
    },
    error = function(e) {
      write_log(paste0("ERROR saving timeseries plot: ", conditionMessage(e)))
    }
  )

  # 1b) Plot topological charge vs PTBC-normalized steps (HMC step * n_replicas)
  if (!is.na(n_replicas)) {
    out_ts_ptbc_pdf <- file.path(output_dir, paste0("topological_charge_timeseries_ptbc", param_suffix, ".pdf"))
    write_log(paste0("analyze_topological_charge: saving PTBC-normalized timeseries plot to ", out_ts_ptbc_pdf))

    # Calculate PTBC-normalized steps
    data_ptbc <- data
    data_ptbc$step_ptbc <- data$step * n_replicas

    # Determine x-axis limits (use provided values or defaults to data range)
    xlim_min_ptbc <- if (is.null(xlim_min)) min(data_ptbc$step_ptbc, na.rm = TRUE) else xlim_min
    xlim_max_ptbc <- if (is.null(xlim_max)) max(data_ptbc$step_ptbc, na.rm = TRUE) else xlim_max

    tryCatch(
      {
        pdf(out_ts_ptbc_pdf, width = 8, height = 4)
        par(bty = "o") # Complete frame around plot
        plot(data_ptbc$step_ptbc, data_ptbc$topo,
          type = "l", col = "steelblue", lwd = 1.5,
          xlab = "PTBC Step (HMC Step × n_replicas)",
          ylab = "Topological Charge Q",
          xlim = c(xlim_min_ptbc, xlim_max_ptbc)
        )
        points(data_ptbc$step_ptbc, data_ptbc$topo, pch = 19, cex = 0.6, col = rgb(0, 0, 0, alpha = 0.7))
        skip_ptbc <- skip_initial * n_replicas
        abline(v = skip_ptbc, lty = 2, col = "red", lwd = 0.8)
        y_range <- max(data_ptbc$topo, na.rm = TRUE) - min(data_ptbc$topo, na.rm = TRUE)
        text(skip_ptbc, max(data_ptbc$topo, na.rm = TRUE) - 0.05 * y_range,
          labels = sprintf("skip=%d", skip_ptbc),
          pos = 4, col = "red", cex = 0.9
        )
        dev.off()
        write_log("analyze_topological_charge: PTBC-normalized timeseries saved")
      },
      error = function(e) {
        write_log(paste0("ERROR saving PTBC-normalized timeseries plot: ", conditionMessage(e)))
      }
    )
  }

  # --- Bootstrap analysis of the topological charge (mean ± error) ---
  # Filter data first (bootstrap.analysis doesn't have a skip parameter)
  data_for_boot <- data$topo[data$step > skip_initial]

  out_boot_pdf <- file.path(output_dir, paste0("topological_charge_bootstrap", param_suffix, ".pdf"))
  out_boot_txt <- file.path(output_dir, paste0("topological_charge_bootstrap", param_suffix, ".txt"))
  write_log(paste0("analyze_topological_charge: computing bootstrap.analysis on ", length(data_for_boot), " points and saving to ", out_boot_pdf))

  boot_original <- NULL
  tryCatch(
    {
      # Open PDF device so bootstrap.analysis can plot
      pdf(out_boot_pdf, width = 8, height = 6)
      # capture printed output from bootstrap.analysis into a character vector
      boot_output_text <- capture.output({
        boot_original <- hadron::bootstrap.analysis(data_for_boot, pl = TRUE)
      })
      dev.off()
      write_log("analyze_topological_charge: bootstrap plot output saved")
    },
    error = function(e) {
      write_log(paste0("analyze_topological_charge: ERROR running bootstrap.analysis: ", conditionMessage(e)))
    }
  )

  # 2) Compute autocorrelation time using hadron::uwerr
  # Skip initial configs if requested
  if (skip_initial < 0) skip_initial <- 0

  # Filter data first, then check if we have enough points
  ac_data <- data$topo[data$step > skip_initial]
  ac_data <- as.numeric(ac_data)
  ac_data <- ac_data[!is.na(ac_data)]

  write_log(paste0("analyze_topological_charge: computing uwerr on ", length(ac_data), " samples after skipping steps <= ", skip_initial))

  if (length(ac_data) < 2) {
    write_log(paste0("WARNING: Not enough data points (", length(ac_data), ") to compute autocorrelation after skipping"))
    return(list(timeseries = timeseries_plot, autocorr = NULL))
  }

  # uwerrprimary handles primary (1D) observables
  uw <- tryCatch(
    {
      lreps <- as.integer(length(ac_data) / r_q)
      diff <- length(ac_data) - (r_q * lreps)
      n_reps <- c(rep(lreps, r_q))
      n_reps[1] <- n_reps[1] + diff # add remainder to last rep
      hadron::uwerrprimary(ac_data, n_reps, pl = TRUE, S = s_q)
    },
    error = function(e) {
      write_log(paste0("ERROR computing uwerrprimary: ", conditionMessage(e)))
      return(NULL)
    }
  )

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

  # Save the uwerr plot to PDF (plot method for uwerr)autocorr.pdf
  autocorr_filename <- paste0("topological_charge_autocorr_Q", param_suffix, ".pdf")
  out_ac_pdf <- file.path(output_dir, autocorr_filename)
  write_log(paste0("analyze_topological_charge: saving autocorr plot to ", out_ac_pdf))
  tryCatch(
    {
      pdf(out_ac_pdf, width = 8, height = 6)
      plot(uw)
      dev.off()
      write_log("analyze_topological_charge: autocorr plot saved")
    },
    error = function(e) {
      write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
    }
  )

  # Save summary to a small text file
  summary_file <- file.path(output_dir, paste0("topological_charge_autocorr_summary", param_suffix, ".txt"))
  # Note: autocorr_filename already set above with param_suffix
  write_log(paste0("analyze_topological_charge: writing summary to ", summary_file))

  # Compute ACF using computeacf (with safe wrapper for cluster compatibility)
  acf_result <- safe_computeacf(ac_data, floor(length(ac_data) / 5), "original data")

  acf_tau <- if (!is.null(acf_result)) acf_result$tau else NA
  acf_dtau <- if (!is.null(acf_result)) acf_result$dtau else NA

  # Extract tau and its error from bootstrap.analysis
  # Use the row with maximum Tauint and propagate error from DError
  boot_tau <- NA
  boot_dtau <- NA
  if (!is.null(boot_original)) {
    max_idx <- which.max(boot_original$Tauint)
    boot_tau <- boot_original$Tauint[max_idx]
    # Propagate error: Tauint = Error^2 / (error.naive^2 * 2)
    # dTauint/dError = 2*Error / (error.naive^2 * 2) = Error / error.naive^2
    error_naive <- sd(ac_data) / sqrt(length(ac_data))
    boot_dtau <- 2 * boot_original$Error[max_idx] * boot_original$DError[max_idx] / (error_naive^2)
  }

  if (!is.null(acf_result)) {
    write_log(paste0("analyze_topological_charge: computeacf results - tau=", acf_tau, ", dtau=", acf_dtau))
  }

  # Write summary for original data
  tryCatch(
    {
      cat(sprintf("n_points_after_skip: %d\n", length(ac_data)), file = summary_file)
      cat("\n--- Original Data ---\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint: %s\n", as.character(tauint)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint: %s\n", as.character(dtauint)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau: %s\n", as.character(acf_tau)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau: %s\n", as.character(acf_dtau)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau: %s\n", as.character(boot_tau)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau: %s\n", as.character(boot_dtau)), file = summary_file, append = TRUE)

      # HMC normalized (multiply by measurement_interval)
      cat("\nHMC normalized:\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint_hmc: %s\n", as.character(if (!is.na(tauint)) tauint * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint_hmc: %s\n", as.character(if (!is.na(dtauint)) dtauint * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau_hmc: %s\n", as.character(if (!is.na(acf_tau)) acf_tau * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau_hmc: %s\n", as.character(if (!is.na(acf_dtau)) acf_dtau * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau_hmc: %s\n", as.character(if (!is.na(boot_tau)) boot_tau * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau_hmc: %s\n", as.character(if (!is.na(boot_dtau)) boot_dtau * measurement_interval else NA)), file = summary_file, append = TRUE)

      # PTBC normalized (multiply by measurement_interval and n_replicas if available)
      if (!is.na(n_replicas)) {
        cat("\nPTBC normalized:\n", file = summary_file, append = TRUE)
        cat(sprintf("uwerr_tauint_ptbc: %s\n", as.character(if (!is.na(tauint)) tauint * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("uwerr_dtauint_ptbc: %s\n", as.character(if (!is.na(dtauint)) dtauint * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_tau_ptbc: %s\n", as.character(if (!is.na(acf_tau)) acf_tau * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_dtau_ptbc: %s\n", as.character(if (!is.na(acf_dtau)) acf_dtau * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis tau_ptbc: %s\n", as.character(if (!is.na(boot_tau)) boot_tau * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis dtau_ptbc: %s\n", as.character(if (!is.na(boot_dtau)) boot_dtau * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
      }

      write_log("analyze_topological_charge: summary for original data written")
    },
    error = function(e) {
      write_log(paste0("ERROR writing summary file for original data: ", conditionMessage(e)))
    }
  )

  # --- Analysis on rounded data ---
  write_log("analyze_topological_charge: starting analysis on rounded data")
  ac_data_rounded <- round(ac_data)

  boot_rounded <- tryCatch(
    {
      hadron::bootstrap.analysis(ac_data_rounded, pl = FALSE)
    },
    error = function(e) {
      return(NULL)
    }
  )
  # uwerr on rounded data
  uw_rounded <- tryCatch(
    {
      lreps <- as.integer(length(ac_data_rounded) / r_q_rounded)
      diff <- length(ac_data_rounded) - (r_q_rounded * lreps)
      n_reps <- c(rep(lreps, r_q_rounded))
      n_reps[1] <- n_reps[1] + diff # add remainder to last rep
      hadron::uwerrprimary(ac_data_rounded, n_reps, pl = FALSE, S = s_q_rounded)
    },
    error = function(e) {
      write_log(paste0("ERROR computing uwerrprimary on rounded data: ", conditionMessage(e)))
      return(NULL)
    }
  )


  # Save the uwerr plot to PDF (plot method for uwerr)autocorr.pdf
  autocorr_filename <- paste0("topological_charge_autocorr_Q_rounded", param_suffix, ".pdf")
  out_ac_pdf <- file.path(output_dir, autocorr_filename)
  write_log(paste0("analyze_topological_charge: saving autocorr plot to ", out_ac_pdf))
  tryCatch(
    {
      pdf(out_ac_pdf, width = 8, height = 6)
      plot(uw_rounded)
      dev.off()
      write_log("analyze_topological_charge Q rounded: autocorr plot saved")
    },
    error = function(e) {
      write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
    }
  )

  tauint_rounded <- if (!is.null(uw_rounded)) uw_rounded$tauint else NA
  dtauint_rounded <- if (!is.null(uw_rounded)) uw_rounded$dtauint else NA
  write_log(paste0("analyze_topological_charge: rounded uwerr results - tauint=", tauint_rounded, ", dtauint=", dtauint_rounded))

  # computeacf on rounded data
  acf_result_rounded <- safe_computeacf(ac_data_rounded, floor(length(ac_data_rounded) / 5), "rounded data")

  acf_tau_rounded <- if (!is.null(acf_result_rounded)) acf_result_rounded$tau else NA
  acf_dtau_rounded <- if (!is.null(acf_result_rounded)) acf_result_rounded$dtau else NA

  # Extract tau and its error from bootstrap.analysis for rounded data
  boot_tau_rounded <- NA
  boot_dtau_rounded <- NA
  if (!is.null(boot_rounded)) {
    max_idx <- which.max(boot_rounded$Tauint)
    boot_tau_rounded <- boot_rounded$Tauint[max_idx]
    error_naive_rounded <- sd(ac_data_rounded) / sqrt(length(ac_data_rounded))
    boot_dtau_rounded <- 2 * boot_rounded$Error[max_idx] * boot_rounded$DError[max_idx] / (error_naive_rounded^2)
  }

  if (!is.null(acf_result_rounded)) {
    write_log(paste0("analyze_topological_charge: rounded computeacf results - tau=", acf_tau_rounded, ", dtau=", acf_dtau_rounded))
  }

  # Append rounded results to summary file
  tryCatch(
    {
      cat("\n\n--- Rounded Data ---\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint: %s\n", as.character(tauint_rounded)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint: %s\n", as.character(dtauint_rounded)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau: %s\n", as.character(acf_tau_rounded)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau: %s\n", as.character(acf_dtau_rounded)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau: %s\n", as.character(boot_tau_rounded)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau: %s\n", as.character(boot_dtau_rounded)), file = summary_file, append = TRUE)

      # HMC normalized
      cat("\nHMC normalized:\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint_hmc: %s\n", as.character(if (!is.na(tauint_rounded)) tauint_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint_hmc: %s\n", as.character(if (!is.na(dtauint_rounded)) dtauint_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau_hmc: %s\n", as.character(if (!is.na(acf_tau_rounded)) acf_tau_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau_hmc: %s\n", as.character(if (!is.na(acf_dtau_rounded)) acf_dtau_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau_hmc: %s\n", as.character(if (!is.na(boot_tau_rounded)) boot_tau_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau_hmc: %s\n", as.character(if (!is.na(boot_dtau_rounded)) boot_dtau_rounded * measurement_interval else NA)), file = summary_file, append = TRUE)

      # PTBC normalized
      if (!is.na(n_replicas)) {
        cat("\nPTBC normalized:\n", file = summary_file, append = TRUE)
        cat(sprintf("uwerr_tauint_ptbc: %s\n", as.character(if (!is.na(tauint_rounded)) tauint_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("uwerr_dtauint_ptbc: %s\n", as.character(if (!is.na(dtauint_rounded)) dtauint_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_tau_ptbc: %s\n", as.character(if (!is.na(acf_tau_rounded)) acf_tau_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_dtau_ptbc: %s\n", as.character(if (!is.na(acf_dtau_rounded)) acf_dtau_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis tau_ptbc: %s\n", as.character(if (!is.na(boot_tau_rounded)) boot_tau_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis dtau_ptbc: %s\n", as.character(if (!is.na(boot_dtau_rounded)) boot_dtau_rounded * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
      }

      write_log("analyze_topological_charge: summary for rounded data written")
    },
    error = function(e) {
      write_log(paste0("ERROR writing summary file for rounded data: ", conditionMessage(e)))
    }
  )

  # --- Analysis on Q^2 data ---
  write_log("analyze_topological_charge: starting analysis on Q^2 data")
  ac_data_q_squared <- ac_data^2

  boot_q_squared <- tryCatch(
    {
      hadron::bootstrap.analysis(ac_data_q_squared, pl = FALSE)
    },
    error = function(e) {
      return(NULL)
    }
  )

  # uwerr on Q^2 data
  uw_q_squared <- tryCatch(
    {
      lreps <- as.integer(length(ac_data_q_squared) / r_q_squared)
      diff <- length(ac_data_q_squared) - (r_q_squared * lreps)
      n_reps <- c(rep(lreps, r_q_squared))
      n_reps[1] <- n_reps[1] + diff # add remainder to last rep
      hadron::uwerrprimary(ac_data_q_squared, n_reps, pl = FALSE, S = s_q_squared)
    },
    error = function(e) {
      write_log(paste0("ERROR computing uwerrprimary on Q^2 data: ", conditionMessage(e)))
      return(NULL)
    }
  )

  # Save the uwerr plot to PDF (plot method for uwerr)autocorr.pdf
  autocorr_filename <- paste0("topological_charge_autocorr_Q_squared", param_suffix, ".pdf")
  out_ac_pdf <- file.path(output_dir, autocorr_filename)
  write_log(paste0("analyze_topological_charge: saving autocorr plot to ", out_ac_pdf))
  tryCatch(
    {
      pdf(out_ac_pdf, width = 8, height = 6)
      plot(uw_q_squared)
      dev.off()
      write_log("analyze_topological_charge Q squared: autocorr plot saved")
    },
    error = function(e) {
      write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
    }
  )

  tauint_q_squared <- if (!is.null(uw_q_squared)) uw_q_squared$tauint else NA
  dtauint_q_squared <- if (!is.null(uw_q_squared)) uw_q_squared$dtauint else NA
  write_log(paste0("analyze_topological_charge: Q^2 uwerr results - tauint=", tauint_q_squared, ", dtauint=", dtauint_q_squared))

  # computeacf on Q^2 data
  acf_result_q_squared <- safe_computeacf(ac_data_q_squared, floor(length(ac_data_q_squared) / 5), "Q^2 data")

  acf_tau_q_squared <- if (!is.null(acf_result_q_squared)) acf_result_q_squared$tau else NA
  acf_dtau_q_squared <- if (!is.null(acf_result_q_squared)) acf_result_q_squared$dtau else NA

  # Extract tau and its error from bootstrap.analysis for Q^2 data
  boot_tau_q_squared <- NA
  boot_dtau_q_squared <- NA
  if (!is.null(boot_q_squared)) {
    max_idx <- which.max(boot_q_squared$Tauint)
    boot_tau_q_squared <- boot_q_squared$Tauint[max_idx]
    error_naive_q_squared <- sd(ac_data_q_squared) / sqrt(length(ac_data_q_squared))
    boot_dtau_q_squared <- 2 * boot_q_squared$Error[max_idx] * boot_q_squared$DError[max_idx] / (error_naive_q_squared^2)
  }

  if (!is.null(acf_result_q_squared)) {
    write_log(paste0("analyze_topological_charge: Q^2 computeacf results - tau=", acf_tau_q_squared, ", dtau=", acf_dtau_q_squared))
  }

  # Append Q^2 results to summary file
  tryCatch(
    {
      cat("\n\n--- Q^2 Data ---\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint: %s\n", as.character(tauint_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint: %s\n", as.character(dtauint_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau: %s\n", as.character(acf_tau_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau: %s\n", as.character(acf_dtau_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau: %s\n", as.character(boot_tau_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau: %s\n", as.character(boot_dtau_q_squared)), file = summary_file, append = TRUE)

      # HMC normalized
      cat("\nHMC normalized:\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint_hmc: %s\n", as.character(if (!is.na(tauint_q_squared)) tauint_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint_hmc: %s\n", as.character(if (!is.na(dtauint_q_squared)) dtauint_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau_hmc: %s\n", as.character(if (!is.na(acf_tau_q_squared)) acf_tau_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau_hmc: %s\n", as.character(if (!is.na(acf_dtau_q_squared)) acf_dtau_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau_hmc: %s\n", as.character(if (!is.na(boot_tau_q_squared)) boot_tau_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau_hmc: %s\n", as.character(if (!is.na(boot_dtau_q_squared)) boot_dtau_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)

      # PTBC normalized
      if (!is.na(n_replicas)) {
        cat("\nPTBC normalized:\n", file = summary_file, append = TRUE)
        cat(sprintf("uwerr_tauint_ptbc: %s\n", as.character(if (!is.na(tauint_q_squared)) tauint_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("uwerr_dtauint_ptbc: %s\n", as.character(if (!is.na(dtauint_q_squared)) dtauint_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_tau_ptbc: %s\n", as.character(if (!is.na(acf_tau_q_squared)) acf_tau_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_dtau_ptbc: %s\n", as.character(if (!is.na(acf_dtau_q_squared)) acf_dtau_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis tau_ptbc: %s\n", as.character(if (!is.na(boot_tau_q_squared)) boot_tau_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis dtau_ptbc: %s\n", as.character(if (!is.na(boot_dtau_q_squared)) boot_dtau_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
      }

      write_log("analyze_topological_charge: summary for Q^2 data written")
    },
    error = function(e) {
      write_log(paste0("ERROR writing summary file for Q^2 data: ", conditionMessage(e)))
    }
  )

  # --- Analysis on Q^2 from rounded data ---
  write_log("analyze_topological_charge: starting analysis on Q^2 from rounded data")
  ac_data_rounded_q_squared <- ac_data_rounded^2

  boot_rounded_q_squared <- tryCatch(
    {
      hadron::bootstrap.analysis(ac_data_rounded_q_squared, pl = FALSE)
    },
    error = function(e) {
      return(NULL)
    }
  )

  # uwerr on Q^2 from rounded data
  uw_rounded_q_squared <- tryCatch(
    {
      lreps <- as.integer(length(ac_data_rounded_q_squared) / r_q_squared_rounded)
      diff <- length(ac_data_rounded_q_squared) - (r_q_squared_rounded * lreps)
      n_reps <- c(rep(lreps, r_q_squared_rounded))
      n_reps[1] <- n_reps[1] + diff # add remainder to last rep
      hadron::uwerrprimary(ac_data_rounded_q_squared, n_reps, pl = FALSE, S = s_q_squared_rounded)
    },
    error = function(e) {
      write_log(paste0("ERROR computing uwerrprimary on Q^2 from rounded data: ", conditionMessage(e)))
      return(NULL)
    }
  )

  # Save the uwerr plot to PDF (plot method for uwerr)autocorr.pdf
  autocorr_filename <- paste0("topological_charge_autocorr_rounded_Q_squared", param_suffix, ".pdf")
  out_ac_pdf <- file.path(output_dir, autocorr_filename)
  write_log(paste0("analyze_topological_charge: saving autocorr plot to ", out_ac_pdf))
  tryCatch(
    {
      pdf(out_ac_pdf, width = 8, height = 6)
      plot(uw_rounded_q_squared)
      dev.off()
      write_log("analyze_topological_charge rounded Q squared: autocorr plot saved")
    },
    error = function(e) {
      write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
    }
  )


  tauint_rounded_q_squared <- if (!is.null(uw_rounded_q_squared)) uw_rounded_q_squared$tauint else NA
  dtauint_rounded_q_squared <- if (!is.null(uw_rounded_q_squared)) uw_rounded_q_squared$dtauint else NA
  write_log(paste0("analyze_topological_charge: Q^2 from rounded uwerr results - tauint=", tauint_rounded_q_squared, ", dtauint=", dtauint_rounded_q_squared))

  # computeacf on Q^2 from rounded data
  acf_result_rounded_q_squared <- safe_computeacf(ac_data_rounded_q_squared, floor(length(ac_data_rounded_q_squared) / 5), "rounded Q^2 data")

  acf_tau_rounded_q_squared <- if (!is.null(acf_result_rounded_q_squared)) acf_result_rounded_q_squared$tau else NA
  acf_dtau_rounded_q_squared <- if (!is.null(acf_result_rounded_q_squared)) acf_result_rounded_q_squared$dtau else NA

  # Extract tau and its error from bootstrap.analysis for rounded Q^2 data
  boot_tau_rounded_q_squared <- NA
  boot_dtau_rounded_q_squared <- NA
  if (!is.null(boot_rounded_q_squared)) {
    max_idx <- which.max(boot_rounded_q_squared$Tauint)
    boot_tau_rounded_q_squared <- boot_rounded_q_squared$Tauint[max_idx]
    error_naive_rounded_q_squared <- sd(ac_data_rounded_q_squared) / sqrt(length(ac_data_rounded_q_squared))
    boot_dtau_rounded_q_squared <- 2 * boot_rounded_q_squared$Error[max_idx] * boot_rounded_q_squared$DError[max_idx] / (error_naive_rounded_q_squared^2)
  }

  if (!is.null(acf_result_rounded_q_squared)) {
    write_log(paste0("analyze_topological_charge: Q^2 from rounded computeacf results - tau=", acf_tau_rounded_q_squared, ", dtau=", acf_dtau_rounded_q_squared))
  }

  # Append Q^2 from rounded results to summary file
  tryCatch(
    {
      cat("\n\n--- Rounded Q^2 Data ---\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint: %s\n", as.character(tauint_rounded_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint: %s\n", as.character(dtauint_rounded_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau: %s\n", as.character(acf_tau_rounded_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau: %s\n", as.character(acf_dtau_rounded_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau: %s\n", as.character(boot_tau_rounded_q_squared)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau: %s\n", as.character(boot_dtau_rounded_q_squared)), file = summary_file, append = TRUE)

      # HMC normalized
      cat("\nHMC normalized:\n", file = summary_file, append = TRUE)
      cat(sprintf("uwerr_tauint_hmc: %s\n", as.character(if (!is.na(tauint_rounded_q_squared)) tauint_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("uwerr_dtauint_hmc: %s\n", as.character(if (!is.na(dtauint_rounded_q_squared)) dtauint_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_tau_hmc: %s\n", as.character(if (!is.na(acf_tau_rounded_q_squared)) acf_tau_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("computeacf_dtau_hmc: %s\n", as.character(if (!is.na(acf_dtau_rounded_q_squared)) acf_dtau_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis tau_hmc: %s\n", as.character(if (!is.na(boot_tau_rounded_q_squared)) boot_tau_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)
      cat(sprintf("bootstrap.analysis dtau_hmc: %s\n", as.character(if (!is.na(boot_dtau_rounded_q_squared)) boot_dtau_rounded_q_squared * measurement_interval else NA)), file = summary_file, append = TRUE)

      # PTBC normalized
      if (!is.na(n_replicas)) {
        cat("\nPTBC normalized:\n", file = summary_file, append = TRUE)
        cat(sprintf("uwerr_tauint_ptbc: %s\n", as.character(if (!is.na(tauint_rounded_q_squared)) tauint_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("uwerr_dtauint_ptbc: %s\n", as.character(if (!is.na(dtauint_rounded_q_squared)) dtauint_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_tau_ptbc: %s\n", as.character(if (!is.na(acf_tau_rounded_q_squared)) acf_tau_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("computeacf_dtau_ptbc: %s\n", as.character(if (!is.na(acf_dtau_rounded_q_squared)) acf_dtau_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis tau_ptbc: %s\n", as.character(if (!is.na(boot_tau_rounded_q_squared)) boot_tau_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
        cat(sprintf("bootstrap.analysis dtau_ptbc: %s\n", as.character(if (!is.na(boot_dtau_rounded_q_squared)) boot_dtau_rounded_q_squared * measurement_interval * n_replicas else NA)), file = summary_file, append = TRUE)
      }

      write_log("analyze_topological_charge: summary for Q^2 from rounded data written")
    },
    error = function(e) {
      write_log(paste0("ERROR writing summary file for Q^2 from rounded data: ", conditionMessage(e)))
    }
  )

  return(list(timeseries = timeseries_plot, autocorr = uw, tauint = tauint, dtauint = dtauint))
}

# CLI entrypoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analysis_topological_charge.R <directory> [skip_initial] [s] [s_squared] [r] [xlim_min] [xlim_max]")
}

directory <- args[1]
logs_dir <- file.path(directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
assign("WF_LOG_FILE", file.path(logs_dir, "analysis_topological_charge.log"), envir = .GlobalEnv)
skip_initial <- if (length(args) >= 2) as.integer(args[2]) else 0
s <- if (length(args) >= 3) as.numeric(args[3]) else 2.5
s_squared <- if (length(args) >= 4) as.numeric(args[4]) else 5.0
r <- if (length(args) >= 5) as.integer(args[5]) else 1
xlim_min <- if (length(args) >= 6) as.numeric(args[6]) else NULL
xlim_max <- if (length(args) >= 7) as.numeric(args[7]) else NULL
analyze_topological_charge(directory,
  skip_initial = skip_initial, s = s, s_squared = s_squared, r = r,
  xlim_min = xlim_min, xlim_max = xlim_max
)
write_log("=== Analysis completed successfully ===")
