library(ggplot2)
library(dplyr)
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory



analyze_acceptance <- function(directory, skip_steps = 0) {
    # Read the simulation log data
    config <- read_yaml_config(directory)
    filename <- config$SimulationLoggingParams$log_filename

    # Build pattern to match simulation log files (not analysis outputs)
    base_prefix <- sub("\\.txt$", "", filename)
    pattern <- paste0("^", base_prefix, "(\\.(rank[0-9]+))?\\.txt$")

    files <- list.files(directory, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
        stop(paste0("No simulation log files found matching pattern '", pattern, "' in directory: ", directory))
    }

    # Process the first (or main) file
    log_file <- files[1]
    message(paste0("Analyzing acceptance from: ", basename(log_file)))

    log_data <- read_data_file(log_file)
    print(head(log_data))

    # Check if acceptance column exists
    if (!"acceptance" %in% names(log_data)) {
        stop("Error: 'acceptance' column not found in simulation log")
    }

    # Check if accept column exists (for bootstrap later)
    has_accept <- "accept" %in% names(log_data)
    if (has_accept) {
        message("Found 'accept' column - will use for bootstrap analysis")
    } else {
        message("No 'accept' column found - using acceptance values directly")
    }

    # Prepare timeseries data
    step_col <- if ("step" %in% names(log_data)) "step" else names(log_data)[1]
    timeseries_dat <- data.frame(
        y = log_data$acceptance,
        t = log_data[[step_col]]
    )

    # Get total number of measurements
    n_total <- length(log_data$acceptance)

    # ===== STEP 1: DETECT THERMALIZATION FIRST =====
    # Determine the measurement interval (spacing between HMC steps in data)
    if (nrow(log_data) > 1) {
        step_intervals <- diff(log_data[[step_col]])
        measurement_interval <- as.integer(median(step_intervals))
    } else {
        measurement_interval <- 1 # Default if only one measurement
    }

    # Detect plateau using running mean and moving window statistics
    window_size <- min(50, floor(n_total / 4)) # Use 50 points or 1/4 of data, whichever is smaller
    thermalization_idx <- 1 # Default: start from beginning

    if (n_total > window_size * 2) {
        # Calculate rolling standard deviation over windows
        rolling_std <- sapply((window_size + 1):n_total, function(i) {
            sd(log_data$acceptance[(i - window_size + 1):i])
        })

        # Calculate rolling mean over windows
        rolling_mean <- sapply((window_size + 1):n_total, function(i) {
            mean(log_data$acceptance[(i - window_size + 1):i])
        })

        # Find where the rolling mean stabilizes (derivative becomes small)
        if (length(rolling_mean) > 10) {
            mean_derivative <- abs(diff(rolling_mean))
            # Normalize by the overall std
            overall_std <- sd(log_data$acceptance)
            normalized_derivative <- mean_derivative / overall_std

            # Find first point where derivative stays below threshold for sustained period
            threshold <- 0.02 # 2% of std per measurement
            sustained_length <- min(20, floor(length(normalized_derivative) / 5))

            for (i in 1:(length(normalized_derivative) - sustained_length)) {
                if (all(normalized_derivative[i:(i + sustained_length)] < threshold)) {
                    thermalization_idx <- i + window_size
                    break
                }
            }
        }
    }

    # Convert to HMC steps
    thermalization_hmc_step <- log_data[[step_col]][thermalization_idx]

    # Use the maximum of provided skip_steps, detected thermalization, or minimum threshold
    minimum_skip <- max(100, measurement_interval * 10) # At least 100 steps or 10 measurements
    recommended_skip <- max(skip_steps, thermalization_hmc_step, minimum_skip)

    # Find the index corresponding to the recommended skip
    stat_start <- max(thermalization_idx, skip_steps + 1)
    stat_end <- n_total

    # ===== STEP 2: PERFORM STATISTICAL ANALYSIS ON THERMALIZED DATA =====

    # If we have the accept column, bootstrap those binary values
    # Otherwise, use uwerr on the acceptance timeseries
    if (has_accept) {
        # Extract thermalized accept values (binary 0/1)
        analysis_data <- log_data$accept[stat_start:stat_end]

        # Bootstrap analysis on binary accept values
        boot_result <- tryCatch(
            {
                hadron::bootstrap.analysis(analysis_data, boot.R = 400, boot.l = 2, pl = FALSE)
            },
            error = function(e) {
                warning(paste("bootstrap.analysis failed:", conditionMessage(e)))
                return(NULL)
            }
        )

        if (!is.null(boot_result)) {
            # Extract result at optimal block size
            max_idx <- which.max(boot_result$Tauint)
            mean_acceptance <- boot_result$Mean[max_idx]
            error_acceptance <- boot_result$Error[max_idx]
            tau_int <- boot_result$Tauint[max_idx]
            tau_error <- boot_result$DError[max_idx]
        } else {
            mean_acceptance <- mean(analysis_data)
            error_acceptance <- sd(analysis_data) / sqrt(length(analysis_data))
            tau_int <- NA
            tau_error <- NA
        }
    } else {
        # Use acceptance timeseries with uwerr
        analysis_data <- log_data$acceptance[stat_start:stat_end]

        uw_result <- tryCatch(
            {
                hadron::uwerrprimary(analysis_data, pl = FALSE)
            },
            error = function(e) {
                warning(paste("uwerr analysis failed:", conditionMessage(e)))
                return(NULL)
            }
        )

        if (!is.null(uw_result)) {
            mean_acceptance <- uw_result$value
            error_acceptance <- uw_result$dvalue
            tau_int <- uw_result$tauint
            tau_error <- uw_result$dtauint
        } else {
            mean_acceptance <- mean(analysis_data)
            error_acceptance <- sd(analysis_data) / sqrt(length(analysis_data))
            tau_int <- NA
            tau_error <- NA
        }
        boot_result <- NULL
    }

    # ===== STEP 3: CREATE VISUALIZATIONS =====
    # Create PDF with multiple plots
    pdf_filename <- if (has_accept) "acceptance_bootstrap_analysis.pdf" else "acceptance_timeseries_analysis.pdf"
    pdf(file.path(directory, pdf_filename), width = 12, height = 10)

    # Set up layout: Two timeseries on top (full width), histogram and analysis below
    layout(matrix(c(1, 1, 2, 2, 3, 4), nrow = 3, byrow = TRUE), heights = c(1.5, 1.5, 1))

    # Plot 1: Full timeseries (takes full width)
    par(mar = c(4, 4, 3, 2))
    plot(timeseries_dat$t, timeseries_dat$y,
        type = "l", col = "steelblue",
        xlab = "HMC Step", ylab = "Acceptance Rate",
        main = "Acceptance Rate Timeseries (Full)",
        ylim = c(0, 1)
    )

    # Mark detected thermalization point
    abline(v = log_data[[step_col]][thermalization_idx], col = "orange", lty = 2, lwd = 2)

    # Add error band
    rect(
        xleft = log_data[[step_col]][stat_start],
        xright = log_data[[step_col]][stat_end],
        ytop = mean_acceptance + error_acceptance,
        ybottom = mean_acceptance - error_acceptance,
        border = NA, col = rgb(0.6, 0, 0, 0.3)
    )

    # Add legend with result
    legend_text <- c(
        sprintf("Mean = %.4f +/- %.4f", mean_acceptance, error_acceptance),
        sprintf("Detected therm: step %d", log_data[[step_col]][thermalization_idx])
    )

    if (!is.na(tau_int)) {
        legend_text <- c(legend_text, sprintf("tau_int = %.2f +/- %.2f", tau_int, tau_error))
    }

    legend("topright",
        legend = legend_text,
        col = c("red", "orange", if (!is.na(tau_int)) "red" else NULL),
        lty = c(NA, 2, if (!is.na(tau_int)) NA else NULL),
        lwd = c(NA, 2, if (!is.na(tau_int)) NA else NULL),
        bty = "n"
    )

    # Plot 2: Timeseries after thermalization (takes full width)
    par(mar = c(4, 4, 3, 2))

    # Prepare data after detected thermalization
    thermalized_data <- log_data[thermalization_idx:nrow(log_data), ]

    plot(thermalized_data[[step_col]], thermalized_data$acceptance,
        type = "l", col = "steelblue",
        xlab = "HMC Step", ylab = "Acceptance Rate",
        main = sprintf(
            "Acceptance Rate (After Thermalization, from step %d)",
            log_data[[step_col]][thermalization_idx]
        ),
        ylim = c(0, 1)
    )

    # Show the analysis region used (if different from detected thermalization)
    if (skip_steps > 0 && stat_start > thermalization_idx) {
        abline(v = log_data[[step_col]][stat_start], col = "red", lty = 2)
    }

    # Add horizontal line at mean
    abline(h = mean_acceptance, col = "red", lwd = 2)

    # Add error band
    rect(
        xleft = par("usr")[1],
        xright = par("usr")[2],
        ytop = mean_acceptance + error_acceptance,
        ybottom = mean_acceptance - error_acceptance,
        border = NA, col = rgb(0.6, 0, 0, 0.2)
    )

    # Add legend
    legend("topright",
        legend = c(
            sprintf("Mean = %.4f +/- %.4f", mean_acceptance, error_acceptance),
            sprintf("N_measurements = %d", length(analysis_data))
        ),
        col = c("red", NA),
        lty = c(1, NA),
        lwd = c(2, NA),
        bty = "n"
    )

    # Plot 3: Histogram of analysis region
    par(mar = c(4, 4, 3, 2))

    # For accept column: show binary histogram
    # For acceptance column: show continuous histogram
    if (has_accept) {
        hist(analysis_data,
            breaks = c(-0.5, 0.5, 1.5), col = "lightblue",
            xlab = "Accept (0=reject, 1=accept)",
            main = sprintf("Accept Distribution (Mean=%.4f)", mean_acceptance),
            freq = TRUE,
            xaxt = "n"
        )
        axis(1, at = c(0, 1))
    } else {
        # Calculate bin width to get exactly 30 bins for the thermalized data
        data_range <- range(thermalized_data$acceptance)
        bin_width <- diff(data_range) / 30
        breaks <- seq(data_range[1], data_range[2], by = bin_width)

        hist(thermalized_data$acceptance,
            breaks = breaks, col = "lightblue",
            xlab = "Acceptance Rate",
            main = "Histogram (after thermalization, 30 bins)",
            freq = TRUE
        )
        abline(v = mean_acceptance, col = "red", lwd = 2)
        abline(
            v = c(mean_acceptance - error_acceptance, mean_acceptance + error_acceptance),
            col = "red", lty = 2
        )
    }

    # Plot 4: Analysis summary plot
    par(mar = c(4, 4, 3, 2))
    if (has_accept && !is.null(boot_result)) {
        # Plot bootstrap analysis results
        plot(boot_result$Blocksize, boot_result$Tauint,
            type = "b", col = "blue", pch = 19,
            xlab = "Block Size", ylab = "Tau_int",
            main = "Bootstrap Analysis (Integrated Autocorrelation Time)"
        )
        grid()

        # Mark optimal block size
        max_idx <- which.max(boot_result$Tauint)
        points(boot_result$Blocksize[max_idx], boot_result$Tauint[max_idx],
            col = "red", pch = 19, cex = 1.5
        )

        legend("topright",
            legend = c(
                sprintf("Optimal block: %d", boot_result$Blocksize[max_idx]),
                sprintf("tau_int = %.2f", boot_result$Tauint[max_idx])
            ),
            bty = "n"
        )
    } else if (!is.null(uw_result)) {
        plot(uw_result, main = "Autocorrelation Analysis (uwerr)")
    } else {
        plot.new()
        text(0.5, 0.5, "Statistical analysis failed", cex = 1.5)
    }

    dev.off()

    # ===== STEP 4: SAVE SUMMARY =====
    summary_file <- file.path(directory, "acceptance_summary.txt")
    cat("Acceptance Rate Analysis Summary\n", file = summary_file)
    cat("================================\n\n", file = summary_file, append = TRUE)
    cat(sprintf("Total measurements: %d\n", n_total), file = summary_file, append = TRUE)
    cat(sprintf("Used for analysis: %d\n\n", length(analysis_data)), file = summary_file, append = TRUE)

    # Write thermalization detection info to summary
    cat(sprintf("--- Thermalization Detection ---\n"), file = summary_file, append = TRUE)
    cat(sprintf("Measurement interval: %d HMC steps\n", measurement_interval), file = summary_file, append = TRUE)
    cat(sprintf(
        "Detected plateau at measurement: %d (HMC step %d)\n",
        thermalization_idx, thermalization_hmc_step
    ), file = summary_file, append = TRUE)
    cat(sprintf("Recommended thermalization skip: %d HMC steps\n", recommended_skip), file = summary_file, append = TRUE)
    cat(sprintf("(Detection method: Running mean stabilization)\n\n"), file = summary_file, append = TRUE)

    # Write statistical results
    cat(sprintf("--- Statistical Analysis (on thermalized data) ---\n"), file = summary_file, append = TRUE)
    if (has_accept) {
        cat(sprintf("Analysis method: Bootstrap on binary accept values\n"), file = summary_file, append = TRUE)
    } else {
        cat(sprintf("Analysis method: Gamma method (uwerr) on acceptance timeseries\n"), file = summary_file, append = TRUE)
    }

    cat(sprintf("Mean acceptance rate: %.6f +/- %.6f\n", mean_acceptance, error_acceptance), file = summary_file, append = TRUE)

    if (!is.na(tau_int)) {
        cat(sprintf(
            "Integrated autocorrelation time: %.4f +/- %.4f measurements\n",
            tau_int, tau_error
        ), file = summary_file, append = TRUE)
        cat(sprintf(
            "                                 = %.0f +/- %.0f HMC steps\n",
            tau_int * measurement_interval,
            tau_error * measurement_interval
        ), file = summary_file, append = TRUE)
    }

    message(sprintf("Acceptance analysis complete. Results saved to %s", directory))

    # Return the recommended skip value (in HMC steps)
    return(invisible(list(
        recommended_skip = recommended_skip,
        thermalization_idx = thermalization_idx,
        mean_acceptance = mean_acceptance,
        error_acceptance = error_acceptance,
        tau_int = tau_int,
        tau_error = tau_error
    )))
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
    stop("Usage: Rscript analysis_acceptance.R <directory> [skip_steps]")
}
directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_steps <- if (length(args) >= 2) as.integer(args[2]) else 0
result <- analyze_acceptance(directory, skip_steps = skip_steps)
cat(sprintf("Recommended skip steps: %d\n", result$recommended_skip))
cat(sprintf("Mean acceptance: %.4f +/- %.4f\n", result$mean_acceptance, result$error_acceptance))
