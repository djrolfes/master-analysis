library(ggplot2)
library(dplyr)
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory



analyze_plaquette <- function(directory, skip_steps = 0) {
  # Read the plaquette data
  plaquette_data <- read_data_plaquette_filename(directory)
  print(head(plaquette_data))

  # Timeseries plot using plot_timeseries
  # Prepare data frame with y (plaquette values) and t (HMC steps)
  timeseries_dat <- data.frame(
    y = plaquette_data$plaquette,
    t = plaquette_data$step
  )

  # Get total number of measurements
  n_total <- length(plaquette_data$plaquette)

  # ===== STEP 1: DETECT THERMALIZATION FIRST =====
  # Determine the measurement interval (spacing between HMC steps in data)
  if (nrow(plaquette_data) > 1) {
    step_intervals <- diff(plaquette_data$step)
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
      sd(plaquette_data$plaquette[(i - window_size + 1):i])
    })

    # Calculate rolling mean over windows
    rolling_mean <- sapply((window_size + 1):n_total, function(i) {
      mean(plaquette_data$plaquette[(i - window_size + 1):i])
    })

    # Find where the rolling mean stabilizes (derivative becomes small)
    if (length(rolling_mean) > 10) {
      mean_derivative <- abs(diff(rolling_mean))
      # Normalize by the overall std
      overall_std <- sd(plaquette_data$plaquette)
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
  thermalization_hmc_step <- plaquette_data$step[thermalization_idx]

  # Use the maximum of provided skip_steps, detected thermalization, or minimum threshold
  minimum_skip <- max(100, measurement_interval * 10) # At least 100 steps or 10 measurements
  recommended_skip <- max(skip_steps, thermalization_hmc_step, minimum_skip)

  # Find the index corresponding to the recommended skip
  stat_start <- max(thermalization_idx, skip_steps + 1)
  stat_end <- n_total

  # ===== STEP 2: PERFORM UWERR ANALYSIS ON THERMALIZED DATA =====
  # Extract thermalized data for analysis
  analysis_data <- plaquette_data$plaquette[stat_start:stat_end]

  # Perform uwerr analysis on thermalized data
  uw_result <- tryCatch(
    {
      hadron::uwerrprimary(analysis_data, pl = FALSE)
    },
    error = function(e) {
      warning(paste("uwerr analysis failed:", conditionMessage(e)))
      return(NULL)
    }
  )

  # Create PDF with multiple plots
  pdf(file.path(directory, "timeseries_plaquette.pdf"), width = 12, height = 10)

  # Set up layout: Two timeseries on top (full width), histogram and uwerr below
  # Layout matrix: 1 and 2 are timeseries plots, 3 is histogram, 4 is uwerr
  layout(matrix(c(1, 1, 2, 2, 3, 4), nrow = 3, byrow = TRUE), heights = c(1.5, 1.5, 1))

  # Plot 1: Full timeseries (takes full width)
  par(mar = c(4, 4, 3, 2))
  plot(timeseries_dat$t, timeseries_dat$y,
    type = "l", col = "steelblue",
    xlab = "HMC Step", ylab = "Plaquette",
    main = "Plaquette Timeseries (Full)"
  )

  # Mark detected thermalization point
  abline(v = plaquette_data$step[thermalization_idx], col = "orange", lty = 2, lwd = 2)

  # Add error band if uwerr succeeded
  if (!is.null(uw_result)) {
    # Highlight the analysis region
    if (skip_steps > 0) {
      abline(v = plaquette_data$step[stat_start], col = "red", lty = 2)
    }
    rect(
      xleft = plaquette_data$step[stat_start],
      xright = plaquette_data$step[stat_end],
      ytop = uw_result$value + uw_result$dvalue,
      ybottom = uw_result$value - uw_result$dvalue,
      border = NA, col = rgb(0.6, 0, 0, 0.3)
    )
    # abline(h = uw_result$value, col = "red", lwd = 2)

    # Add legend with result
    legend("topright",
      legend = c(
        sprintf("Mean = %.6f +/- %.6f", uw_result$value, uw_result$dvalue),
        sprintf("tau_int = %.2f +/- %.2f", uw_result$tauint, uw_result$dtauint),
        sprintf("Detected therm: step %d", plaquette_data$step[thermalization_idx])
      ),
      col = c("red", "red", "orange"),
      lty = c(NA, NA, 2),
      lwd = c(NA, NA, 2),
      bty = "n"
    )
  } else {
    # Even without uwerr, show thermalization
    legend("topright",
      legend = sprintf("Detected therm: step %d", plaquette_data$step[thermalization_idx]),
      col = "orange",
      lty = 2,
      lwd = 2,
      bty = "n"
    )
  }

  # Plot 2: Timeseries after thermalization (takes full width)
  par(mar = c(4, 4, 3, 2))

  # Prepare data after detected thermalization
  thermalized_data <- plaquette_data[thermalization_idx:nrow(plaquette_data), ]

  plot(thermalized_data$step, thermalized_data$plaquette,
    type = "l", col = "steelblue",
    xlab = "HMC Step", ylab = "Plaquette",
    main = sprintf(
      "Plaquette Timeseries (After Thermalization, from step %d)",
      plaquette_data$step[thermalization_idx]
    )
  )

  # Add error band if uwerr succeeded
  if (!is.null(uw_result)) {
    # Show the analysis region used (if different from detected thermalization)
    if (skip_steps > 0 && stat_start > thermalization_idx) {
      abline(v = plaquette_data$step[stat_start], col = "red", lty = 2)
    }

    # Add horizontal line at mean
    abline(h = uw_result$value, col = "red", lwd = 2)

    # Add error band
    rect(
      xleft = par("usr")[1],
      xright = par("usr")[2],
      ytop = uw_result$value + uw_result$dvalue,
      ybottom = uw_result$value - uw_result$dvalue,
      border = NA, col = rgb(0.6, 0, 0, 0.2)
    )

    # Add legend
    legend("topright",
      legend = c(
        sprintf("Mean = %.6f +/- %.6f", uw_result$value, uw_result$dvalue),
        sprintf("N_measurements = %d", length(analysis_data))
      ),
      col = c("red", NA),
      lty = c(1, NA),
      lwd = c(2, NA),
      bty = "n"
    )
  }

  # Plot 3: Histogram of analysis region (30 bins based on thermalized data)
  par(mar = c(4, 4, 3, 2))
  if (!is.null(uw_result)) {
    # Calculate bin width to get exactly 30 bins for the thermalized data
    data_range <- range(thermalized_data$plaquette)
    bin_width <- diff(data_range) / 30
    breaks <- seq(data_range[1], data_range[2], by = bin_width)

    hist(thermalized_data$plaquette,
      breaks = breaks, col = "lightblue",
      xlab = "Plaquette", main = "Histogram (after thermalization, 30 bins)",
      freq = TRUE
    )
    abline(v = uw_result$value, col = "red", lwd = 2)
    abline(
      v = c(uw_result$value - uw_result$dvalue, uw_result$value + uw_result$dvalue),
      col = "red", lty = 2
    )
  } else {
    # Calculate bin width to get exactly 30 bins
    data_range <- range(thermalized_data$plaquette)
    bin_width <- diff(data_range) / 30
    breaks <- seq(data_range[1], data_range[2], by = bin_width)

    hist(thermalized_data$plaquette,
      breaks = breaks, col = "lightblue",
      xlab = "Plaquette", main = "Histogram (after thermalization, 30 bins)",
      freq = TRUE
    )
  }

  # Plot 4: uwerr summary plot
  par(mar = c(4, 4, 3, 2))
  if (!is.null(uw_result)) {
    plot(uw_result, main = "Autocorrelation Analysis")
  } else {
    plot.new()
    text(0.5, 0.5, "uwerr analysis failed", cex = 1.5)
  }

  dev.off()

  # Save summary to text file
  summary_file <- file.path(directory, "plaquette_summary.txt")
  cat("Plaquette Analysis Summary\n", file = summary_file)
  cat("==========================\n\n", file = summary_file, append = TRUE)
  cat(sprintf("Total measurements: %d\n", n_total), file = summary_file, append = TRUE)
  cat(sprintf("Used for analysis: %d\n\n", length(analysis_data)), file = summary_file, append = TRUE)

  # Write thermalization detection info to summary
  cat(sprintf("\n--- Thermalization Detection ---\n"), file = summary_file, append = TRUE)
  cat(sprintf("Measurement interval: %d HMC steps\n", measurement_interval), file = summary_file, append = TRUE)
  cat(sprintf(
    "Detected plateau at measurement: %d (HMC step %d)\n",
    thermalization_idx, thermalization_hmc_step
  ), file = summary_file, append = TRUE)
  cat(sprintf("Recommended thermalization skip: %d HMC steps\n", recommended_skip), file = summary_file, append = TRUE)
  cat(sprintf("(Detection method: Running mean stabilization)\n"), file = summary_file, append = TRUE)

  if (!is.null(uw_result)) {
    cat(sprintf("\nMean value: %.8f +/- %.8f\n", uw_result$value, uw_result$dvalue), file = summary_file, append = TRUE)
    cat(sprintf(
      "Integrated autocorrelation time: %.4f +/- %.4f measurements\n",
      uw_result$tauint, uw_result$dtauint
    ), file = summary_file, append = TRUE)
    cat(sprintf(
      "                                 = %.0f +/- %.0f HMC steps\n",
      uw_result$tauint * measurement_interval,
      uw_result$dtauint * measurement_interval
    ), file = summary_file, append = TRUE)
    cat(sprintf("Optimal window length: %d\n", uw_result$Wopt), file = summary_file, append = TRUE)
  } else {
    cat("\nuwerr analysis failed\n", file = summary_file, append = TRUE)
  }

  # Write to a separate file for easy parsing by analysis.py
  skip_file <- file.path(directory, "recommended_skip.txt")
  cat(sprintf("%d\n", recommended_skip), file = skip_file)

  message(sprintf("Plaquette analysis complete. Results saved to %s", directory))

  # Return the recommended skip value (in HMC steps)
  # This is already calculated above in the thermalization detection section
  return(recommended_skip)
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript analysis_plaquette.R <directory> [skip_steps]")
}
directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_steps <- if (length(args) >= 2) as.integer(args[2]) else 0
recommended_skip <- analyze_plaquette(directory, skip_steps = skip_steps)
cat(sprintf("Recommended skip steps: %d\n", recommended_skip))
