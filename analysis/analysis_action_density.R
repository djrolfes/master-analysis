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

analyze_action_density <- function(directory, skip_initial = 0) {
    # Configure logfile for this run
    assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
    write_log(paste0("analyze_action_density: start directory=", directory, " skip_initial=", skip_initial))

    # Read YAML config to find filename and Wilson flow time
    cfg <- tryCatch(
        {
            read_yaml_config(directory)
        },
        error = function(e) {
            write_log(paste0("ERROR reading YAML config: ", conditionMessage(e)))
            stop(e)
        }
    )

    # Get action density filename
    action_density_filename <- cfg$GaugeObservableParams$action_density_filename
    if (is.null(action_density_filename) || action_density_filename == "") {
        write_log("ERROR: action_density_filename not found in YAML under GaugeObservableParams")
        stop("action_density_filename not found in YAML")
    }

    # Get Wilson flow time (tau) from WilsonFlowParams
    if (is.null(cfg$GaugeObservableParams$WilsonFlowParams)) {
        write_log("ERROR: WilsonFlowParams not found in YAML under GaugeObservableParams")
        stop("WilsonFlowParams not found in YAML")
    }

    t0 <- cfg$GaugeObservableParams$WilsonFlowParams$tau
    if (is.null(t0)) {
        write_log("ERROR: tau not found in WilsonFlowParams")
        stop("tau not found in WilsonFlowParams")
    }

    write_log(paste0("analyze_action_density: Wilson flow time t0 = ", t0))

    action_density_path <- file.path(directory, action_density_filename)
    if (!file.exists(action_density_path)) {
        write_log(paste0("ERROR: action density file not found: ", action_density_path))
        stop(paste0("Action density file not found: ", action_density_path))
    }

    # Read the data file using the shared data_io reader
    df <- tryCatch(
        {
            read_data_file(action_density_path)
        },
        error = function(e) {
            write_log(paste0("ERROR reading action density file: ", conditionMessage(e)))
            stop(e)
        }
    )

    write_log(paste0("analyze_action_density: read data with ", nrow(df), " rows and columns: ", paste(names(df), collapse = ", ")))

    # Detect column names: step and action_density
    col_names <- tolower(names(df))
    step_col <- names(df)[which(grepl("step", col_names) | grepl("hmc_step", col_names))][1]

    # Look for action_density column
    action_col_candidates <- names(df)[which(grepl("action", col_names) | grepl("density", col_names))]
    action_col <- if (length(action_col_candidates) >= 1) action_col_candidates[1] else names(df)[2]

    if (is.na(step_col) || length(step_col) == 0) {
        write_log("WARNING: could not find 'step' column; using first column as step")
        step_col <- names(df)[1]
    }
    if (is.na(action_col) || length(action_col) == 0) {
        write_log("ERROR: could not find action density column")
        stop("Could not find action density column in file")
    }

    write_log(paste0("analyze_action_density: using step column '", step_col, "' and action density column '", action_col, "'"))

    # Prepare a tidy data.frame with numeric values
    data <- df %>%
        transmute(
            step = as.integer(as.numeric(.data[[step_col]])),
            action_density = as.numeric(.data[[action_col]])
        ) %>%
        filter(!is.na(step) & !is.na(action_density))

    write_log(paste0("analyze_action_density: prepared data rows = ", nrow(data)))

    # Skip initial steps if requested
    if (skip_initial < 0) skip_initial <- 0
    if (skip_initial >= nrow(data)) {
        write_log(paste0("WARNING: skip_initial >= number of rows (", nrow(data), "); nothing to analyze"))
        stop("skip_initial exceeds available data")
    }

    if (skip_initial > 0) {
        data <- data[(skip_initial + 1):nrow(data), ]
        write_log(paste0("analyze_action_density: after skipping ", skip_initial, " steps, ", nrow(data), " rows remain"))
    }

    # 1) Plot action density vs HMC step
    timeseries_plot <- ggplot(data, aes(x = step, y = action_density)) +
        geom_line(color = "steelblue") +
        geom_point(size = 0.8, alpha = 0.7) +
        labs(
            title = "Action Density vs HMC Step",
            x = "HMC Step",
            y = "Action Density E(t)"
        ) +
        theme_minimal()

    out_ts_pdf <- file.path(directory, "action_density_timeseries.pdf")
    write_log(paste0("analyze_action_density: saving timeseries plot to ", out_ts_pdf))
    tryCatch(
        {
            ggsave(out_ts_pdf, plot = timeseries_plot, width = 8, height = 4)
            write_log("analyze_action_density: timeseries plot saved")
        },
        error = function(e) {
            write_log(paste0("ERROR saving timeseries plot: ", conditionMessage(e)))
        }
    )

    # 2) Autocorrelation-corrected error analysis using uwerr
    write_log("analyze_action_density: computing autocorrelation-corrected errors for <E>")

    E_data <- data$action_density

    # Use uwerrprimary to get proper error analysis accounting for autocorrelation
    uw_E <- tryCatch(
        {
            hadron::uwerrprimary(E_data, pl = FALSE)
        },
        error = function(e) {
            write_log(paste0("ERROR computing uwerrprimary for E: ", conditionMessage(e)))
            stop(e)
        }
    )

    mean_E <- uw_E$value
    error_E <- uw_E$dvalue
    tauint_E <- uw_E$tauint

    write_log(paste0("analyze_action_density: <E> = ", mean_E, " ± ", error_E, " (tau_int = ", tauint_E, ")"))

    # 3) Error propagation for t0^2 * <E>
    # For simple multiplication f(x) = c*x, error propagates as: error_f = c * error_x
    # The uwerr error already accounts for autocorrelation, so we can use simple propagation
    write_log(paste0("analyze_action_density: computing t0^2 * <E> where t0 = ", t0))

    t0_squared <- t0^2
    mean_t0_squared_E <- t0_squared * mean_E
    error_t0_squared_E <- t0_squared * error_E

    write_log(paste0("analyze_action_density: t0^2 * <E> = ", mean_t0_squared_E, " ± ", error_t0_squared_E))

    # 4) Generate bootstrap samples for visualization (using parametric bootstrap from uwerr result)
    # Note: These are for histogram visualization only; the error is already computed correctly above
    write_log("analyze_action_density: generating bootstrap samples for histogram visualization")

    set.seed(1234) # For reproducibility
    # Generate samples from normal distribution using uwerr mean and error
    # This gives a visual representation of the uncertainty
    E_boot_samples <- rnorm(1500, mean = mean_E, sd = error_E)
    t0_squared_E_samples <- t0_squared * E_boot_samples

    write_log(paste0("analyze_action_density: histogram samples generated (", length(t0_squared_E_samples), " samples)"))

    # 5) Create histogram of t0^2 * <E>
    write_log("analyze_action_density: creating histogram of t0^2 * <E>")

    hist_data <- data.frame(t0_squared_E = t0_squared_E_samples)

    histogram_plot <- ggplot(hist_data, aes(x = t0_squared_E)) +
        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
        geom_vline(xintercept = mean_t0_squared_E, color = "red", linetype = "dashed", size = 1) +
        geom_vline(xintercept = mean_t0_squared_E - error_t0_squared_E, color = "orange", linetype = "dotted", size = 0.8) +
        geom_vline(xintercept = mean_t0_squared_E + error_t0_squared_E, color = "orange", linetype = "dotted", size = 0.8) +
        labs(
            title = expression(paste("Bootstrap Distribution of ", t[0]^2, " ", symbol("\341"), "E", symbol("\361"))),
            x = expression(paste(t[0]^2, " ", symbol("\341"), "E", symbol("\361"))),
            y = "Count"
        ) +
        annotate("text",
            x = mean_t0_squared_E,
            y = Inf,
            label = sprintf("%.4f ± %.4f", mean_t0_squared_E, error_t0_squared_E),
            vjust = 2,
            color = "red",
            size = 4
        ) +
        theme_minimal()

    out_hist_pdf <- file.path(directory, "action_density_t0_squared_E_histogram.pdf")
    write_log(paste0("analyze_action_density: saving histogram to ", out_hist_pdf))
    tryCatch(
        {
            ggsave(out_hist_pdf, plot = histogram_plot, width = 8, height = 6)
            write_log("analyze_action_density: histogram saved")
        },
        error = function(e) {
            write_log(paste0("ERROR saving histogram: ", conditionMessage(e)))
        }
    )

    # 6) Save results to text file
    out_txt <- file.path(directory, "action_density_t0_squared_E_results.txt")
    write_log(paste0("analyze_action_density: saving results to ", out_txt))

    tryCatch(
        {
            cat("Action Density Bootstrap Analysis Results\n", file = out_txt)
            cat("==========================================\n\n", file = out_txt, append = TRUE)
            cat(sprintf("Wilson flow time t0: %.6f\n", t0), file = out_txt, append = TRUE)
            cat(sprintf("t0^2: %.6f\n", t0_squared), file = out_txt, append = TRUE)
            cat(sprintf("Number of configurations (after skip): %d\n", nrow(data)), file = out_txt, append = TRUE)
            cat(sprintf("Number of bootstrap samples (for histogram): %d\n", length(t0_squared_E_samples)), file = out_txt, append = TRUE)
            cat(sprintf("Integrated autocorrelation time (tau_int): %.6f\n", tauint_E), file = out_txt, append = TRUE)
            cat("\n", file = out_txt, append = TRUE)
            cat(sprintf("<E> = %.6f ± %.6f (error includes autocorrelation)\n", mean_E, error_E), file = out_txt, append = TRUE)
            cat(sprintf("t0^2 * <E> = %.6f ± %.6f\n", mean_t0_squared_E, error_t0_squared_E), file = out_txt, append = TRUE)
            write_log("analyze_action_density: results saved to text file")
        },
        error = function(e) {
            write_log(paste0("ERROR writing results file: ", conditionMessage(e)))
        }
    )

    # 7) Save autocorrelation plot and results
    write_log("analyze_action_density: saving autocorrelation plot and results")

    uw <- uw_E # Use the uwerr result from earlier

    if (!is.null(uw)) {
        tauint <- if (!is.null(uw$tauint)) uw$tauint else NA
        dtauint <- if (!is.null(uw$dtauint)) uw$dtauint else NA

        n_eff <- nrow(data) / (2 * tauint) # Effective number of independent measurements
        write_log(paste0("analyze_action_density: uwerr results - tauint=", tauint, ", dtauint=", dtauint, ", n_eff=", n_eff))

        # Save the uwerr plot to PDF
        out_ac_pdf <- file.path(directory, "action_density_autocorr.pdf")
        write_log(paste0("analyze_action_density: saving autocorr plot to ", out_ac_pdf))
        tryCatch(
            {
                pdf(out_ac_pdf, width = 8, height = 6)
                plot(uw)
                dev.off()
                write_log("analyze_action_density: autocorr plot saved")
            },
            error = function(e) {
                write_log(paste0("ERROR saving autocorr plot: ", conditionMessage(e)))
            }
        )

        # Append autocorrelation results to text file
        tryCatch(
            {
                cat("\nDetailed Autocorrelation Analysis\n", file = out_txt, append = TRUE)
                cat("==================================\n", file = out_txt, append = TRUE)
                cat(sprintf("Integrated autocorrelation time: %.6f ± %.6f\n", tauint, dtauint), file = out_txt, append = TRUE)
                cat(sprintf("Effective number of independent measurements: %.1f\n", n_eff), file = out_txt, append = TRUE)
                cat(sprintf("Error enhancement factor: sqrt(2*tau_int) = %.3f\n", sqrt(2 * tauint)), file = out_txt, append = TRUE)
                cat("\nNote: Errors above already include autocorrelation correction via uwerr.\n", file = out_txt, append = TRUE)
                write_log("analyze_action_density: autocorrelation results written to text file")
            },
            error = function(e) {
                write_log(paste0("ERROR appending autocorrelation results: ", conditionMessage(e)))
            }
        )
    }

    write_log("analyze_action_density: analysis complete")

    return(list(
        timeseries = timeseries_plot,
        histogram = histogram_plot,
        mean_E = mean_E,
        error_E = error_E,
        mean_t0_squared_E = mean_t0_squared_E,
        error_t0_squared_E = error_t0_squared_E,
        t0 = t0,
        autocorr = uw
    ))
}

# CLI entrypoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript analysis_action_density.R <directory> [skip_initial]")
}

directory <- args[1]
logs_dir <- file.path(directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
assign("WF_LOG_FILE", file.path(logs_dir, "analysis_action_density.log"), envir = .GlobalEnv)
skip_initial <- if (length(args) >= 2) as.integer(args[2]) else 0
analyze_action_density(directory, skip_initial = skip_initial)
