# Set up library paths to include project-local R library
# This ensures packages installed by check_and_install_hadron.sh are found
script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
)
project_root <- dirname(script_dir)
local_r_lib <- file.path(project_root, ".R", "library")

if (dir.exists(local_r_lib)) {
    .libPaths(c(local_r_lib, .libPaths()))
}

library(hadron)
library(dplyr)
library(ggplot2)
library(tidyr)
source("data_io.R")

# Simple logging helper
write_log <- function(msg) {
    logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
    if (is.na(logfile)) logfile <- "wilson_temp_analysis.log"
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- sprintf("[%s] %s\n", timestamp, msg)
    cat(entry, file = logfile, append = TRUE)
}

# Automatic plateau detection
# Strategy: prefer longer plateaus with chi^2/dof close to 1 and good Q-value

find_best_plateau <- function(effmass_obj, min_window = 5, max_t2_frac = 0.8,
                              chi2_target = 1.0, chi2_tolerance = 0.5,
                              min_qval = 0.1, max_qval = 0.9) {
    tmax <- length(effmass_obj$effMass)
    max_t2 <- floor(tmax * max_t2_frac)

    best_score <- -Inf
    best_t1 <- NA
    best_t2 <- NA
    best_chisqr_dof <- NA
    best_qval <- NA
    best_fit <- NULL

    results <- data.frame()

    # Try different fit ranges
    # Iterate t2 from high to low to find longer plateaus first
    for (t1 in 1:(max_t2 - min_window)) {
        for (t2 in max_t2:(t1 + min_window)) {
            tryCatch(
                {
                    fit <- fit.effectivemass(effmass_obj, t1 = t1, t2 = t2, useCov = TRUE)
                    chisqr_dof <- fit$chisqr / fit$dof
                    qval <- fit$Qval
                    m <- fit$effmassfit$t0[[1]]
                    dm <- fit$effmassfit$se

                    # Only consider fits with reasonable chi^2/dof and Q-value
                    if (chisqr_dof < (chi2_target + 1.5) && qval >= min_qval && qval <= max_qval && (m + dm) > 0) {
                        window_size <- (t2 - t1 + 1) / (max_t2 - 1)
                        window_score <- window_size * 0.5

                        # Score: reward long plateaus, penalize deviation from chi2_target and bad Q-values
                        chi2_penalty <- abs(chisqr_dof - chi2_target) / chi2_tolerance

                        # Q-value penalty: penalize values too close to 0 or 1
                        # Optimal Q-value is around 0.5
                        qval_penalty <- abs(1 - qval) # testing

                        # error penalty: penalize fits with large standard error
                        error_penalty <- abs(dm / m)

                        na_count <- sum(is.na(effmass_obj$effMass.tsboot[, t1:t2]))
                        total_elements <- length(effmass_obj$effMass.tsboot[, t1:t2])
                        na_fraction <- na_count / total_elements

                        score <- window_score - 0.5 * chi2_penalty - 0.5 * qval_penalty - 2 * error_penalty - 1000 * na_fraction

                        results <- rbind(results, data.frame(
                            t1 = t1, t2 = t2, window = window_score,
                            chisqr_dof = chisqr_dof, qval = qval, error_penalty = error_penalty, score = score
                        ))

                        if (score > best_score) {
                            best_score <- score
                            best_t1 <- t1
                            best_t2 <- t2
                            best_chisqr_dof <- chisqr_dof
                            best_qval <- qval
                            best_error_penalty <- error_penalty
                            best_fit <- fit
                        }
                    }
                },
                error = function(e) {
                    # Skip fits that fail
                }
            )
        }
    }

    # Print top 5 candidates
    cat("\nTop 5 plateau candidates:\n")
    top_results <- head(results[order(-results$score), ], 5)
    print(top_results)

    cat(sprintf(
        "\nSelected: t1 = %d, t2 = %d (window = %d, chi^2/dof = %.2f, Q = %.3f, error penalty = %.3f)\n",
        best_t1, best_t2, best_t2 - best_t1 + 1, best_chisqr_dof, best_qval, best_error_penalty
    ))

    return(best_fit)
}

# Function to fit the static potential to temporal Wilson loops
# Model: <W(L,T)> = a * exp(-V(L)*T)
# where V(L) = A + B/L + sigma*L (static potential)
fit_static_potential <- function(directory, skip_steps = 0) {
    write_log(paste0("fit_static_potential: start for directory=", directory))

    # Read the W_temp data
    w_temp_data <- read_data_W_temp_filename(directory)
    write_log(paste0("fit_static_potential: read W_temp data with ", nrow(w_temp_data), " rows"))

    # Data is always in long format with columns: "step" (or "hmc_step"), "L", "T", "W_temp"
    col_names <- colnames(w_temp_data)
    write_log(paste0("fit_static_potential: columns: ", paste(col_names, collapse = ", ")))

    # Standardize column names
    if ("W.temp" %in% col_names) {
        colnames(w_temp_data)[colnames(w_temp_data) == "W.temp"] <- "W_temp"
    }
    if ("step" %in% col_names) {
        colnames(w_temp_data)[colnames(w_temp_data) == "step"] <- "hmc_step"
    }

    # Skip thermalization steps
    if (skip_steps > 0 && "hmc_step" %in% colnames(w_temp_data)) {
        w_temp_data <- w_temp_data %>% filter(hmc_step > skip_steps)
        write_log(paste0("fit_static_potential: after skipping ", skip_steps, " steps, ", nrow(w_temp_data), " rows remain"))
    }

    # Create plots directory
    plots_dir <- file.path(directory, "static_potential_plots")
    if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
        write_log(paste0("fit_static_potential: created plots directory: ", plots_dir))
    }

    unique_L <- unique(w_temp_data$L)
    unique_T <- unique(w_temp_data$T)
    unique_steps <- unique(w_temp_data$hmc_step)

    w_matrices <- lapply(unique_L, function(L_val) {
        # Filter data for this L value
        data_L <- w_temp_data[w_temp_data$L == L_val, ]

        # Reshape to matrix: rows = steps, columns = T values
        data_wide <- pivot_wider(data_L,
            id_cols = hmc_step,
            names_from = T,
            values_from = W_temp,
            names_prefix = "T_"
        )

        # Convert to matrix (exclude step column)
        mat <- as.matrix(data_wide[, -1])
        rownames(mat) <- data_wide$hmc_step

        return(mat)
    })

    # Name the list elements by L value
    names(w_matrices) <- paste0("L_", unique_L)

    v_of_L_data <- data.frame()
    wilson_fits <- list() # Store fit objects for plotting

    # Load or create plateau override file
    plateau_override_file <- file.path(plots_dir, "plateau_override.txt")
    plateau_overrides <- data.frame(L = integer(), t1 = integer(), t2 = integer(), skip = integer())

    if (file.exists(plateau_override_file)) {
        write_log(paste0("fit_static_potential: loading plateau overrides from ", plateau_override_file))
        plateau_overrides <- read.table(plateau_override_file, header = TRUE, comment.char = "#")

        # If skip column is missing, this is an old format file - recreate it with skip column
        if (!"skip" %in% colnames(plateau_overrides)) {
            write_log(paste0("fit_static_potential: old format detected (no skip column), recreating file with new format"))

            # Backup old file
            backup_file <- paste0(plateau_override_file, ".old")
            file.copy(plateau_override_file, backup_file, overwrite = TRUE)
            write_log(paste0("fit_static_potential: backed up old file to ", backup_file))

            # Add skip=0 column
            plateau_overrides$skip <- 0

            # Rewrite file with new format
            sink(plateau_override_file)
            cat("# Plateau Override File\n")
            cat("# Edit the t1 and t2 values below to manually set plateau ranges for each L value\n")
            cat("# Set skip=1 to exclude a point from the V(L) potential fit (will show with reduced opacity)\n")
            cat("# Lines starting with # are comments and will be ignored\n")
            cat("# Format: L  t1  t2  skip\n")
            cat("#\n")
            write.table(plateau_overrides, row.names = FALSE, col.names = TRUE, quote = FALSE)
            sink()

            write_log(paste0("fit_static_potential: recreated plateau override file with skip column"))
        }
    } else {
        write_log(paste0("fit_static_potential: no plateau override file found, will create one after first run"))
    }

    # Store automatically detected plateaus for generating override file
    detected_plateaus <- data.frame(L = integer(), t1 = integer(), t2 = integer(), skip = integer())

    for (i in seq_along(unique_L)) {
        L_val <- unique_L[i]
        L_matrix <- w_matrices[[i]]

        if (nrow(L_matrix) < 2 || ncol(L_matrix) < 2) {
            write_log(paste0("fit_static_potential: skipping L=", L_val, " (insufficient data points)"))
            next
        }

        write_log(paste0("fit_static_potential: fitting L=", L_val, " with ", nrow(L_matrix), " measurements and ", ncol(L_matrix), " T values"))
        write_log(paste0("fit_static_potential: L=", L_val, " T range: [1, ", ncol(L_matrix), "]"))

        # Check for numerical issues
        if (any(!is.finite(L_matrix))) {
            write_log(paste0("fit_static_potential: L=", L_val, " has non-finite values, skipping"))
            next
        }

        tryCatch(
            {
                Time_extent <- ncol(L_matrix)
                wloop <- cf_orig(cf = L_matrix)
                wloop <- cf_meta(wloop, nrObs = 1, Time = Time_extent, symmetrise = FALSE)

                boot.R <- 300
                dbboot.R <- 200
                wloop.boot <- bootstrap.cf(wloop, boot.R = boot.R, boot.l = 2)
                wloop.boot <- double_bootstrap.cf(wloop.boot, dbboot.R = dbboot.R)

                # Calculate bootstrap effective mass
                wloop.efm <- bootstrap.effectivemass(wloop.boot, type = "log")

                # Check if there's an override for this L value
                override_idx <- which(plateau_overrides$L == L_val)

                if (length(override_idx) > 0) {
                    # Use overridden plateau range
                    t1_override <- plateau_overrides$t1[override_idx]
                    t2_override <- plateau_overrides$t2[override_idx]
                    skip_val <- plateau_overrides$skip[override_idx]
                    write_log(paste0("fit_static_potential: L=", L_val, " using OVERRIDE plateau: t1=", t1_override, ", t2=", t2_override, ", skip=", skip_val))

                    wloop.efm.fit <- fit.effectivemass(wloop.efm, t1 = t1_override, t2 = t2_override, useCov = TRUE)

                    # Store the override in detected_plateaus
                    detected_plateaus <- rbind(detected_plateaus, data.frame(L = L_val, t1 = t1_override, t2 = t2_override, skip = skip_val))
                } else {
                    # Find best plateau automatically
                    wloop.efm.fit <- find_best_plateau(wloop.efm,
                        min_window = 3, max_t2_frac = 0.7,
                        chi2_target = 1.0, chi2_tolerance = 0.9,
                        min_qval = 0.1, max_qval = 0.9
                    )
                    # Store detected plateau for override file generation
                    detected_plateaus <- rbind(detected_plateaus, data.frame(
                        L = L_val,
                        t1 = wloop.efm.fit$t1,
                        t2 = wloop.efm.fit$t2,
                        skip = 0
                    ))
                }

                # Calculate y-limits based on fit result with margin
                V_fit <- wloop.efm.fit$effmassfit$t0[[1]]
                V_fit_error <- wloop.efm.fit$effmassfit$se
                y_margin <- 0.1 * V_fit # 3-sigma margin
                ylim_linear <- c(V_fit - y_margin, V_fit + y_margin)

                # Save fitted effective mass plot (linear scale)
                efm_fit_plot_file <- file.path(plots_dir, sprintf("L_%d_effective_mass_fit.pdf", L_val))
                pdf(efm_fit_plot_file, width = 8, height = 6)
                plot(wloop.efm.fit, ylab = "V_eff", xlab = "t/a", xlim = c(1, Time_extent))
                dev.off()

                # Save fitted effective mass plot (log scale)
                efm_fit_plot_log_file <- file.path(plots_dir, sprintf("L_%d_effective_mass_fit_log.pdf", L_val))
                pdf(efm_fit_plot_log_file, width = 8, height = 6)
                plot(wloop.efm.fit, ylab = "V_eff", xlab = "t/a", xlim = c(1, Time_extent), ylim = ylim_linear, log = "y")
                dev.off()

                # Save summary to text file
                summary_file <- file.path(plots_dir, sprintf("L_%d_fit_summary.txt", L_val))
                sink(summary_file)
                cat(sprintf("=== Wilson Loop Effective Mass Fit Summary for L = %d ===\n\n", L_val))
                print(summary(wloop.efm.fit))
                sink()

                V <- wloop.efm.fit$effmassfit$t0[[1]]
                dV <- wloop.efm.fit$effmassfit$se
                chisqr_per_dof <- wloop.efm.fit$chisqr / wloop.efm.fit$dof

                # Extract bootstrap samples for V(L) from the fit
                # wloop.efm.fit$effmassfit$t is a matrix: rows=bootstrap samples, cols=[parameter, chisqr]
                # We only want the first column (the fitted parameter V)
                V_bootstrap_samples <- wloop.efm.fit$effmassfit$t[, 1]

                write_log(paste0(
                    "fit_static_potential: L=", L_val, " fit results: V=",
                    round(V, 6), " ± ", round(dV, 6),
                    ", chi2/dof=", round(chisqr_per_dof, 3),
                    " (chi2=", round(wloop.efm.fit$chisqr, 2), ", dof=", wloop.efm.fit$dof, ")",
                    ", bootstrap samples: ", length(V_bootstrap_samples)
                ))
                # Check if this L value should be skipped in the fit
                skip_in_fit <- 0
                if (length(override_idx) > 0 && "skip" %in% colnames(plateau_overrides)) {
                    skip_in_fit <- plateau_overrides$skip[override_idx]
                }

                v_of_L_data <- rbind(v_of_L_data, data.frame(
                    L = L_val,
                    V = V,
                    V_error = dV,
                    skip = skip_in_fit,
                    boot_samples = I(list(V_bootstrap_samples))
                ))

                # Store fit object for later use
                wilson_fits[[as.character(L_val)]] <- wloop.efm.fit
            },
            error = function(e) {
                write_log(paste0("fit_static_potential: ERROR fitting L=", L_val, ": ", conditionMessage(e)))
            }
        )
    }

    if (nrow(v_of_L_data) == 0) {
        write_log("fit_static_potential: ERROR - no successful V(L) fits")
        return(NULL)
    }

    write_log(paste0("fit_static_potential: successfully fit V(L) for ", nrow(v_of_L_data), " L values"))

    # Generate/update plateau override file if it didn't exist
    if (!file.exists(plateau_override_file) && nrow(detected_plateaus) > 0) {
        write_log(paste0("fit_static_potential: creating plateau override file: ", plateau_override_file))

        sink(plateau_override_file)
        cat("# Plateau Override File\n")
        cat("# Edit the t1 and t2 values below to manually set plateau ranges for each L value\n")
        cat("# Set skip=1 to exclude a point from the V(L) potential fit (will show with reduced opacity)\n")
        cat("# Lines starting with # are comments and will be ignored\n")
        cat("# Format: L  t1  t2  skip\n")
        cat("#\n")
        write.table(detected_plateaus, row.names = FALSE, col.names = TRUE, quote = FALSE)
        sink()

        write_log(paste0("fit_static_potential: wrote ", nrow(detected_plateaus), " plateau ranges to override file"))
        write_log("fit_static_potential: edit plateau_override.txt and rerun to use custom plateau ranges")
    }
    # Save combined V(L) data to text file
    v_of_L_file <- file.path(plots_dir, "V_of_L_data.txt")
    write_log(paste0("fit_static_potential: writing V(L) data to ", v_of_L_file))
    sink(v_of_L_file)
    cat("=== Extracted Static Potential V(L) for all L values ===\n\n")
    cat(sprintf("%-10s %-20s %-20s %-15s %-10s\n", "L", "V(L)", "V_error", "N_boot", "skip"))
    cat(strrep("-", 75), "\n")
    for (i in 1:nrow(v_of_L_data)) {
        cat(sprintf(
            "%-10d %-20.10f %-20.10f %-15d %-10d\n",
            v_of_L_data$L[i],
            v_of_L_data$V[i],
            v_of_L_data$V_error[i],
            length(v_of_L_data$boot_samples[[i]]),
            v_of_L_data$skip[i]
        ))
    }
    sink()
    sink()

    # Save raw bootstrap samples for V(L) to a file
    boot_samples_file <- file.path(plots_dir, "V_of_L_bootstrap_samples.txt")
    write_log(paste0("fit_static_potential: writing V(L) bootstrap samples to ", boot_samples_file))

    # Create matrix of bootstrap samples (rows = bootstrap replicas, cols = L values)
    n_boot <- length(v_of_L_data$boot_samples[[1]])
    boot_matrix <- matrix(NA, nrow = n_boot, ncol = nrow(v_of_L_data))
    colnames(boot_matrix) <- paste0("L_", v_of_L_data$L)

    for (i in 1:nrow(v_of_L_data)) {
        boot_matrix[, i] <- v_of_L_data$boot_samples[[i]]
    }

    write.table(boot_matrix, boot_samples_file,
        row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
    )

    # Now fit V(L) = A - B/L + sigma*L to the extracted V(L) values
    write_log("fit_static_potential: fitting V(L) = A - B/L + sigma*L")

    tryCatch(
        {
            # Initial parameter guesses
            # V(L) = A + B/L + sigma*L (renamed: A=constant, B=Coulomb, sigma=string tension)
            A_init <- mean(v_of_L_data$V)
            B_init <- -0.1
            sigma_init <- 0.01

            write_log(paste0(
                "fit_static_potential: V(L) initial params: A=", round(A_init, 6),
                ", B=", round(B_init, 6), ", sigma=", round(sigma_init, 6)
            ))

            # Define fit function for V(L) with required boot.r parameter
            fn_potential <- function(par, x, boot.r, ...) {
                par[1] + par[2] / x + par[3] * x
            }
            # Filter out skipped points for the fit
            v_of_L_fit <- v_of_L_data[v_of_L_data$skip == 0, ]
            v_of_L_skipped <- v_of_L_data[v_of_L_data$skip == 1, ]

            if (nrow(v_of_L_skipped) > 0) {
                write_log(paste0("fit_static_potential: skipping ", nrow(v_of_L_skipped), " L values in fit: ", paste(v_of_L_skipped$L, collapse = ", ")))
            }

            # Prepare bootstrap samples matrix (only for non-skipped points)
            # Each column is V(L) for a specific L, each row is a bootstrap replica
            n_boot <- length(v_of_L_fit$boot_samples[[1]])
            bsamples <- matrix(NA, nrow = n_boot, ncol = nrow(v_of_L_fit))
            for (i in 1:nrow(v_of_L_fit)) {
                bsamples[, i] <- v_of_L_fit$boot_samples[[i]]
            }

            write_log(paste0("fit_static_potential: using ", n_boot, " bootstrap samples for ", nrow(v_of_L_fit), " L values (", nrow(v_of_L_data), " total, ", nrow(v_of_L_skipped), " skipped)"))

            # Bootstrap fit using bootstrap.nlsfit with pre-generated samples
            potential_fit <- bootstrap.nlsfit(
                fn = fn_potential,
                par.guess = c(A_init, B_init, sigma_init),
                y = v_of_L_fit$V,
                x = v_of_L_fit$L,
                bsamples = bsamples
            )

            # Extract parameters from bootstrap fit
            # potential_fit$t0 contains [A, B, sigma]
            # potential_fit$se contains bootstrap standard errors
            params <- potential_fit$t0
            names(params) <- c("A", "B", "sigma")
            param_errors <- potential_fit$se
            names(param_errors) <- c("A", "B", "sigma")

            write_log(paste0("fit_static_potential: V(L) fit results:"))
            write_log(paste0("  A = ", round(params["A"], 6), " ± ", round(param_errors["A"], 6)))
            write_log(paste0("  B = ", round(params["B"], 6), " ± ", round(param_errors["B"], 6)))
            write_log(paste0("  sigma = ", round(params["sigma"], 6), " ± ", round(param_errors["sigma"], 6)))
            write_log(paste0("  chi2/dof = ", round(potential_fit$chisqr / potential_fit$dof, 3)))
            write_log(paste0("  Qval (p-value) = ", round(potential_fit$Qval, 4)))

            # Calculate r_0 (Sommer scale) and lattice spacing
            # r_0 = sqrt((1.65 + B) / sigma)
            r_0 <- sqrt((1.65 + params["B"]) / params["sigma"])

            # Error propagation for r_0 using bootstrap samples
            # Extract bootstrap samples for B and sigma
            B_samples <- potential_fit$t[, 2]
            sigma_samples <- potential_fit$t[, 3]
            r_0_samples <- sqrt((1.65 + B_samples) / sigma_samples)
            r_0_error <- sd(r_0_samples, na.rm = TRUE)

            # Lattice spacing: a = 0.5 fm / r_0
            a_fm <- 0.5 / r_0
            a_fm_error <- a_fm * (r_0_error / r_0) # Relative error propagation

            write_log(paste0("  r_0 = ", round(r_0, 6), " ± ", round(r_0_error, 6)))
            write_log(paste0("  lattice spacing a = ", round(a_fm, 6), " ± ", round(a_fm_error, 6), " fm"))

            # Write scale setting to file
            scale_file <- file.path(directory, "scale_setting_sommer.txt")
            write_log(paste0("fit_static_potential: writing scale setting to ", scale_file))

            scale_content <- sprintf(
                "# Static Potential Fit Results (Sommer Scale Setting)\n# V(L) = A + B/L + sigma*L\n\n# Fit Parameters:\nA = %.6f ± %.6f\nB = %.6f ± %.6f\nsigma = %.6f ± %.6f\n\n# Fit Quality:\nchi2/dof = %.3f\nQval (p-value) = %.4f\n\n# Scale Setting:\nr_0 = %.6f ± %.6f (lattice units)\nlattice_spacing_a = %.6f ± %.6f fm\n",
                params["A"], param_errors["A"],
                params["B"], param_errors["B"],
                params["sigma"], param_errors["sigma"],
                potential_fit$chisqr / potential_fit$dof,
                potential_fit$Qval,
                r_0, r_0_error,
                a_fm, a_fm_error
            )

            writeLines(scale_content, scale_file)
            write_log(paste0("fit_static_potential: successfully wrote ", scale_file))

            # Save V(L) fit plot to plots directory
            out_pdf <- file.path(plots_dir, "V_of_L_fit.pdf")
            write_log(paste0("fit_static_potential: saving V(L) fit plot to ", out_pdf))

            pdf(out_pdf, width = 10, height = 6)

            # Plot 1: V(L) fit using hadron's plot (better for small errors)
            # Extend plot range to include r_0 (whether r_0 > L_max or r_0 < L_max)
            L_range <- range(v_of_L_data$L)
            # Ensure we extend to at least r_0, and beyond if we have data beyond r_0
            plot_range <- c(L_range[1], max(L_range[2], r_0 * 1.1))

            # Calculate y-limits that include ALL data points (fitted + skipped) with error bars
            # Use all v_of_L_data (not just v_of_L_fit) to ensure skipped points are visible
            y_min <- min(v_of_L_data$V - v_of_L_data$V_error)
            y_max <- max(v_of_L_data$V + v_of_L_data$V_error)
            # Add margin
            y_margin <- 0.02 * (y_max - y_min)
            y_limits <- c(y_min - y_margin, y_max + y_margin)

            write_log(paste0("fit_static_potential: plot_range = [", plot_range[1], ", ", plot_range[2], "], r_0 = ", round(r_0, 4)))
            write_log(paste0("fit_static_potential: y_limits = [", round(y_limits[1], 4), ", ", round(y_limits[2], 4), "]"))

            # Create the plot
            plot(potential_fit,
                xlab = "Spatial Separation L (lattice units)",
                ylab = "Potential V(L)",
                main = sprintf(
                    "Static Potential V(L) = A + B/L + sigma*L\nA=%.4f±%.4f, B=%.4f±%.4f, sigma=%.4f±%.4f\nr_0=%.4f±%.4f, a=%.4f±%.4f fm",
                    params["A"], param_errors["A"],
                    params["B"], param_errors["B"],
                    params["sigma"], param_errors["sigma"],
                    r_0, r_0_error,
                    a_fm, a_fm_error
                ),
                col.line = "black",
                col.band = "black",
                opacity.band = 0.5,
                lwd = 2,
                plot.range = plot_range,
                xlim = plot_range,
                ylim = y_limits,
                rep = FALSE
            )

            # Overlay skipped points with reduced opacity (if any)
            if (nrow(v_of_L_skipped) > 0) {
                points(v_of_L_skipped$L, v_of_L_skipped$V,
                    col = rgb(1, 0, 0, alpha = 0.3), pch = 19, cex = 1.2
                )
                arrows(v_of_L_skipped$L, v_of_L_skipped$V - v_of_L_skipped$V_error,
                    v_of_L_skipped$L, v_of_L_skipped$V + v_of_L_skipped$V_error,
                    angle = 90, code = 3, length = 0.05,
                    col = rgb(1, 0, 0, alpha = 0.3), lwd = 1.5
                )
            }

            # Add error band for r_0 (shaded region)
            rect(r_0 - r_0_error, par("usr")[3], r_0 + r_0_error, par("usr")[4],
                col = rgb(0, 0, 1, alpha = 0.15), border = NA
            )

            # Add vertical line at r_0
            abline(v = r_0, col = "blue", lwd = 1, lty = 1)

            # Add text label for r_0
            text(r_0, par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
                sprintf("r_0 = %.3f ± %.3f", r_0, r_0_error),
                col = "blue", pos = 4, cex = 0.9
            )

            dev.off()

            write_log(paste0("fit_static_potential: successfully saved ", out_pdf))

            return(list(
                v_of_L_data = v_of_L_data,
                fit_params = params,
                fit_errors = param_errors,
                potential_fit = potential_fit,
                wilson_fits = wilson_fits,
                r_0 = r_0,
                r_0_error = r_0_error,
                a_fm = a_fm,
                a_fm_error = a_fm_error
            ))
        },
        error = function(e) {
            write_log(paste0("fit_static_potential: ERROR fitting V(L) model: ", conditionMessage(e)))
            return(list(v_of_L_data = v_of_L_data))
        }
    )
}

# Main execution
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 1 || length(args) > 2) {
        stop("Usage: Rscript analysis_wilson_temp.R <directory> [skip_steps]")
    }

    directory <- args[1]
    skip_steps <- if (length(args) == 2) as.integer(args[2]) else 0

    logs_dir <- file.path(directory, "logs")
    if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
    assign("WF_LOG_FILE", file.path(logs_dir, "analysis_wilson_temp.log"), envir = .GlobalEnv)

    result <- fit_static_potential(directory, skip_steps)
    write_log("=== Analysis completed successfully ===")

    if (!is.null(result) && !is.null(result$fit_params)) {
        cat("\n=== Static Potential Fit Results ===\n")
        cat(sprintf("A     = %.6f ± %.6f\n", result$fit_params["A"], result$fit_errors["A"]))
        cat(sprintf("B     = %.6f ± %.6f\n", result$fit_params["B"], result$fit_errors["B"]))
        cat(sprintf("sigma = %.6f ± %.6f\n", result$fit_params["sigma"], result$fit_errors["sigma"]))
        cat("\n=== Scale Setting (Sommer) ===\n")
        cat(sprintf("r_0   = %.6f ± %.6f (lattice units)\n", result$r_0, result$r_0_error))
    }
}
