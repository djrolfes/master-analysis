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
source("data_io.R")

# Simple logging helper
write_log <- function(msg) {
    logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
    if (is.na(logfile)) logfile <- "wilson_temp_analysis.log"
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- sprintf("[%s] %s\n", timestamp, msg)
    cat(entry, file = logfile, append = TRUE)
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

    # Group by L and T, calculate mean and error using uwerr (Gamma method)
    # This properly accounts for autocorrelation (e.g., from topological freezing)
    fit_data <- w_temp_data %>%
        group_by(L, T) %>%
        summarise(
            mean = mean(W_temp, na.rm = TRUE),
            error = tryCatch(
                {
                    # Use uwerr (Gamma method) for proper autocorrelation analysis
                    W_clean <- W_temp[!is.na(W_temp)]
                    if (length(W_clean) > 10) {
                        # nrep should be the total length for a single replica
                        uwerr_result <- uwerr(data = W_clean, nrep = length(W_clean), S = 1.5, pl = FALSE)
                        # Extract integrated autocorrelation time
                        tau_int <- if (!is.null(uwerr_result$tauint)) uwerr_result$tauint else 0.5

                        # Log autocorrelation information (especially important with topological freezing)
                        write_log(paste0(
                            "fit_static_potential: L=", L[1], ", T=", T[1],
                            " tau_int=", round(tau_int, 2),
                            ", n_eff=", round(length(W_clean) / (2 * tau_int), 1),
                            ", error=", format(uwerr_result$dvalue, scientific = TRUE, digits = 3)
                        ))

                        # Return the uwerr error which includes autocorrelation correction
                        uwerr_result$dvalue
                    } else {
                        sd(W_clean) / sqrt(length(W_clean))
                    }
                },
                error = function(e) {
                    write_log(paste0(
                        "fit_static_potential: ERROR in uwerr for L=", L[1], ", T=", T[1],
                        ": ", conditionMessage(e)
                    ))
                    sd(W_temp, na.rm = TRUE) / sqrt(sum(!is.na(W_temp)))
                }
            ),
            n_points = n(),
            .groups = "drop"
        )

    write_log(paste0("fit_static_potential: prepared ", nrow(fit_data), " data points for fitting"))

    if (nrow(fit_data) == 0) {
        write_log("fit_static_potential: ERROR - no valid data points found")
        return(NULL)
    }

    # Print summary of L and T values
    write_log(paste0("fit_static_potential: L values: ", paste(sort(unique(fit_data$L)), collapse = ", ")))
    write_log(paste0("fit_static_potential: T values: ", paste(sort(unique(fit_data$T)), collapse = ", ")))

    # For each L value, fit the exponential decay <W(L,T)> = a * exp(-V(L)*T)
    # This gives us V(L) for each L

    unique_L <- sort(unique(fit_data$L))
    v_of_L_data <- data.frame()
    v_of_L_data_nls <- data.frame() # For comparison with plain nls
    wilson_fits <- list() # Store fit objects for plotting (hadron)
    wilson_fits_nls <- list() # Store fit objects for plotting (plain nls)

    for (L_val in unique_L) {
        L_data <- fit_data %>%
            filter(L == L_val, T > 1) %>% # Exclude T=1 from fitting
            arrange(T)

        if (nrow(L_data) < 2) {
            write_log(paste0("fit_static_potential: skipping L=", L_val, " (insufficient data points)"))
            next
        }

        write_log(paste0("fit_static_potential: fitting L=", L_val, " with ", nrow(L_data), " T values (T>1)"))
        write_log(paste0("fit_static_potential: L=", L_val, " T range: [", min(L_data$T), ", ", max(L_data$T), "]"))

        # Perform exponential fit: W(T) = a * exp(-V*T)
        # Check for numerical issues first
        if (any(!is.finite(L_data$mean)) || any(!is.finite(L_data$error))) {
            write_log(paste0("fit_static_potential: L=", L_val, " has non-finite values, skipping"))
            next
        }

        tryCatch(
            {
                # Initial parameter guesses
                # a ~ W(T=1), V ~ -log(W(T_max)/W(T_min)) / (T_max - T_min)
                T_range <- range(L_data$T)
                W_at_Tmin <- L_data$mean[which.min(L_data$T)]
                W_at_Tmax <- L_data$mean[which.max(L_data$T)]

                if (W_at_Tmax > 0 && W_at_Tmin > 0) {
                    V_init <- -log(W_at_Tmax / W_at_Tmin) / (T_range[2] - T_range[1])
                    a_init <- W_at_Tmin * exp(V_init * T_range[1])
                } else {
                    V_init <- 0.1
                    a_init <- mean(L_data$mean)
                }

                # Ensure positive starting values
                if (V_init <= 0) V_init <- 0.1
                if (a_init <= 0) a_init <- abs(mean(L_data$mean))

                write_log(paste0(
                    "fit_static_potential: L=", L_val, " initial params: a=",
                    round(a_init, 6), ", V=", round(V_init, 6)
                ))

                # Use parametric.nlsfit with exponential model
                # W(T) = a * exp(-V*T)
                fn_exp <- function(par, x, boot.r, ...) {
                    par[1] * exp(-par[2] * x)
                }

                fit.result <- parametric.nlsfit(
                    fn = fn_exp,
                    par.guess = c(a_init, V_init),
                    y = L_data$mean,
                    dy = L_data$error,
                    x = L_data$T,
                    boot.R = 200
                )

                # Extract V(L) and its error from bootstrap fit
                # fit.result$t0 contains [a, V]
                a_fit <- fit.result$t0[1]
                V_fit <- fit.result$t0[2]
                a_error <- fit.result$se[1]
                V_error <- fit.result$se[2]

                # Get chi²/dof for this fit
                chisqr_per_dof <- if (!is.null(fit.result$chisqr) && !is.null(fit.result$dof)) {
                    fit.result$chisqr / fit.result$dof
                } else {
                    NA
                }

                write_log(paste0(
                    "fit_static_potential: L=", L_val, " fit results: V=",
                    round(V_fit, 6), " ± ", round(V_error, 6),
                    ", chi2/dof=", round(chisqr_per_dof, 3),
                    " (chi2=", round(fit.result$chisqr, 2), ", dof=", fit.result$dof, ")"
                ))

                v_of_L_data <- rbind(v_of_L_data, data.frame(
                    L = L_val,
                    V = V_fit,
                    V_error = V_error,
                    a = a_fit,
                    a_error = a_error
                ))

                # Store fit parameters and their errors for plotting
                wilson_fits[[as.character(L_val)]] <- list(
                    a = a_fit,
                    V = V_fit,
                    a_error = a_error,
                    V_error = V_error,
                    T_range = range(L_data$T)
                )

                # Also fit with plain nls for comparison
                tryCatch(
                    {
                        nls_fit <- nls(mean ~ a * exp(-V * T),
                            data = L_data,
                            start = list(a = a_init, V = V_init),
                            weights = 1 / (L_data$error^2)
                        )

                        nls_summary <- summary(nls_fit)
                        V_fit_nls <- coef(nls_fit)["V"]
                        a_fit_nls <- coef(nls_fit)["a"]
                        V_error_nls <- nls_summary$coefficients["V", "Std. Error"]
                        a_error_nls <- nls_summary$coefficients["a", "Std. Error"]

                        write_log(paste0(
                            "fit_static_potential: L=", L_val, " plain nls results: V=",
                            round(V_fit_nls, 6), " ± ", round(V_error_nls, 6)
                        ))

                        v_of_L_data_nls <- rbind(v_of_L_data_nls, data.frame(
                            L = L_val,
                            V = V_fit_nls,
                            V_error = V_error_nls,
                            a = a_fit_nls,
                            a_error = a_error_nls
                        ))

                        wilson_fits_nls[[as.character(L_val)]] <- list(
                            a = a_fit_nls,
                            V = V_fit_nls,
                            a_error = a_error_nls,
                            V_error = V_error_nls,
                            T_range = range(L_data$T)
                        )
                    },
                    error = function(e) {
                        write_log(paste0("fit_static_potential: ERROR in plain nls for L=", L_val, ": ", conditionMessage(e)))
                    }
                )
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

            # Bootstrap fit using parametric.nlsfit
            potential_fit <- parametric.nlsfit(
                fn = fn_potential,
                par.guess = c(A_init, B_init, sigma_init),
                boot.R = 200,
                y = v_of_L_data$V,
                dy = v_of_L_data$V_error,
                x = v_of_L_data$L
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

            # Now fit V(L) using plain nls for comparison
            params_nls <- NULL
            param_errors_nls <- NULL
            r_0_nls <- NA
            r_0_error_nls <- NA
            a_fm_nls <- NA
            a_fm_error_nls <- NA

            if (nrow(v_of_L_data_nls) > 0) {
                write_log("fit_static_potential: fitting V(L) = A + B/L + sigma*L with plain nls")

                tryCatch(
                    {
                        nls_potential_fit <- nls(V ~ A + B / L + sigma * L,
                            data = v_of_L_data_nls,
                            start = list(A = A_init, B = B_init, sigma = sigma_init),
                            weights = 1 / (v_of_L_data_nls$V_error^2)
                        )

                        nls_pot_summary <- summary(nls_potential_fit)
                        params_nls <- coef(nls_potential_fit)
                        names(params_nls) <- c("A", "B", "sigma")
                        param_errors_nls <- nls_pot_summary$coefficients[, "Std. Error"]
                        names(param_errors_nls) <- c("A", "B", "sigma")

                        # Calculate chi2/dof for nls fit
                        residuals_nls <- v_of_L_data_nls$V - predict(nls_potential_fit)
                        weights_nls <- 1 / (v_of_L_data_nls$V_error^2)
                        chisq_nls <- sum(weights_nls * residuals_nls^2)
                        dof_nls <- nrow(v_of_L_data_nls) - 3

                        write_log(paste0("fit_static_potential: V(L) plain nls results:"))
                        write_log(paste0("  A = ", round(params_nls["A"], 6), " ± ", round(param_errors_nls["A"], 6)))
                        write_log(paste0("  B = ", round(params_nls["B"], 6), " ± ", round(param_errors_nls["B"], 6)))
                        write_log(paste0("  sigma = ", round(params_nls["sigma"], 6), " ± ", round(param_errors_nls["sigma"], 6)))
                        write_log(paste0("  chi2/dof = ", round(chisq_nls / dof_nls, 3)))

                        # Calculate r_0 for nls
                        r_0_nls <- sqrt((1.65 + params_nls["B"]) / params_nls["sigma"])

                        # Error propagation for r_0: r_0 = sqrt((1.65 + B) / sigma)
                        # dr_0/dB = 1/(2*sqrt((1.65+B)/sigma)) * 1/sigma
                        # dr_0/dsigma = 1/(2*sqrt((1.65+B)/sigma)) * (-(1.65+B)/sigma^2)
                        dr0_dB <- 1 / (2 * r_0_nls * params_nls["sigma"])
                        dr0_dsigma <- -(1.65 + params_nls["B"]) / (2 * r_0_nls * params_nls["sigma"]^2)
                        r_0_error_nls <- sqrt((dr0_dB * param_errors_nls["B"])^2 + (dr0_dsigma * param_errors_nls["sigma"])^2)

                        a_fm_nls <- 0.5 / r_0_nls
                        a_fm_error_nls <- a_fm_nls * (r_0_error_nls / r_0_nls)

                        write_log(paste0("  r_0 = ", round(r_0_nls, 6), " ± ", round(r_0_error_nls, 6)))
                        write_log(paste0("  lattice spacing a = ", round(a_fm_nls, 6), " ± ", round(a_fm_error_nls, 6), " fm"))

                        # Append nls results to scale setting file
                        scale_content_nls <- sprintf(
                            "\n\n# Plain NLS Fit Results (for comparison)\n# V(L) = A + B/L + sigma*L\n\n# Fit Parameters (NLS):\nA_nls = %.6f ± %.6f\nB_nls = %.6f ± %.6f\nsigma_nls = %.6f ± %.6f\n\n# Fit Quality (NLS):\nchi2/dof_nls = %.3f\n\n# Scale Setting (NLS):\nr_0_nls = %.6f ± %.6f (lattice units)\nlattice_spacing_a_nls = %.6f ± %.6f fm\n",
                            params_nls["A"], param_errors_nls["A"],
                            params_nls["B"], param_errors_nls["B"],
                            params_nls["sigma"], param_errors_nls["sigma"],
                            chisq_nls / dof_nls,
                            r_0_nls, r_0_error_nls,
                            a_fm_nls, a_fm_error_nls
                        )

                        write(scale_content_nls, file = scale_file, append = TRUE)
                        write_log(paste0("fit_static_potential: appended nls results to ", scale_file))
                    },
                    error = function(e) {
                        write_log(paste0("fit_static_potential: ERROR in plain nls V(L) fit: ", conditionMessage(e)))
                    }
                )
            }

            # Create plots
            # Plot 2: All Wilson loops on one plot with fits (using ggplot2)
            # Create fit curves with error bands for each L
            fit_curves_list <- lapply(names(wilson_fits), function(L_val) {
                fit_params <- wilson_fits[[L_val]]
                T_range <- seq(fit_params$T_range[1], fit_params$T_range[2], length.out = 50)

                # Central fit: W = a * exp(-V * T)
                W_fit <- fit_params$a * exp(-fit_params$V * T_range)

                # Error propagation: W = a * exp(-V * T)
                # dW/da = exp(-V * T)
                # dW/dV = -a * T * exp(-V * T)
                # sigma_W^2 = (dW/da)^2 * sigma_a^2 + (dW/dV)^2 * sigma_V^2
                dW_da <- exp(-fit_params$V * T_range)
                dW_dV <- -fit_params$a * T_range * exp(-fit_params$V * T_range)
                W_error <- sqrt((dW_da * fit_params$a_error)^2 + (dW_dV * fit_params$V_error)^2)

                data.frame(
                    L = as.numeric(L_val),
                    T = T_range,
                    W_fit = W_fit,
                    W_lower = W_fit - W_error,
                    W_upper = W_fit + W_error
                )
            })
            fit_curves_data <- do.call(rbind, fit_curves_list)

            # Filter out T=1 for plotting
            fit_data_plot <- fit_data %>% filter(T > 1)

            wilson_plot <- ggplot() +
                geom_ribbon(
                    data = fit_curves_data,
                    aes(x = T, ymin = W_lower, ymax = W_upper, fill = factor(L), group = L),
                    alpha = 0.2
                ) +
                geom_line(
                    data = fit_curves_data,
                    aes(x = T, y = W_fit, color = factor(L), group = L),
                    linewidth = 1, alpha = 0.9
                ) +
                geom_point(
                    data = fit_data_plot,
                    aes(x = T, y = mean, color = factor(L)),
                    size = 2.5, alpha = 1
                ) +
                geom_errorbar(
                    data = fit_data_plot,
                    aes(x = T, y = mean, ymin = mean - error, ymax = mean + error, color = factor(L)),
                    width = 0.15, alpha = 0.9
                ) +
                scale_y_log10() +
                scale_color_brewer(palette = "Set1") +
                scale_fill_brewer(palette = "Set1") +
                labs(
                    title = "Temporal Wilson Loops <W(L,T)> vs T (T > 1)",
                    x = "Temporal Extent T",
                    y = "<W(L,T)> (log scale)",
                    color = "Spatial\nSeparation L",
                    fill = "Spatial\nSeparation L"
                ) +
                theme_minimal(base_size = 14) +
                theme(legend.position = "right")

            # Save plots to PDF
            out_pdf <- file.path(directory, "static_potential_fit.pdf")
            write_log(paste0("fit_static_potential: saving PDF to ", out_pdf))

            pdf(out_pdf, width = 10, height = 6)

            # Plot 1: V(L) fit using hadron's plot (better for small errors)
            # Extend plot range to include r_0 (whether r_0 > L_max or r_0 < L_max)
            L_range <- range(v_of_L_data$L)
            # Ensure we extend to at least r_0, and beyond if we have data beyond r_0
            plot_range <- c(L_range[1], max(L_range[2], r_0 * 1.1))

            # Calculate y-range that includes both data points and fit curve at extended range
            # V(L) = A - B/L + sigma*L
            L_extended <- seq(plot_range[1], plot_range[2], length.out = 100)
            V_extended <- params["A"] - params["B"] / L_extended + params["sigma"] * L_extended
            y_data_range <- range(c(v_of_L_data$V - v_of_L_data$V_error, v_of_L_data$V + v_of_L_data$V_error))
            y_fit_range <- range(V_extended)
            y_range <- range(c(y_data_range, y_fit_range))
            # Add 5% padding to y-range
            y_padding <- 0.05 * diff(y_range)
            y_range <- y_range + c(-y_padding, y_padding)

            write_log(paste0("fit_static_potential: plot_range = [", plot_range[1], ", ", plot_range[2], "], r_0 = ", round(r_0, 4)))
            write_log(paste0("fit_static_potential: ylim = [", round(y_range[1], 4), ", ", round(y_range[2], 4), "]"))

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
                col.band = "gray",
                opacity.band = 0.5,
                lwd = 2,
                plot.range = plot_range,
                xlim = plot_range,
                ylim = y_range
            )

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

            # Plot 2: Combined Wilson loops
            print(wilson_plot)

            dev.off()

            write_log(paste0("fit_static_potential: successfully saved ", out_pdf))

            # Create separate PDF for plain nls fits if available
            if (!is.null(params_nls) && nrow(v_of_L_data_nls) > 0) {
                out_pdf_nls <- file.path(directory, "static_potential_fit_nonhadron.pdf")
                write_log(paste0("fit_static_potential: saving nls PDF to ", out_pdf_nls))

                pdf(out_pdf_nls, width = 10, height = 6)

                # Plot 1: V(L) fit using plain nls (match hadron plot style)
                # Use same plot range as hadron plot
                L_range_nls <- range(v_of_L_data_nls$L)
                plot_range_nls <- c(L_range_nls[1], max(L_range_nls[2], r_0_nls * 1.1))

                # Calculate fit curve over extended range
                L_extended_nls <- seq(plot_range_nls[1], plot_range_nls[2], length.out = 100)
                V_extended_nls <- params_nls["A"] + params_nls["B"] / L_extended_nls + params_nls["sigma"] * L_extended_nls

                # Calculate y-range to match hadron plot style
                y_data_range_nls <- range(c(
                    v_of_L_data_nls$V - v_of_L_data_nls$V_error,
                    v_of_L_data_nls$V + v_of_L_data_nls$V_error
                ))
                y_fit_range_nls <- range(V_extended_nls)
                y_range_nls <- range(c(y_data_range_nls, y_fit_range_nls))
                y_padding_nls <- 0.05 * diff(y_range_nls)
                y_range_nls <- y_range_nls + c(-y_padding_nls, y_padding_nls)

                # Plot data points with error bars
                plot(v_of_L_data_nls$L, v_of_L_data_nls$V,
                    xlim = plot_range_nls,
                    ylim = y_range_nls,
                    xlab = "Spatial Separation L (lattice units)",
                    ylab = "Potential V(L)",
                    main = sprintf(
                        "Static Potential V(L) = A + B/L + sigma*L (Plain NLS)\nA=%.4f±%.4f, B=%.4f±%.4f, sigma=%.4f±%.4f\nr_0=%.4f±%.4f, a=%.4f±%.4f fm",
                        params_nls["A"], param_errors_nls["A"],
                        params_nls["B"], param_errors_nls["B"],
                        params_nls["sigma"], param_errors_nls["sigma"],
                        r_0_nls, r_0_error_nls,
                        a_fm_nls, a_fm_error_nls
                    ),
                    pch = 19,
                    col = "black",
                    cex = 1.2
                )

                # Add error bars
                arrows(v_of_L_data_nls$L, v_of_L_data_nls$V - v_of_L_data_nls$V_error,
                    v_of_L_data_nls$L, v_of_L_data_nls$V + v_of_L_data_nls$V_error,
                    length = 0.05, angle = 90, code = 3, col = "black", lwd = 1.5
                )

                # Add fit curve with gray confidence band (approximate)
                # Calculate error band using parameter uncertainties
                V_error_band <- sqrt(
                    param_errors_nls["A"]^2 +
                        (param_errors_nls["B"] / L_extended_nls)^2 +
                        (param_errors_nls["sigma"] * L_extended_nls)^2
                )

                polygon(c(L_extended_nls, rev(L_extended_nls)),
                    c(V_extended_nls - V_error_band, rev(V_extended_nls + V_error_band)),
                    col = rgb(0.5, 0.5, 0.5, alpha = 0.5), border = NA
                )

                # Add fit curve on top
                lines(L_extended_nls, V_extended_nls, col = "black", lwd = 2)

                # Add r_0 marker if valid
                if (!is.na(r_0_nls)) {
                    rect(r_0_nls - r_0_error_nls, par("usr")[3],
                        r_0_nls + r_0_error_nls, par("usr")[4],
                        col = rgb(0, 0, 1, alpha = 0.15), border = NA
                    )
                    abline(v = r_0_nls, col = "blue", lwd = 1, lty = 1)
                    text(r_0_nls, par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
                        sprintf("r_0 = %.3f ± %.3f", r_0_nls, r_0_error_nls),
                        col = "blue", pos = 4, cex = 0.9
                    )
                }

                # Plot 2: Wilson loops with nls fits
                if (length(wilson_fits_nls) > 0) {
                    fit_curves_list_nls <- lapply(names(wilson_fits_nls), function(L_val) {
                        fit_params <- wilson_fits_nls[[L_val]]
                        T_range <- seq(fit_params$T_range[1], fit_params$T_range[2], length.out = 50)

                        W_fit <- fit_params$a * exp(-fit_params$V * T_range)
                        dW_da <- exp(-fit_params$V * T_range)
                        dW_dV <- -fit_params$a * T_range * exp(-fit_params$V * T_range)
                        W_error <- sqrt((dW_da * fit_params$a_error)^2 + (dW_dV * fit_params$V_error)^2)

                        data.frame(
                            L = as.numeric(L_val),
                            T = T_range,
                            W_fit = W_fit,
                            W_lower = W_fit - W_error,
                            W_upper = W_fit + W_error
                        )
                    })
                    fit_curves_data_nls <- do.call(rbind, fit_curves_list_nls)

                    wilson_plot_nls <- ggplot() +
                        geom_ribbon(
                            data = fit_curves_data_nls,
                            aes(x = T, ymin = W_lower, ymax = W_upper, fill = factor(L), group = L),
                            alpha = 0.2
                        ) +
                        geom_line(
                            data = fit_curves_data_nls,
                            aes(x = T, y = W_fit, color = factor(L), group = L),
                            linewidth = 1, alpha = 0.9
                        ) +
                        geom_point(
                            data = fit_data_plot,
                            aes(x = T, y = mean, color = factor(L)),
                            size = 2.5, alpha = 1
                        ) +
                        geom_errorbar(
                            data = fit_data_plot,
                            aes(x = T, y = mean, ymin = mean - error, ymax = mean + error, color = factor(L)),
                            width = 0.15, alpha = 0.9
                        ) +
                        scale_y_log10() +
                        scale_color_brewer(palette = "Set1") +
                        scale_fill_brewer(palette = "Set1") +
                        labs(
                            title = "Temporal Wilson Loops <W(L,T)> vs T (Plain NLS Fits)",
                            x = "Temporal Extent T",
                            y = "<W(L,T)> (log scale)",
                            color = "Spatial\nSeparation L",
                            fill = "Spatial\nSeparation L"
                        ) +
                        theme_minimal(base_size = 14) +
                        theme(legend.position = "right")

                    print(wilson_plot_nls)
                }

                dev.off()
                write_log(paste0("fit_static_potential: successfully saved ", out_pdf_nls))
            }

            return(list(
                v_of_L_data = v_of_L_data,
                fit_params = params,
                fit_errors = param_errors,
                potential_fit = potential_fit,
                wilson_fits = wilson_fits,
                fit_data = fit_data,
                r_0 = r_0,
                r_0_error = r_0_error,
                a_fm = a_fm,
                a_fm_error = a_fm_error
            ))
        },
        error = function(e) {
            write_log(paste0("fit_static_potential: ERROR fitting V(L) model: ", conditionMessage(e)))
            return(list(v_of_L_data = v_of_L_data, fit_data = fit_data))
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
        cat(sprintf("a     = %.6f ± %.6f fm\n", result$a_fm, result$a_fm_error))
    }
}
if (interactive()) {
    directory <- "../example_data" # Replace with your data directory
    skip_steps <- if (length(args) == 2) as.integer(args[2]) else 0

    assign("WF_LOG_FILE", file.path(directory, "wilson_temp_analysis.log"), envir = .GlobalEnv)
    result <- fit_static_potential(directory, skip_steps)

    if (!is.null(result) && !is.null(result$fit_params)) {
        cat("\n=== Static Potential Fit Results ===\n")
        cat(sprintf("A     = %.6f ± %.6f\n", result$fit_params["A"], result$fit_errors["A"]))
        cat(sprintf("B     = %.6f ± %.6f\n", result$fit_params["B"], result$fit_errors["B"]))
        cat(sprintf("sigma = %.6f ± %.6f\n", result$fit_params["sigma"], result$fit_errors["sigma"]))
        cat("\n=== Scale Setting (Sommer) ===\n")
        cat(sprintf("r_0   = %.6f ± %.6f (lattice units)\n", result$r_0, result$r_0_error))
        cat(sprintf("a     = %.6f ± %.6f fm\n", result$a_fm, result$a_fm_error))
    }
}
