library(ggplot2)
library(dplyr)
library(tidyr) # Added to use pivot_longer
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory

# Simple logging helper that writes to a per-run logfile set in analyze_wilsonflow
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) {
    # fallback to current working directory
    logfile <- "analysis_debug.log"
  }
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

bootstrap_analysis <- function(data, n_boot = 400) {
  write_log(paste0("bootstrap_analysis: start with n_boot=", n_boot))
  results <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    group_by(flow_time) %>%
    summarise(
      mean = mean(observable),
      error = hadron::bootstrap.meanerror(observable, n_boot),
      .groups = "drop"
    )

  write_log("bootstrap_analysis: completed aggregation")
  # Note: make.names() prepends "X" to numeric column names, so we need to remove it
  # Base R plot with hadron style
  flow_times <- as.numeric(gsub("^X", "", results$flow_time))
  par(bty = "o") # Complete frame
  plot(flow_times, results$mean,
    type = "n",
    xlab = "Wilson Flow Time t", ylab = "Mean Observable",
    main = "Bootstrap Analysis of Wilson Flow Observable"
  )
  # Error band
  polygon(c(flow_times, rev(flow_times)),
    c(results$mean + results$error, rev(results$mean - results$error)),
    col = rgb(128 / 255, 128 / 255, 128 / 255, alpha = 0.65), border = NA
  )
  lines(flow_times, results$mean, col = "black", lwd = 2)
}

plot_heatmap <- function(data) {
  write_log("plot_heatmap: start")
  # Note: make.names() prepends "X" to numeric column names, so we need to remove it
  data_long <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    mutate(flow_time = as.numeric(gsub("^X", "", flow_time)))

  p <- ggplot(data_long, aes(x = flow_time, y = hmc_step, fill = observable)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1) + # Flip the color scale: yellow low, blue high
    labs(title = "Heatmap of Wilson Flow Observables", x = "Wilson Flow Time", y = "HMC Step") +
    theme_minimal()
  write_log("plot_heatmap: completed")
  return(p)
}


# Function to compute average distance to closest integer for each flow_time and plot
compute_avg_dist_to_integer <- function(data, skip_steps = 0, n_boot = 400) {
  write_log(paste0("compute_avg_dist_to_integer: start (skip_steps=", skip_steps, ", n_boot=", n_boot, ")"))
  # Skip initial thermalization steps if requested
  if (skip_steps > 0) {
    data <- data %>% filter(hmc_step > skip_steps)
  }
  write_log(paste0("compute_avg_dist_to_integer: rows after skip = ", nrow(data)))
  avg_dist <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "dist") %>%
    group_by(flow_time) %>%
    summarise(
      mean = mean(dist),
      error = hadron::bootstrap.meanerror(dist, n_boot),
      .groups = "drop"
    )
  write_log("compute_avg_dist_to_integer: completed aggregation")
  print(head(avg_dist))
  # Note: make.names() prepends "X" to numeric column names, so we need to remove it
  # Base R plot with hadron style
  flow_times <- as.numeric(gsub("^X", "", avg_dist$flow_time))
  par(bty = "o") # Complete frame
  plot(flow_times, avg_dist$mean,
    type = "n",
    xlab = "Wilson Flow Time t",
    ylab = "Mean Distance ± Error [dimensionless]",
    main = "Average Distance to Closest Integer (Topological Charge)"
  )
  # Error band
  polygon(c(flow_times, rev(flow_times)),
    c(avg_dist$mean + avg_dist$error, rev(avg_dist$mean - avg_dist$error)),
    col = rgb(128 / 255, 128 / 255, 128 / 255, alpha = 0.65), border = NA
  )
  lines(flow_times, avg_dist$mean, col = "black", lwd = 2)
  avg_dist_plot <- recordPlot()
  write_log("compute_avg_dist_to_integer: plot created")
  return(list(data = avg_dist, plot = avg_dist_plot))
}

analyze_combined_action_density <- function(directory, skip_steps = 200, n_boot = 400, target_ad_ft2 = 0.1, target_w = 0.1) {
  write_log(paste0("analyze_combined_action_density: start for directory=", directory))
  # Read action density data using data_io.R convention
  action_density_collection <- read_action_densities(directory)

  write_log(paste0("analyze_combined_action_density: found ", length(action_density_collection), " collections"))
  print(summary(action_density_collection))

  # Loop over each element in the collection and perform analysis
  results_list <- lapply(names(action_density_collection), function(name) {
    write_log(paste0("analyze_combined_action_density: processing ", name))
    action_density_data <- action_density_collection[[name]]
    # Skip initial thermalization steps if requested
    if (skip_steps > 0) {
      action_density_data <- action_density_data %>% filter(hmc_step > skip_steps)
    }

    write_log(paste0("analyze_combined_action_density: rows after skip for ", name, " = ", nrow(action_density_data)))

    # ===== ANALYSIS 1: Bootstrap analysis for action density =====
    results <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      group_by(flow_time) %>%
      summarise(
        mean = mean(action_density),
        error = hadron::bootstrap.meanerror(action_density, n_boot)
      )

    results <- results %>%
      mutate(flow_time_num = as.numeric(gsub("^X", "", flow_time)))

    # ===== ANALYSIS 2: Bootstrap analysis for t^2 * E(t) =====
    boot_results <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      mutate(
        flow_time_num = as.numeric(gsub("^X", "", flow_time)),
        ad_ft2 = action_density * flow_time_num**2
      ) %>%
      group_by(flow_time) %>%
      summarise(
        mean_ad_ft2 = mean(ad_ft2),
        error_ad_ft2 = hadron::bootstrap.meanerror(ad_ft2, n_boot)
      ) %>%
      mutate(flow_time_num = as.numeric(gsub("^X", "", flow_time)))

    # ===== ANALYSIS 3: W(t) = t * d/dt(t^2 <E(t)>) using central differences =====
    # Using product rule: W(t) = 2t^2 <E(t)> + t^3 * d<E(t)>/dt
    data_long <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      mutate(flow_time_num = as.numeric(gsub("^X", "", flow_time))) %>%
      arrange(hmc_step, flow_time_num)

    w_data <- data_long %>%
      group_by(hmc_step) %>%
      arrange(flow_time_num) %>%
      mutate(
        dE_dt = case_when(
          row_number() == 1 ~ (lead(action_density, 1) - action_density) / (lead(flow_time_num, 1) - flow_time_num),
          row_number() == n() ~ (action_density - lag(action_density, 1)) / (flow_time_num - lag(flow_time_num, 1)),
          TRUE ~ (lead(action_density, 1) - lag(action_density, 1)) / (lead(flow_time_num, 1) - lag(flow_time_num, 1))
        ),
        w_t = 2 * flow_time_num^2 * action_density + flow_time_num^3 * dE_dt
      ) %>%
      ungroup()

    w_results <- w_data %>%
      group_by(flow_time_num) %>%
      summarise(
        mean_w = mean(w_t, na.rm = TRUE),
        error_w = hadron::bootstrap.meanerror(w_t[!is.na(w_t)], n_boot),
        .groups = "drop"
      )

    write_log(paste0("analyze_combined_action_density: computed all analyses for ", nrow(w_results), " flow times"))

    # ===== Find crossing points for t^2 * E(t) =====
    target_ad_ft2_reached <- any(boot_results$mean_ad_ft2 >= target_ad_ft2)
    target_ad_ft2_info <- NULL

    if (target_ad_ft2_reached) {
      below_target <- boot_results %>% filter(mean_ad_ft2 < target_ad_ft2)
      above_target <- boot_results %>% filter(mean_ad_ft2 >= target_ad_ft2)

      if (nrow(below_target) > 0 && nrow(above_target) > 0) {
        pt_below <- below_target %>% slice_tail(n = 1)
        pt_above <- above_target %>% slice_head(n = 1)

        t0 <- pt_below$flow_time_num
        t1 <- pt_above$flow_time_num
        y0 <- pt_below$mean_ad_ft2
        y1 <- pt_above$mean_ad_ft2

        flow_time_at_target <- t0 + (target_ad_ft2 - y0) * (t1 - t0) / (y1 - y0)

        err0 <- pt_below$error_ad_ft2
        err1 <- pt_above$error_ad_ft2
        error_y_at_target <- err0 + (target_ad_ft2 - y0) * (err1 - err0) / (y1 - y0)

        slope <- (y1 - y0) / (t1 - t0)
        error_at_target <- error_y_at_target / abs(slope)

        target_ad_ft2_info <- list(
          flow_time = flow_time_at_target,
          error = error_at_target,
          value = target_ad_ft2
        )

        write_log(paste0(
          "analyze_combined_action_density: target t^2*E=", target_ad_ft2, " reached at flow_time = ",
          round(flow_time_at_target, 4), " ± ", round(error_at_target, 4)
        ))
      } else if (nrow(above_target) > 0) {
        pt <- above_target %>% slice_head(n = 1)
        if (nrow(above_target) > 1) {
          pt_next <- above_target %>% slice(2)
          slope <- (pt_next$mean_ad_ft2 - pt$mean_ad_ft2) / (pt_next$flow_time_num - pt$flow_time_num)
          error_at_target <- pt$error_ad_ft2 / abs(slope)
        } else {
          error_at_target <- 0
        }

        target_ad_ft2_info <- list(
          flow_time = pt$flow_time_num,
          error = error_at_target,
          value = pt$mean_ad_ft2
        )
        write_log(paste0(
          "analyze_combined_action_density: target t^2*E=", target_ad_ft2, " reached at first point"
        ))
      }
    } else {
      write_log(paste0("analyze_combined_action_density: target t^2*E=", target_ad_ft2, " not reached"))
    }

    # ===== Find crossing points for W(t) =====
    target_w_reached <- any(w_results$mean_w >= target_w, na.rm = TRUE)
    target_w_info <- NULL

    if (target_w_reached) {
      below_target <- w_results %>% filter(mean_w < target_w)
      above_target <- w_results %>% filter(mean_w >= target_w)

      if (nrow(below_target) > 0 && nrow(above_target) > 0) {
        pt_below <- below_target %>% slice_tail(n = 1)
        pt_above <- above_target %>% slice_head(n = 1)

        t0 <- pt_below$flow_time_num
        t1 <- pt_above$flow_time_num
        y0 <- pt_below$mean_w
        y1 <- pt_above$mean_w

        flow_time_at_target <- t0 + (target_w - y0) * (t1 - t0) / (y1 - y0)

        err0 <- pt_below$error_w
        err1 <- pt_above$error_w
        error_y_at_target <- err0 + (target_w - y0) * (err1 - err0) / (y1 - y0)

        slope <- (y1 - y0) / (t1 - t0)
        error_at_target <- error_y_at_target / abs(slope)

        target_w_info <- list(
          flow_time = flow_time_at_target,
          error = error_at_target,
          value = target_w
        )

        write_log(paste0(
          "analyze_combined_action_density: target W=", target_w, " reached at flow_time = ",
          round(flow_time_at_target, 4), " ± ", round(error_at_target, 4)
        ))
      } else if (nrow(above_target) > 0) {
        pt <- above_target %>% slice_head(n = 1)
        if (nrow(above_target) > 1) {
          pt_next <- above_target %>% slice(2)
          slope <- (pt_next$mean_w - pt$mean_w) / (pt_next$flow_time_num - pt$flow_time_num)
          error_at_target <- pt$error_w / abs(slope)
        } else {
          error_at_target <- 0
        }

        target_w_info <- list(
          flow_time = pt$flow_time_num,
          error = error_at_target,
          value = pt$mean_w
        )
        write_log(paste0(
          "analyze_combined_action_density: target W=", target_w, " reached at first point"
        ))
      }
    } else {
      write_log(paste0("analyze_combined_action_density: target W=", target_w, " not reached"))
    }

    # ===== Find intersection of t^2*E and W curves =====
    intersection_info <- NULL
    # Merge the two datasets
    combined_data <- inner_join(
      boot_results %>% select(flow_time_num, mean_ad_ft2),
      w_results %>% select(flow_time_num, mean_w),
      by = "flow_time_num"
    )

    # Find where they cross (where sign of difference changes)
    combined_data <- combined_data %>%
      mutate(diff = mean_ad_ft2 - mean_w)

    # Look for all sign changes and keep the largest flow_time
    all_intersections <- list()
    for (i in 2:nrow(combined_data)) {
      # Check for NA values before comparing signs
      if (!is.na(combined_data$diff[i - 1]) && !is.na(combined_data$diff[i]) &&
        sign(combined_data$diff[i - 1]) != sign(combined_data$diff[i])) {
        # Found a crossing
        t0 <- combined_data$flow_time_num[i - 1]
        t1 <- combined_data$flow_time_num[i]
        y0_ad <- combined_data$mean_ad_ft2[i - 1]
        y1_ad <- combined_data$mean_ad_ft2[i]
        y0_w <- combined_data$mean_w[i - 1]
        y1_w <- combined_data$mean_w[i]

        # Linear interpolation to find crossing point
        # At crossing: y_ad = y_w, so find t where (y0_ad + slope_ad*(t-t0)) = (y0_w + slope_w*(t-t0))
        slope_ad <- (y1_ad - y0_ad) / (t1 - t0)
        slope_w <- (y1_w - y0_w) / (t1 - t0)

        if (abs(slope_ad - slope_w) > 1e-10) {
          t_intersect <- t0 + (y0_w - y0_ad) / (slope_ad - slope_w)
          y_intersect <- y0_ad + slope_ad * (t_intersect - t0)

          all_intersections <- c(all_intersections, list(list(
            flow_time = t_intersect,
            value = y_intersect
          )))
        }
      }
    }

    # Use the intersection with the largest flow_time
    if (length(all_intersections) > 0) {
      # Find the intersection with maximum flow_time
      max_idx <- which.max(sapply(all_intersections, function(x) x$flow_time))
      intersection_info <- all_intersections[[max_idx]]

      write_log(paste0(
        "analyze_combined_action_density: found ", length(all_intersections), " intersections, using largest at flow_time = ",
        round(intersection_info$flow_time, 4), ", value = ", round(intersection_info$value, 6)
      ))
    } else {
      write_log("analyze_combined_action_density: no intersection found between t^2*E and W curves")
    }

    # ===== CREATE PLOTS =====

    # Plot 1: Action Density (base R with hadron style)
    par(bty = "o") # Complete frame
    plot(results$flow_time_num, results$mean,
      type = "n",
      xlab = "Wilson Flow Time t",
      ylab = "Mean Action Density E(t) ± Error",
      main = paste("Action Density <E(t)>:", name)
    )
    # Error band
    polygon(c(results$flow_time_num, rev(results$flow_time_num)),
      c(results$mean + results$error, rev(results$mean - results$error)),
      col = rgb(128 / 255, 128 / 255, 128 / 255, alpha = 0.65), border = NA
    )
    lines(results$flow_time_num, results$mean, col = "black", lwd = 2)
    action_density_plot <- recordPlot()

    # Plot 2: t^2 * E(t) with target (base R with hadron style)
    plot_title_ad_ft2 <- paste("t² × <E(t)>:", name)

    if (!is.null(target_ad_ft2_info)) {
      plot_title_ad_ft2 <- paste0(plot_title_ad_ft2, sprintf(
        "\nt0(t^2*E=%.3f) = %.4f +/- %.4f",
        target_ad_ft2, target_ad_ft2_info$flow_time, target_ad_ft2_info$error
      ))
    } else {
      plot_title_ad_ft2 <- paste0(plot_title_ad_ft2, sprintf("\n(Target t^2*E=%.3f not reached)", target_ad_ft2))
    }

    par(bty = "o") # Complete frame
    plot(boot_results$flow_time_num, boot_results$mean_ad_ft2,
      type = "n",
      xlab = "Wilson Flow Time t",
      ylab = "t² <E(t)> ± Error",
      main = plot_title_ad_ft2
    )
    # Error band in red
    polygon(c(boot_results$flow_time_num, rev(boot_results$flow_time_num)),
      c(
        boot_results$mean_ad_ft2 + boot_results$error_ad_ft2,
        rev(boot_results$mean_ad_ft2 - boot_results$error_ad_ft2)
      ),
      col = rgb(1, 0, 0, alpha = 0.2), border = NA
    )
    lines(boot_results$flow_time_num, boot_results$mean_ad_ft2, col = "darkred", lwd = 2)

    if (!is.null(target_ad_ft2_info)) {
      abline(h = target_ad_ft2, lty = 2, col = "blue", lwd = 0.8)
      abline(v = target_ad_ft2_info$flow_time, lty = 2, col = "blue", lwd = 0.8)
      points(target_ad_ft2_info$flow_time, target_ad_ft2, col = "blue", pch = 19, cex = 1.5)
      # Horizontal error bar
      arrows(target_ad_ft2_info$flow_time - target_ad_ft2_info$error, target_ad_ft2,
        target_ad_ft2_info$flow_time + target_ad_ft2_info$error, target_ad_ft2,
        angle = 90, code = 3, length = 0.05, col = "blue", lwd = 0.8
      )
    }
    ad_ft2_plot <- recordPlot()

    # Plot 3: W(t) with target (base R with hadron style)
    plot_title_w <- paste("W(t) = t × d/dt(t² <E(t)>):", name)

    if (!is.null(target_w_info)) {
      plot_title_w <- paste0(plot_title_w, sprintf(
        "\nt0(W=%.3f) = %.4f +/- %.4f",
        target_w, target_w_info$flow_time, target_w_info$error
      ))
    } else {
      plot_title_w <- paste0(plot_title_w, sprintf("\n(Target W=%.3f not reached)", target_w))
    }

    par(bty = "o") # Complete frame
    plot(w_results$flow_time_num, w_results$mean_w,
      type = "n",
      xlab = "Wilson Flow Time t",
      ylab = "W(t) ± Error",
      main = plot_title_w
    )
    # Error band in green
    polygon(c(w_results$flow_time_num, rev(w_results$flow_time_num)),
      c(
        w_results$mean_w + w_results$error_w,
        rev(w_results$mean_w - w_results$error_w)
      ),
      col = rgb(0, 1, 0, alpha = 0.2), border = NA
    )
    lines(w_results$flow_time_num, w_results$mean_w, col = "darkgreen", lwd = 2)

    if (!is.null(target_w_info)) {
      abline(h = target_w, lty = 2, col = "blue", lwd = 0.8)
      abline(v = target_w_info$flow_time, lty = 2, col = "blue", lwd = 0.8)
      points(target_w_info$flow_time, target_w, col = "blue", pch = 19, cex = 1.5)
      # Horizontal error bar
      arrows(target_w_info$flow_time - target_w_info$error, target_w,
        target_w_info$flow_time + target_w_info$error, target_w,
        angle = 90, code = 3, length = 0.05, col = "blue", lwd = 0.8
      )
    }
    w_plot <- recordPlot()

    # Plot 4: Overlay of t^2*E and W with intersection
    plot_title_overlay <- paste("Overlay: t² <E(t)> and W(t):", name)
    overlay_plot <- ggplot() +
      geom_line(data = boot_results, aes(x = flow_time_num, y = mean_ad_ft2, color = "t² × E(t)"), linewidth = 1) +
      geom_ribbon(data = boot_results, aes(
        x = flow_time_num, ymin = mean_ad_ft2 - error_ad_ft2,
        ymax = mean_ad_ft2 + error_ad_ft2
      ), fill = "red", alpha = 0.1) +
      geom_line(data = w_results, aes(x = flow_time_num, y = mean_w, color = "W(t)"), linewidth = 1) +
      geom_ribbon(data = w_results, aes(
        x = flow_time_num, ymin = mean_w - error_w,
        ymax = mean_w + error_w
      ), fill = "green", alpha = 0.1) +
      scale_color_manual(values = c("t² × E(t)" = "darkred", "W(t)" = "darkgreen"))

    if (!is.null(intersection_info)) {
      overlay_plot <- overlay_plot +
        geom_vline(xintercept = intersection_info$flow_time, linetype = "dashed", color = "purple", linewidth = 0.8) +
        annotate("point",
          x = intersection_info$flow_time, y = intersection_info$value,
          color = "purple", size = 3, shape = 18
        )
      plot_title_overlay <- paste0(plot_title_overlay, sprintf(
        "\nIntersection at t = %.4f",
        intersection_info$flow_time
      ))
    }

    overlay_plot <- overlay_plot +
      labs(
        title = plot_title_overlay,
        x = "Wilson Flow Time t",
        y = "Observable Value",
        color = "Observable"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

    # Create clean action density plot without title, with t0 annotation
    clean_ad_plot <- ggplot(boot_results, aes(x = flow_time_num, y = mean_ad_ft2)) +
      geom_line(color = "darkred", linewidth = 1) +
      geom_ribbon(aes(ymin = mean_ad_ft2 - error_ad_ft2, ymax = mean_ad_ft2 + error_ad_ft2), fill = "red", alpha = 0.2)

    if (!is.null(target_ad_ft2_info)) {
      # Calculate text position near bottom of plot
      y_min <- min(boot_results$mean_ad_ft2 - boot_results$error_ad_ft2, na.rm = TRUE)
      y_max <- max(boot_results$mean_ad_ft2 + boot_results$error_ad_ft2, na.rm = TRUE)
      text_y_pos <- y_min + 0.05 * (y_max - y_min)

      clean_ad_plot <- clean_ad_plot +
        geom_hline(yintercept = target_ad_ft2, linetype = "dashed", color = "blue", linewidth = 0.4) +
        geom_vline(xintercept = target_ad_ft2_info$flow_time, linetype = "dashed", color = "blue", linewidth = 0.4) +
        annotate("errorbarh",
          y = target_ad_ft2,
          xmin = target_ad_ft2_info$flow_time - target_ad_ft2_info$error,
          xmax = target_ad_ft2_info$flow_time + target_ad_ft2_info$error,
          color = "blue", height = target_ad_ft2 * 0.05, linewidth = 0.6
        ) +
        annotate("point", x = target_ad_ft2_info$flow_time, y = target_ad_ft2, color = "blue", size = 1.5) +
        annotate("text",
          x = target_ad_ft2_info$flow_time + 0.1,
          y = text_y_pos,
          label = {
            # Format as value(error) with proper significant digits
            # Find the order of magnitude of the error
            error_magnitude <- floor(log10(target_ad_ft2_info$error))
            # Determine decimal places to show (2 significant digits in error)
            decimal_places <- max(0, -error_magnitude)
            # Round value and error to appropriate precision
            value_rounded <- round(target_ad_ft2_info$flow_time, decimal_places)
            error_rounded <- round(target_ad_ft2_info$error * 10^decimal_places)
            sprintf("t0 = %.*f(%d)", decimal_places, value_rounded, as.integer(error_rounded))
          },
          color = "blue", size = 5, hjust = 0, vjust = 0
        )
    }

    clean_ad_plot <- clean_ad_plot +
      labs(
        x = "Wilson Flow Time t",
        y = "t² <E(t)>"
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )

    # Save all 4 plots to a single PDF
    out_pdf <- file.path(directory, paste0("combined_analysis_", sub("\\.txt$", "", name), ".pdf"))
    write_log(paste0("analyze_combined_action_density: saving PDF to ", out_pdf))
    tryCatch(
      {
        pdf(out_pdf, width = 10, height = 12)
        print(ad_ft2_plot)
        print(w_plot)
        print(overlay_plot)
        print(action_density_plot)
        dev.off()
        write_log(paste0("analyze_combined_action_density: successfully saved ", out_pdf))
      },
      error = function(e) {
        write_log(paste0("analyze_combined_action_density: ERROR saving PDF for ", name, ": ", conditionMessage(e)))
      }
    )

    # Save clean action density plot to separate single-page PDF
    out_clean_pdf <- file.path(directory, paste0("action_density_clean_", sub("\\.txt$", "", name), ".pdf"))
    write_log(paste0("analyze_combined_action_density: saving clean PDF to ", out_clean_pdf))
    tryCatch(
      {
        pdf(out_clean_pdf, width = 8, height = 6)
        print(clean_ad_plot)
        dev.off()
        write_log(paste0("analyze_combined_action_density: successfully saved ", out_clean_pdf))
      },
      error = function(e) {
        write_log(paste0("analyze_combined_action_density: ERROR saving clean PDF for ", name, ": ", conditionMessage(e)))
      }
    )

    return(list(
      action_density_data = results,
      ad_ft2_data = boot_results,
      w_data = w_results,
      target_ad_ft2_info = target_ad_ft2_info,
      target_w_info = target_w_info,
      intersection_info = intersection_info,
      name = name
    ))
  })
  write_log("analyze_combined_action_density: completed all collections")
  return(results_list)
}

analyze_action_density <- function(directory, skip_steps = 200, n_boot = 200, target_ad_ft2 = 0.1) {
  write_log(paste0("analyze_action_density: start for directory=", directory, " target_ad_ft2=", target_ad_ft2))
  # Read action density data using data_io.R convention
  action_density_collection <- read_action_densities(directory)

  write_log(paste0("analyze_action_density: found ", length(action_density_collection), " collections"))
  print(summary(action_density_collection))

  # Loop over each element in the collection and perform analysis
  results_list <- lapply(names(action_density_collection), function(name) {
    write_log(paste0("analyze_action_density: processing ", name))
    action_density_data <- action_density_collection[[name]]
    # Skip initial thermalization steps if requested
    if (skip_steps > 0) {
      action_density_data <- action_density_data %>% filter(hmc_step > skip_steps)
    }

    write_log(paste0("analyze_action_density: rows after skip for ", name, " = ", nrow(action_density_data)))

    # Bootstrap analysis for each flow time
    results <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      group_by(flow_time) %>%
      summarise(
        mean = mean(action_density),
        error = hadron::bootstrap.meanerror(action_density, n_boot)
      )

    # Compute mean and error for action_density * flow_time^2 using bootstrap
    # Note: make.names() prepends "X" to numeric column names, so we need to remove it
    results <- results %>%
      mutate(
        flow_time_num = as.numeric(gsub("^X", "", flow_time))
      )

    # For bootstrap, need to resample action_density * flow_time^2
    # Note: make.names() prepends "X" to numeric column names, so we need to remove it
    boot_results <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      mutate(
        flow_time_num = as.numeric(gsub("^X", "", flow_time)),
        ad_ft2 = action_density * flow_time_num**2
      ) %>%
      group_by(flow_time) %>%
      summarise(
        mean_ad_ft2 = mean(ad_ft2),
        error_ad_ft2 = hadron::bootstrap.meanerror(ad_ft2, n_boot)
      )

    # Find flow time where ad_ft2 crosses target value (e.g., 0.1)
    boot_results_num <- boot_results %>%
      mutate(flow_time_num = as.numeric(gsub("^X", "", flow_time)))

    # Check if target value is reached
    target_reached <- any(boot_results_num$mean_ad_ft2 >= target_ad_ft2)
    target_info <- NULL

    if (target_reached) {
      # Find the crossing point by linear interpolation
      # Find the two points that bracket the target
      below_target <- boot_results_num %>% filter(mean_ad_ft2 < target_ad_ft2)
      above_target <- boot_results_num %>% filter(mean_ad_ft2 >= target_ad_ft2)

      if (nrow(below_target) > 0 && nrow(above_target) > 0) {
        # Get the last point below and first point above
        pt_below <- below_target %>% slice_tail(n = 1)
        pt_above <- above_target %>% slice_head(n = 1)

        # Linear interpolation for flow time
        t0 <- pt_below$flow_time_num
        t1 <- pt_above$flow_time_num
        y0 <- pt_below$mean_ad_ft2
        y1 <- pt_above$mean_ad_ft2

        flow_time_at_target <- t0 + (target_ad_ft2 - y0) * (t1 - t0) / (y1 - y0)

        # Error propagation: error on t from error on y
        # At the crossing point, interpolate the error on y
        err0 <- pt_below$error_ad_ft2
        err1 <- pt_above$error_ad_ft2
        error_y_at_target <- err0 + (target_ad_ft2 - y0) * (err1 - err0) / (y1 - y0)

        # Slope at the crossing point (derivative dy/dt)
        slope <- (y1 - y0) / (t1 - t0)

        # Error on t from error on y: Δt = Δy / |slope|
        error_at_target <- error_y_at_target / abs(slope)

        target_info <- list(
          flow_time = flow_time_at_target,
          error = error_at_target,
          value = target_ad_ft2
        )

        write_log(paste0(
          "analyze_action_density: target ", target_ad_ft2, " reached at flow_time = ",
          round(flow_time_at_target, 4), " ± ", round(error_at_target, 4),
          " (slope=", round(slope, 6), ", Δy=", round(error_y_at_target, 6), ")"
        ))
      } else if (nrow(above_target) > 0) {
        # Target is reached at first point
        pt <- above_target %>% slice_head(n = 1)
        # If we have a next point, estimate slope
        if (nrow(above_target) > 1) {
          pt_next <- above_target %>% slice(2)
          slope <- (pt_next$mean_ad_ft2 - pt$mean_ad_ft2) / (pt_next$flow_time_num - pt$flow_time_num)
          error_at_target <- pt$error_ad_ft2 / abs(slope)
        } else {
          error_at_target <- 0 # Cannot estimate without slope
        }

        target_info <- list(
          flow_time = pt$flow_time_num,
          error = error_at_target,
          value = pt$mean_ad_ft2
        )
        write_log(paste0(
          "analyze_action_density: target ", target_ad_ft2, " reached at first point: flow_time = ",
          round(pt$flow_time_num, 4), " ± ", round(error_at_target, 4)
        ))
      }
    } else {
      write_log(paste0("analyze_action_density: target ", target_ad_ft2, " not reached in data"))
    }

    # Plot mean and error vs flow time for action_density
    action_density_plot <- ggplot(results, aes(x = flow_time_num, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.2) +
      # scale_y_log10() +
      labs(
        title = paste("Bootstrap Analysis of Action Density:", name),
        x = "Wilson Flow Time t",
        y = "Mean Action Density E(t) ± Error [lattice units]"
      ) +
      theme_minimal()

    # Plot mean and error vs flow time for action_density * flow_time^2
    ad_ft2_plot <- ggplot(boot_results, aes(x = as.numeric(gsub("^X", "", flow_time)), y = mean_ad_ft2)) +
      geom_line(color = "darkred") +
      geom_ribbon(aes(ymin = mean_ad_ft2 - error_ad_ft2, ymax = mean_ad_ft2 + error_ad_ft2), fill = "red", alpha = 0.2)

    # Add target value indicator if it was reached
    plot_title <- paste("Bootstrap Analysis of Action Density × flow_time²:", name)
    if (!is.null(target_info)) {
      # Add horizontal line at target value
      ad_ft2_plot <- ad_ft2_plot +
        geom_hline(yintercept = target_ad_ft2, linetype = "dashed", color = "blue", linewidth = 0.8) +
        geom_vline(xintercept = target_info$flow_time, linetype = "dashed", color = "blue", linewidth = 0.8) +
        annotate("point", x = target_info$flow_time, y = target_ad_ft2, color = "blue", size = 1) +
        annotate("errorbar",
          y = target_ad_ft2, xmin = target_info$flow_time - target_info$error,
          xmax = target_info$flow_time + target_info$error,
          orientation = "y", color = "blue", width = target_ad_ft2 * 0.1, linewidth = 0.8
        )

      plot_title <- paste0(plot_title, sprintf(
        "\nt0(ad*t^2=%.3f) = %.4f +/- %.4f",
        target_ad_ft2, target_info$flow_time, target_info$error
      ))
    } else {
      plot_title <- paste0(plot_title, sprintf("\n(Target ad*t^2=%.3f not reached)", target_ad_ft2))
    }

    ad_ft2_plot <- ad_ft2_plot +
      labs(
        title = plot_title,
        x = "Wilson Flow Time t",
        y = "t² E(t) ± Error [lattice units]"
      ) +
      theme_minimal()

    # Save both plots to the same PDF
    out_pdf <- file.path(directory, paste0("bootstrap_", sub("\\.txt$", "", name), ".pdf"))
    write_log(paste0("analyze_action_density: saving PDF to ", out_pdf))
    tryCatch(
      {
        pdf(out_pdf, width = 8, height = 10)
        print(action_density_plot)
        print(ad_ft2_plot)
        dev.off()
        write_log(paste0("analyze_action_density: successfully saved ", out_pdf))
      },
      error = function(e) {
        write_log(paste0("analyze_action_density: ERROR saving PDF for ", name, ": ", conditionMessage(e)))
      }
    )

    return(list(data = results, plot = action_density_plot, ad_ft2_data = boot_results, ad_ft2_plot = ad_ft2_plot, name = name))
  })
  write_log("analyze_action_density: completed all collections")
  return(results_list)
}

plot_topological_charge_samples <- function(data, directory, n_samples = 400, skip_initial = 100) {
  write_log(paste0("plot_topological_charge_samples: start (n_samples=", n_samples, ", skip_initial=", skip_initial, ")"))
  # Skip initial M configurations
  if (nrow(data) > skip_initial) {
    data <- data %>% slice((skip_initial + 1):n())
  } else {
    # Not enough data to skip, plot nothing or return a message
    message("Not enough data to skip the specified number of initial configurations. No sample plot generated.")
    write_log("plot_topological_charge_samples: not enough rows after skip; aborting sample plot")
    return(NULL)
  }

  # Handle if total data is lower than M+N
  if (nrow(data) < n_samples) {
    sampled_data <- data
    write_log(paste0("plot_topological_charge_samples: available < n_samples; using all ", nrow(sampled_data), " rows"))
  } else {
    sampled_data <- data %>% sample_n(n_samples)
    write_log(paste0("plot_topological_charge_samples: sampled ", n_samples, " rows"))
  }

  # Pivot to long format for plotting
  # Note: make.names() prepends "X" to numeric column names, so we need to remove it
  data_long <- sampled_data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "topological_charge") %>%
    mutate(flow_time = as.numeric(gsub("^X", "", flow_time)))

  # Determine integer y-values for the dashed lines
  y_range <- range(data_long$topological_charge, na.rm = TRUE)
  y_intercepts <- floor(y_range[1]):ceiling(y_range[2])

  # Create the plot
  sample_plot <- ggplot(data_long, aes(x = flow_time, y = topological_charge, group = hmc_step)) +
    geom_hline(yintercept = y_intercepts, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(alpha = 0.5) +
    labs(
      title = paste("Sample of", nrow(sampled_data), "Topological Charge Configurations"),
      x = "Wilson Flow Time t",
      y = "Topological Charge Q [dimensionless]"
    ) +
    theme_minimal()

  out_file <- file.path(directory, "topological_charge_samples.pdf")
  write_log(paste0("plot_topological_charge_samples: saving to ", out_file))
  tryCatch(
    {
      ggsave(out_file, plot = sample_plot)
      write_log(paste0("plot_topological_charge_samples: saved ", out_file))
    },
    error = function(e) {
      write_log(paste0("plot_topological_charge_samples: ERROR saving plot: ", conditionMessage(e)))
    }
  )

  return(sample_plot)
}

analyze_wilsonflow <- function(directory, skip_steps = 200, target_ad_ft2 = 0.1, target_w = 0.1) {
  # set up logfile for this run
  logs_dir <- file.path(directory, "logs")
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
  assign("WF_LOG_FILE", file.path(logs_dir, "analysis_wilsonflow_tests.log"), envir = .GlobalEnv)
  write_log(paste0("analyze_wilsonflow: start for directory=", directory, " skip_steps=", skip_steps, " target_ad_ft2=", target_ad_ft2, " target_w=", target_w))

  # Wrap entire analysis in tryCatch so errors are logged
  tryCatch(
    {
      # Read the Wilson flow data
      write_log("analyze_wilsonflow: reading topological_charge data")
      topological_charge_data <- read_wilsonflow_data(directory, "topological_charge_cumulative.txt")
      write_log(paste0("analyze_wilsonflow: topological_charge rows = ", nrow(topological_charge_data)))

      # Plot a sample of the topological charge data
      plot_topological_charge_samples(topological_charge_data, directory, n_samples = 400, skip_initial = skip_steps)

      # Compute distance to closest integer for each value
      topological_charge_dist <- topological_charge_data %>%
        mutate(across(-hmc_step, ~ abs(. - round(.))))

      # do a bootstrap
      result <- compute_avg_dist_to_integer(topological_charge_dist, skip_steps = skip_steps)
      avg_dist <- result$data
      avg_dist_plot <- result$plot
      out_avg_file <- file.path(directory, "topological_charge_avg_dist_bootstrap.pdf")
      write_log(paste0("analyze_wilsonflow: saving avg dist plot to ", out_avg_file))
      tryCatch(
        {
          pdf(out_avg_file, width = 8, height = 6)
          replayPlot(avg_dist_plot)
          dev.off()
          write_log(paste0("analyze_wilsonflow: saved avg dist plot to ", out_avg_file))
        },
        error = function(e) {
          write_log(paste0("analyze_wilsonflow: ERROR saving avg dist plot: ", conditionMessage(e)))
        }
      )

      # Combined analysis of action density, t^2*E, and W(t)
      combined_result <- analyze_combined_action_density(directory,
        skip_steps = skip_steps,
        n_boot = 200, target_ad_ft2 = target_ad_ft2,
        target_w = target_w
      )

      write_log("analyze_wilsonflow: completed successfully")
    },
    error = function(e) {
      write_log(paste0("analyze_wilsonflow: ERROR: ", conditionMessage(e)))
      # attempt to capture a traceback
      tb <- tryCatch(
        {
          paste(capture.output(traceback()), collapse = "\n")
        },
        error = function(x) "<traceback failed>"
      )
      write_log(paste0("analyze_wilsonflow: TRACEBACK:\n", tb))
      stop(e)
    }
  )
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: Rscript analysis_wilsonflow_tests.R <directory> <skip_steps> [target_ad_ft2] [target_w]")
}
directory <- args[1]
logs_dir <- file.path(directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
assign("WF_LOG_FILE", file.path(logs_dir, "analysis_wilsonflow_tests.log"), envir = .GlobalEnv)
skip_steps <- as.integer(args[2])
target_ad_ft2 <- if (length(args) >= 3) as.numeric(args[3]) else 0.1
target_w <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
analyze_wilsonflow(directory, skip_steps = skip_steps, target_ad_ft2 = target_ad_ft2, target_w = target_w)
write_log("=== Analysis completed successfully ===")
