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

bootstrap_analysis <- function(data, n_boot = 200) {
  write_log("bootstrap_analysis: start")
  results <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    group_by(flow_time) %>%
    summarise(
      boot_result = hadron::bootstrap.meanerror(observable, n_boot),
      mean = boot_result["mean"],
      error = boot_result["error"]
    )

  write_log("bootstrap_analysis: completed aggregation")
  ggplot(results, aes(x = as.numeric(flow_time), y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.2) +
    labs(title = "Bootstrap Analysis of Wilson Flow Observable", x = "Wilson Flow Time", y = "Mean Observable") +
    theme_minimal()
}

plot_heatmap <- function(data) {
  write_log("plot_heatmap: start")
  data_long <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    mutate(flow_time = as.numeric(flow_time))

  p <- ggplot(data_long, aes(x = flow_time, y = hmc_step, fill = observable)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1) + # Flip the color scale: yellow low, blue high
    labs(title = "Heatmap of Wilson Flow Observables", x = "Wilson Flow Time", y = "HMC Step") +
    theme_minimal()
  write_log("plot_heatmap: completed")
  return(p)
}


# Function to compute average distance to closest integer for each flow_time and plot
compute_avg_dist_to_integer <- function(data, skip_steps = 0, n_boot = 200) {
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
      error = hadron::bootstrap.meanerror(dist, n_boot)
    )
  write_log("compute_avg_dist_to_integer: completed aggregation")
  print(head(avg_dist))
  avg_dist_plot <- ggplot(avg_dist, aes(x = as.numeric(flow_time), y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.2) +
    labs(
      title = "Average Distance to Closest Integer (Topological Charge)",
      x = "Wilson Flow Time",
      y = "Mean Distance ± Error"
    ) +
    theme_minimal()
  write_log("compute_avg_dist_to_integer: plot created")
  return(list(data = avg_dist, plot = avg_dist_plot))
}

analyze_action_density <- function(directory, skip_steps = 200, n_boot = 200) {
  write_log(paste0("analyze_action_density: start for directory=", directory))
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
    results <- results %>%
      mutate(
        flow_time_num = as.numeric(flow_time)
      )

    # For bootstrap, need to resample action_density * flow_time^2
    boot_results <- action_density_data %>%
      pivot_longer(-hmc_step, names_to = "flow_time", values_to = "action_density") %>%
      mutate(
        flow_time_num = as.numeric(flow_time),
        ad_ft2 = action_density * flow_time_num**2
      ) %>%
      group_by(flow_time) %>%
      summarise(
        mean_ad_ft2 = mean(ad_ft2),
        error_ad_ft2 = hadron::bootstrap.meanerror(ad_ft2, n_boot)
      )

    # Plot mean and error vs flow time for action_density
    action_density_plot <- ggplot(results, aes(x = flow_time_num, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.2) +
      # scale_y_log10() +
      labs(
        title = paste("Bootstrap Analysis of Action Density:", name),
        x = "Wilson Flow Time",
        y = "Mean Action Density ± Error (log scale)"
      ) +
      theme_minimal()

    # Plot mean and error vs flow time for action_density * flow_time^2
    ad_ft2_plot <- ggplot(boot_results, aes(x = as.numeric(flow_time), y = mean_ad_ft2)) +
      geom_line(color = "darkred") +
      geom_ribbon(aes(ymin = mean_ad_ft2 - error_ad_ft2, ymax = mean_ad_ft2 + error_ad_ft2), fill = "red", alpha = 0.2) +
      labs(
        title = paste("Bootstrap Analysis of Action Density × flow_time²:", name),
        x = "Wilson Flow Time",
        y = "Mean (Action Density × flow_time²) ± Error"
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
  data_long <- sampled_data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "topological_charge") %>%
    mutate(flow_time = as.numeric(flow_time))

  # Determine integer y-values for the dashed lines
  y_range <- range(data_long$topological_charge, na.rm = TRUE)
  y_intercepts <- floor(y_range[1]):ceiling(y_range[2])

  # Create the plot
  sample_plot <- ggplot(data_long, aes(x = flow_time, y = topological_charge, group = hmc_step)) +
    geom_hline(yintercept = y_intercepts, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(alpha = 0.5) +
    labs(
      title = paste("Sample of", nrow(sampled_data), "Topological Charge Configurations"),
      x = "Wilson Flow Time",
      y = "Topological Charge"
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

analyze_wilsonflow <- function(directory, skip_steps = 200) {
  # set up logfile for this run
  assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
  write_log(paste0("analyze_wilsonflow: start for directory=", directory, " skip_steps=", skip_steps))

  # Wrap entire analysis in tryCatch so errors are logged
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
          ggsave(out_avg_file, plot = avg_dist_plot)
          write_log(paste0("analyze_wilsonflow: saved avg dist plot to ", out_avg_file))
        },
        error = function(e) {
          write_log(paste0("analyze_wilsonflow: ERROR saving avg dist plot: ", conditionMessage(e)))
        }
      )

      action_density_result <- analyze_action_density(directory, skip_steps = skip_steps)

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
if (length(args) != 2) {
  stop("Usage: Rscript analysis_wilsonflow_tests.R <directory> <skip_steps>")
}
directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_steps <- args[2]
analyze_wilsonflow(directory, skip_steps = as.integer(skip_steps))
