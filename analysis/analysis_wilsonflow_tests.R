library(ggplot2)
library(dplyr)
library(tidyr) # Added to use pivot_longer
library(hadron) # Assuming hadron package is installed and available
source("data_io.R") # Assuming data_io.R is in the working directory

bootstrap_analysis <- function(data, n_boot = 200) {
  results <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    group_by(flow_time) %>%
    summarise(
      boot_result = hadron::bootstrap.meanerror(observable, n_boot),
      mean = boot_result["mean"],
      error = boot_result["error"]
    )

  ggplot(results, aes(x = as.numeric(flow_time), y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.2) +
    labs(title = "Bootstrap Analysis of Wilson Flow Observable", x = "Wilson Flow Time", y = "Mean Observable") +
    theme_minimal()
}

plot_heatmap <- function(data) {
  data_long <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "observable") %>%
    mutate(flow_time = as.numeric(flow_time))

  ggplot(data_long, aes(x = flow_time, y = hmc_step, fill = observable)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1) + # Flip the color scale: yellow low, blue high
    labs(title = "Heatmap of Wilson Flow Observables", x = "Wilson Flow Time", y = "HMC Step") +
    theme_minimal()
}


# Function to compute average distance to closest integer for each flow_time and plot
compute_avg_dist_to_integer <- function(data, skip_steps = 0, n_boot = 200) {
  # Skip initial thermalization steps if requested
  if (skip_steps > 0) {
    data <- data %>% filter(hmc_step > skip_steps)
  }
  avg_dist <- data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "dist") %>%
    group_by(flow_time) %>%
    summarise(
      mean = mean(dist),
      error = hadron::bootstrap.meanerror(dist, n_boot)
    )
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
  return(list(data = avg_dist, plot = avg_dist_plot))
}

analyze_action_density <- function(directory, skip_steps = 200, n_boot = 200) {
  # Read action density data using data_io.R convention
  action_density_collection <- read_action_densities(directory)

  print(summary(action_density_collection))

  # Loop over each element in the collection and perform analysis
  results_list <- lapply(names(action_density_collection), function(name) {
    action_density_data <- action_density_collection[[name]]
    # Skip initial thermalization steps if requested
    if (skip_steps > 0) {
      action_density_data <- action_density_data %>% filter(hmc_step > skip_steps)
    }

    # if (grepl("clover", name, ignore.case = TRUE)) {
    #   action_density_data <- action_density_data %>%
    #     mutate(across(-hmc_step, ~ . - 0.5))
    # }

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
    pdf(file.path(directory, paste0("bootstrap_", sub("\\.txt$", "", name), ".pdf")), width = 8, height = 10)
    print(action_density_plot)
    print(ad_ft2_plot)
    dev.off()

    return(list(data = results, plot = action_density_plot, ad_ft2_data = boot_results, ad_ft2_plot = ad_ft2_plot, name = name))
  })
  return(results_list)
}

plot_topological_charge_samples <- function(data, directory, n_samples = 400, skip_initial = 100) {
  # Skip initial M configurations
  if (nrow(data) > skip_initial) {
    data <- data %>% slice((skip_initial + 1):n())
  } else {
    # Not enough data to skip, plot nothing or return a message
    message("Not enough data to skip the specified number of initial configurations. No sample plot generated.")
    return(NULL)
  }

  # Handle if total data is lower than M+N
  if (nrow(data) < n_samples) {
    sampled_data <- data
    message(paste("Number of available configurations after skipping is less than n_samples. Plotting all", nrow(sampled_data), "available configurations."))
  } else {
    sampled_data <- data %>% sample_n(n_samples)
  }

  # Pivot to long format for plotting
  data_long <- sampled_data %>%
    pivot_longer(-hmc_step, names_to = "flow_time", values_to = "topological_charge") %>%
    mutate(flow_time = as.numeric(flow_time))

  # Create the plot
  sample_plot <- ggplot(data_long, aes(x = flow_time, y = topological_charge, group = hmc_step)) +
    geom_line(alpha = 0.5) +
    labs(
      title = paste("Sample of", nrow(sampled_data), "Topological Charge Configurations"),
      x = "Wilson Flow Time",
      y = "Topological Charge"
    ) +
    theme_minimal()

  # Save the plot
  ggsave(file.path(directory, "topological_charge_samples.pdf"), plot = sample_plot)
  
  return(sample_plot)
}

analyze_wilsonflow <- function(directory, skip_steps = 200) {
  # Read the Wilson flow data
  topological_charge_data <- read_wilsonflow_data(directory, "topological_charge_cumulative.txt")

  # Plot a sample of the topological charge data
  plot_topological_charge_samples(topological_charge_data, directory, n_samples = 400, skip_initial = 100)

  # Compute distance to closest integer for each value
  topological_charge_dist <- topological_charge_data %>%
    mutate(across(-hmc_step, ~ abs(. - round(.))))

  # do a bootstrap
  result <- compute_avg_dist_to_integer(topological_charge_dist, skip_steps = skip_steps)
  avg_dist <- result$data
  avg_dist_plot <- result$plot
  ggsave(file.path(directory, "topological_charge_avg_dist_bootstrap.pdf"), plot = avg_dist_plot)


  action_density_result <- analyze_action_density(directory, skip_steps = skip_steps)
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analysis_wilsonflow_tests.R <directory> <skip_steps>")
}
directory <- args[1]
skip_steps <- args[2]
analyze_wilsonflow(directory, skip_steps = as.integer(skip_steps))
