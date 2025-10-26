library(ggplot2)
library(dplyr)
library(hadron)
library(zoo)
source("data_io.R") # relies on analysis/ being the working dir when called

# Simple per-run logger (compatible with other analysis_* scripts)
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}


generate_plot <- function(df_topo, df_wilsonflow, output_path, skip_initial) {
  # Merge dataframes on the HMC step
  df_merged <- df_topo %>%
    inner_join(df_wilsonflow, by = "step") %>%
    filter(step > skip_initial) # Skip initial configurations

  if (nrow(df_merged) == 0) {
    write_log("No data to plot after filtering and merging.")
    return()
  }

  df_merged$topological_charge_rounded <- round(df_merged$topological_charge)
  # Ensure sp_max_deriv is numeric before calculations
  df_merged$sp_max_deriv <- as.numeric(df_merged$sp_max_deriv)
  df_merged$sp_max_init <- as.numeric(df_merged$sp_max_init)

  df_merged$sp_deriv_mut <- abs(df_merged$sp_max_deriv - mean(df_merged$sp_max_deriv, na.rm = TRUE))
  df_merged$sp_max_deriv_diff <- (df_merged$sp_max_init - df_merged$sp_max_deriv) / df_merged$flow_step

  coeff <- max(df_merged$topological_charge_rounded^2, na.rm = TRUE) / max(df_merged$sp_max_deriv_diff, na.rm = TRUE)
  # coeff <- 1 / max(df_merged$sp_max_deriv_diff, na.rm = TRUE)
  # Calculate running average for sp_max_deriv
  df_merged <- df_merged %>%
    mutate(sp_max_deriv_avg = zoo::rollmean(sp_max_deriv, k = 5, fill = NA, align = "left"))

  # Create the plot
  p <- ggplot(df_merged, aes(x = step)) +
    # geom_line(aes(y = sp_max_deriv_avg * coeff, color = "sp_max_deriv_avg"), linetype = "dashed") +
    geom_line(aes(y = topological_charge_rounded^2, color = "Topological Charge^2")) +
    geom_line(aes(y = sp_max_deriv_diff * coeff, color = "sp_max_deriv_diff"), linetype = "dashed") +
    # geom_line(aes(y = sp_max_deriv * coeff, color = "sp_max_deriv")) +
    scale_y_continuous(
      name = "Topological Charge^2",
      sec.axis = sec_axis(~ . / coeff, name = "sp_max_deriv")
    ) +
    labs(x = "step", color = "Legend") +
    theme_minimal()

  # Save the plot
  ggsave(output_path, plot = p)
  write_log(paste0("Plot saved to ", output_path))
}

analyze_wilsonflow_details <- function(directory, skip_initial = 0) {
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

  topo_path <- file.path(directory, topo_filename)
  if (!file.exists(topo_path)) {
    write_log(paste0("ERROR: topological charge file not found: ", topo_path))
    stop(paste0("Topological charge file not found: ", topo_path))
  }

  df_topo <- tryCatch(
    {
      read_data_file(topo_path)
    },
    error = function(e) {
      write_log(paste0("ERROR reading topological charge file: ", conditionMessage(e)))
      stop(e)
    }
  )

  df_wilsonflow <- tryCatch(
    {
      read_wilsonflow_details(directory)
    },
    error = function(e) {
      write_log(paste0("ERROR reading Wilson flow details: ", conditionMessage(e)))
      stop(e)
    }
  )

  # Generate the plot
  output_path <- file.path(directory, "topological_charge_sp_max_deriv_plot.pdf")
  generate_plot(df_topo, df_wilsonflow, output_path, skip_initial)
}

# CLI entrypoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analysis_wilsonflow_details.R <directory> [skip_initial]")
}

directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_initial <- if (length(args) >= 2) as.integer(args[2]) else 0
analyze_wilsonflow_details(directory, skip_initial = skip_initial)
