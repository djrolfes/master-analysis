library(ggplot2)
library(dplyr)
library(hadron)
library(zoo)
source("data_io.R")  # relies on analysis/ being the working dir when called

# Simple per-run logger (compatible with other analysis_* scripts)
write_log <- function(msg) {
  logfile <- get0("WF_LOG_FILE", ifnotfound = NA)
  if (is.na(logfile)) logfile <- "analysis_debug.log"
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s\n", timestamp, msg)
  cat(entry, file = logfile, append = TRUE)
}

# Function to analyze acceptance rates
analyze_acceptance <- function(data) {
  if (nrow(data) == 0) {
    write_log("Acceptance analysis skipped: no data.")
    return(NULL)
  }
  
  # The 'accepts' column is a list of vectors; transform into a matrix
  accepts_matrix <- do.call(rbind, lapply(data$accepts, function(x) x[[1]]))
  
  # Use bootstrap.analysis on each column (rank)
  bootstrap_results <- apply(accepts_matrix, 2, function(col) {
    # hadron::bootstrap.analysis expects a function that returns a single value
    boot_func <- function(data) {
      return(mean(data))
    }
    result <- hadron::bootstrap.analysis(data = col, boot.R = 100, boot.fun = boot_func)
    # The result from bootstrap.analysis is a 'cf' object, extract mean and se
    c(mean = result$t0, error = result$se)
  })
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Rank = 0:(ncol(accepts_matrix) - 1),
    Mean = bootstrap_results["mean", ],
    Error = bootstrap_results["error", ]
  )
  
  # Generate and save the plot
  p <- ggplot(plot_data, aes(x = Rank, y = Mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), width = 0.2) +
    labs(
      title = "Acceptance Rate per Rank",
      x = "Rank",
      y = "Acceptance Rate"
    ) +
    theme_minimal()
  
  ggsave("acceptance_by_rank.pdf", plot = p, width = 8, height = 6)
  write_log("Acceptance analysis plot saved to acceptance_by_rank.pdf")
  
  return(plot_data)
}

# Function to track defect swaps
track_defect_swaps <- function(data, config) {
  if (is.null(config$PTBCSimulationLoggingParams$initial_defects)) {
    write_log("Defect tracking skipped: initial_defects not found in YAML.")
    return(NULL)
  }
  
  initial_defects <- config$PTBCSimulationLoggingParams$initial_defects
  if (nrow(data) == 0) {
    write_log("Defect tracking skipped: no data.")
    return(NULL)
  }
  
  num_defects <- length(initial_defects)
  current_defects <- initial_defects
  swap_details <- list()
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    swap_start <- row$swap_start
    accepts <- row$accepts[[1]]
    delta_H_swaps <- row$delta_H_swap[[1]]
    
    temp_defects <- current_defects
    
    for (j in 1:length(accepts)) {
      idx1 <- (swap_start + j - 1) %% num_defects
      idx2 <- (idx1 + 1) %% num_defects
      
      # R uses 1-based indexing
      r_idx1 <- idx1 + 1
      r_idx2 <- idx2 + 1
      
      defect1 <- temp_defects[r_idx1]
      defect2 <- temp_defects[r_idx2]
      accepted <- accepts[j]
      delta_H <- delta_H_swaps[j]
      
      swap_details[[length(swap_details) + 1]] <- list(
        step = row$step,
        swap_attempt = j,
        defect1 = defect1,
        defect2 = defect2,
        accepted = accepted,
        delta_H_swap = delta_H
      )
      
      if (accepted == 1) {
        # Perform the swap
        temp_defects[r_idx1] <- defect2
        temp_defects[r_idx2] <- defect1
      }
    }
    current_defects <- temp_defects
  }
  
  swap_df <- do.call(rbind, lapply(swap_details, as.data.frame))
  write.csv(swap_df, "defect_swap_details.csv", row.names = FALSE)
  write_log("Defect swap details saved to defect_swap_details.csv")
  
  return(swap_df)
}

# Main analysis function for PTBC log
analyze_ptbc_log <- function(directory, skip_initial = 0) {
  write_log(paste("Starting PTBC log analysis for directory:", directory))
  
  # Read config and data
  config <- read_yaml_config(directory)
  ptbc_data <- read_ptbc_simulation_log(directory)
  
  if (skip_initial > 0) {
    ptbc_data <- ptbc_data %>% filter(step > skip_initial)
  }
  
  # Analyze acceptance
  acceptance_results <- analyze_acceptance(ptbc_data)
  
  if (!is.null(acceptance_results)) {
    write_log("Acceptance analysis results:")
    write_log(paste(capture.output(print(acceptance_results)), collapse = "\n"))
  }
  
  # Track defects
  defect_swap_results <- track_defect_swaps(ptbc_data, config)
  if (!is.null(defect_swap_results)) {
    write_log("Defect tracking analysis complete.")
  }
  
  write_log("PTBC log analysis finished.")
}

# CLI entrypoint
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analysis_ptbc_log.R <directory> [skip_initial]")
}

directory <- args[1]
assign("WF_LOG_FILE", file.path(directory, "analysis_debug.log"), envir = .GlobalEnv)
skip_initial <- if (length(args) >= 2) as.integer(args[2]) else 0
analyze_ptbc_log(directory, skip_initial = skip_initial)