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

# Function to analyze acceptance rates
analyze_acceptance <- function(directory, data) {
  if (nrow(data) == 0) {
    write_log("Acceptance analysis skipped: no data.")
    return(NULL)
  }
  num_defects <- length(data$defects[[1]])
  expected_len <- max(0, num_defects - 1)

  # Store the normalized accepts back into a local variable for the rest of the function
  # (the later code will use data$accepts and expects the canonical ordering)
  swap_attempts_per_replica <- list()
  total_accepts_per_replica <- vector("numeric", length = length(data$defects[[1]]))
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    ascending <- as.logical(row$ascending)
    accepts <- unlist(row$accepts)
    sorted_defects <- sort(data$defects[[i]], !ascending)

    if (!is.na(ascending) && !ascending) {
      accepts <- rev(accepts)
    }


    for (j in seq_along(accepts)) {
      replica_index <- j
      partner_index <- j + 1

      if (length(swap_attempts_per_replica) < replica_index) {
        swap_attempts_per_replica[[replica_index]] <- c(as.integer(accepts[j]))
      } else {
        swap_attempts_per_replica[[replica_index]] <- append(swap_attempts_per_replica[[replica_index]], as.integer(accepts[j]))
      }

      if (length(swap_attempts_per_replica) < partner_index) {
        swap_attempts_per_replica[[partner_index]] <- c(as.integer(accepts[j]))
      } else {
        swap_attempts_per_replica[[partner_index]] <- append(swap_attempts_per_replica[[partner_index]], as.integer(accepts[j]))
      }
      total_accepts_per_replica[replica_index] <- total_accepts_per_replica[replica_index] + accepts[j]
      total_accepts_per_replica[partner_index] <- total_accepts_per_replica[partner_index] + accepts[j]
    }
  }
  print(length(swap_attempts_per_replica))
  mean_attempts_per_replica <- sapply(swap_attempts_per_replica, mean)
  err_attempts_per_replica <- sapply(swap_attempts_per_replica, function(x) {
    hadron::bootstrap.meanerror(x, R = 400)
  })
  # print(head(swap_attempts_per_replica))
  # acceptance_rates <- total_accepts_per_replica / swap_attempts_per_replica

  # TODO: save the swaps of the replicas into lists and do a bootstrap
  # accepts_matrix <- do.call(rbind, data$accepts)
  # Use bootstrap.meanError on each column (replica)
  # Increased R from 100 to 400 for better error estimates (best practices)
  # bootstrap_results <- apply(accepts_matrix, 2, function(col) {
  # bootstrap.meanError is suitable for bootstrapping the mean of a vector
  #   result <- hadron::bootstrap.meanerror(col, R = 400)
  # The result is a vector with mean and error
  #  c(mean = mean(col), error = result)
  # })

  # Prepare data for plotting
  plot_data <- data.frame(
    Replica = sort(data$defects[[1]]),
    Mean = mean_attempts_per_replica,
    Error = err_attempts_per_replica
  )

  weighted_results <- weighted_mean_errors(plot_data$Mean, plot_data$Error)

  # Generate and save the plot
  p <- ggplot(plot_data, aes(x = Replica, y = Mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), width = 0.05) +
    geom_hline(yintercept = weighted_results["mean"], linetype = "dashed", color = "red") +
    geom_ribbon(aes(ymin = weighted_results["mean"] - weighted_results["error"], ymax = weighted_results["mean"] + weighted_results["error"]), fill = "red", alpha = 0.1) +
    labs(
      title = "Acceptance Rate per Replica",
      x = "Replica / Defect",
      y = "Acceptance Rate"
    )

  ggsave(file.path(directory, "acceptance_by_replica.pdf"), plot = p, width = 8, height = 6)
  write_log("Acceptance analysis plot saved to acceptance_by_replica.pdf")

  return(plot_data)
}

weighted_mean_errors <- function(values, errors) {
  weights <- 1 / (errors^2)
  normalized_weights <- weights / sum(weights)
  weighted_mean <- sum(values * normalized_weights)
  weighted_error <- errors %*% normalized_weights
  return(c(mean = weighted_mean, error = weighted_error))
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
  number_defects <- length(ptbc_data$defects[[1]])
  # Analyze acceptance
  acceptance_results <- analyze_acceptance(directory, ptbc_data)

  if (!is.null(acceptance_results)) {
    write_log("Acceptance analysis results:")
    write_log(paste(capture.output(print(acceptance_results)), collapse = "\n"))
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
