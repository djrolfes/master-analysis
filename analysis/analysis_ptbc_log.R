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
  swap_attempts_per_replica_asc <- list()
  swap_attempts_per_replica_desc <- list()
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

      if (ascending) {
        if (length(swap_attempts_per_replica_asc) < replica_index) {
          swap_attempts_per_replica_asc[[replica_index]] <- c(as.integer(accepts[j]))
        } else {
          swap_attempts_per_replica_asc[[replica_index]] <- append(swap_attempts_per_replica_asc[[replica_index]], as.integer(accepts[j]))
        }
      } else {
        if (length(swap_attempts_per_replica_desc) < replica_index) {
          swap_attempts_per_replica_desc[[replica_index]] <- c(as.integer(accepts[j]))
        } else {
          swap_attempts_per_replica_desc[[replica_index]] <- append(swap_attempts_per_replica_desc[[replica_index]], as.integer(accepts[j]))
        }
      }


      total_accepts_per_replica[replica_index] <- total_accepts_per_replica[replica_index] + accepts[j]
    }
  }
  mean_attempts_per_replica_asc <- sapply(swap_attempts_per_replica_asc, mean)
  err_attempts_per_replica_asc <- sapply(swap_attempts_per_replica_asc, function(x) {
    hadron::bootstrap.meanerror(x, R = 400)
  })

  swap_attempts_per_replica_desc <- rev(swap_attempts_per_replica_desc)
  mean_attempts_per_replica_desc <- sapply(swap_attempts_per_replica_desc, mean)
  err_attempts_per_replica_desc <- sapply(swap_attempts_per_replica_desc, function(x) {
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
  # Prepare data for plotting: use only the first (N-1) defects to match swap pairs
  defects_full <- sort(data$defects[[1]])
  num_defects <- length(defects_full)
  Replica_asc_vec <- defects_full[seq_len(num_defects - 1)]
  Replica_desc_vec <- rev(rev(defects_full)[seq_len(num_defects - 1)])
  rep_sorted <- defects_full[seq_len(num_defects)]

  n <- length(rep_sorted) - 1

  # Ensure mean/error vectors are numeric and handle differing lengths
  n_asc <- length(mean_attempts_per_replica_asc)
  n_desc <- length(mean_attempts_per_replica_desc)

  Mean_asc_vec <- rep(NA_real_, n)
  Error_asc_vec <- rep(NA_real_, n)
  Mean_desc_vec <- rep(NA_real_, n)
  Error_desc_vec <- rep(NA_real_, n)

  if (n_asc > 0 && n > 0) {
    idx_asc <- seq_len(n_asc)
    Mean_asc_vec[idx_asc] <- as.numeric(mean_attempts_per_replica_asc[seq_len(n_asc)])
    Error_asc_vec[idx_asc] <- as.numeric(err_attempts_per_replica_asc[seq_len(n_asc)])
  }

  if (n_desc > 0 && n > 0) {
    idx_desc <- seq_len(n_desc)
    Mean_desc_vec[idx_desc] <- as.numeric(mean_attempts_per_replica_desc[seq_len(n_desc)])
    Error_desc_vec[idx_desc] <- as.numeric(err_attempts_per_replica_desc[seq_len(n_desc)])
  }

  plot_data <- data.frame(
    Replica_asc = Replica_asc_vec,
    Replica_desc = Replica_desc_vec,
    Mean_asc = Mean_asc_vec,
    Error_asc = Error_asc_vec,
    Mean_desc = rev(Mean_desc_vec),
    Error_desc = rev(Error_desc_vec),
    stringsAsFactors = FALSE
  )

  weighted_results_asc <- weighted_mean_errors(plot_data$Mean_asc, plot_data$Error_asc)
  weighted_results_desc <- weighted_mean_errors(plot_data$Mean_desc, plot_data$Error_desc)

  # Generate and save the plot
  eps <- 0.005
  x_min <- min(rep_sorted, na.rm = TRUE)
  x_max <- max(rep_sorted, na.rm = TRUE)

  asc_mean <- as.numeric(weighted_results_asc["mean"])
  asc_err <- as.numeric(weighted_results_asc["error"])
  desc_mean <- as.numeric(weighted_results_desc["mean"])
  desc_err <- as.numeric(weighted_results_desc["error"])

  p <- ggplot(plot_data, aes(x = Replica_asc, y = Mean_asc)) +
    # ascending points and errors (shifted left)
    geom_point(aes(x = Replica_asc - eps, y = Mean_asc), size = 1) +
    geom_errorbar(aes(x = Replica_asc - eps, ymin = Mean_asc - Error_asc, ymax = Mean_asc + Error_asc), width = 0.03) +

    # descending points and errors (shifted right)
    geom_point(aes(x = Replica_desc + eps, y = Mean_desc), size = 1, color = "blue") +
    geom_errorbar(aes(x = Replica_desc + eps, ymin = Mean_desc - Error_desc, ymax = Mean_desc + Error_desc), width = 0.03, color = "blue") +

    # full-span shaded bands using annotate (covers the full Replica x-range)
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = asc_mean - asc_err, ymax = asc_mean + asc_err,
      fill = "red", alpha = 0.1
    ) +
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = desc_mean - desc_err, ymax = desc_mean + desc_err,
      fill = "blue", alpha = 0.1
    ) +

    # horizontal lines for the weighted means (span full width by using yintercept param)
    geom_hline(yintercept = asc_mean, linetype = "dashed", color = "red") +
    geom_hline(yintercept = desc_mean, linetype = "dashed", color = "blue") +

    # labels with numeric formatting
    # annotate("text",
    #  x = x_max - 0.1 * (x_max - x_min), y = asc_mean + 0.05,
    #  label = sprintf("Weighted Mean ascending: %.3f ± %.3f", asc_mean, asc_err), color = "red", hjust = 1
    # ) +
    # annotate("text",
    #  x = x_max - 0.1 * (x_max - x_min), y = desc_mean + 0.05,
    #  label = sprintf("Weighted Mean descending: %.3f ± %.3f", desc_mean, desc_err), color = "blue", hjust = 1
    # ) +
    labs(
      title = "Acceptance Rate per Replica",
      x = "Replica / Defect",
      y = "Acceptance Rate"
    ) +
    theme_minimal()

  write_log(sprintf(
    "Weighted Mean Acceptance Rates: ascending = %.3f ± %.3f, descending = %.3f ± %.3f",
    asc_mean, asc_err, desc_mean, desc_err
  ))

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
