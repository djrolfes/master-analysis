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

# Function to compute per-defect acceptance from rank log files
compute_per_defect_acceptance <- function(directory, data) {
  write_log("Computing per-defect acceptance from rank log files...")

  # Get unique defect values (sorted)
  defects_full <- sort(data$defects[[1]])
  num_defects <- length(defects_full)
  num_steps <- nrow(data)

  write_log(sprintf("Processing %d steps with %d defects/ranks", num_steps, num_defects))

  # Pre-read all rank files to avoid repeated file I/O
  write_log("Pre-reading all rank log files...")
  rank_data_cache <- list()
  for (rank in seq_len(num_defects) - 1) {
    rank_file <- file.path(directory, sprintf("simulation_log.rank%d.txt", rank))
    if (file.exists(rank_file)) {
      rank_data_cache[[as.character(rank)]] <- read_data_file(rank_file)
      write_log(sprintf("Loaded rank %d (%d rows)", rank, nrow(rank_data_cache[[as.character(rank)]])))
    } else {
      write_log(sprintf("Warning: Rank file not found: %s", rank_file))
    }
  }
  write_log("Finished pre-reading rank files")

  # Initialize list to store accepts for each defect
  accepts_per_defect <- vector("list", length = num_defects)
  names(accepts_per_defect) <- as.character(defects_full)

  # For each step in ptbc log
  progress_interval <- max(1, floor(num_steps / 10))
  for (i in seq_len(num_steps)) {
    if (i %% progress_interval == 0) {
      write_log(sprintf("Processing step %d/%d (%.1f%%)", i, num_steps, 100 * i / num_steps))
    }

    step <- data$step[i]
    prev_defects <- unlist(data$prev_defects[i])

    # Perform argsort to determine which rank had which defect value
    # order() gives the indices that would sort the array
    rank_order <- order(prev_defects)

    # For each defect value, find which rank had it
    for (sorted_idx in seq_along(prev_defects)) {
      # rank_order[sorted_idx] gives the original position (1-indexed)
      # The rank is this position minus 1 (0-indexed)
      rank <- rank_order[sorted_idx] - 1
      defect_value <- prev_defects[rank_order[sorted_idx]]

      # Get cached rank data
      rank_key <- as.character(rank)
      if (!rank_key %in% names(rank_data_cache)) {
        next
      }

      rank_data <- rank_data_cache[[rank_key]]
      matching_row <- rank_data[rank_data$step == step, ]

      if (nrow(matching_row) == 1) {
        accept_value <- matching_row$accept
        defect_key <- as.character(defect_value)

        if (is.null(accepts_per_defect[[defect_key]])) {
          accepts_per_defect[[defect_key]] <- c(accept_value)
        } else {
          accepts_per_defect[[defect_key]] <- c(accepts_per_defect[[defect_key]], accept_value)
        }
      }
    }
  }

  write_log("Finished processing all steps")

  # Compute mean and bootstrap errors for each defect
  mean_per_defect <- numeric(num_defects)
  err_per_defect <- numeric(num_defects)

  for (i in seq_len(num_defects)) {
    defect_key <- as.character(defects_full[i])
    accepts <- accepts_per_defect[[defect_key]]

    if (!is.null(accepts) && length(accepts) > 0) {
      mean_per_defect[i] <- mean(accepts)
      err_per_defect[i] <- hadron::bootstrap.meanerror(accepts, R = 400)
      write_log(sprintf(
        "Defect %.3f: %d samples, mean = %.4f ± %.4f",
        defects_full[i], length(accepts), mean_per_defect[i], err_per_defect[i]
      ))
    } else {
      mean_per_defect[i] <- NA
      err_per_defect[i] <- NA
      write_log(sprintf("Defect %.3f: No data found", defects_full[i]))
    }
  }

  return(list(
    defects = defects_full,
    mean = mean_per_defect,
    error = err_per_defect
  ))
}

# Function to analyze acceptance rates
analyze_acceptance <- function(directory, data) {
  if (nrow(data) == 0) {
    write_log("Acceptance analysis skipped: no data.")
    return(NULL)
  }
  num_defects <- length(data$defects[[1]])
  expected_len <- max(0, num_defects - 1)

  # Compute per-defect acceptance from rank logs
  per_defect_results <- compute_per_defect_acceptance(directory, data)

  # Collect swap acceptance data per replica (normalized to ascending order)
  swap_attempts_per_replica <- list()
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    ascending <- as.logical(row$ascending)
    accepts <- unlist(row$accepts)

    # Normalize to ascending order
    if (!is.na(ascending) && !ascending) {
      accepts <- rev(accepts)
    }

    for (j in seq_along(accepts)) {
      replica_index <- j
      if (length(swap_attempts_per_replica) < replica_index) {
        swap_attempts_per_replica[[replica_index]] <- c(as.integer(accepts[j]))
      } else {
        swap_attempts_per_replica[[replica_index]] <- append(swap_attempts_per_replica[[replica_index]], as.integer(accepts[j]))
      }
    }
  }

  # Compute mean and bootstrap errors for swap acceptance
  mean_swap_acceptance <- sapply(swap_attempts_per_replica, mean)
  err_swap_acceptance <- sapply(swap_attempts_per_replica, function(x) {
    hadron::bootstrap.meanerror(x, R = 400)
  })

  # Prepare data for plotting: use only the first (N-1) defects to match swap pairs
  defects_full <- sort(data$defects[[1]])
  num_defects <- length(defects_full)
  rep_sorted <- defects_full[seq_len(num_defects)]

  # Replica positions are midpoints between adjacent defects for swap acceptance
  n <- length(rep_sorted) - 1
  Replica_vec <- numeric(n)
  for (i in seq_len(n)) {
    Replica_vec[i] <- (rep_sorted[i] + rep_sorted[i + 1]) / 2
  }

  # Ensure mean/error vectors are numeric and handle differing lengths
  n_swap <- length(mean_swap_acceptance)
  Mean_swap_vec <- rep(NA_real_, n)
  Error_swap_vec <- rep(NA_real_, n)

  if (n_swap > 0 && n > 0) {
    idx_swap <- seq_len(min(n_swap, n))
    Mean_swap_vec[idx_swap] <- as.numeric(mean_swap_acceptance[idx_swap])
    Error_swap_vec[idx_swap] <- as.numeric(err_swap_acceptance[idx_swap])
  }

  plot_data <- data.frame(
    Replica = Replica_vec,
    Mean_swap = Mean_swap_vec,
    Error_swap = Error_swap_vec,
    stringsAsFactors = FALSE
  )

  # Add per-defect data to plot_data
  plot_data_defects <- data.frame(
    Defect = per_defect_results$defects,
    Mean_defect = per_defect_results$mean,
    Error_defect = per_defect_results$error,
    stringsAsFactors = FALSE
  )

  weighted_results_swap <- weighted_mean_errors(plot_data$Mean_swap, plot_data$Error_swap)
  weighted_results_defect <- weighted_mean_errors(plot_data_defects$Mean_defect, plot_data_defects$Error_defect)

  # Generate and save the plot
  x_min <- min(rep_sorted, na.rm = TRUE)
  x_max <- max(rep_sorted, na.rm = TRUE)

  swap_mean <- as.numeric(weighted_results_swap["mean"])
  swap_err <- as.numeric(weighted_results_swap["error"])
  defect_mean <- as.numeric(weighted_results_defect["mean"])
  defect_err <- as.numeric(weighted_results_defect["error"])

  eps <- 0.005 # Small shift to separate overlapping points

  # Combined plot (both swap and per-defect)
  p_combined <- ggplot() +
    # Swap acceptance points and errors (shifted left slightly)
    geom_point(data = plot_data, aes(x = Replica - eps, y = Mean_swap), size = 1.5, color = "red") +
    geom_errorbar(data = plot_data, aes(x = Replica - eps, ymin = Mean_swap - Error_swap, ymax = Mean_swap + Error_swap), width = 0.03, color = "red") +

    # Per-defect acceptance points and errors (shifted right slightly)
    geom_point(data = plot_data_defects, aes(x = Defect + eps, y = Mean_defect), size = 1.5, color = "blue") +
    geom_errorbar(data = plot_data_defects, aes(x = Defect + eps, ymin = Mean_defect - Error_defect, ymax = Mean_defect + Error_defect), width = 0.03, color = "blue") +

    # Full-span shaded bands for weighted means
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = swap_mean - swap_err, ymax = swap_mean + swap_err,
      fill = "red", alpha = 0.1
    ) +
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = defect_mean - defect_err, ymax = defect_mean + defect_err,
      fill = "blue", alpha = 0.1
    ) +

    # Horizontal lines for the weighted means
    geom_hline(yintercept = swap_mean, linetype = "dashed", color = "red") +
    geom_hline(yintercept = defect_mean, linetype = "dashed", color = "blue") +
    labs(
      title = "Acceptance Rates: Swap (red) vs Per-Defect (blue)",
      x = "Replica / Defect",
      y = "Acceptance Rate"
    ) +
    # Set x ticks at defects_full positions
    scale_x_continuous(breaks = rep_sorted, labels = rep_sorted) +
    theme_minimal() +
    theme(
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black")
    )

  # Swap-only plot
  p_swap <- ggplot(plot_data, aes(x = Replica, y = Mean_swap)) +
    geom_point(size = 1.5, color = "red") +
    geom_errorbar(aes(ymin = Mean_swap - Error_swap, ymax = Mean_swap + Error_swap), width = 0.03, color = "red") +
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = swap_mean - swap_err, ymax = swap_mean + swap_err,
      fill = "red", alpha = 0.1
    ) +
    geom_hline(yintercept = swap_mean, linetype = "dashed", color = "red") +
    labs(
      title = "Swap Acceptance Rate per Replica",
      x = "Replica / Defect",
      y = "Acceptance Rate"
    ) +
    scale_x_continuous(breaks = rep_sorted, labels = rep_sorted) +
    theme_minimal() +
    theme(
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black")
    )

  # Per-defect-only plot
  p_defect <- ggplot(plot_data_defects, aes(x = Defect, y = Mean_defect)) +
    geom_point(size = 1.5, color = "blue") +
    geom_errorbar(aes(ymin = Mean_defect - Error_defect, ymax = Mean_defect + Error_defect), width = 0.03, color = "blue") +
    annotate("rect",
      xmin = x_min, xmax = x_max,
      ymin = defect_mean - defect_err, ymax = defect_mean + defect_err,
      fill = "blue", alpha = 0.1
    ) +
    geom_hline(yintercept = defect_mean, linetype = "dashed", color = "blue") +
    labs(
      title = "Per-Defect Acceptance Rate",
      x = "Defect",
      y = "Acceptance Rate"
    ) +
    scale_x_continuous(breaks = rep_sorted, labels = rep_sorted) +
    theme_minimal() +
    theme(
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black")
    )

  write_log(sprintf(
    "Weighted Mean Swap Acceptance Rate: %.3f ± %.3f",
    swap_mean, swap_err
  ))
  write_log(sprintf(
    "Weighted Mean Per-Defect Acceptance Rate: %.3f ± %.3f",
    defect_mean, defect_err
  ))

  ggsave(file.path(directory, "acceptance_by_replica.pdf"), plot = p_combined, width = 8, height = 6)
  write_log("Combined acceptance plot saved to acceptance_by_replica.pdf")

  ggsave(file.path(directory, "acceptance_swap_only.pdf"), plot = p_swap, width = 8, height = 6)
  write_log("Swap acceptance plot saved to acceptance_swap_only.pdf")

  ggsave(file.path(directory, "acceptance_per_defect_only.pdf"), plot = p_defect, width = 8, height = 6)
  write_log("Per-defect acceptance plot saved to acceptance_per_defect_only.pdf")

  return(plot_data)
}

weighted_mean_errors <- function(values, errors) {
  weights <- 1 / (errors^2)
  weighted_error <- sqrt(1 / sum(weights))
  weighted_mean <- sum(values * weights) * weighted_error^2
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
