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
  swap_attempts_per_replica <- vector("numeric", length = length(data$defects[[1]]))
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

      swap_attempts_per_replica[replica_index] <- swap_attempts_per_replica[replica_index] + 1
      total_accepts_per_replica[replica_index] <- total_accepts_per_replica[replica_index] + accepts[j]
      swap_attempts_per_replica[partner_index] <- swap_attempts_per_replica[partner_index] + 1
      total_accepts_per_replica[partner_index] <- total_accepts_per_replica[partner_index] + accepts[j]
    }
  }
  acceptance_rates <- total_accepts_per_replica / swap_attempts_per_replica

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
    Replica = 0:(length(data$defects[[1]]) - 1),
    Mean = acceptance_rates
    # Error = bootstrap_results["error", ]
  )

  # Generate and save the plot
  p <- ggplot(plot_data, aes(x = Replica, y = Mean)) +
    geom_point(size = 3) +
    # geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), width = 0.2) +
    labs(
      title = "Acceptance Rate per Relplica",
      x = "Replica",
      y = "Acceptance Rate"
    ) +
    theme_minimal()

  ggsave(file.path(directory, "acceptance_by_rank.pdf"), plot = p, width = 8, height = 6)
  write_log("Acceptance analysis plot saved to acceptance_by_rank.pdf")

  return(plot_data)
}

# Function to process one row for defect swaps, designed for parallel execution
get_swap_details_for_row <- function(i, data) {
  # For a given row 'i', get the state of defects from the previous row 'i-1'
  current_defects <- data$defects[[i - 1]]
  num_defects <- length(current_defects)

  row <- data[i, ]
  swap_start <- row$swap_start
  accepts <- row$accepts[[1]]
  delta_H_swaps <- row$delta_H_swap[[1]]

  row_swap_details <- list()

  for (j in seq_along(accepts)) {
    idx1 <- (swap_start + j - 1) %% num_defects
    idx2 <- (idx1 + 1) %% num_defects

    r_idx1 <- idx1 + 1
    r_idx2 <- idx2 + 1

    defect1 <- current_defects[r_idx1]
    defect2 <- current_defects[r_idx2]
    accepted <- accepts[j]
    delta_H <- delta_H_swaps[j]

    # Note: This logic assumes the 'defects' in row 'i-1' are the starting point
    # for the swaps in row 'i'. The swaps are calculated but not chained within this function.
    row_swap_details[[length(row_swap_details) + 1]] <- list(
      step = row$step,
      swap_attempt = j,
      defect1 = defect1,
      defect2 = defect2,
      accepted = accepted,
      delta_H_swap = delta_H
    )

    if (accepted == 1) {
      # Swap for the next attempt in the same step
      current_defects[r_idx1] <- defect2
      current_defects[r_idx2] <- defect1
    }
  }

  return(row_swap_details)
}

# Function to track defect swaps across all steps
track_defect_swaps <- function(directory, data) {
  if (nrow(data) < 2) {
    write_log("Defect tracking skipped: not enough data rows.")
    return(NULL)
  }

  # Use lapply to process rows 2 to N. This can be swapped for a parallel version.
  list_of_swap_lists <- lapply(2:nrow(data), get_swap_details_for_row, data = data)

  # Combine the list of lists into a single data frame
  all_swap_details <- unlist(list_of_swap_lists, recursive = FALSE)

  swap_df <- dplyr::bind_rows(all_swap_details)
  #   output_path <- file.path(directory, "defect_swap_details.csv")
  #   write.csv(swap_df, output_path, row.names = FALSE)
  #   write_log(paste("Defect swap details saved to", output_path))

  return(swap_df)
}

# Function to analyze acceptance based on defect distance
analyze_acceptance_by_defect_distance <- function(directory, swap_df, num_bins = 10) {
  if (is.null(swap_df) || nrow(swap_df) == 0) {
    write_log("Acceptance by defect distance analysis skipped: no swap data.")
    return(NULL)
  }

  # Calculate defect distance and create bins for it
  binned_data <- swap_df %>%
    mutate(defect_dist = abs(defect1 - defect2)) %>%
    filter(!is.na(defect_dist)) %>%
    mutate(dist_bin = cut(defect_dist, breaks = num_bins))

  # Group by the new bins and calculate the mean acceptance rate and its error
  # Increased R from 100 to 400 for better error estimates (best practices)
  acceptance_by_bin <- binned_data %>%
    group_by(dist_bin) %>%
    summarise(
      acceptance_rate = mean(accepted),
      acceptance_error = hadron::bootstrap.meanerror(accepted, R = 400),
      .groups = "drop"
    )

  # Create and save the histogram using the binned data
  p <- ggplot(acceptance_by_bin, aes(x = dist_bin, y = acceptance_rate)) +
    geom_col(fill = "steelblue") +
    geom_errorbar(
      aes(ymin = acceptance_rate - acceptance_error, ymax = acceptance_rate + acceptance_error),
      width = 0.25
    ) +
    labs(
      title = "Swap Acceptance Rate by Binned Defect Distance",
      x = "Defect Distance Bin",
      y = "Swap Acceptance Rate"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_path <- file.path(directory, "acceptance_by_defect_distance.pdf")
  ggsave(output_path, plot = p, width = 10, height = 6)
  write_log(paste("Acceptance by defect distance plot saved to", output_path))

  return(acceptance_by_bin)
}

# Function to analyze delta_H distributions based on defect distance
analyze_delta_h_by_defect_distance <- function(directory, swap_df, num_bins = 30) {
  if (is.null(swap_df) || nrow(swap_df) == 0) {
    write_log("Delta_H by defect distance analysis skipped: no swap data.")
    return(NULL)
  }

  # Calculate defect distance
  plot_data <- swap_df %>%
    mutate(defect_dist = abs(defect1 - defect2)) %>%
    filter(!is.na(defect_dist))

  # Create bins with a fixed border at 0
  max_dist <- max(plot_data$defect_dist, na.rm = TRUE)
  breaks <- seq(0, max_dist, length.out = num_bins + 1)
  breaks <- unique(round(breaks, 5))

  plot_data <- plot_data %>%
    mutate(dist_bin = cut(defect_dist, breaks = breaks, include.lowest = TRUE))

  # Define histogram breaks to ensure 0 is an edge
  min_h <- min(plot_data$delta_H_swap, na.rm = TRUE)
  max_h <- max(plot_data$delta_H_swap, na.rm = TRUE)
  binwidth <- (max_h - min_h) / num_bins # bins, as before

  # Create breaks from 0 outwards
  breaks_neg <- seq(0, min_h - binwidth, by = -binwidth)
  breaks_pos <- seq(0, max_h + binwidth, by = binwidth)
  hist_breaks <- sort(unique(c(breaks_neg, breaks_pos)))

  # Create and save the faceted histogram
  p <- ggplot(plot_data, aes(x = delta_H_swap)) +
    geom_histogram(breaks = hist_breaks, fill = "blue", alpha = 0.7) +
    facet_wrap(~dist_bin, scales = "free_y") +
    labs(
      title = "Distribution of delta_H_swap by Binned Defect Distance",
      x = "delta_H_swap",
      y = "Count"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_path <- file.path(directory, "delta_h_by_defect_distance.pdf")
  ggsave(output_path, plot = p, width = 12, height = 8)
  write_log(paste("Delta_H by defect distance plot saved to", output_path))

  return(NULL)
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

  return()

  # Track defects
  defect_swap_results <- track_defect_swaps(directory, ptbc_data)
  if (!is.null(defect_swap_results)) {
    write_log("Defect tracking analysis complete.")
    # New analysis based on defect swaps
    analyze_acceptance_by_defect_distance(directory, defect_swap_results, num_bins = number_defects - 1)
    analyze_delta_h_by_defect_distance(directory, defect_swap_results, num_bins = 30)
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
