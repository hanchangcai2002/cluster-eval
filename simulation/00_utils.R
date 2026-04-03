# =============================================================================
# 00_utils.R
# Core utility functions shared across the pipeline
#
# Contents:
#   1. Data reshaping utilities
#   2. Metric computation helpers
#
# Removed (no longer used):
#   - trim01, draw_correlated_U, InverseCDF.Sim, decompose_var,
#     gen_correlated_e  (May09 math helpers — only used by the deprecated
#     cross-sectional simulate_ecdf_copula)
#   - long_to_wide_matrix  (unused)
#   - estimate_icc         (unused)
#   - simulate_ecdf_copula (cross-sectional / May09 version — superseded by
#     simulate_longitudinal_ecdf_copula in 02_simulate.R)
# =============================================================================

suppressPackageStartupMessages({
  library(mclust)       # adjustedRandIndex
  library(clusterCrit)  # intCriteria (CHI)
  library(tidyverse)
})


# =============================================================================
# 1. Data reshaping utilities
# =============================================================================

# Convert long-format data to a subject x (feature * time) flat matrix.
# Used primarily for computing CHI.
#
#   dat          long-format data frame
#   subject_var  column name for subject ID
#   time_var     column name for time/visit
#   feature_vars character vector of feature column names
#   time_levels  ordered time levels (inferred from data if NULL)
#
# Returns a list:
#   $matrix      numeric matrix, n_subjects x (n_features * n_times)
#   $subjects    subject IDs in row order
make_flat_matrix <- function(dat, subject_var, time_var, feature_vars,
                              time_levels = NULL) {
  dat <- as.data.frame(dat)
  if (is.null(time_levels))
    time_levels <- sort(unique(dat[[time_var]]))

  subjects   <- unique(dat[[subject_var]])
  n_subjects <- length(subjects)

  dat <- dat[order(match(dat[[subject_var]], subjects),
                   match(dat[[time_var]], time_levels)), ]

  mat <- matrix(NA_real_, nrow = n_subjects,
                ncol = length(time_levels) * length(feature_vars))

  for (i in seq_along(subjects)) {
    subj_rows <- dat[[subject_var]] == subjects[i]
    subj_dat  <- dat[subj_rows, ]
    idx <- match(time_levels, subj_dat[[time_var]])
    mat[i, ] <- unlist(lapply(feature_vars, function(v) subj_dat[[v]][idx]))
  }
  list(matrix = mat, subjects = subjects)
}


# Convert wide simulation output back to long format.
#
#   sim_result   list from simulate_longitudinal_ecdf_copula():
#                $sim_data, $col_names, $n_sim
#   group_label  group name string (e.g., "Control")
#   start_id     integer: starting subject ID
#   subject_var, time_var, group_var  column names for output data frame
#   time_levels  ordered character vector of time levels
wide_to_long <- function(sim_result, group_label, start_id = 1,
                          subject_var = "Subject", time_var = "Time",
                          group_var = "Group", time_levels) {
  sim_data <- sim_result$sim_data
  colnames(sim_data) <- sim_result$col_names
  n_sim <- sim_result$n_sim

  sim_df <- as.data.frame(sim_data)
  sim_df[[subject_var]] <- seq(start_id, start_id + n_sim - 1L)

  # Build a regex that matches any time level (escape special characters)
  esc <- function(x) gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", x)
  time_pattern <- paste(esc(time_levels), collapse = "|")
  pattern <- paste0("^(.*?)_(", time_pattern, ")$")

  long <- sim_df %>%
    pivot_longer(
      cols          = -all_of(subject_var),
      names_to      = c("variable", time_var),
      names_pattern = pattern,
      values_to     = "value"
    ) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    mutate(
      !!group_var := factor(group_label),
      !!time_var  := factor(.data[[time_var]], levels = time_levels)
    ) %>%
    arrange(.data[[subject_var]], .data[[time_var]])

  long
}


# =============================================================================
# 2. Metric computation helpers
# =============================================================================

# Compute ARI (Adjusted Rand Index) with error handling
compute_ari <- function(pred, true) {
  tryCatch(
    adjustedRandIndex(pred, true),
    error = function(e) { warning("ARI failed: ", e$message); NA_real_ }
  )
}

# Compute Calinski-Harabasz Index
compute_chi <- function(flat_mat, pred_labels) {
  tryCatch(
    intCriteria(traj  = flat_mat,
                part  = as.integer(pred_labels),
                crit  = "Calinski_Harabasz")$calinski_harabasz,
    error = function(e) { warning("CHI failed: ", e$message); NA_real_ }
  )
}
