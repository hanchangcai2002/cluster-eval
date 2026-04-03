# =============================================================================
# 03_prepare.R
# Data preparation: standardization and dimensionality reduction
#
# Two preparation strategies (matching the paper's SCE2 and SCE3):
#
#   SCE2 — use summary/composite feature scores directly (standardized)
#   SCE3 — apply PCA across all item-level features × time points, then use
#           the leading PC scores as clustering inputs
#
# Both functions add a numeric time column ("time_num") required by GBMT,
# flexmix, and multlcmm. The mapping from visit labels to numbers is
# controlled via make_time_map().
#
# Output format for clustering (04_cluster.R):
#   Long-format data frame with columns:
#     config$subject_var   subject ID
#     config$time_var      original time labels (factor)
#     time_num             numeric time (for model-based methods)
#     config$group_var     true group label (for ARI evaluation only)
#     <feature columns>    standardized clustering features
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))


# =============================================================================
# Time mapping helper
# =============================================================================

#' Build a named numeric vector mapping time labels to numeric values
#'
#' @param time_levels  Ordered character vector of time labels (from config)
#' @param time_values  Numeric values to assign. Defaults to 0-indexed integers.
#'                     Pass a named or positional vector to use actual time
#'                     units (e.g., months).
#'
#' @examples
#' # Default: 0, 1, 2, ..., 13
#' make_time_map(config$time_levels)
#'
#' # Actual month values for the original dataset
#' make_time_map(config$time_levels,
#'               c(-1, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24))
make_time_map <- function(time_levels, time_values = NULL) {
  if (is.null(time_values))
    time_values <- seq(0, length(time_levels) - 1L)
  stopifnot(length(time_values) == length(time_levels))
  setNames(as.numeric(time_values), time_levels)
}


# =============================================================================
# SCE2 preparation — standardize composite/summary features
# =============================================================================

#' Prepare SCE2: select and standardize summary-level features
#'
#' Selects the specified feature columns from the simulated long-format data,
#' standardizes them (z-score across all observations), and adds a numeric
#' time column for use with model-based clustering methods.
#'
#' @param dat          Long-format simulated data frame
#' @param config       Config list from make_config()
#' @param feature_vars Character vector: which columns to use as clustering
#'                     features (a subset of config$outcome_vars, e.g. the
#'                     composite scores rather than item-level scores)
#' @param time_map     Named numeric vector from make_time_map(). If NULL,
#'                     defaults to 0-indexed integers.
#'
#' @return Long-format data frame with columns:
#'   subject_var, time_var (factor), time_num (numeric), group_var,
#'   and one standardized column per feature_var
prepare_sce2 <- function(dat, config, feature_vars, time_map = NULL) {
  dat <- as.data.frame(dat)
  sv  <- config$subject_var
  tv  <- config$time_var
  gv  <- config$group_var

  if (is.null(time_map))
    time_map <- make_time_map(config$time_levels)

  stopifnot(all(feature_vars %in% names(dat)))
  stopifnot(all(as.character(unique(dat[[tv]])) %in% names(time_map)))

  dat %>%
    select(all_of(c(sv, tv, gv, feature_vars))) %>%
    mutate(
      across(all_of(feature_vars), ~ as.vector(scale(.))),
      time_num = time_map[as.character(.data[[tv]])]
    ) %>%
    arrange(.data[[sv]], .data[[tv]])
}


# =============================================================================
# SCE3 preparation — PCA across item-level features × time points
# =============================================================================

#' Prepare SCE3: PCA dimensionality reduction on item-level longitudinal features
#'
#' Builds a subject × (feature × time) matrix, runs PCA with standardization,
#' and attaches the leading PC scores as time-invariant clustering features.
#' The PC scores are the same for every time point of a subject (they summarise
#' the full longitudinal profile).
#'
#' @param dat          Long-format simulated data frame
#' @param config       Config list from make_config()
#' @param feature_vars Character vector: item-level variables to use for PCA
#'                     (typically all item-level columns, e.g. ma_1:10 + hars_1:14_sev)
#' @param n_pc         Number of PCs to retain (default 10)
#' @param time_map     Named numeric vector from make_time_map(). If NULL,
#'                     defaults to 0-indexed integers.
#'
#' @return A list with:
#'   $data      Long-format data frame (same columns as prepare_sce2 but features
#'              are PC1 … PCn_pc instead of the original item columns)
#'   $pca       prcomp object (for diagnostics, elbow plots, etc.)
#'   $pc_names  Character vector of PC column names used
#'   $var_exp   Numeric vector: proportion of variance explained per PC
prepare_sce3 <- function(dat, config, feature_vars, n_pc = 10, time_map = NULL) {
  dat <- as.data.frame(dat)
  sv  <- config$subject_var
  tv  <- config$time_var
  gv  <- config$group_var

  if (is.null(time_map))
    time_map <- make_time_map(config$time_levels)

  stopifnot(all(feature_vars %in% names(dat)))

  time_levels <- config$time_levels

  # Step 1: standardize item-level features globally
  dat[, feature_vars] <- scale(dat[, feature_vars])

  # Step 2: build n_subjects × (n_features * n_times) matrix for PCA
  #         column order: (feat1_t1, feat2_t1, …, featP_t1, feat1_t2, …)
  subjects   <- unique(dat[[sv]])
  n_subjects <- length(subjects)

  pca_mat <- matrix(NA_real_, nrow = n_subjects,
                    ncol = length(feature_vars) * length(time_levels))

  for (i in seq_along(subjects)) {
    subj_rows <- dat[[sv]] == subjects[i]
    subj_dat  <- dat[subj_rows, ]
    # Align rows to time_levels order
    idx       <- match(time_levels, as.character(subj_dat[[tv]]))
    subj_dat  <- subj_dat[idx, ]
    pca_mat[i, ] <- unlist(lapply(feature_vars, function(v) subj_dat[[v]]))
  }

  # Step 3: PCA (center + scale the flattened matrix columns)
  pca_result  <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
  n_pc_actual <- min(n_pc, ncol(pca_result$x))
  pc_names    <- paste0("PC", seq_len(n_pc_actual))
  var_exp     <- (pca_result$sdev^2) / sum(pca_result$sdev^2)

  # Step 4: attach PC scores to the long-format data
  #         PCs are time-invariant — the same value is repeated for every
  #         time point of a subject (the full longitudinal profile is already
  #         encoded in the PC scores)
  pc_df        <- as.data.frame(pca_result$x[, seq_len(n_pc_actual), drop = FALSE])
  colnames(pc_df) <- pc_names
  pc_df[[sv]]  <- subjects

  out_data <- dat %>%
    select(all_of(c(sv, tv, gv))) %>%
    left_join(pc_df, by = sv) %>%
    mutate(time_num = time_map[as.character(.data[[tv]])]) %>%
    arrange(.data[[sv]], .data[[tv]])

  list(
    data     = out_data,
    pca      = pca_result,
    pc_names = pc_names,
    var_exp  = var_exp,
    n_pc     = n_pc_actual
  )
}
