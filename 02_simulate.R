# =============================================================================
# 02_simulate.R
# eCDF-Copula resampling framework for multivariate longitudinal data
#
# Implements the longitudinal eCDF-Copula simulation described in the paper
# (Algorithm 1 / sim6). This is the CORRECT version for longitudinal data.
#
# Key difference from the May09 / cross-sectional approach:
#   Cross-sectional (WRONG): data is first flattened to a wide matrix; a single
#     exchangeable scalar drives the copula; ICC decomposition then patches in
#     temporal structure. This treats all (feature × time) dimensions as
#     exchangeable and discards the true temporal correlation structure.
#
#   Longitudinal (THIS FILE): within-visit Spearman correlations (p×p at each
#     time point) and between-visit cross-correlations (p×p for each pair of
#     time points) are computed directly from the data and assembled into a
#     full (t·p)×(t·p) block correlation matrix. The Gaussian copula uses this
#     block matrix, and the inverse-CDF step operates per time point using
#     that time point's empirical marginal distribution. This preserves both
#     marginal distributions and the temporal dependence structure.
#
# Workflow:
#   1. Impute data using 01_impute.R  →  named list of completed datasets
#   2. Call run_all_simulations()    →  named list of simulated datasets
#      - For each imputation method (CC, MICE-L, …): split by group,
#        simulate n_per_group subjects per group, reassemble to long format
#      - For eCDF-Copula: use the original (possibly incomplete) data directly;
#        the simulation handles NAs within each time point's empirical CDF
#
# Main functions (public):
#   simulate_longitudinal_ecdf_copula()   core simulation for one group
#   run_all_simulations()                 wrapper over all methods
# =============================================================================

suppressPackageStartupMessages({
  library(copula)
  library(Matrix)
  library(tidyverse)
})

# 00_utils.R must be sourced before this file (run_all.R handles this)


# =============================================================================
# Step 1: Within-visit correlation matrices
# =============================================================================

# Compute a p×p Spearman correlation matrix for each time point.
# If fewer than min_n subjects have data at a time point, fall back to a weak
# positive-correlation default (avoids singular copula).
#
#   dat          long-format data for one group
#   feature_vars character vector of feature column names
#   time_var     time/visit column name
#   min_n        minimum subjects required to estimate correlation (default 5)
.get_within_visit_cor <- function(dat, feature_vars, time_var, min_n = 5) {
  t_vals <- sort(unique(dat[[time_var]]))
  p      <- length(feature_vars)

  setNames(lapply(t_vals, function(tv) {
    sub <- dat[dat[[time_var]] == tv, feature_vars, drop = FALSE]
    sub <- sub[complete.cases(sub), , drop = FALSE]

    if (nrow(sub) < min_n) {
      message(sprintf("  within-visit [%s]: only %d complete rows, using default cor",
                      tv, nrow(sub)))
      m <- matrix(0.1, p, p); diag(m) <- 1
    } else {
      m <- cor(sub, method = "spearman", use = "pairwise.complete.obs")
      m[is.na(m)] <- 0
      diag(m) <- 1
    }
    m
  }), as.character(t_vals))
}


# =============================================================================
# Step 2: Between-visit correlation matrices
# =============================================================================

# Compute a p×p cross-correlation matrix (Cor(features at t1, features at t2))
# for every ordered pair of time points, using only subjects observed at both.
#
#   subject_var  subject ID column name
.get_between_visit_cor <- function(dat, feature_vars, time_var, subject_var,
                                    min_n = 5) {
  t_vals <- sort(unique(dat[[time_var]]))
  p      <- length(feature_vars)
  pairs  <- combn(seq_along(t_vals), 2, simplify = FALSE)

  cors <- list()
  for (pair in pairs) {
    i <- pair[1]; j <- pair[2]
    tv1 <- t_vals[i]; tv2 <- t_vals[j]
    key <- paste0(as.character(tv1), "__", as.character(tv2))

    sub1 <- dat[dat[[time_var]] == tv1, ]
    sub2 <- dat[dat[[time_var]] == tv2, ]
    common <- intersect(sub1[[subject_var]], sub2[[subject_var]])

    if (length(common) < min_n) {
      message(sprintf("  between-visit [%s — %s]: only %d common subjects, using default",
                      tv1, tv2, length(common)))
      cors[[key]] <- matrix(0.2, p, p)
    } else {
      d1 <- sub1[sub1[[subject_var]] %in% common, feature_vars, drop = FALSE]
      d2 <- sub2[sub2[[subject_var]] %in% common, feature_vars, drop = FALSE]
      m  <- cor(d1, d2, method = "spearman", use = "pairwise.complete.obs")
      m[is.na(m)] <- 0.1
      cors[[key]] <- m
    }
  }
  cors
}


# =============================================================================
# Step 3: Assemble (t·p) × (t·p) block correlation matrix
# =============================================================================

# Fills diagonal blocks from within_cors and off-diagonal blocks from
# between_cors, then ensures positive definiteness via nearPD.
#
# Block layout:  rows/cols are ordered (t1·p features, t2·p features, …)
# Key format for between_cors: "{tv1}__{tv2}" (double underscore)
.assemble_block_cor <- function(within_cors, between_cors, t_vals, p) {
  t         <- length(t_vals)
  total_dim <- t * p
  big       <- matrix(0, total_dim, total_dim)

  # Diagonal blocks
  for (i in seq_along(t_vals)) {
    key <- as.character(t_vals[i])
    rs  <- (i - 1) * p + 1;  re <- i * p
    if (key %in% names(within_cors)) {
      big[rs:re, rs:re] <- within_cors[[key]]
    } else {
      diag(big[rs:re, rs:re]) <- 1
    }
  }

  # Off-diagonal blocks (symmetric)
  for (i in 1:(t - 1)) {
    for (j in (i + 1):t) {
      key <- paste0(as.character(t_vals[i]), "__", as.character(t_vals[j]))
      rs  <- (i - 1) * p + 1;  re <- i * p
      cs  <- (j - 1) * p + 1;  ce <- j * p
      m   <- if (key %in% names(between_cors)) between_cors[[key]] else matrix(0.1, p, p)
      big[rs:re, cs:ce] <- m
      big[cs:ce, rs:re] <- t(m)
    }
  }

  # Ensure positive definiteness (nearPD is more robust than diagonal shift)
  ev <- min(eigen(big, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 0) {
    message(sprintf("  Block cor matrix not PD (min eigenvalue = %.4f), applying nearPD", ev))
    big <- as.matrix(Matrix::nearPD(big, corr = TRUE, keepDiag = TRUE)$mat)
  }

  big
}


# =============================================================================
# Step 4: Inverse-CDF transformation (longitudinal, per time point)
# =============================================================================

# For each (time point, feature) combination, inverts the copula-sampled
# uniform values through the empirical CDF of the observed data at that
# specific time point.  NAs in the original data are automatically excluded
# when building the empirical CDF, so this works on incomplete datasets.
#
#   Fx_matrix    n_sim × (t*p) matrix of copula-sampled uniform values
#                Column order: (all features at t1, all features at t2, …)
#   dat          long-format data (one group, possibly with NAs)
#   feature_vars feature column names
#   time_var     time column name
#   t_vals       ordered time levels used for the simulation
.invcdf_longitudinal <- function(Fx_matrix, dat, feature_vars, time_var, t_vals) {
  n_sim   <- nrow(Fx_matrix)
  p       <- length(feature_vars)
  Sim     <- matrix(0, n_sim, length(t_vals) * p)
  col_idx <- 1L

  for (tv in t_vals) {
    time_data <- dat[dat[[time_var]] == tv, , drop = FALSE]

    for (fv in feature_vars) {
      x_all <- time_data[[fv]]
      x_obs <- x_all[!is.na(x_all)]

      Fx_col <- Fx_matrix[, col_idx]

      if (length(x_obs) == 0) {
        # No data at this time point for this feature: sample from U(0,1)
        Sim[, col_idx] <- runif(n_sim)
      } else {
        sc_ranks <- rank(x_obs, ties.method = "first") / length(x_obs)

        for (si in seq_len(n_sim)) {
          if (is.na(Fx_col[si])) {
            Sim[si, col_idx] <- sample(x_obs, 1)
          } else {
            valid <- which(Fx_col[si] <= sc_ranks)
            if (!length(valid)) {
              Sim[si, col_idx] <- max(x_obs)
            } else {
              Sim[si, col_idx] <- x_obs[valid[which.min(abs(sc_ranks[valid] - Fx_col[si]))]]
            }
          }
        }
      }
      col_idx <- col_idx + 1L
    }
  }
  Sim
}


# =============================================================================
# Main simulation function (one group)
# =============================================================================

#' Simulate n_sim longitudinal subjects from one group using eCDF-Copula
#'
#' @param dat_group   Long-format data for a SINGLE group. May contain NAs
#'                    (they are excluded per time point in the eCDF step).
#' @param feature_vars Character vector: variables to simulate
#' @param time_var    Column name for time/visit
#' @param subject_var Column name for subject ID (used for between-visit matching)
#' @param time_levels Ordered character vector of time levels. If NULL, inferred
#'                    from data (but using config$time_levels is preferred).
#' @param n_sim       Number of subjects to generate
#' @param seed        Random seed
#'
#' @return List with:
#'   $sim_data   n_sim × (t*p) numeric matrix
#'   $col_names  column names in format {feature}_{time_level}
#'   $n_sim
#'   $cor_matrix the (t·p)×(t·p) block correlation matrix (for diagnostics)
simulate_longitudinal_ecdf_copula <- function(dat_group, feature_vars,
                                               time_var, subject_var,
                                               time_levels = NULL,
                                               n_sim = 200, seed = 123) {
  set.seed(seed)
  dat_group <- as.data.frame(dat_group)

  if (is.null(time_levels))
    time_levels <- sort(unique(as.character(dat_group[[time_var]])))

  t_vals    <- time_levels
  t         <- length(t_vals)
  p         <- length(feature_vars)
  total_dim <- t * p

  cat(sprintf("  Time points: %d | Features: %d | Dimensions: %d\n",
              t, p, total_dim))

  # Step 1: correlation matrices
  cat("  Step 1: computing correlations\n")
  within_cors  <- .get_within_visit_cor(dat_group, feature_vars, time_var)
  between_cors <- .get_between_visit_cor(dat_group, feature_vars, time_var, subject_var)

  # Step 2: assemble block matrix and extract parameter vector
  cat("  Step 2: assembling block correlation matrix\n")
  big_cor   <- .assemble_block_cor(within_cors, between_cors, t_vals, p)
  param_vec <- big_cor[lower.tri(big_cor)]

  # Step 3: sample correlated uniforms from the Gaussian copula
  cat("  Step 3: sampling from copula\n")
  Mv_unif <- tryCatch({
    mycop  <- normalCopula(param = param_vec, dim = total_dim, dispstr = "un")
    mydist <- mvdc(mycop,
                   rep("unif", total_dim),
                   rep(list(list(min = 0, max = 1)), total_dim))
    rMvdc(n_sim, mydist)
  }, error = function(e) {
    warning("Copula sampling failed (", conditionMessage(e),
            "); falling back to independent uniforms.")
    matrix(runif(n_sim * total_dim), n_sim, total_dim)
  })

  # Step 4: inverse-CDF transformation per time point
  cat("  Step 4: inverse-CDF transformation\n")
  sim_data <- .invcdf_longitudinal(Mv_unif, dat_group, feature_vars,
                                    time_var, t_vals)

  # Column names: {feature}_{time_level}  (compatible with wide_to_long())
  col_names <- unlist(lapply(t_vals, function(tv)
    paste0(feature_vars, "_", as.character(tv))))
  colnames(sim_data) <- col_names

  cat("  Done.\n")
  list(sim_data   = sim_data,
       col_names  = col_names,
       n_sim      = n_sim,
       cor_matrix = big_cor)
}


# =============================================================================
# Wrapper: simulate for all imputation methods
# =============================================================================

#' Run eCDF-Copula simulation for all imputed datasets (and optionally raw data)
#'
#' For each dataset in imputed_list, the function:
#'   1. Splits subjects by group (config$group_var)
#'   2. Calls simulate_longitudinal_ecdf_copula() per group
#'   3. Converts each group's wide matrix to long format (via wide_to_long())
#'   4. Binds groups and returns one long-format data frame
#'
#' The eCDF-Copula method (sim6) uses the raw, possibly incomplete data directly;
#' this is passed via raw_data. If raw_data = NULL, sim6 is skipped.
#'
#' @param imputed_list Named list from run_all_imputations(); keys are method
#'                     abbreviations ("CC", "MICE-L", etc.)
#' @param raw_data     Original long-format data (with NAs) for the eCDF-Copula
#'                     method. Pass NULL to skip.
#' @param config       Config list from make_config() in 01_impute.R
#' @param n_per_group  Subjects to simulate per group (default 200)
#' @param seed         Base random seed. Group g gets seed + (g-1)*1000.
#' @param save         If TRUE, save each result as sim_{method}.csv in
#'                     config$dir_output
#'
#' @return Named list of long-format data frames, keyed by method name.
#'         For imputation-based methods the keys match imputed_list.
#'         The eCDF-Copula method is keyed "eCDF-Copula".
run_all_simulations <- function(imputed_list, raw_data = NULL, config,
                                n_per_group = 200, seed = 123, save = FALSE) {

  if (save)
    dir.create(config$dir_output, recursive = TRUE, showWarnings = FALSE)

  # Combine: imputed datasets + optional raw data for eCDF-Copula
  data_list <- imputed_list
  if (!is.null(raw_data))
    data_list[["eCDF-Copula"]] <- raw_data

  results <- list()

  for (method_name in names(data_list)) {
    dat    <- as.data.frame(data_list[[method_name]])
    groups <- unique(dat[[config$group_var]])

    cat(sprintf("\n=== Simulating: %s (%d groups, %d per group) ===\n",
                method_name, length(groups), n_per_group))

    group_dfs <- vector("list", length(groups))

    for (g_idx in seq_along(groups)) {
      g_label  <- as.character(groups[g_idx])
      dat_g    <- dat[as.character(dat[[config$group_var]]) == g_label, ]
      start_id <- (g_idx - 1L) * n_per_group + 1L
      g_seed   <- seed + (g_idx - 1L) * 1000L

      cat(sprintf("  Group: %s (n_obs = %d, seed = %d)\n",
                  g_label, nrow(dat_g), g_seed))

      sim_result <- simulate_longitudinal_ecdf_copula(
        dat_group    = dat_g,
        feature_vars = config$outcome_vars,
        time_var     = config$time_var,
        subject_var  = config$subject_var,
        time_levels  = config$time_levels,
        n_sim        = n_per_group,
        seed         = g_seed
      )

      group_dfs[[g_idx]] <- wide_to_long(
        sim_result  = sim_result,
        group_label = g_label,
        start_id    = start_id,
        subject_var = config$subject_var,
        time_var    = config$time_var,
        group_var   = config$group_var,
        time_levels = config$time_levels
      )
    }

    result_df          <- do.call(rbind, group_dfs)
    results[[method_name]] <- result_df

    if (save) {
      fname <- file.path(
        config$dir_output,
        paste0("sim_", gsub("[^A-Za-z0-9]", "_", method_name), ".csv")
      )
      write_csv(result_df, fname)
      cat(sprintf("  Saved: %s\n", fname))
    }
  }

  results
}
