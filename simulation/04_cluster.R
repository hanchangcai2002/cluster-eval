# =============================================================================
# 04_cluster.R
# Five longitudinal clustering methods
#
# Methods (named as in the manuscript):
#   GBMT    Group-based multivariate trajectory model   cluster_gbmt()
#   KML3D   Multivariate longitudinal k-means           cluster_kml3d()
#   flexmix Finite mixture model (GBTM family)          cluster_flexmix()
#   LCMM    Latent class mixed-effects model            cluster_lcmm()
#   mixAK   Bayesian mixture model (MCMC)               cluster_mixak()
#
# Input data format (output of 03_prepare.R):
#   Long-format data frame with columns:
#     config$subject_var   subject ID
#     config$time_var      original time labels (factor, for KML3D ordering)
#     time_num             numeric time (required by GBMT, flexmix, LCMM, mixAK)
#     config$group_var     true group label (used ONLY for ARI, not for clustering)
#     <feature_vars>       standardized clustering features
#
# All functions return a named list keyed by "k{K}" (e.g. "k2", "k3", ...),
# each element being a list with:
#   $method   character string
#   $k        integer
#   $pred     named integer vector: cluster assignment per subject (1..k)
#   $metrics  data.frame(k, ARI, CHI, BIC, PED, Runtime_seconds)
#             BIC/PED are NA when unavailable for a method
#
# Convenience wrapper:
#   run_all_clustering() runs any subset of methods and returns a nested list
#   plus a flat metrics data.frame.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(gbmt)
  library(kml3d)
  library(flexmix)
  library(lcmm)
  library(mixAK)
  library(mclust)      # adjustedRandIndex
  library(clusterCrit) # intCriteria (CHI)
})

# 00_utils.R must be sourced before this file (run_all.R handles this)


# =============================================================================
# Internal helper: extract true labels and flat matrix
# =============================================================================

# Returns a named integer vector of true group labels (1, 2, …) per subject,
# in the same order as unique(dat[[sv]]).
.true_labels <- function(dat, sv, gv) {
  subj_order <- unique(dat[[sv]])
  grp        <- dat[[gv]][match(subj_order, dat[[sv]])]
  as.integer(as.factor(grp))
}

# Wrapper around make_flat_matrix() from 00_utils.R.
.flat_mat <- function(dat, config, feature_vars) {
  make_flat_matrix(dat,
                   subject_var  = config$subject_var,
                   time_var     = config$time_var,
                   feature_vars = feature_vars,
                   time_levels  = config$time_levels)$matrix
}


# =============================================================================
# GBMT — Group-Based Multivariate Trajectory Model
# =============================================================================

#' Fit GBMT for k in k_range and return metrics
#'
#' @param dat          Prepared long-format data (from prepare_sce2/sce3)
#' @param config       Config list
#' @param feature_vars Clustering feature column names
#' @param k_range      Integer vector of cluster numbers to try (default 2:5)
#' @param d            Polynomial degree for trajectory (default 2 = quadratic)
#' @param nstart       Number of random EM restarts (default 3)
#'
#' NOTE: scaling = 0 because features are already standardized in 03_prepare.R
cluster_gbmt <- function(dat, config, feature_vars, k_range = 2:5,
                          d = 2, nstart = 3) {
  dat  <- as.data.frame(dat)
  sv   <- config$subject_var
  gv   <- config$group_var
  true <- .true_labels(dat, sv, gv)
  flat <- .flat_mat(dat, config, feature_vars)

  # GBMT needs Subject as factor and numeric time
  dat[[sv]] <- as.factor(dat[[sv]])

  results <- list()
  for (k in k_range) {
    cat(sprintf("GBMT k=%d\n", k))
    t0  <- Sys.time()
    fit <- gbmt(
      x.names = feature_vars,
      unit    = sv,
      time    = "time_num",
      ng      = k,
      d       = d,
      data    = dat,
      scaling = 0,        # already standardized
      nstart  = nstart,
      maxit   = 1000,
      quiet   = TRUE
    )
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    pred <- fit$assign          # integer vector, one entry per subject
    ari  <- compute_ari(pred, true)
    chi  <- compute_chi(flat, pred)
    bic  <- unname(fit$ic["bic"])

    results[[paste0("k", k)]] <- list(
      method  = "GBMT",
      k       = k,
      pred    = setNames(as.integer(pred), levels(dat[[sv]])),
      metrics = data.frame(k = k, ARI = ari, CHI = chi,
                            BIC = bic, PED = NA_real_,
                            Runtime_seconds = runtime)
    )
  }
  results
}


# =============================================================================
# KML3D — Multivariate Longitudinal k-Means
# =============================================================================

#' Fit KML3D for k in k_range and return metrics
#'
#' @param nb_redrawing Number of random restarts (default 5)
#'
#' NOTE: KML3D does not produce BIC; CHI is computed from the flat matrix.
cluster_kml3d <- function(dat, config, feature_vars, k_range = 2:5,
                           nb_redrawing = 5) {
  dat  <- as.data.frame(dat)
  sv   <- config$subject_var
  tv   <- config$time_var
  gv   <- config$group_var
  true <- .true_labels(dat, sv, gv)

  # Build n × t × p array (required format for kml3d)
  subjects    <- unique(dat[[sv]])
  time_levels <- config$time_levels
  n_subjects  <- length(subjects)
  n_times     <- length(time_levels)
  n_vars      <- length(feature_vars)

  dat <- dat[order(match(dat[[sv]], subjects),
                   match(as.character(dat[[tv]]), time_levels)), ]

  data_array <- array(NA_real_, dim = c(n_subjects, n_times, n_vars))
  for (v in seq_along(feature_vars)) {
    mat <- matrix(dat[[feature_vars[v]]], nrow = n_times, ncol = n_subjects)
    data_array[, , v] <- t(mat)
  }

  flat  <- matrix(NA_real_, n_subjects, n_times * n_vars)
  for (i in seq_len(n_subjects))
    flat[i, ] <- as.vector(data_array[i, , ])

  cld3d <- clusterLongData3d(data_array)

  results <- list()
  for (k in k_range) {
    cat(sprintf("KML3D k=%d\n", k))
    t0 <- Sys.time()
    kml3d(cld3d, nbClusters = k, nbRedrawing = nb_redrawing, toPlot = "none")
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    pred_char <- getClusters(cld3d, k)   # character labels A, B, C, …
    pred      <- as.integer(as.factor(pred_char))
    ari       <- compute_ari(pred, true)
    chi       <- compute_chi(flat, pred)

    results[[paste0("k", k)]] <- list(
      method  = "KML3D",
      k       = k,
      pred    = setNames(pred, subjects),
      metrics = data.frame(k = k, ARI = ari, CHI = chi,
                            BIC = NA_real_, PED = NA_real_,
                            Runtime_seconds = runtime)
    )
  }
  results
}


# =============================================================================
# flexmix — Finite Mixture Model (GBTM family)
# =============================================================================

#' Fit flexmix for k in k_range and return metrics
#'
#' Each feature is modelled with FLXMRglm(feature ~ time_num).
#' Subject-level cluster assignment is determined by majority vote across
#' all observations for that subject.
cluster_flexmix <- function(dat, config, feature_vars, k_range = 2:5) {
  dat  <- as.data.frame(dat)
  sv   <- config$subject_var
  gv   <- config$group_var
  true <- .true_labels(dat, sv, gv)
  flat <- .flat_mat(dat, config, feature_vars)

  # One FLXMRglm component per feature (linear in numeric time)
  model_list <- lapply(feature_vars, function(v)
    FLXMRglm(as.formula(paste(v, "~ time_num"))))

  subjects <- unique(dat[[sv]])

  results <- list()
  for (k in k_range) {
    cat(sprintf("flexmix k=%d\n", k))
    t0  <- Sys.time()
    fit <- flexmix(as.formula(paste0("~ 1 | ", sv)),
                   data = dat, k = k, model = model_list)
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # Majority-vote cluster per subject across all their observations
    cluster_obs  <- clusters(fit)
    cluster_subj <- tapply(cluster_obs, dat[[sv]], function(x)
      as.integer(names(which.max(table(x)))))
    pred <- as.integer(cluster_subj[as.character(subjects)])

    ari  <- compute_ari(pred, true)
    chi  <- compute_chi(flat, pred)
    bic  <- BIC(fit)

    results[[paste0("k", k)]] <- list(
      method  = "flexmix",
      k       = k,
      pred    = setNames(pred, subjects),
      metrics = data.frame(k = k, ARI = ari, CHI = chi,
                            BIC = bic, PED = NA_real_,
                            Runtime_seconds = runtime)
    )
  }
  results
}


# =============================================================================
# LCMM — Latent Class Mixed-Effects Model (multivariate, multlcmm)
# =============================================================================

#' Fit multlcmm for k in k_range and return metrics
#'
#' A single ng=1 model is fit first and used to initialize all ng>1 models
#' (B = m1). Only a random intercept per feature is included to avoid
#' over-parameterization; class-specific intercepts are estimated via
#' mixture = ~1.
#'
#' @param maxiter  Maximum EM iterations (default 200)
#'
#' NOTE: CHI is computed here (was NA in the original code; now fixed).
cluster_lcmm <- function(dat, config, feature_vars, k_range = 2:5,
                          maxiter = 200) {
  dat  <- as.data.frame(dat)
  sv   <- config$subject_var
  gv   <- config$group_var
  true <- .true_labels(dat, sv, gv)
  flat <- .flat_mat(dat, config, feature_vars)

  formula_fixed <- as.formula(
    paste(paste(feature_vars, collapse = " + "), "~ time_num"))

  cat("LCMM: fitting initial ng=1 model\n")
  m1 <- multlcmm(
    fixed   = formula_fixed,
    random  = ~ 1,
    subject = sv,
    data    = dat,
    link    = "linear",
    maxiter = maxiter,
    verbose = FALSE
  )

  subjects <- unique(dat[[sv]])
  results  <- list()

  for (k in k_range) {
    cat(sprintf("LCMM k=%d\n", k))
    t0  <- Sys.time()
    fit <- multlcmm(
      fixed   = formula_fixed,
      mixture = ~ 1,      # class-specific intercepts only
      random  = ~ 1,      # random intercept per subject-feature
      subject = sv,
      data    = dat,
      link    = "linear",
      ng      = k,
      B       = m1,       # initialize from ng=1 solution
      maxiter = maxiter,
      verbose = FALSE
    )
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # pprob rows correspond to subjects in the order they appear in fit$pprob
    pred_df <- fit$pprob
    pred    <- setNames(as.integer(pred_df$class),
                        as.character(pred_df[[sv]]))
    pred_ordered <- pred[as.character(subjects)]

    ari  <- compute_ari(pred_ordered, true)
    chi  <- compute_chi(flat, pred_ordered)   # fixed: was NA in original
    bic  <- fit$BIC

    results[[paste0("k", k)]] <- list(
      method  = "LCMM",
      k       = k,
      pred    = pred_ordered,
      metrics = data.frame(k = k, ARI = ari, CHI = chi,
                            BIC = bic, PED = NA_real_,
                            Runtime_seconds = runtime)
    )
  }
  results
}


# =============================================================================
# mixAK — Bayesian Mixture Model (MCMC via GLMM_MCMC)
# =============================================================================

#' Fit mixAK / GLMM_MCMC for k in k_range and return metrics
#'
#' Model structure (as in the paper):
#'   - No fixed effects (x = "empty")
#'   - Random slope for time (z = time_num column, random.intercept = FALSE)
#'   - Mixture distribution on the random effects → subjects clustered by
#'     their time trajectory shape
#'
#' Two MCMC chains are run automatically by GLMM_MCMC; posterior component
#' probabilities are averaged across chains for more stable assignments.
#'
#' @param nMCMC  Named integer vector c(burn, keep, thin, info)
cluster_mixak <- function(dat, config, feature_vars, k_range = 2:5,
                           nMCMC = c(burn = 400, keep = 800,
                                     thin = 3, info = 100)) {
  dat  <- as.data.frame(dat)
  sv   <- config$subject_var
  gv   <- config$group_var

  # Ensure consistent subject ordering
  subjects <- unique(dat[[sv]])
  dat      <- dat[order(match(dat[[sv]], subjects)), ]

  true     <- .true_labels(dat, sv, gv)
  flat     <- .flat_mat(dat, config, feature_vars)
  n_vars   <- length(feature_vars)

  # Response matrix (observations × features, continuous outcomes first)
  y  <- as.matrix(dat[, feature_vars])
  id <- as.integer(as.factor(dat[[sv]]))

  # No fixed effects; random slope for time (no random intercept)
  X_list <- rep(list("empty"), n_vars)
  Z_list <- replicate(n_vars,
                       matrix(dat[["time_num"]], ncol = 1),
                       simplify = FALSE)
  rand_int <- rep(FALSE, n_vars)   # random slope, no separate intercept

  results <- list()
  for (k in k_range) {
    cat(sprintf("mixAK k=%d\n", k))
    t0  <- Sys.time()
    fit <- GLMM_MCMC(
      y               = y,
      dist            = rep("gaussian", n_vars),
      id              = id,
      x               = X_list,
      z               = Z_list,
      random.intercept = rand_int,
      nMCMC           = nMCMC,
      prior.b         = list(Kmax = k),
      init.b          = list(K   = k),
      PED             = TRUE
    )
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # Average posterior component probabilities across the two MCMC chains
    # (chain 1 = fit[[1]], chain 2 = fit[[2]]; averaging gives more stable estimates)
    prob_mat     <- (fit[[1]]$poster.comp.prob_b +
                     fit[[2]]$poster.comp.prob_b) / 2
    pred         <- as.integer(apply(prob_mat, 1, which.max))
    names(pred)  <- as.character(subjects)

    ped <- tryCatch(unname(fit$PED["PED"]), error = function(e) NA_real_)
    ari <- compute_ari(pred, true)
    chi <- compute_chi(flat, pred)

    results[[paste0("k", k)]] <- list(
      method  = "mixAK",
      k       = k,
      pred    = pred,
      metrics = data.frame(k = k, ARI = ari, CHI = chi,
                            BIC = NA_real_, PED = ped,
                            Runtime_seconds = runtime)
    )
  }
  results
}


# =============================================================================
# Convenience wrapper
# =============================================================================

#' Run a subset of clustering methods and return combined results
#'
#' @param dat          Prepared long-format data (from 03_prepare.R)
#' @param config       Config list from make_config()
#' @param feature_vars Clustering feature column names
#' @param methods      Character vector, any subset of:
#'                     c("GBMT", "KML3D", "flexmix", "LCMM", "mixAK")
#' @param k_range      Integer vector of cluster numbers to evaluate (default 2:5)
#' @param ...          Additional arguments forwarded to individual cluster_*()
#'                     functions (e.g. d=1 for GBMT, nb_redrawing=10 for KML3D)
#'
#' @return A list with:
#'   $results   Nested list: results[[method]][[paste0("k",k)]] → see above
#'   $metrics   Flat data.frame with columns:
#'              method, k, ARI, CHI, BIC, PED, Runtime_seconds
run_all_clustering <- function(dat, config, feature_vars,
                                methods = c("GBMT", "KML3D", "flexmix",
                                            "LCMM", "mixAK"),
                                k_range = 2:5, ...) {
  valid <- c("GBMT", "KML3D", "flexmix", "LCMM", "mixAK")
  unknown <- setdiff(methods, valid)
  if (length(unknown))
    stop("Unknown method(s): ", paste(unknown, collapse = ", "))

  dots    <- list(...)
  results <- list()

  for (m in methods) {
    cat(sprintf("\n=== %s ===\n", m))
    results[[m]] <- switch(m,
      "GBMT"    = do.call(cluster_gbmt,
                          c(list(dat, config, feature_vars, k_range),
                            dots[intersect(names(dots), c("d", "nstart"))])),
      "KML3D"   = do.call(cluster_kml3d,
                          c(list(dat, config, feature_vars, k_range),
                            dots[intersect(names(dots), "nb_redrawing")])),
      "flexmix" = cluster_flexmix(dat, config, feature_vars, k_range),
      "LCMM"    = do.call(cluster_lcmm,
                          c(list(dat, config, feature_vars, k_range),
                            dots[intersect(names(dots), "maxiter")])),
      "mixAK"   = do.call(cluster_mixak,
                          c(list(dat, config, feature_vars, k_range),
                            dots[intersect(names(dots), "nMCMC")]))
    )
  }

  # Build flat metrics table
  metrics_rows <- lapply(names(results), function(m) {
    lapply(names(results[[m]]), function(kk) {
      cbind(method = m, results[[m]][[kk]]$metrics)
    })
  })
  metrics_df <- do.call(rbind, unlist(metrics_rows, recursive = FALSE))
  rownames(metrics_df) <- NULL

  list(results = results, metrics = metrics_df)
}
