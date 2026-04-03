# =============================================================================
# 01_impute.R
# Five imputation strategies for multivariate longitudinal data
#
# Methods (named as in the manuscript):
#   CC      Complete Case          impute_complete_case()
#   MICE-L  MICE-Longitudinal      impute_mice_longitudinal()
#   LME-F   LME-Full               impute_lme(..., selective = FALSE)
#   LME-S   LME-Selective          impute_lme(..., selective = TRUE)
#   MICE-CS MICE-Cross-sectional   impute_mice_crosssectional()
#
# All functions share the interface:
#   impute_*(data, config, ...)
#
#   data    long-format data frame
#   config  list produced by make_config()
#
# All functions return a long-format data frame with the same column structure
# as the input, with NAs in outcome_vars filled in.
#
# NOTE on MICE m parameter
# ─────────────────────────
# All MICE functions use m = 1 (single imputation). This is intentional:
# each method produces ONE completed dataset that feeds into the eCDF-Copula
# simulation (02_simulate.R). There is no standard Rubin's-rules pooling for
# cluster assignments, and the 50-seed replication in the simulation already
# covers downstream uncertainty. Using m > 1 without pooling (e.g. taking
# complete(imp, 1)) is equivalent to m = 1 but wastes computation.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mice)
  library(nlme)
})


# =============================================================================
# Config
# =============================================================================

#' Build an imputation config list
#'
#' @param subject_var      Column name for subject ID
#' @param time_var         Column name for time/visit
#' @param group_var        Column name for group label (kept but never imputed)
#' @param outcome_vars     Character vector: variables to impute
#' @param demographic_vars Character vector: time-invariant covariates used as
#'                         predictors in MICE methods (never imputed themselves)
#' @param time_levels      Ordered character vector of ALL expected time levels
#' @param lme_knot_index   Integer (0-indexed): position of the broken-stick
#'                         spline knot. E.g. 3 means timeplus = pmax(0, time_num - 3).
#'                         Set NULL for a plain linear trend (no spline).
#' @param dir_output       Directory for CSV output when save = TRUE
#'
#' @examples
#' # Original SCE2/SCE3 dataset
#' cfg <- make_config(
#'   subject_var      = "subj_num",
#'   time_var         = "visit",
#'   group_var        = "condition",
#'   outcome_vars     = c("ma_tot", "hars_score",
#'                        paste0("ma_", 1:10),
#'                        paste0("hars_", 1:14, "_sev")),
#'   demographic_vars = c("age", "sex", "sex_xcount", "race", "ethnicity"),
#'   time_levels      = c("screening", "baseline",
#'                        paste0("month_", sprintf("%02d", seq(2, 24, by = 2)))),
#'   lme_knot_index   = 3
#' )
make_config <- function(subject_var      = "Subject",
                        time_var         = "Time",
                        group_var        = "Group",
                        outcome_vars,
                        demographic_vars = character(0),
                        time_levels,
                        lme_knot_index   = NULL,
                        dir_output       = "output/imputed") {
  stopifnot(
    is.character(outcome_vars), length(outcome_vars) >= 1,
    is.character(time_levels),  length(time_levels) >= 2
  )
  list(
    subject_var      = subject_var,
    time_var         = time_var,
    group_var        = group_var,
    outcome_vars     = outcome_vars,
    demographic_vars = demographic_vars,
    time_levels      = time_levels,
    n_timepoints     = length(time_levels),
    lme_knot_index   = lme_knot_index,
    dir_output       = dir_output
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# Expand to a complete subject × time grid and left-join observed data.
# Subjects missing some time points get NA for those rows.
.expand_grid <- function(data, config) {
  sv <- config$subject_var
  tv <- config$time_var
  grid <- expand.grid(
    setNames(list(unique(data[[sv]]), config$time_levels), c(sv, tv)),
    stringsAsFactors = FALSE
  )
  grid[[tv]] <- factor(grid[[tv]], levels = config$time_levels)
  left_join(grid, data, by = c(sv, tv)) %>%
    arrange(.data[[sv]], .data[[tv]])
}

# Extract one row per subject keeping time-invariant columns (group + demographics).
.get_demo <- function(full_dat, config) {
  sv   <- config$subject_var
  keep <- intersect(c(sv, config$group_var, config$demographic_vars),
                    names(full_dat))
  full_dat %>%
    group_by(.data[[sv]]) %>%
    slice(1) %>%
    ungroup() %>%
    select(all_of(keep))
}


# =============================================================================
# CC — Complete Case
# =============================================================================

#' CC: retain only subjects with complete records across all time points
#'
#' No imputation is performed. A subject is kept if and only if they have
#' exactly n_timepoints rows AND no NAs in any of the checked variables.
#'
#' @param key_vars  Subset of outcome_vars to check for completeness.
#'                  Defaults to all outcome_vars.
impute_complete_case <- function(data, config, key_vars = NULL) {
  cat("Imputation: CC (Complete Case)\n")
  if (is.null(key_vars)) key_vars <- config$outcome_vars
  sv <- config$subject_var

  complete_ids <- data %>%
    group_by(.data[[sv]]) %>%
    filter(
      n() == config$n_timepoints,
      all(!is.na(across(all_of(key_vars))))
    ) %>%
    ungroup() %>%
    pull(.data[[sv]]) %>%
    unique()

  n_kept  <- length(complete_ids)
  n_total <- length(unique(data[[sv]]))
  cat(sprintf("  Retained %d / %d subjects with complete data\n", n_kept, n_total))

  data[data[[sv]] %in% complete_ids, ]
}


# =============================================================================
# MICE-L — MICE-Longitudinal
# =============================================================================

#' MICE-L: MICE exploiting within-subject temporal correlation
#'
#' Each outcome variable is pivoted wide (one column per time point). MICE
#' then uses the subject's own trajectory across other time points (plus
#' demographics) to impute missing values via predictive mean matching (PMM).
#'
#' Uses m = 1 (single imputation); see file-level note for rationale.
#'
#' @param seed  Random seed (default 123)
impute_mice_longitudinal <- function(data, config, seed = 123) {
  cat("Imputation: MICE-L (MICE-Longitudinal)\n")
  sv <- config$subject_var
  tv <- config$time_var

  full_dat <- .expand_grid(data, config) %>%
    group_by(.data[[sv]]) %>%
    fill(all_of(c(config$group_var, config$demographic_vars)), .direction = "updown") %>%
    ungroup()
  demo_dat <- .get_demo(full_dat, config)

  impute_one_var <- function(varname) {
    wide <- full_dat %>%
      select(all_of(c(sv, tv, varname))) %>%
      pivot_wider(
        names_from   = all_of(tv),
        values_from  = all_of(varname),
        names_prefix = paste0(varname, "_")
      ) %>%
      left_join(demo_dat, by = sv) %>%
      mutate(across(where(is.character), as.factor)) %>%
      as.data.frame()

    outcome_cols <- grep(paste0("^", varname, "_"), colnames(wide), value = TRUE)

    ini  <- mice(wide, maxit = 0, printFlag = FALSE)
    pred <- ini$predictorMatrix
    meth <- ini$method

    # Subject ID: not a predictor and not imputed
    pred[sv, ] <- 0; pred[, sv] <- 0; meth[sv] <- ""

    # Non-outcome columns (group, demographics): predictor-only, never imputed
    non_outcome <- setdiff(colnames(wide), c(sv, outcome_cols))
    if (length(non_outcome)) {
      pred[non_outcome, ] <- 0
      meth[non_outcome]   <- ""
    }

    # All-NA outcome columns: no donors for PMM, skip imputation
    all_na_cols <- outcome_cols[vapply(wide[outcome_cols],
                                       function(x) all(is.na(x)), logical(1))]
    if (length(all_na_cols)) {
      pred[all_na_cols, ] <- 0; pred[, all_na_cols] <- 0
      meth[all_na_cols]   <- ""
    }

    diag(pred) <- 0

    imp       <- mice(wide, method = meth, predictorMatrix = pred,
                      m = 1, maxit = 10, seed = seed, printFlag = FALSE)
    completed <- complete(imp, 1)

    completed %>%
      select(all_of(sv), starts_with(paste0(varname, "_"))) %>%
      pivot_longer(
        cols         = starts_with(paste0(varname, "_")),
        names_to     = tv,
        names_prefix = paste0(varname, "_"),
        values_to    = varname
      ) %>%
      mutate(!!tv := factor(.data[[tv]], levels = config$time_levels))
  }

  imputed_vars <- map(config$outcome_vars, impute_one_var) %>%
    reduce(full_join, by = c(sv, tv))

  # Replace outcome columns in full_dat with imputed versions
  full_dat %>%
    select(-all_of(config$outcome_vars)) %>%
    left_join(imputed_vars, by = c(sv, tv)) %>%
    arrange(.data[[sv]], .data[[tv]])
}


# =============================================================================
# LME-F / LME-S — LME-Full and LME-Selective
# =============================================================================

#' LME-F / LME-S: linear mixed-effects model imputation
#'
#' Fits one LME per outcome variable:
#'   outcome ~ group * time_num  [ + group * timeplus if lme_knot_index is set ]
#' with a subject-specific random intercept and random slope for time.
#' Predictions are floored at 0 (appropriate for non-negative clinical scores).
#'
#' LME-F (selective = FALSE): replace ALL values with model predictions
#' LME-S (selective = TRUE):  replace only NAs; keep observed values as-is
#'
#' If the full broken-stick model fails to converge, a plain linear model is
#' used as fallback.
#'
#' @param selective  FALSE → LME-F (Method 3); TRUE → LME-S (Method 4)
impute_lme <- function(data, config, selective = FALSE) {
  label <- if (selective) "LME-S (LME-Selective)" else "LME-F (LME-Full)"
  cat("Imputation:", label, "\n")

  sv <- config$subject_var
  tv <- config$time_var
  gv <- config$group_var

  full_dat <- .expand_grid(data, config) %>%
    group_by(.data[[sv]]) %>%
    fill(all_of(c(config$group_var, config$demographic_vars)), .direction = "updown") %>%
    ungroup() %>%
    mutate(
      time_num = as.integer(factor(.data[[tv]], levels = config$time_levels)) - 1L
    )

  if (!is.null(config$lme_knot_index))
    full_dat$timeplus <- pmax(0L, full_dat$time_num - as.integer(config$lme_knot_index))

  has_knot  <- !is.null(config$lme_knot_index)
  rand_form <- as.formula(sprintf("~ time_num | %s", sv))

  predictions <- setNames(vector("list", length(config$outcome_vars)),
                           config$outcome_vars)

  for (var in config$outcome_vars) {
    cat(sprintf("  Fitting: %s\n", var))

    f_complex <- as.formula(
      if (has_knot)
        sprintf("%s ~ %s * (time_num + timeplus)", var, gv)
      else
        sprintf("%s ~ %s * time_num", var, gv)
    )
    f_simple <- as.formula(sprintf("%s ~ %s * time_num", var, gv))

    fit <- tryCatch(
      lme(f_complex, random = rand_form, data = full_dat,
          method = "REML", na.action = na.omit,
          control = lmeControl(returnObject = TRUE)),
      error = function(e) {
        message("  -> fallback to simple model for ", var)
        lme(f_simple, random = rand_form, data = full_dat,
            method = "REML", na.action = na.omit,
            control = lmeControl(returnObject = TRUE))
      }
    )
    predictions[[var]] <- predict(fit, newdata = full_dat, level = 1)
  }

  result <- full_dat %>% select(all_of(c(sv, tv, gv)))

  for (var in config$outcome_vars) {
    result[[var]] <- if (selective) {
      ifelse(is.na(full_dat[[var]]), predictions[[var]], full_dat[[var]])
    } else {
      predictions[[var]]
    }
  }

  # Floor predictions at 0 (LME can produce negative values)
  result[config$outcome_vars] <- lapply(result[config$outcome_vars],
                                         function(x) pmax(x, 0))

  result %>% arrange(.data[[sv]], .data[[tv]])
}


# =============================================================================
# MICE-CS — MICE-Cross-sectional
# =============================================================================

#' MICE-CS: MICE using within-visit cross-sectional information
#'
#' Missing values at each time point are imputed using demographic covariates
#' and other outcome variables observed at the same visit. Time/visit and group
#' columns serve as predictors but are never imputed.
#'
#' Uses m = 1 (single imputation); see file-level note for rationale.
#' The original code used m = 5 but discarded 4 of the 5 datasets via
#' complete(imp, 1), which was inconsistent. This is corrected here.
#'
#' @param seed  Random seed (default 123)
impute_mice_crosssectional <- function(data, config, seed = 123) {
  cat("Imputation: MICE-CS (MICE-Cross-sectional)\n")
  sv <- config$subject_var
  tv <- config$time_var
  gv <- config$group_var

  fixed_vars <- c(gv, tv, config$demographic_vars)
  keep_cols  <- intersect(c(fixed_vars, config$outcome_vars), names(data))
  data_mice  <- as.data.frame(data)[, keep_cols, drop = FALSE]

  ini  <- mice(data_mice, maxit = 0, printFlag = FALSE)
  pred <- ini$predictorMatrix

  # Fixed vars: not imputed (row = 0) but used as predictors (col = 1)
  for (v in intersect(fixed_vars, rownames(pred))) pred[v, ] <- 0
  for (v in intersect(fixed_vars, colnames(pred))) pred[, v] <- 1
  diag(pred) <- 0

  imp       <- mice(data_mice, method = "pmm", predictorMatrix = pred,
                    m = 1, maxit = 10, seed = seed, printFlag = FALSE)
  completed <- complete(imp, 1)

  result <- as.data.frame(data)
  result[, config$outcome_vars] <- completed[, config$outcome_vars]
  result %>% arrange(.data[[sv]], .data[[tv]])
}


# =============================================================================
# Convenience wrapper
# =============================================================================

#' Run a subset of imputation methods and return a named list
#'
#' @param methods  Character vector of method names, any subset of:
#'                 c("CC", "MICE-L", "LME-F", "LME-S", "MICE-CS")
#'                 Defaults to all five.
#' @param seed     Random seed for stochastic methods (MICE-L, MICE-CS)
#' @param save     If TRUE, save each result as impt_{method}.csv in config$dir_output
#'
#' @return Named list of imputed data frames, keyed by method abbreviation
run_all_imputations <- function(data, config,
                                methods = c("CC", "MICE-L", "LME-F", "LME-S", "MICE-CS"),
                                seed = 123, save = FALSE) {
  valid_methods <- c("CC", "MICE-L", "LME-F", "LME-S", "MICE-CS")
  unknown <- setdiff(methods, valid_methods)
  if (length(unknown))
    stop("Unknown method(s): ", paste(unknown, collapse = ", "),
         "\nValid options: ", paste(valid_methods, collapse = ", "))

  if (save)
    dir.create(config$dir_output, recursive = TRUE, showWarnings = FALSE)

  results <- list()

  for (m in methods) {
    cat(sprintf("\n=== %s ===\n", m))

    res <- switch(m,
      "CC"      = impute_complete_case(data, config),
      "MICE-L"  = impute_mice_longitudinal(data, config, seed = seed),
      "LME-F"   = impute_lme(data, config, selective = FALSE),
      "LME-S"   = impute_lme(data, config, selective = TRUE),
      "MICE-CS" = impute_mice_crosssectional(data, config, seed = seed)
    )

    results[[m]] <- res

    if (save) {
      fname <- file.path(config$dir_output,
                         paste0("impt_", gsub("-", "_", m), ".csv"))
      write_csv(res, fname)
      cat(sprintf("  Saved: %s\n", fname))
    }
  }

  results
}

