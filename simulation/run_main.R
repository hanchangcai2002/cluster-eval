# =============================================================================
# run_main.R
# Three independent phases — each reads from disk, writes to disk.
#
# Phase A  prepare   simulate + prepare SCE2/SCE3 → save to output/prepared/
# Phase B  cluster   non-mixAK clustering on prepared data → output/cluster/
# Phase C  mixak     mixAK clustering on prepared data    → output/cluster/
# Phase D  evaluate  aggregate all seed results           → output/results/
#
# Usage (SLURM — run A then B then C in separate array jobs):
#   Rscript run_main.R --phase prepare  --seed $SLURM_ARRAY_TASK_ID
#   Rscript run_main.R --phase cluster  --seed $SLURM_ARRAY_TASK_ID
#   Rscript run_main.R --phase mixak    --seed $SLURM_ARRAY_TASK_ID
#   Rscript run_main.R --phase evaluate
#
# Local sequential test (all phases, seeds 1-5):
#   Rscript run_main.R --phase all --seeds 1:5
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- grep("--file=", args, value = TRUE)
  if (length(f)) {
    normalizePath(dirname(sub("--file=", "", f[1])))
  } else {
    tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
             error = function(e) getwd())
  }
})
for (f in c("00_utils.R", "01_impute.R", "02_simulate.R", "03_prepare.R",
             "04_cluster.R", "05_evaluate.R")) source(file.path(script_dir, f))


# =============================================================================
# USER CONFIGURATION
# =============================================================================

DIR_BASE      <- "output"
DIR_SIM       <- file.path(DIR_BASE, "sim")
DIR_PREPARED  <- file.path(DIR_BASE, "prepared")
K_RANGE       <- 2:5
N_PER_GROUP   <- 200

DATA_PATH     <- file.path(script_dir, "dat_sce2.csv")

IMPUTED_FILES <- list(
  "CC"      = file.path(DIR_BASE, "imputed", "impt_CC.csv"),
  "MICE-L"  = file.path(DIR_BASE, "imputed", "impt_MICE_L.csv"),
  "LME-F"   = file.path(DIR_BASE, "imputed", "impt_LME_F.csv"),
  "LME-S"   = file.path(DIR_BASE, "imputed", "impt_LME_S.csv"),
  "MICE-CS" = file.path(DIR_BASE, "imputed", "impt_MICE_CS.csv")
)

CFG <- make_config(
  subject_var      = "subj_num",
  time_var         = "visit",
  group_var        = "condition",
  outcome_vars     = c("ma_tot", "hars_score",
                       paste0("ma_", 1:10),
                       paste0("hars_", 1:14, "_sev")),
  demographic_vars = c("age", "sex", "sex_xcount", "race", "ethnicity"),
  time_levels      = c("screening", "baseline",
                       paste0("month_", sprintf("%02d", seq(2, 24, by = 2)))),
  lme_knot_index   = 3,
  dir_output       = file.path(DIR_BASE, "imputed")
)

TIME_MAP <- make_time_map(
  CFG$time_levels,
  c(-1, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
)

FEAT_SCE2       <- c("ma_tot", "hars_score")
FEAT_SCE3       <- c(paste0("ma_", 1:10), paste0("hars_", 1:14, "_sev"))
N_PC_SCE3       <- 10

CLUSTER_METHODS <- c("GBMT", "KML3D", "flexmix", "LCMM")
CLUSTER_METHODS_MIXAK <- c("mixAK")

CLUSTER_ARGS    <- list(d = 2, nstart = 3, nb_redrawing = 5, maxiter = 200,
                        nMCMC = c(burn = 400, keep = 800, thin = 2, info = 100),
                        lcmm_timeout = 7200L)   # max elapsed seconds for entire LCMM fit


# =============================================================================
# File path helpers
# =============================================================================

# Hyphens are valid in filenames — only replace truly problematic characters
.safe <- function(x) gsub("[^A-Za-z0-9-]", "_", x)

.sce2_path <- function(seed, sim_method)
  file.path(DIR_PREPARED, sprintf("sce2_seed%03d_%s.csv", seed, .safe(sim_method)))

.sce3_path <- function(seed, sim_method)
  file.path(DIR_PREPARED, sprintf("sce3_seed%03d_%s.csv", seed, .safe(sim_method)))

# Reads back the sim method names saved during prepare phase
.load_sim_methods <- function(seed) {
  p <- file.path(DIR_PREPARED, sprintf("sim_methods_seed%03d.rds", seed))
  if (!file.exists(p))
    stop("Prepared files not found for seed ", seed, ". Run --phase prepare first.")
  readRDS(p)
}

# Columns that are not clustering features in the prepared CSVs
.meta_vars <- function() c(CFG$subject_var, CFG$time_var, CFG$group_var, "time_num")


# =============================================================================
# Load imputed datasets
# =============================================================================

load_imputed <- function() {
  missing <- Filter(function(p) !file.exists(p), IMPUTED_FILES)
  if (length(missing))
    stop("Imputed file(s) not found:\n",
         paste(" ", names(missing), "->", unlist(missing), collapse = "\n"),
         "\nRun run_impute.R first, or update IMPUTED_FILES paths in run_main.R.")
  lapply(IMPUTED_FILES, read_csv, show_col_types = FALSE)
}


# =============================================================================
# Phase A: simulate + prepare SCE2 and SCE3, save to disk
# =============================================================================

run_one_seed_prepare <- function(seed, impt_list = NULL) {
  cat(sprintf("\n========== Seed %d | Phase A: prepare ==========\n", seed))

  if (is.null(impt_list)) impt_list <- load_imputed()
  raw_data <- read_csv(DATA_PATH, show_col_types = FALSE)

  cat("Simulating...\n")
  sim_list <- run_all_simulations(
    imputed_list = impt_list,
    raw_data     = raw_data,
    config       = CFG,
    n_per_group  = N_PER_GROUP,
    seed         = seed,
    save         = FALSE
  )

  dir.create(DIR_SIM,      recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_PREPARED, recursive = TRUE, showWarnings = FALSE)

  for (m in names(sim_list)) {
    sim_dat <- sim_list[[m]]

    # Raw sim (saved for reference / debugging)
    write_csv(sim_dat,
              file.path(DIR_SIM,
                        sprintf("sim_seed%03d_%s.csv", seed, .safe(m))))

    # SCE2 prepared
    sce2_dat <- prepare_sce2(sim_dat, CFG, FEAT_SCE2, TIME_MAP)
    write_csv(sce2_dat, .sce2_path(seed, m))

    # SCE3 prepared
    sce3_obj <- prepare_sce3(sim_dat, CFG, FEAT_SCE3, N_PC_SCE3, TIME_MAP)
    write_csv(sce3_obj$data, .sce3_path(seed, m))

    cat(sprintf("  Prepared: %s\n", m))
  }

  # Save original method names so cluster phases can recover them
  saveRDS(names(sim_list),
          file.path(DIR_PREPARED, sprintf("sim_methods_seed%03d.rds", seed)))

  invisible(NULL)
}


# =============================================================================
# Phase B: non-mixAK clustering on prepared data
# =============================================================================

run_one_seed_cluster <- function(seed) {
  cat(sprintf("\n========== Seed %d | Phase B: cluster ==========\n", seed))

  sim_methods <- .load_sim_methods(seed)
  meta_vars   <- .meta_vars()

  out_dir <- file.path(DIR_BASE, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Helper: partial-file paths per (seed, sim_method, scenario)
  .partial_path <- function(m, sce)
    file.path(out_dir, sprintf("cluster_seed%03d_%s_%s.csv", seed, .safe(m), sce))

  for (m in sim_methods) {
    cat(sprintf("\n  Sim method: %s\n", m))

    # --- SCE2 ---
    sce2_f <- .partial_path(m, "sce2")
    if (file.exists(sce2_f)) {
      cat(sprintf("  [skip] SCE2 already saved: %s\n", basename(sce2_f)))
    } else {
      sce2_dat <- read_csv(.sce2_path(seed, m), show_col_types = FALSE)
      cl2 <- run_all_clustering(
        dat          = sce2_dat,
        config       = CFG,
        feature_vars = FEAT_SCE2,
        methods      = CLUSTER_METHODS,
        k_range      = K_RANGE,
        d            = CLUSTER_ARGS$d,
        nstart       = CLUSTER_ARGS$nstart,
        nb_redrawing = CLUSTER_ARGS$nb_redrawing,
        maxiter      = CLUSTER_ARGS$maxiter,
        nMCMC        = CLUSTER_ARGS$nMCMC,
        lcmm_timeout = CLUSTER_ARGS$lcmm_timeout
      )
      sce2_metrics <- cl2$metrics %>%
        rename(cluster_method = method) %>%
        mutate(scenario = "SCE2", sim_method = m)
      write_csv(sce2_metrics, sce2_f)
      cat(sprintf("  Saved SCE2 partial: %s\n", basename(sce2_f)))
    }

    # --- SCE3 (wrapped: LCMM can fail/hang; save whatever completes) ---
    sce3_f <- .partial_path(m, "sce3")
    if (file.exists(sce3_f)) {
      cat(sprintf("  [skip] SCE3 already saved: %s\n", basename(sce3_f)))
    } else {
      tryCatch({
        sce3_dat <- read_csv(.sce3_path(seed, m), show_col_types = FALSE)
        pc_names <- setdiff(names(sce3_dat), meta_vars)
        cl3 <- run_all_clustering(
          dat          = sce3_dat,
          config       = CFG,
          feature_vars = pc_names,
          methods      = CLUSTER_METHODS,
          k_range      = K_RANGE,
          d            = CLUSTER_ARGS$d,
          nstart       = CLUSTER_ARGS$nstart,
          nb_redrawing = CLUSTER_ARGS$nb_redrawing,
          maxiter      = CLUSTER_ARGS$maxiter,
          nMCMC        = CLUSTER_ARGS$nMCMC,
          lcmm_timeout = CLUSTER_ARGS$lcmm_timeout
        )
        sce3_metrics <- cl3$metrics %>%
          rename(cluster_method = method) %>%
          mutate(scenario = "SCE3", sim_method = m)
        write_csv(sce3_metrics, sce3_f)
        cat(sprintf("  Saved SCE3 partial: %s\n", basename(sce3_f)))
      }, error = function(e) {
        cat(sprintf("  WARNING: SCE3 failed for sim method '%s': %s\n",
                    m, conditionMessage(e)))
      })
    }
  }

  # Combine all available partial results into the final seed file
  partial_list <- list()
  for (m in sim_methods) {
    for (sce in c("sce2", "sce3")) {
      f <- .partial_path(m, sce)
      if (file.exists(f))
        partial_list[[paste0(m, "_", sce)]] <-
          read_csv(f, show_col_types = FALSE)
    }
  }

  if (length(partial_list) == 0) {
    warning("No partial results available for seed ", seed,
            " — nothing written.")
    return(invisible(NULL))
  }

  metrics_df <- do.call(rbind, partial_list)
  rownames(metrics_df) <- NULL
  write_csv(metrics_df,
            file.path(out_dir, sprintf("cluster_seed%03d.csv", seed)))
  cat(sprintf("  Saved: cluster_seed%03d.csv  (%d rows from %d partial file(s))\n",
              seed, nrow(metrics_df), length(partial_list)))

  invisible(metrics_df)
}


# =============================================================================
# Phase C: mixAK clustering on prepared data
# =============================================================================

run_one_seed_mixak <- function(seed) {
  cat(sprintf("\n========== Seed %d | Phase C: mixAK ==========\n", seed))

  sim_methods <- .load_sim_methods(seed)
  meta_vars   <- .meta_vars()

  out_dir   <- file.path(DIR_BASE, "cluster")
  trace_dir <- file.path(out_dir, "traceplots")
  dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(trace_dir, recursive = TRUE, showWarnings = FALSE)

  .mixak_partial <- function(m, sce)
    file.path(out_dir, sprintf("mixak_seed%03d_%s_%s.csv", seed, .safe(m), sce))

  .trace_prefix <- function(m, sce)
    file.path(trace_dir, sprintf("mixak_trace_seed%03d_%s_%s", seed, .safe(m), sce))

  MIXAK_SCE3_TIMEOUT <- 4L * 3600L   # 4 hours per sim_method SCE3 section

  for (m in sim_methods) {
    cat(sprintf("\n  Sim method: %s\n", m))

    # --- SCE2 ---
    sce2_f <- .mixak_partial(m, "sce2")
    if (file.exists(sce2_f)) {
      cat(sprintf("  [skip] SCE2 already saved: %s\n", basename(sce2_f)))
    } else {
      sce2_dat <- read_csv(.sce2_path(seed, m), show_col_types = FALSE)
      cl2 <- run_all_clustering(
        dat          = sce2_dat,
        config       = CFG,
        feature_vars = FEAT_SCE2,
        methods      = CLUSTER_METHODS_MIXAK,
        k_range      = K_RANGE,
        d            = CLUSTER_ARGS$d,
        nstart       = CLUSTER_ARGS$nstart,
        nb_redrawing = CLUSTER_ARGS$nb_redrawing,
        maxiter      = CLUSTER_ARGS$maxiter,
        trace_prefix = .trace_prefix(m, "sce2"),
        scenario     = "SCE2"
      )
      sce2_metrics <- cl2$metrics %>%
        rename(cluster_method = method) %>%
        mutate(scenario = "SCE2", sim_method = m)
      write_csv(sce2_metrics, sce2_f)
      cat(sprintf("  Saved SCE2 partial: %s\n", basename(sce2_f)))
    }

    # --- SCE3 (wrapped: may be very slow or crash; 4h timelimit per method) ---
    sce3_f <- .mixak_partial(m, "sce3")
    if (file.exists(sce3_f)) {
      cat(sprintf("  [skip] SCE3 already saved: %s\n", basename(sce3_f)))
    } else {
      tryCatch({
        setTimeLimit(elapsed = MIXAK_SCE3_TIMEOUT, transient = FALSE)
        sce3_dat <- read_csv(.sce3_path(seed, m), show_col_types = FALSE)
        pc_names <- setdiff(names(sce3_dat), meta_vars)
        cl3 <- run_all_clustering(
          dat          = sce3_dat,
          config       = CFG,
          feature_vars = pc_names,
          methods      = CLUSTER_METHODS_MIXAK,
          k_range      = K_RANGE,
          d            = CLUSTER_ARGS$d,
          nstart       = CLUSTER_ARGS$nstart,
          nb_redrawing = CLUSTER_ARGS$nb_redrawing,
          maxiter      = CLUSTER_ARGS$maxiter,
          trace_prefix = .trace_prefix(m, "sce3"),
          scenario     = "SCE3"
        )
        setTimeLimit(elapsed = Inf)
        sce3_metrics <- cl3$metrics %>%
          rename(cluster_method = method) %>%
          mutate(scenario = "SCE3", sim_method = m)
        write_csv(sce3_metrics, sce3_f)
        cat(sprintf("  Saved SCE3 partial: %s\n", basename(sce3_f)))
      }, error = function(e) {
        setTimeLimit(elapsed = Inf)
        cat(sprintf("  WARNING: SCE3 skipped for '%s': %s\n", m, conditionMessage(e)))
      })
    }
  }

  # Combine all available partials into final mixak seed file
  partial_list <- list()
  for (m in sim_methods) {
    for (sce in c("sce2", "sce3")) {
      f <- .mixak_partial(m, sce)
      if (file.exists(f))
        partial_list[[paste0(m, "_", sce)]] <- read_csv(f, show_col_types = FALSE)
    }
  }

  if (length(partial_list) == 0) {
    warning("No mixAK partial results for seed ", seed, " — nothing written.")
    return(invisible(NULL))
  }

  metrics_df <- do.call(rbind, partial_list)
  rownames(metrics_df) <- NULL
  write_csv(metrics_df,
            file.path(out_dir, sprintf("mixak_seed%03d.csv", seed)))
  cat(sprintf("  Saved: mixak_seed%03d.csv  (%d rows from %d partial file(s))\n",
              seed, nrow(metrics_df), length(partial_list)))

  invisible(metrics_df)
}


# =============================================================================
# Phase D: Evaluation — loads cluster + mixak (if present), combines, evaluates
# =============================================================================

run_evaluation_phase <- function(seeds = 1:50) {
  cat("=== Phase D: Evaluation ===\n")

  cluster_dir  <- file.path(DIR_BASE, "cluster")
  seed_results <- setNames(
    lapply(seeds, function(s) {
      main_f  <- file.path(cluster_dir, sprintf("cluster_seed%03d.csv", s))
      mixak_f <- file.path(cluster_dir, sprintf("mixak_seed%03d.csv",   s))

      if (!file.exists(main_f)) {
        warning("Missing cluster file: ", main_f)
        return(NULL)
      }

      metrics <- read_csv(main_f, show_col_types = FALSE)
      if (file.exists(mixak_f))
        metrics <- rbind(metrics, read_csv(mixak_f, show_col_types = FALSE))

      list(metrics = metrics, seed = s)
    }),
    as.character(seeds)
  )
  seed_results <- Filter(Negate(is.null), seed_results)
  cat(sprintf("Loaded %d / %d seed results\n", length(seed_results), length(seeds)))

  run_evaluation(seed_results, dir_output = file.path(DIR_BASE, "results"))
}


# =============================================================================
# Command-line dispatch
# =============================================================================

.parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- list(phase = "prepare", seed = 1L, seeds = 1:50)
  i <- 1L
  while (i <= length(args)) {
    if      (args[i] == "--phase") { out$phase <- args[i + 1L]; i <- i + 2L }
    else if (args[i] == "--seed")  { out$seed  <- as.integer(args[i + 1L]); i <- i + 2L }
    else if (args[i] == "--seeds") { out$seeds <- eval(parse(text = args[i + 1L])); i <- i + 2L }
    else i <- i + 1L
  }
  out
}

if (!interactive()) {
  args <- .parse_args()

  dir.create(DIR_BASE, recursive = TRUE, showWarnings = FALSE)
  sink(file.path(DIR_BASE, "session_info.txt")); sessionInfo(); sink()

  switch(args$phase,
    "prepare"  = run_one_seed_prepare(args$seed),
    "cluster"  = run_one_seed_cluster(args$seed),
    "mixak"    = run_one_seed_mixak(args$seed),
    "evaluate" = run_evaluation_phase(args$seeds),
    "all"      = {
      impt <- load_imputed()
      for (s in args$seeds) {
        run_one_seed_prepare(s, impt_list = impt)
        run_one_seed_cluster(s)
        run_one_seed_mixak(s)
      }
      seed_results <- setNames(
        lapply(args$seeds, function(s) {
          cluster_dir <- file.path(DIR_BASE, "cluster")
          main_f  <- file.path(cluster_dir, sprintf("cluster_seed%03d.csv", s))
          mixak_f <- file.path(cluster_dir, sprintf("mixak_seed%03d.csv",   s))
          metrics <- read_csv(main_f, show_col_types = FALSE)
          if (file.exists(mixak_f))
            metrics <- rbind(metrics, read_csv(mixak_f, show_col_types = FALSE))
          list(metrics = metrics, seed = s)
        }),
        as.character(args$seeds)
      )
      run_evaluation(seed_results, dir_output = file.path(DIR_BASE, "results"))
    },
    stop("Unknown --phase: ", args$phase,
         "\nValid: prepare | cluster | mixak | evaluate | all")
  )
}
