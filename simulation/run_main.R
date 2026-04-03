# =============================================================================
# run_main.R
# Phase 2+3 — Simulation + Clustering + Evaluation
#
# Reads pre-imputed CSVs directly; does NOT require run_impute.R to have been
# run in the same session.
#
# Usage:
#   Phase 2 (one seed):
#     Rscript run_main.R --phase run --seed 1
#
#   Phase 3 (evaluate after all seeds done):
#     Rscript run_main.R --phase evaluate
#     Rscript run_main.R --phase evaluate --seeds 1:50
#
#   Local sequential test:
#     Rscript run_main.R --phase all --seeds 1:5
#
#   SLURM array (Phase 2):
#     #SBATCH --array=1-50
#     Rscript run_main.R --phase run --seed $SLURM_ARRAY_TASK_ID
# =============================================================================


# dat=prepare_sce2(sim_dat, CFG, FEAT_SCE2, TIME_MAP)
# config       <- CFG
# feature_vars <- FEAT_SCE2
# methods      <- CLUSTER_METHODS
# k_range      <- K_RANGE
# dots <- list(
#   d            = CLUSTER_ARGS$d,
#   nstart       = CLUSTER_ARGS$nstart,
#   nb_redrawing = CLUSTER_ARGS$nb_redrawing,
#   maxiter      = CLUSTER_ARGS$maxiter,
#   nMCMC        = CLUSTER_ARGS$nMCMC
# )

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

DIR_BASE    <- "output"
DIR_SIM     <- file.path(DIR_BASE, "sim")   # where per-seed sim datasets are saved
K_RANGE     <- 2:5
N_PER_GROUP <- 200

DATA_PATH   <- file.path(script_dir, "dat_sce2.csv")   # original data (for eCDF-Copula)

# --- Imputed dataset paths ---
# Edit these to point at your pre-imputed CSVs, e.g.:
#   "CC" = "/Users/hazel/Dropbox/jinyuan/Project2/0pipeline/output/imputed/impt_CC.csv"
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

FEAT_SCE2 <- c("ma_tot", "hars_score")
FEAT_SCE3 <- c(paste0("ma_", 1:10), paste0("hars_", 1:14, "_sev"))
N_PC_SCE3 <- 10

# mixAK is excluded here — run it separately with run_mixAK.R
CLUSTER_METHODS <- c("GBMT", "KML3D", "flexmix", "LCMM")
CLUSTER_ARGS    <- list(d = 2, nstart = 3, nb_redrawing = 5, maxiter = 200,
                         nMCMC = c(burn = 400, keep = 800, thin = 3, info = 100))


# =============================================================================
# Load imputed datasets from IMPUTED_FILES
# =============================================================================

load_imputed <- function() {
  missing <- Filter(function(p) !file.exists(p), IMPUTED_FILES)
  if (length(missing)) {
    stop("Imputed file(s) not found:\n",
         paste(" ", names(missing), "->", unlist(missing), collapse = "\n"),
         "\nRun run_impute.R first, or update IMPUTED_FILES paths in run_main.R.")
  }
  lapply(IMPUTED_FILES, read_csv, show_col_types = FALSE)
}


# =============================================================================
# Phase 2: One seed
# =============================================================================

run_one_seed <- function(seed, impt_list = NULL) {
  cat(sprintf("\n========== Seed %d ==========\n", seed))

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

  # Save per-seed sim datasets so run_mixAK.R can use them later
  dir.create(DIR_SIM, recursive = TRUE, showWarnings = FALSE)
  for (m in names(sim_list)) {
    fname <- file.path(DIR_SIM,
                       sprintf("sim_seed%03d_%s.csv", seed,
                               gsub("[^A-Za-z0-9]", "_", m)))
    write_csv(sim_list[[m]], fname)
    cat(sprintf("  Saved sim: %s\n", fname))
  }

  all_metrics <- list()

  for (sim_method in names(sim_list)) {
    sim_dat <- sim_list[[sim_method]]
    cat(sprintf("\n  Sim method: %s\n", sim_method))

    sce2_dat <- prepare_sce2(sim_dat, CFG, FEAT_SCE2, TIME_MAP)
    cl2      <- run_all_clustering(
      dat          = sce2_dat,
      config       = CFG,
      feature_vars = FEAT_SCE2,
      methods      = CLUSTER_METHODS,
      k_range      = K_RANGE,
      d            = CLUSTER_ARGS$d,
      nstart       = CLUSTER_ARGS$nstart,
      nb_redrawing = CLUSTER_ARGS$nb_redrawing,
      maxiter      = CLUSTER_ARGS$maxiter,
      nMCMC        = CLUSTER_ARGS$nMCMC
    )

    sce2_metrics <- cl2$metrics %>%
      rename(cluster_method = method) %>%
      mutate(scenario = "SCE2", sim_method = sim_method)

    sce3_obj <- prepare_sce3(sim_dat, CFG, FEAT_SCE3, N_PC_SCE3, TIME_MAP)
    cl3      <- run_all_clustering(sce3_obj$data, CFG, sce3_obj$pc_names,
                                    CLUSTER_METHODS, K_RANGE,
                                    d            = CLUSTER_ARGS$d,
                                    nstart       = CLUSTER_ARGS$nstart,
                                    nb_redrawing = CLUSTER_ARGS$nb_redrawing,
                                    maxiter      = CLUSTER_ARGS$maxiter,
                                    nMCMC        = CLUSTER_ARGS$nMCMC)

    sce3_metrics <- cl3$metrics %>%
      rename(cluster_method = method) %>%
      mutate(scenario = "SCE3", sim_method = sim_method)

    all_metrics[[sim_method]] <- rbind(sce2_metrics, sce3_metrics)
  }

  metrics_df <- do.call(rbind, all_metrics)
  rownames(metrics_df) <- NULL

  out_dir <- file.path(DIR_BASE, "cluster")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(metrics_df,
            file.path(out_dir, sprintf("cluster_seed%03d.csv", seed)))

  list(metrics = metrics_df, seed = seed)
}


# =============================================================================
# Phase 3: Evaluation
# =============================================================================

run_evaluation_phase <- function(seeds = 1:50) {
  cat("=== Phase 3: Evaluation ===\n")

  cluster_dir  <- file.path(DIR_BASE, "cluster")
  seed_results <- setNames(
    lapply(seeds, function(s) {
      fname <- file.path(cluster_dir, sprintf("cluster_seed%03d.csv", s))
      if (!file.exists(fname)) { warning("Missing: ", fname); return(NULL) }
      list(metrics = read_csv(fname, show_col_types = FALSE), seed = s)
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
  out  <- list(phase = "run", seed = 1L, seeds = 1:50)
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
    "run"      = run_one_seed(args$seed),
    "evaluate" = run_evaluation_phase(args$seeds),
    "all"      = {
      impt <- load_imputed()
      seed_results <- setNames(
        lapply(args$seeds, function(s) run_one_seed(s, impt_list = impt)),
        as.character(args$seeds)
      )
      run_evaluation(seed_results, dir_output = file.path(DIR_BASE, "results"))
    },
    stop("Unknown --phase: ", args$phase, "\nValid: run | evaluate | all")
  )
}
