# =============================================================================
# run_mixAK.R
# Run mixAK clustering on pre-saved simulated datasets and merge results
# into the existing cluster CSV produced by run_main.R.
#
# Prerequisites:
#   run_main.R --phase run --seed N   must have completed first
#   (creates output/sim/sim_seed{NNN}_{method}.csv and
#           output/cluster/cluster_seed{NNN}.csv)
#
# Usage:
#   Single seed:
#     Rscript run_mixAK.R --seed 1
#
#   Multiple seeds (sequential):
#     Rscript run_mixAK.R --seeds 1:10
#
# Output:
#   Appends mixAK rows into output/cluster/cluster_seed{NNN}.csv
#   (any pre-existing mixAK rows are replaced, so re-runs are safe)
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
for (f in c("00_utils.R", "03_prepare.R", "04_cluster.R")) source(file.path(script_dir, f))


# =============================================================================
# USER CONFIGURATION  (keep in sync with run_main.R)
# =============================================================================

DIR_BASE    <- "output"
DIR_SIM     <- file.path(DIR_BASE, "sim")
DIR_CLUSTER <- file.path(DIR_BASE, "cluster")

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
K_RANGE   <- 2:5

# --- mixAK MCMC parameters ---
# Tune these independently of run_main.R
NMCMC <- c(burn = 1000, keep = 2000, thin = 5, info = 500)


# =============================================================================
# Core: run mixAK for one seed
# =============================================================================

run_mixAK_one_seed <- function(seed) {
  cat(sprintf("\n========== mixAK | Seed %d ==========\n", seed))

  # Discover sim files for this seed
  sim_files <- list.files(DIR_SIM,
                           pattern  = sprintf("^sim_seed%03d_.*\\.csv$", seed),
                           full.names = TRUE)
  if (!length(sim_files)) {
    stop(sprintf(
      "No sim files found for seed %d in %s\n  Run: Rscript run_main.R --phase run --seed %d",
      seed, DIR_SIM, seed))
  }

  # Recover method names from filenames:
  #   sim_seed001_CC.csv       → "CC"
  #   sim_seed001_MICE_L.csv   → "MICE-L"   (underscores back to hyphens for 2+ consecutive)
  # Strategy: strip prefix "sim_seed{NNN}_" and suffix ".csv", then reverse the
  # gsub("[^A-Za-z0-9]", "_", m) by matching known method names.
  known_methods <- c("CC", "MICE-L", "LME-F", "LME-S", "MICE-CS", "eCDF-Copula")
  method_keys   <- setNames(
    gsub("[^A-Za-z0-9]", "_", known_methods),
    known_methods
  )  # e.g. "MICE_L" → "MICE-L"

  sim_list <- list()
  for (fp in sim_files) {
    raw_key <- sub(sprintf("^sim_seed%03d_", seed), "",
                   sub("\\.csv$", "", basename(fp)))
    method <- names(method_keys)[method_keys == raw_key]
    if (!length(method)) method <- raw_key   # fallback: use as-is
    sim_list[[method]] <- read_csv(fp, show_col_types = FALSE)
    cat(sprintf("  Loaded: %s  (%s)\n", basename(fp), method))
  }

  # Run mixAK for each sim method
  all_metrics <- list()

  for (sim_method in names(sim_list)) {
    sim_dat <- sim_list[[sim_method]]
    cat(sprintf("\n  Sim method: %s\n", sim_method))

    # SCE2
    cat("    SCE2...\n")
    sce2_dat     <- prepare_sce2(sim_dat, CFG, FEAT_SCE2, TIME_MAP)
    cl2          <- cluster_mixak(sce2_dat, CFG, FEAT_SCE2, K_RANGE, nMCMC = NMCMC)
    sce2_metrics <- do.call(rbind, lapply(cl2, `[[`, "metrics")) %>%
      mutate(method = "mixAK") %>%
      rename(cluster_method = method) %>%
      mutate(scenario = "SCE2", sim_method = sim_method)

    # SCE3
    cat("    SCE3...\n")
    sce3_obj     <- prepare_sce3(sim_dat, CFG, FEAT_SCE3, N_PC_SCE3, TIME_MAP)
    cl3          <- cluster_mixak(sce3_obj$data, CFG, sce3_obj$pc_names, K_RANGE, nMCMC = NMCMC)
    sce3_metrics <- do.call(rbind, lapply(cl3, `[[`, "metrics")) %>%
      mutate(method = "mixAK") %>%
      rename(cluster_method = method) %>%
      mutate(scenario = "SCE3", sim_method = sim_method)

    all_metrics[[sim_method]] <- rbind(sce2_metrics, sce3_metrics)
  }

  new_rows <- do.call(rbind, all_metrics)
  rownames(new_rows) <- NULL

  # Merge into existing cluster CSV
  cluster_file <- file.path(DIR_CLUSTER, sprintf("cluster_seed%03d.csv", seed))

  if (file.exists(cluster_file)) {
    existing <- read_csv(cluster_file, show_col_types = FALSE)
    # Drop any old mixAK rows (safe to re-run)
    existing <- existing %>% filter(cluster_method != "mixAK")
    combined <- bind_rows(existing, new_rows)
    cat(sprintf("\n  Merged with existing: %s\n", basename(cluster_file)))
  } else {
    warning(sprintf(
      "No existing cluster file for seed %d (%s). Saving mixAK-only results.",
      seed, cluster_file))
    combined <- new_rows
  }

  dir.create(DIR_CLUSTER, recursive = TRUE, showWarnings = FALSE)
  write_csv(combined, cluster_file)
  cat(sprintf("  Saved: %s\n", cluster_file))

  invisible(new_rows)
}


# =============================================================================
# Command-line dispatch
# =============================================================================

.parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- list(seed = NULL, seeds = NULL)
  i <- 1L
  while (i <= length(args)) {
    if      (args[i] == "--seed")  { out$seed  <- as.integer(args[i + 1L]); i <- i + 2L }
    else if (args[i] == "--seeds") { out$seeds <- eval(parse(text = args[i + 1L])); i <- i + 2L }
    else i <- i + 1L
  }
  if (is.null(out$seed) && is.null(out$seeds))
    stop("Provide --seed N  or  --seeds N:M")
  if (!is.null(out$seed))
    out$seeds <- out$seed
  out
}

if (!interactive()) {
  a <- .parse_args()
  for (s in a$seeds) run_mixAK_one_seed(s)
  cat("\nAll done.\n")
}
