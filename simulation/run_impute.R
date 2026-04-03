# =============================================================================
# run_impute.R
# Phase 1 — Imputation only (run ONCE, seed-independent)
#
# Usage:
#   Rscript run_impute.R
#   Rscript run_impute.R --data /path/to/data.csv --out output/imputed
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
for (f in c("00_utils.R", "01_impute.R")) source(file.path(script_dir, f))


# =============================================================================
# USER CONFIGURATION
# =============================================================================

DATA_PATH  <- file.path(script_dir, "dat_sce2.csv")
DIR_IMPUTED <- "output/imputed"

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
  dir_output       = DIR_IMPUTED
)


# =============================================================================
# Run
# =============================================================================

.parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- list(data = DATA_PATH, out = DIR_IMPUTED)
  i <- 1L
  while (i <= length(args)) {
    if      (args[i] == "--data") { out$data <- args[i + 1L]; i <- i + 2L }
    else if (args[i] == "--out")  { out$out  <- args[i + 1L]; i <- i + 2L }
    else i <- i + 1L
  }
  out
}

main <- function() {
  if (!interactive()) {
    a <- .parse_args()
    DATA_PATH         <<- a$data
    CFG$dir_output    <<- a$out
  }

  cat("=== Phase 1: Imputation ===\n")
  cat(sprintf("Input:  %s\n", DATA_PATH))
  cat(sprintf("Output: %s\n", CFG$dir_output))

  data <- read_csv(DATA_PATH, show_col_types = FALSE)

  run_all_imputations(
    data    = data,
    config  = CFG,
    methods = c("CC", "MICE-L", "LME-F", "LME-S", "MICE-CS"),
    seed    = 111,
    save    = TRUE
  )

  cat(sprintf("\nDone. Imputed files saved to: %s\n", CFG$dir_output))
}

main()
