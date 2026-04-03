# =============================================================================
# 05_evaluate.R
# Aggregation and evaluation across 50 simulation replications
#
# Workflow:
#   Each seed produces one run_output list (from run_one_seed() in run_all.R):
#     run_output$metrics  ← flat data.frame for that seed
#
#   This script collects those outputs and produces:
#     1. A long-format table of all metrics (all seeds × methods × k)
#     2. Summary statistics: mean / SD / SE per (scenario, sim_method,
#        cluster_method, k)
#     3. Best-k selection per (scenario, sim_method, cluster_method) via CHI
#     4. ARI at the best k, aggregated across seeds
#     5. CSV exports
#
# Main functions:
#   bind_seed_results()       combine per-seed metric tables into one
#   summarize_metrics()       mean / SD / SE of ARI, CHI, BIC, PED, Runtime
#   select_best_k()           pick k with highest mean CHI per method × sim
#   ari_at_best_k()           extract ARI rows for the CHI-selected k
#   write_results_tables()    save summary and best-k CSVs
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))


# =============================================================================
# 1. Collect per-seed outputs
# =============================================================================

#' Extract the flat metrics table from one seed's output
#'
#' @param seed_output  A single element of the list produced by run_one_seed()
#'                     Must contain $metrics (a data.frame with columns
#'                     scenario, sim_method, cluster_method, k, ARI, CHI,
#'                     BIC, PED, Runtime_seconds)
#' @param seed         Integer seed value (added as a column)
#' @return data.frame with an additional "seed" column
extract_seed_metrics <- function(seed_output, seed) {
  df       <- seed_output$metrics
  df$seed  <- as.integer(seed)
  df
}


#' Combine metric tables from all seeds into one long data.frame
#'
#' @param seed_results Named list: names are seed integers (or coercible to
#'                     integer), values are outputs of run_one_seed()
#' @return data.frame with columns:
#'   seed, scenario, sim_method, cluster_method, k,
#'   ARI, CHI, BIC, PED, Runtime_seconds
bind_seed_results <- function(seed_results) {
  rows <- mapply(
    extract_seed_metrics,
    seed_output = seed_results,
    seed        = as.integer(names(seed_results)),
    SIMPLIFY    = FALSE
  )
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  # Enforce column types
  df$seed           <- as.integer(df$seed)
  df$k              <- as.integer(df$k)
  df$ARI            <- as.numeric(df$ARI)
  df$CHI            <- as.numeric(df$CHI)
  df$BIC            <- as.numeric(df$BIC)
  df$PED            <- as.numeric(df$PED)
  df$Runtime_seconds <- as.numeric(df$Runtime_seconds)
  df
}


# =============================================================================
# 2. Summary statistics
# =============================================================================

#' Compute mean / SD / SE across seeds for every (scenario, sim_method,
#' cluster_method, k) combination
#'
#' @param metrics_df  Output of bind_seed_results()
#' @return data.frame with one row per (scenario, sim_method, cluster_method, k)
#'   and columns: n_seeds, ARI_mean, ARI_sd, ARI_se,
#'                CHI_mean, CHI_sd, CHI_se,
#'                BIC_mean, BIC_sd,
#'                PED_mean, PED_sd,
#'                Runtime_mean
summarize_metrics <- function(metrics_df) {
  metrics_df %>%
    group_by(scenario, sim_method, cluster_method, k) %>%
    summarise(
      n_seeds      = sum(!is.na(ARI)),
      ARI_mean     = mean(ARI,             na.rm = TRUE),
      ARI_sd       = sd(ARI,               na.rm = TRUE),
      ARI_se       = ARI_sd / sqrt(n_seeds),
      CHI_mean     = mean(CHI,             na.rm = TRUE),
      CHI_sd       = sd(CHI,               na.rm = TRUE),
      CHI_se       = CHI_sd / sqrt(n_seeds),
      BIC_mean     = mean(BIC,             na.rm = TRUE),
      BIC_sd       = sd(BIC,               na.rm = TRUE),
      PED_mean     = mean(PED,             na.rm = TRUE),
      PED_sd       = sd(PED,               na.rm = TRUE),
      Runtime_mean = mean(Runtime_seconds, na.rm = TRUE),
      .groups = "drop"
    )
}


# =============================================================================
# 3. Best-k selection via CHI
# =============================================================================

#' For each (scenario, sim_method, cluster_method), select the k with the
#' highest mean CHI across seeds
#'
#' @param summary_df  Output of summarize_metrics()
#' @return data.frame: same grouping columns + best_k, CHI_mean at best_k
select_best_k <- function(summary_df) {
  summary_df %>%
    group_by(scenario, sim_method, cluster_method) %>%
    slice_max(CHI_mean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(scenario, sim_method, cluster_method,
           best_k = k, CHI_mean_at_best_k = CHI_mean)
}


#' Extract per-seed ARI values at the CHI-selected best k
#'
#' @param metrics_df  Output of bind_seed_results()
#' @param best_k_df   Output of select_best_k()
#' @return Long data.frame with columns:
#'   seed, scenario, sim_method, cluster_method, best_k, ARI
ari_at_best_k <- function(metrics_df, best_k_df) {
  metrics_df %>%
    inner_join(
      best_k_df %>%
        select(scenario, sim_method, cluster_method, best_k),
      by = c("scenario", "sim_method", "cluster_method",
             "k" = "best_k")
    ) %>%
    rename(best_k = k) %>%
    select(seed, scenario, sim_method, cluster_method, best_k, ARI)
}


#' Summarize ARI at best k across seeds
#'
#' @param ari_bestk_df  Output of ari_at_best_k()
#' @return data.frame: (scenario, sim_method, cluster_method, best_k,
#'                      ARI_mean, ARI_sd, ARI_se)
summarize_ari_best_k <- function(ari_bestk_df) {
  ari_bestk_df %>%
    group_by(scenario, sim_method, cluster_method, best_k) %>%
    summarise(
      n_seeds  = sum(!is.na(ARI)),
      ARI_mean = mean(ARI, na.rm = TRUE),
      ARI_sd   = sd(ARI,   na.rm = TRUE),
      ARI_se   = ARI_sd / sqrt(n_seeds),
      .groups = "drop"
    )
}


# =============================================================================
# 4. Export
# =============================================================================

#' Write all result tables to CSV files
#'
#' @param metrics_df    Output of bind_seed_results()
#' @param summary_df    Output of summarize_metrics()
#' @param best_k_df     Output of select_best_k()
#' @param ari_summary   Output of summarize_ari_best_k()
#' @param dir_output    Directory to write files into
write_results_tables <- function(metrics_df, summary_df,
                                  best_k_df, ari_summary,
                                  dir_output = "output/results") {
  dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)

  write.csv(metrics_df,
            file.path(dir_output, "all_metrics_by_seed.csv"),
            row.names = FALSE)

  write.csv(summary_df,
            file.path(dir_output, "summary_by_method_k.csv"),
            row.names = FALSE)

  write.csv(best_k_df,
            file.path(dir_output, "best_k_by_chi.csv"),
            row.names = FALSE)

  write.csv(ari_summary,
            file.path(dir_output, "ari_at_best_k.csv"),
            row.names = FALSE)

  cat("Results written to:", dir_output, "\n")
  invisible(list(
    all_metrics  = metrics_df,
    summary      = summary_df,
    best_k       = best_k_df,
    ari_best_k   = ari_summary
  ))
}


# =============================================================================
# 5. One-call wrapper
# =============================================================================

#' Run the full evaluation pipeline from a named list of seed results
#'
#' @param seed_results  Named list (seed → run_one_seed output)
#' @param dir_output    Directory for CSV output
#' @return List with all four tables
run_evaluation <- function(seed_results, dir_output = "output/results") {
  cat("Binding", length(seed_results), "seed results...\n")
  metrics_df  <- bind_seed_results(seed_results)

  cat("Summarizing metrics...\n")
  summary_df  <- summarize_metrics(metrics_df)

  cat("Selecting best k by CHI...\n")
  best_k_df   <- select_best_k(summary_df)

  cat("Extracting ARI at best k...\n")
  ari_bestk   <- ari_at_best_k(metrics_df, best_k_df)
  ari_summary <- summarize_ari_best_k(ari_bestk)

  write_results_tables(metrics_df, summary_df, best_k_df,
                        ari_summary, dir_output)
}
