# Longitudinal Clustering Pipeline

## Project Structure

```
simulation/       Simulation study pipeline
realworld/        Real-world data analysis
```

---

## Simulation Pipeline

### Module Scripts (core logic — do not run directly)

| File | Role |
|------|------|
| `00_utils.R` | Shared helpers: `make_flat_matrix()`, `wide_to_long()`, `compute_ari()`, `compute_chi()` |
| `01_impute.R` | Five imputation methods (CC, MICE-L, LME-F, LME-S, MICE-CS); `make_config()` |
| `02_simulate.R` | eCDF-Copula resampling preserving within/between-visit correlation |
| `03_prepare.R` | SCE2 (composite scores) and SCE3 (PCA on item-level features) preparation |
| `04_cluster.R` | Five clustering methods: GBMT, KML3D, flexmix, LCMM, mixAK |
| `05_evaluate.R` | Aggregate metrics across seeds; select best-k; write summary tables |

---

### Entry-Point Scripts

```
run_impute.R   Step 1 — imputation (run once, seed-independent)
run_main.R     Steps 2-4 — four phases: prepare → cluster → mixak → evaluate
```

---

### Step 1 — Imputation: `run_impute.R`

Run **once** before any simulation. Produces five imputed datasets.

```bash
Rscript run_impute.R
Rscript run_impute.R --data /path/to/data.csv --out output/imputed
```

| Argument | Default          | Description        |
|----------|------------------|--------------------|
| `--data` | `dat_sce2.csv`   | Path to input data |
| `--out`  | `output/imputed` | Output directory   |

**Output:** `output/imputed/impt_CC.csv`, `impt_MICE_L.csv`, `impt_LME_F.csv`, `impt_LME_S.csv`, `impt_MICE_CS.csv`

---

### Steps 2–4 — `run_main.R`

Four independent phases. Each reads from disk and writes to disk, so any phase can be re-run or submitted as a separate SLURM job.

```bash
Rscript run_main.R --phase <phase> [--seed <int>] [--seeds <expr>]
```

| Argument  | Default   | Description |
|-----------|-----------|-------------|
| `--phase` | `prepare` | One of `prepare`, `cluster`, `mixak`, `evaluate`, `all` |
| `--seed`  | `1`       | Single integer seed (phases A–C) |
| `--seeds` | `1:50`    | R expression for seed vector (phase D / `all`) |

#### Phase A — `prepare`

Simulate data and prepare SCE2 / SCE3 feature matrices. Must run before B or C.

```bash
Rscript run_main.R --phase prepare --seed 1
```

Reads: `output/imputed/`, `dat_sce2.csv`
Writes: `output/sim/sim_seed001_*.csv`, `output/prepared/sce2_seed001_*.csv`, `output/prepared/sce3_seed001_*.csv`, `output/prepared/sim_methods_seed001.rds`

#### Phase B — `cluster`

Non-mixAK clustering (GBMT, KML3D, flexmix, LCMM) on prepared SCE2 and SCE3 data.

```bash
Rscript run_main.R --phase cluster --seed 1
```

Reads: `output/prepared/`
Writes: `output/cluster/cluster_seed001.csv`

#### Phase C — `mixak`

mixAK clustering only. Kept separate due to long MCMC runtime.

```bash
Rscript run_main.R --phase mixak --seed 1
```

Reads: `output/prepared/`
Writes: `output/cluster/mixak_seed001.csv`

#### Phase D — `evaluate`

Merge `cluster_seed*.csv` + `mixak_seed*.csv` (if present) across all seeds, compute summaries.

```bash
Rscript run_main.R --phase evaluate
Rscript run_main.R --phase evaluate --seeds "1:50"
```

Writes: `output/results/`

#### Local sequential test — `all`

Runs A → B → C for each seed in order, then evaluates. Useful for local debugging.

```bash
Rscript run_main.R --phase all --seeds "1:5"
```

---

### SLURM Usage

Submit phases A, B, C as separate array jobs. Phase D runs once after all seeds finish.

```bash
# Phase A — prepare (submit first)
#SBATCH --array=1-50
Rscript run_main.R --phase prepare --seed $SLURM_ARRAY_TASK_ID

# Phase B — cluster (submit after A finishes)
#SBATCH --array=1-50
Rscript run_main.R --phase cluster --seed $SLURM_ARRAY_TASK_ID

# Phase C — mixAK (submit after A finishes; independent of B)
#SBATCH --array=1-50
Rscript run_main.R --phase mixak --seed $SLURM_ARRAY_TASK_ID

# Phase D — evaluate (submit after B and C finish)
Rscript run_main.R --phase evaluate --seeds "1:50"
```

---

### `--seeds` Format

Accepts any valid R expression returning an integer vector:

```
--seeds "1:50"
--seeds "c(1, 5, 10, 20)"
--seeds "seq(1, 100, by=2)"
```

---

### Output Directory Layout

```
output/
  imputed/       impt_CC.csv, impt_MICE_L.csv, ...       (Step 1)
  sim/           sim_seed001_CC.csv, ...                  (Phase A: raw simulated data)
  prepared/      sce2_seed001_CC.csv, ...                 (Phase A: clustering-ready features)
                 sce3_seed001_CC.csv, ...
                 sim_methods_seed001.rds
  cluster/       cluster_seed001.csv, ...                 (Phase B: non-mixAK metrics)
                 mixak_seed001.csv, ...                   (Phase C: mixAK metrics)
  results/       all_metrics_by_seed.csv                  (Phase D)
                 summary_by_method_k.csv
                 best_k_by_chi.csv
                 ari_at_best_k.csv
  session_info.txt
```

---

### Notes

- Phase A must complete before B or C can run (they read `output/prepared/`).
- Phases B and C are independent of each other and can run in parallel.
- Phase D automatically merges B and C results; running only B without C (or vice versa) is safe.
- The `if (!interactive())` guard skips the CLI dispatch block in RStudio.
- `DIR_BASE` (`output/`) is created automatically on first run.
- `session_info.txt` records R and package versions for reproducibility.

---

## Real-World Analysis

```bash
Rscript realworld/main.R
```

Or via SLURM:

```bash
sbatch realworld/submit.sh
```

Runs flexmix, GBMT, LCMM, and mixAK on the real-world dataset. Results saved as `Sep22_realworld_res.RDS`.
