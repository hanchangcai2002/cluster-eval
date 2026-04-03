# Longitudinal Clustering Pipeline

## Project Structure

```
simulation/       Simulation study pipeline
realworld/        Real-world data analysis
```

---

## Simulation Pipeline

Three scripts run in sequence:

```
Step 1: run_impute.R   — imputation (run once)
Step 2: run_main.R     — simulation + clustering (one seed at a time, parallelisable)
Step 3: run_mixAK.R    — mixAK clustering (run after Step 2, one seed at a time)
```

---

### Step 1 — Imputation: `run_impute.R`

Run once. Seed-independent.

```bash
Rscript run_impute.R
Rscript run_impute.R --data /path/to/data.csv --out output/imputed
```

| Argument | Default               | Description              |
|----------|-----------------------|--------------------------|
| `--data` | `dat_sce2.csv`        | Path to input data       |
| `--out`  | `output/imputed`      | Output directory         |

---

### Step 2 — Simulation + Clustering: `run_main.R`

#### Syntax

```bash
Rscript run_main.R --phase <phase> [--seed <int>] [--seeds <expr>]
```

| Argument  | Default | Description                              |
|-----------|---------|------------------------------------------|
| `--phase` | `run`   | One of `run`, `evaluate`, `all`          |
| `--seed`  | `1`     | Single integer seed (for `run`)          |
| `--seeds` | `1:50`  | R expression evaluated via `eval(parse(...))` |

#### Phases

**run** — simulate + cluster for a single seed

```bash
Rscript run_main.R --phase run --seed 42
```

**evaluate** — aggregate results across seeds

```bash
Rscript run_main.R --phase evaluate --seeds "1:50"
```

**all** — run all seeds then evaluate (local sequential test)

```bash
Rscript run_main.R --phase all --seeds "c(1,2,3,4,5)"
```

#### SLURM Array Example

```bash
#!/bin/bash
#SBATCH --array=1-50

Rscript run_main.R --phase run --seed $SLURM_ARRAY_TASK_ID
```

Once all jobs finish:

```bash
Rscript run_main.R --phase evaluate --seeds "1:50"
```

---

### Step 3 — mixAK Clustering: `run_mixAK.R`

Run **after** Step 2 completes. `mixAK` is kept separate due to long MCMC runtime.

```bash
# Single seed
Rscript run_mixAK.R --seed 1

# Multiple seeds (sequential)
Rscript run_mixAK.R --seeds "1:10"
```

| Argument  | Description                                          |
|-----------|------------------------------------------------------|
| `--seed`  | Single integer seed                                  |
| `--seeds` | R expression returning an integer vector             |

Results are merged into the existing `output/cluster/cluster_seed{NNN}.csv`. Re-runs are safe — old mixAK rows are replaced.

---

### `--seeds` Format

Accepts any valid R expression that returns an integer vector:

```
--seeds "1:50"
--seeds "c(1, 5, 10, 20)"
--seeds "seq(1, 100, by=2)"
```

---

### Output Directory Layout

```
output/
  imputed/      impt_CC.csv, impt_MICE_L.csv, ...
  sim/          sim_seed001_CC.csv, ...    (per-seed simulated datasets)
  cluster/      cluster_seed001.csv, ...   (per-seed clustering metrics)
  results/      all_metrics_by_seed.csv, summary_by_method_k.csv, ...
  session_info.txt
```

---

### Notes

- The `if (!interactive())` guard means the CLI block is skipped in RStudio — scripts only dispatch via command line.
- `DIR_BASE` (`output/`) is created automatically on each run.
- `session_info.txt` records R version and package versions for reproducibility.

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
