#!/bin/bash
#SBATCH --job-name=sim_cluster_s10
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --array=10
#SBATCH --mail-user=h6cai@ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH -o /ocean/projects/med220007p/hcai5/cluster-eval/simulation/output/logs/phaseB_%a.out
#SBATCH -e /ocean/projects/med220007p/hcai5/cluster-eval/simulation/output/logs/phaseB_%a.err

module load gcc/13.3.1-p20240614
module load anaconda3/2024.10-1
conda activate bioc_env

SCRIPT_DIR=/ocean/projects/med220007p/hcai5/cluster-eval/simulation

mkdir -p ${SCRIPT_DIR}/output/logs

Rscript ${SCRIPT_DIR}/run_main.R --phase cluster --seed ${SLURM_ARRAY_TASK_ID}
