#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time=24:00:00
#SBATCH -c 60
#SBATCH --mail-user=h6cai@ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH -o realworld-%j

cd /ocean/projects/med220007p/hcai5/jinyuan/out/realworld

Rscript /ocean/projects/med220007p/hcai5/jinyuan/code/realworld_0920.R
