#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=16
#SBATCH --time=20:00:00
#SBATCH --mem=200gb
#SBATCH --export=NONE
export OMP_NUM_THREADS=${SLURM_NTASKS}

module load devel/miniconda/3
eval "$(conda shell.bash hook)"
conda activate sc

python guide_assignment/quantiles.py --config guide_assignment/scripts/config.json
