#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=16
#SBATCH --time=10:00:00
#SBATCH --mem=200gb
#SBATCH --export=NONE
export OMP_NUM_THREADS=${SLURM_NTASKS}

module load devel/miniconda/3
eval "$(conda shell.bash hook)"
conda activate sc

python guide_calling/guide_calling_beta_2.py --config guide_calling/scripts/config_2beta.json
