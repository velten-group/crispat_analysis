#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=16
#SBATCH --time=3:00:00
#SBATCH --mem=100gb
#SBATCH --export=NONE
export OMP_NUM_THREADS=${SLURM_NTASKS}

module load devel/miniconda/3
eval "$(conda shell.bash hook)"
conda activate sc

python guide_calling/guide_calling_UMI_t_pairs.py --config guide_calling/scripts/config_UMI.json
