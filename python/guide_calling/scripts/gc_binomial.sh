#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=16
#SBATCH --time=20:00:00
#SBATCH --mem=100gb
#SBATCH --export=NONE
export OMP_NUM_THREADS=${SLURM_NTASKS}

module load devel/miniconda/3
eval "$(conda shell.bash hook)"
conda activate sc

start_grna="$1"
echo "Guide calling binomial with start gRNA $start_grna"
python guide_calling/guide_calling_binomial.py --config guide_calling/scripts/config_binomial.json --start_gRNA $start_grna
