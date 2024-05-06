#!/bin/bash

# List of .sh scripts
scripts=("scripts/run_power_check_2Beta.sh" "scripts/run_power_check_3Beta.sh" "scripts/run_power_check_SCEPTRE.sh" "scripts/run_power_check_cellranger.sh" "scripts/run_power_check_negative_binomial.sh" "scripts/run_power_check_replogle.sh" "scripts/run_power_check_binomial.sh" "scripts/run_power_check_max.sh")

# Loop through the scripts and submit each one as a job with sbatch
for script in "${scripts[@]}"; do
    sbatch "$script"
    sleep 2
done
