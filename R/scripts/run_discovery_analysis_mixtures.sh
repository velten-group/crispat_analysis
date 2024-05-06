#!/bin/bash

# List of .sh scripts
scripts=("scripts/run_discovery_analysis_2Beta.sh" "scripts/run_discovery_analysis_3Beta.sh" "scripts/run_discovery_analysis_SCEPTRE.sh" "scripts/run_discovery_analysis_cellranger.sh" "scripts/run_discovery_analysis_negative_binomial.sh" "scripts/run_discovery_analysis_replogle.sh" "scripts/run_discovery_analysis_binomial.sh" "scripts/run_discovery_analysis_max.sh")

# Loop through the scripts and submit each one as a job with sbatch
for script in "${scripts[@]}"; do
    sbatch "$script"
    sleep 2
done
