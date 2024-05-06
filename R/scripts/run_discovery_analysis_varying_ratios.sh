#!/bin/bash

# List of ratio thresholds for which to perform the discovery analysis 
# Note: guide calling has to be run previously for all used thresholds
thresholds=("0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9")

# Loop through the scripts and submit each one as a job with sbatch
for param in "${thresholds[@]}"; do
    sbatch scripts/run_discovery_analysis_ratio_t.sh "${param}"
done
