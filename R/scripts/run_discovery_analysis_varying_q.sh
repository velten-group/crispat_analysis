#!/bin/bash

# List of quantile thresholds for which to perform the discovery analysis 
# Note: guide calling has to be run previously for all used thresholds
thresholds=("0.01" "0.025" "0.05" "0.075" "0.1" "0.2" "0.3" "0.4" "0.5")

# Loop through the scripts and submit each one as a job with sbatch
for param in "${thresholds[@]}"; do
    sbatch scripts/run_discovery_analysis_q.sh "${param}"
done
