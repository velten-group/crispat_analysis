#!/bin/bash

# List of UMI thresholds for which to perform the power check 
# Note: guide calling has to be run previously for all used thresholds
thresholds=("2" "5" "10" "15" "20" "30" "40" "50" "60" "80" "100" "120" "140" "160" "180" "200")

# Loop through the scripts and submit each one as a job with sbatch
for param in "${thresholds[@]}"; do
    sbatch scripts/run_power_check_UMI_t.sh "${param}"
done
