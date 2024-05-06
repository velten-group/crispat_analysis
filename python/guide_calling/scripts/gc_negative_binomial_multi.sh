#!/bin/bash

# List of start gRNAs
start_grna=("0" "200" "400" "600" "800" "1000" "1200" "1400" "1600" "1800" "2000" "2200")

# Loop through the scripts and submit each one as a job with sbatch
for param in "${start_grna[@]}"; do
    sbatch guide_calling/scripts/gc_negative_binomial.sh "${param}"
done

