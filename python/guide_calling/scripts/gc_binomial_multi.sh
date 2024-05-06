#!/bin/bash

# List of start gRNAs
start_grna=("0" "400" "800" "1200" "1600" "2000")

# Loop through the scripts and submit each one as a job with sbatch
for param in "${start_grna[@]}"; do
    sbatch guide_calling/scripts/gc_binomial.sh "${param}"
done

