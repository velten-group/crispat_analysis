Analyses in R
----------------

This directory contains the R code to reproduce the analyses for K562 cells of the Replogle CRISPR interference screen. Results for the other data sets can be obtained by adjusting the file paths accordingly. This directory includes the following files:
- requirements.txt: List of all R packages and versions we used. 
- SCEPTRE_power_check.R: R script that runs the SCEPTRE pipeline until the power check step for a given assignment
- SCEPTRE_discovery.R:  R script that runs the full SCEPTRE pipeline on a subset of data (40 selected gRNAs and 5000 genes) for a given assignment
- scripts folder: contains the bash scripts that were used to run the SCEPTRE_power_check.R and SCEPTRE_discovery.R for all assignments. 
- utils_SCEPTRE.R: contains some helper functions to read in our assignments and create a SCEPTRE object which are used by the SCEPTRE_power_check.R and SCEPTRE_discovery.R scripts.
- varying_thresholds.Rmd/.html: This markdown contains our analysis for the methods in which a user-defined threshold has to be chosen. Based on these results, we chose the thresholds that maximise the total discoveries for the comparison with the other methods. 
- discovery_comparison.Rmd/.html: This markdown loads all results obtained from running the SCEPTRE pipeline for every assignment and creates the plots shown in our paper to compare the median log2 fold change, number of significant target genes, as well as the number of total and false discoveries across assignment methods. 