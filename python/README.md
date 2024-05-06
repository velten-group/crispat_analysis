Analyses in python
----------------------

This directory contains all python code used for our analyses and consists of the following files:
- requirements.txt: List of used python packages and versions.
- import_Replogle_K562_essential.py/import_Replogle_RPE.py/import_Schraivogel.py: These scripts were used to create an AnnData object (.h5ad) of the gRNA count matrix which is then used as input for the guide assignment methods. Since we applied the assignment of the Replogle data set on the gRNA pair level, counts of individual gRNAs that belong to the same pair were summarized.
- guide_calling folder: This directory contains the python scripts of the 11 different guide assignment methods analysed in our paper. There are some small differences between this implementation and the one incorporated into our crispat package, which allows more flexibility through incorporating additional parameters. Additionally, we also provide the bash scripts used to run the guide assignment on the Replogle K562 cells as an example of how we ran all assignment methods on a high-performance cluster.
- ga_comparison_Replogle.ipynb/ga_comparison_Schraivogel.ipynb: These notebooks were used to load in the assignment results of all methods and create some initial visualizations such as the number of assigned cells and the overlap of the different assignment methods. 
- utils_guide_calling.py: This python script contains some helper functions which are used in the ga_comparison_Replogle.ipynb and ga_comparison_Schraivogel.ipynb scripts. 
- gRNA_selection.ipynb: For the discovery analyses performed in R, a subset of 40 targeting gRNAs was selected in this notebook.
- sum_Replogle_gRNAs.py: For an additional experiment, we grouped the gRNAs from Replogle et al. into 500 and 86 groups which is shown in this python script. 
