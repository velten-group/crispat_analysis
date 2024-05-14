# crispat analysis repository

This repository contains the code used for the analysis shown in our paper ["Guide assignment in single-cell CRISPR screens using crispat"](https://www.biorxiv.org/content/10.1101/2024.05.06.592692v1). The code, installation instructions and documentation for our package **crispat**, can be found [here](https://github.com/velten-group/crispat).

Pooled single cell CRISPR screens have emerged as a powerful tool in functional genomics to probe the effect of genetic interventions at scale. A crucial step in the analysis of the resulting data is the assignment of cells to gRNAs corresponding to a specific genetic intervention. However, this step is challenging due to a lack of systematic benchmarks and accessible software to apply and compare different guide assignment strategies. To address this, we propose crispat (CRISPR guide assignment tool), a Python package to facilitate the choice of a suitable guide assignment strategy for single cell CRISPR screens.

In this repository, we provide the code we used to demonstrate the package on four single cell CRISPR interference screens from two studies, where crispat identifies strong differences in the number of assigned cells, downregulation of the target genes and number of discoveries across different guide assignment strategies, highlighting the need for a suitable guide assignment strategy to obtain optimal power in single cell CRISPR screens. 
