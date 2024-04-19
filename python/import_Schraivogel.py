import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

def main():
    
    # Parse command line argument
    data_dir = sys.argv[1]
    
    # Load csv files with the gRNA counts
    print('Load data')
    tap = pd.read_csv(data_dir + 'Schraivogel_TAP_gRNA_counts.csv').set_index('gRNA')
    whole = pd.read_csv(data_dir + 'Schraivogel_whole_gRNA_counts.csv').set_index('gRNA')

    # Create anndata object
    print('Create anndata objects')
    tap = ad.AnnData(tap.transpose())
    tap.obs['batch'] = 1
    tap.X = csr_matrix(tap.X)
    whole = ad.AnnData(whole.transpose())
    whole.obs['batch'] = 1
    whole.X = csr_matrix(whole.X)
    
    # Save as h5ad objects
    tap.write(data_dir + 'Schraivogel_TAP_gRNA_counts.h5ad')
    whole.write(data_dir + 'Schraivogel_whole_gRNA_counts.h5ad')
    print('Done: h5ad objects for the Schraivogel data are saved in ' + data_dir)

if __name__ == "__main__":
    main()
