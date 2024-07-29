import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from itertools import chain
from scipy import sparse, stats, special
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate
from torch.distributions.mixture_same_family import MixtureSameFamily
    
    
def get_ratios(adata_crispr):
    '''
    Calculates the proportion of guide counts per cell
    Args:
        adata_crispr: (AnnData) CRISPR gRNA counts per cell
    Returns:
        A long dataframe with the proportion of guide counts per cell
    '''
    # Normalization to get the percentages of the gRNA counts in each cell
    sc.pp.normalize_total(adata_crispr, target_sum = 1)
    
    # Convert adata object to long df
    percent_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    percent_df['cell'] = percent_df.index
    percent_df = pd.melt(percent_df, id_vars=['cell'], var_name='gRNA', value_name='percent_counts')

    # Add batch column
    percent_df['batch'] = [adata_crispr.obs.loc[cell, 'batch'] for cell in percent_df['cell']]

    return percent_df


def main():
    t_start = time.time()
    
    # Parse command line arguments
    argv = sys.argv[1:]
    try:
        options, args = getopt.getopt(argv, "c:", ["config="])
    except:
        print('Incorrect arguments!')
    for name, value in options:
        if name in ('-c', '--config'):
            config_filename = value
    
    # Parse config file
    config = json.load(open(config_filename))
    
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(config['out_dir']):
        os.makedirs(config['out_dir'])
        print("The output directory " + config['out_dir'] +  " was created")
    
    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    # subset to specified batches
    adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(config['batch_list'])]
    
    # Calculate ratios per cell
    data = get_ratios(adata_crispr)
        
    # Get the cells with highest ratios per gRNA
    print('Get cells with highest ratios per gRNA')
    data = data[data['percent_counts'] != 0].sort_values(by=['gRNA', 'percent_counts'], ascending=[True, False])
    
    for quantile in config['quantiles']:
        perturbations = data.groupby('gRNA').apply(lambda x: x.head(int(quantile * len(x))))
        # Save the resulting dataframe
        perturbations.to_csv(config['out_dir'] + 'perturbations_q' + str(quantile) + '.csv', index = False)
        
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['quantiles'], 'time': [t_end - t_start]})
    run_time.to_csv(output_dir + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + config['out_dir'])
    
    
if __name__ == "__main__":
    main()
