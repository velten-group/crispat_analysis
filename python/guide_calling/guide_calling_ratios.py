import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from itertools import chain


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
    
    if not os.path.exists(config['out_dir']):
        os.makedirs(config['out_dir'])
        print("The output directory " + config['out_dir'] +  " was created")
    
    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    # subset to specified batches
    adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(config['batch_list'])]
    
    # Normalization to get the percentages of the gRNA counts in each cell
    sc.pp.normalize_total(adata_crispr, target_sum = 1)
    
    # Convert adata object to long df
    percent_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    percent_df['cell'] = percent_df.index
    percent_df = pd.melt(percent_df, id_vars=['cell'], var_name='gRNA', value_name='percent_counts')
    
    # Get maximum per cell
    max_df = percent_df.groupby('cell').agg({'percent_counts': max})
    max_df.to_csv(config['out_dir'] + 'max_df.csv')
    max_df = max_df.reset_index() 
    max_df = max_df.merge(percent_df, on = ['cell', 'percent_counts'])

    # Save the results
    for t in config['thresholds']: 
        perturbations = max_df[max_df['percent_counts'] > t]
        perturbations.to_csv(config['out_dir'] + 'perturbations_t' + str(t) + '.csv', index = False)
        
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['ratios'], 'time': [t_end - t_start]})
    run_time.to_csv(config['out_dir'] + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + config['out_dir'])
    
    
if __name__ == "__main__":
    main()