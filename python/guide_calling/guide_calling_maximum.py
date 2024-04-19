import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
from tqdm import tqdm
import time


def get_perturbed_cells(adata_crispr):
    '''
    Gets perturbed cells based on maximum assignment 
    Args:
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
    Returns:
        Dataframe containing for each cell the name of the gRNA with highest counts
    '''
   
    # Convert adata object to long df
    count_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    count_df['cell'] = count_df.index
    count_df = pd.melt(count_df, id_vars=['cell'], var_name='gRNA', value_name='UMI_counts')
    
    # Get maximum counts per cell
    count_df['max'] = count_df.groupby('cell')['UMI_counts'].transform('max')
    max_df = count_df[count_df['UMI_counts'] == count_df['max']]
    
    # Remove cells with maximum of 0
    max_df = max_df[max_df['UMI_counts'] != 0]                  
    # Remove cells where multiple gRNAs share the maximum
    occurences = max_df.groupby(['cell']).size().reset_index(name='n_max_gRNAs')
    max_df = max_df.merge(occurences, on = ['cell'])
    max_df = max_df[max_df['n_max_gRNAs'] == 1]
    
    return(max_df)


def main():
    t_start = time.time()
    print('Maximum assignment')
    
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
    output_dir = config['out_dir'] 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created") 

    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    # subset to specified batches
    adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(config['batch_list'])]
    
    gRNA_list = adata_crispr.var_names.tolist()
    
    if config['gRNA_list'] != 'all':
        gRNA_list = gRNA_list[config['start_gRNA']:config['end_gRNA'] + 1]

    # Get maximum assignment of each cell
    print('Get maximum assignment of each cell')
    perturbations = get_perturbed_cells(adata_crispr)
        
    # Save data frames with the results
    perturbations.to_csv(output_dir + 'perturbations.csv', index = False)
    
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['maximum'], 'time': [t_end - t_start]})
    run_time.to_csv(output_dir + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    main()
