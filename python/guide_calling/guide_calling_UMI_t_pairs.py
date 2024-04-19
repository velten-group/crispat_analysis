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


def get_perturbed_cells(gRNA, adata_crispr, threshold):
    '''
    Gets perturbed cells for one gRNA based on a UMI threshold 
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
        threshold: (int) UMI threshold
    Returns:
        List of cells perturbed with the specified gRNA
    '''
    # Get cells with UMI counts higher than the threshold for specified gRNA
    selected_guide = adata_crispr[:,[gRNA]].X
    perturbed_cells = adata_crispr.obs_names[selected_guide.toarray().reshape(-1) >= threshold].tolist()
    return(perturbed_cells)


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

    # Get perturbed cells for each gRNA based on a fixed UMI threshold
    for threshold in config['thresholds']:
        perturbations = pd.DataFrame({'cell': [], 'gRNA': []})
        print('Get perturbed cells for each gRNA with UMI threshold = ' + str(threshold))
        for gRNA in gRNA_list:
            perturbed_cells = get_perturbed_cells(gRNA, adata_crispr, threshold)
            if len(perturbed_cells) != 0:
                df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
                perturbations = pd.concat([perturbations, df], ignore_index = True)

        # Save data frames with the results
        perturbations.to_csv(output_dir + 'perturbations_t' + str(threshold) + '.csv', index = False)
        
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['UMI_t'], 'time': [t_end - t_start]})
    run_time.to_csv(output_dir + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    main()
