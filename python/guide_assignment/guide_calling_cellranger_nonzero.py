import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from scipy import sparse, stats
from sklearn import mixture
import time


def call_presence_with_gmm_ab(umi_counts: np.ndarray, n_components: int = 2) -> np.ndarray:
    '''
    call_presence_with_gmm_ab function from cellranger github repo: https://github.com/10XGenomics/cellranger/blob/cfa9ac1a9a0fcfbd123dc574934a6a72889d2f70/lib/python/cellranger/feature/feature_assigner.py#L213C1-L235C1
    Given the UMI counts for a specific antibody, separate signal from background.
    '''
    if len(umi_counts) == 0:
        return np.repeat(False, len(umi_counts))
    
    if np.max(umi_counts) == 0 or max(umi_counts.shape) < 2:
        # there are no UMIs, or only one UMI, each barcode has 0 count
        return np.repeat(False, len(umi_counts))

    # Turn the numpy array into 2D log scale array
    umi_counts = np.reshape(umi_counts, (len(umi_counts), 1))
    log_ab = np.log10(umi_counts + 1.0)

    # Initialize and fit the gaussian Mixture Model
    gmm = mixture.GaussianMixture(n_components, n_init=10, covariance_type="tied", random_state=0)
    gmm.fit(log_ab)

    # Calculate the threshold
    umi_posterior = gmm.predict_proba(log_ab)
    high_umi_component = np.argmax(gmm.means_)

    in_high_umi_component = umi_posterior[:, high_umi_component] > 0.5

    return in_high_umi_component


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
    
    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    
    for batch in config['batch_list']:
        print('Fit Gaussian Mixture model for batch ' + str(batch))
        # If output_dir doesn't exist, the output folders are created
        output_dir = config['out_dir'] + 'batch' + str(batch) + '/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print("The output directory " + output_dir +  " was created")

        # Subset to selected batch
        adata_crispr_batch = adata_crispr[adata_crispr.obs['batch'] == batch]
        gRNA_list = adata_crispr_batch.var_names.tolist()

        if config['gRNA_list'] != 'all':
            gRNA_list = gRNA_list[config['start_gRNA']:config['end_gRNA'] + 1]

        # Fit Gaussian Mixture Model (GMM) for each gRNA
        perturbations = pd.DataFrame({'cell': [], 'gRNA': []})
        
        for gRNA in gRNA_list:
            # Select data for one gRNA
            selected_guide = adata_crispr_batch[:,[gRNA]].X
            data = selected_guide.toarray() 
            data_nonzero = data[data != 0]
            
            # Fit gaussian mixture model
            perturbed = call_presence_with_gmm_ab(data_nonzero)
            indices_true = np.where(perturbed)[0]
            indices_full_data = np.where(data != 0)[0][indices_true]

            perturbed_cells = adata_crispr_batch.obs_names[indices_full_data].tolist()
            df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
            perturbations = pd.concat([perturbations, df], ignore_index = True)

        # Save data frame with the perturbations assigned to each cell
        perturbations.to_csv(output_dir + 'perturbations.csv', index = False)
    
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['cellranger'], 'time': [t_end - t_start]})
    run_time.to_csv(config['out_dir'] + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)
       
    
if __name__ == "__main__":
    main()