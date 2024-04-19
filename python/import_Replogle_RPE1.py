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

def load_data(batch_list, input_dir):
    '''
    Loads the specified batches from the data set from Reploge et al.
    Args:
        batch_list: (list) list of batch numbers to be included
        input_dir: (str) directory containing the data
    Returns:
        An AnnData object with the gene expression and CRISPR UMI counts
    '''
    batches = []
    for batch_number in tqdm(batch_list):
        time.sleep(1)
        # Load data of one batch 
        batch = sc.read_10x_mtx(input_dir + 'Replogle_RPE1/batch' + str(batch_number)+'/', 
                                var_names='gene_symbols',
                                gex_only = False, 
                                prefix = '') 
        batch.obs_names = [name.split('-')[0] for name in batch.obs_names]
        batch.obs['batch'] = batch_number
        batches.append(batch)
    # Concatenate batches into one AnnData object
    adata = ad.concat(batches, merge = "same", label = "batch", keys = batch_list, index_unique="-")
    return adata


def get_anno(data_dir):
    # create annotation df from the Supplementary Information file of Replogle et al. and deal with the duplicated gRNAs
    anno = pd.read_csv(data_dir + 'RPE1_library_annotation.csv', sep = ';')
    anno = pd.concat([pd.DataFrame({'gRNA': anno['sgID_A'],
                             'gRNA_pair_id': anno['unique sgRNA pair ID'], 
                             'target_gene': anno['gene'], 
                             'ensembl_id': anno['ensembl gene id'],
                             'target_seq': anno['targeting sequence A']}),
                      pd.DataFrame({'gRNA': anno['sgID_B'],
                             'gRNA_pair_id': anno['unique sgRNA pair ID'], 
                             'target_gene': anno['gene'], 
                             'ensembl_id': anno['ensembl gene id'],
                             'target_seq': anno['targeting sequence B']})], ignore_index = True)

    # for duplicated gRNAs (same target sequence) we only keep one (the first) of multiple target gene names
    gene_seq = anno[['target_seq', 'target_gene', 'ensembl_id']].groupby('target_seq').agg(['first'])
    anno['target_gene'] = anno.apply(lambda row: gene_seq.loc[row['target_seq'], 'target_gene'], axis=1)['first']
    anno['ensembl_id'] = anno.apply(lambda row: gene_seq.loc[row['target_seq'], 'ensembl_id'], axis=1)['first']

    # dataset contains gRNAs that in the annotation file contain ',' but in the variable names '-'
    anno['gRNA'] = [gRNA.replace(',', '-') for gRNA in anno['gRNA']]
    return anno


def main():
    
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
    
    # If output_dir doesn't exist, the output folder is created
    output_dir = config['out_dir'] 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created")

    # Load AnnData object of selected batch 
    print('Load data')
    adata = load_data(config['batch_list'], config['input_dir'])

    # Subset to CRISPR gRNA features
    adata_crispr = adata[:, adata.var["feature_types"] == "CRISPR Guide Capture"]
    adata_crispr.var_names = [gRNA[:-5] for gRNA in adata_crispr.var_names]
    
    # Annotation for each gRNA containing the name of the target gene, its ensembl id and the target sequence
    anno = get_anno(config['input_dir'])

    # Sum the gRNA counts over the pairs
    anno_f = anno[anno['gRNA'].isin(adata_crispr.var_names)]
    group_sums = {}
    for group, variables in anno_f.groupby("gRNA_pair_id")["gRNA"]:
        group_sums[group] = adata_crispr[:, variables].X.sum(axis=1)
    merged_data = csr_matrix(np.column_stack(list(group_sums.values())))
    adata_crispr = ad.AnnData(X=merged_data, obs=adata.obs, var=pd.DataFrame(index=list(group_sums.keys())))
    
    # Save as h5ad object
    adata_crispr.write(output_dir + 'Replogle_RPE1_gRNA_counts.h5ad')


if __name__ == "__main__":
    main()
