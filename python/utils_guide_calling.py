import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy
import torch
from upsetplot import from_contents, UpSet
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time

def load_data(batch_list, pre, input_dir):
    '''
    Loads the specified batches from the data set from Reploge et al.
    Args:
        batch_list: (list) list of batch numbers to be included
        pre: (str) prefix of the file names
        input_dir: (str) directory containing the data
    Returns:
        An AnnData object with the gene expression and CRISPR UMI counts
    '''
    batches = []
    for batch_number in tqdm(batch_list):
        # Load data of one batch
        if pre == 'RD7':
            prefix = pre+'_'+str(batch_number)+'_'
        else:
            prefix = pre+'_'+str(batch_number)+'_essential_'
            
        batch = sc.read_10x_mtx(input_dir, 
                                var_names='gene_symbols',
                                gex_only = False, 
                                prefix = prefix) 
        batch.obs_names = [name.split('-')[0] for name in batch.obs_names]
        batches.append(batch)
    # Concatenate batches into one AnnData object
    adata = ad.concat(batches, merge = "same", label = "batch", keys = batch_list, index_unique="-")
    return adata


def get_filtered_data(batch_list, pre, data_dir):
    data = load_data(batch_list, pre, data_dir)
        
    # Preprocessing based on the gene expression counts
    # Get gene expression data (without guide counts)
    ge_data = data[:, data.var["feature_types"] == "Gene Expression"]

    # Filter cells with too many mitochondrial genes
    ge_data.var['mt'] = ge_data.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(ge_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    if pre == 'RD7': 
        ge_data = ge_data[ge_data.obs.pct_counts_mt < 11, :]
    else: 
        ge_data = ge_data[ge_data.obs.pct_counts_mt < 20, :]
    
    # Filter cells with too few total counts
    ge_data = ge_data[ge_data.obs.total_counts > 3000, :]
    
    return data[ge_data.obs_names, :]
    
    
def get_anno(data_dir):
    # create annotation df from the Supplementary Information file of Replogle et al. and deal with the duplicated gRNAs
    anno = pd.read_csv(data_dir + 'K562_essential_library_annotation.csv', sep = ';')
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
    anno.loc[anno['target_gene'] == 'PMF1-BGLAP', 'target_gene'] = 'BGLAP'
    anno.loc[anno['target_gene'] == 'RPS10-NUDT3', 'target_gene'] = 'RPS10'
    anno.loc[anno['target_gene'] == 'NUP62', 'target_gene'] = 'ATF5'
    return anno


def get_thresholds(file_name, name):
    thresholds = pd.read_csv(file_name)
    thresholds['model'] = name
    thresholds['log2t'] = np.log2(thresholds['threshold'])
    return thresholds


def get_perturbations(file_name, anno, filtered_data, all_gRNAs, level):
    perturbations = pd.read_csv(file_name)[['cell', level]]

    if level == 'gRNA':
        perturbations['gRNA'] = [gRNA[:-5] for gRNA in perturbations['gRNA']]
    
    # filter for gRNAs occuring in Replogle et al. (especially to remove non-targeting gRNAs not passing the AD test)
    perturbations = perturbations[[gRNA in all_gRNAs for gRNA in perturbations[level]]]
    
    # filter for valid cells (enough total counts, not too many mitochondrial genes)
    perturbations = perturbations[[cell in filtered_data.obs_names for cell in perturbations['cell']]]
    
    # add annotation info (ensembl id of the target and target sequence)
    perturbations = pd.merge(perturbations, anno, how = 'inner', on = [level])
    
    return perturbations
    

def get_single_assignments(perturbations):
    grouped_df = perturbations[['method', 'cell', 'gRNA']].drop_duplicates()
    grouped_df = grouped_df.groupby(['method', 'cell']).size().reset_index(name='pert_count')
    single_perturbations = perturbations.merge(grouped_df, on=['method', 'cell'])
    single_perturbations = single_perturbations[single_perturbations['pert_count'] == 1]
    return single_perturbations


def get_perturbations_multi_batch(file_dir, batch_list,  anno, filtered_data, all_gRNAs):
    perturbations = pd.DataFrame()
    for batch in batch_list:
        b = get_perturbations(file_dir + "batch" + str(batch) + "/perturbations.csv", anno, filtered_data, all_gRNAs)
        perturbations = pd. concat([perturbations, b])
    return perturbations


def get_perturbations_wo_cell_filtering(file_name, anno, all_gRNAs):
    perturbations = pd.read_csv(file_name)[['cell', 'gRNA']]
    perturbations['gRNA'] = [gRNA[:-5] for gRNA in perturbations['gRNA']]
    
    # filter for gRNAs occuring in Replogle et al. (especially to remove non-targeting gRNAs not passing the AD test)
    perturbations = perturbations[[gRNA in all_gRNAs for gRNA in perturbations['gRNA']]]
    
    # add annotation info (ensembl id of the target and target sequence)
    perturbations = pd.merge(perturbations, anno, how = 'inner', on = ['gRNA'])
    
    # aggregate data per cell
    perturbations = perturbations.groupby('cell').agg({'gRNA':lambda x: x, 
                                                       'target_gene':lambda x: x.unique(),
                                                       'gRNA_pair_id': lambda x: x.unique(), 
                                                       'ensembl_id': lambda x: x.unique(), 
                                                       'target_seq': lambda x: x.unique()})
    
    # filter for cells with one perturbed target
    perturbations['n_targets']  = [len(targets) for targets in perturbations['target_gene']]
    perturbations = perturbations[perturbations['n_targets'] == 1]
    perturbations['target_gene'] = [gene[0] for gene in perturbations['target_gene']]
    return perturbations


def plot_intersections(perturbations):
    upset_dict = {}
    methods = perturbations['method'].unique()
    for method in methods:
        subset = perturbations[perturbations['method'] == method]
        cell_gRNA_list = subset['cell'] + '_' + subset['gRNA'] 
        upset_dict[method] = cell_gRNA_list.unique()
    
    intersections = from_contents(upset_dict)
    upset = UpSet(intersections, subset_size='count', show_counts = False,
                 sort_by="cardinality", min_subset_size = 100).plot()
    

def calculate_intersection_union(pert_dict):
    # calculate pairwise intersection / union of the sets
    matrix = np.zeros((len(pert_dict), len(pert_dict)))
    results = pd.DataFrame(matrix, index=pert_dict.keys(), columns=pert_dict.keys())

    for i in list(pert_dict.keys()):
        for j in list(pert_dict.keys()):
            set1 = set(pert_dict[i])
            set2 = set(pert_dict[j])    
            intersection_size = len(set1.intersection(set2))
            set1_size = len(set1)
            similarity = intersection_size / set1_size
            results.loc[i, j] = similarity
    return results


def plot_intersection_heatmap(perturbations, colors):
    # create dictionary with the cell-perturbation pairs per method
    pert_dict = {}
    methods = perturbations['method'].unique()
    for method in methods:
        subset = perturbations[perturbations['method'] == method]
        cell_gRNA_list = subset['cell'] + '_' + subset['gRNA'] 
        pert_dict[method] = cell_gRNA_list.unique()
    
    # calculate matrix with intersection / total in row
    matrix = calculate_intersection_union(pert_dict)
    intersections = pd.DataFrame(matrix, columns = methods, index = methods)
    
    # calculate total number of assignments per method
    n_assignments_per_method = perturbations[['cell', 'method', 
                                         'gRNA']].drop_duplicates().groupby(['method']).size().reset_index(name='count')
    n_assignments_per_method = n_assignments_per_method.sort_values(by='count')
    # sort the intersections according to the total number of assignments
    order = list(n_assignments_per_method['method'])
    intersections = intersections.loc[order, order]
    
    # plot heatmap
    fig, axs = plt.subplots(figsize = (9, 3), nrows = 1, ncols = 2, sharey = True, width_ratios = (0.6, 0.4))
    cax = axs[0].matshow(intersections, cmap='YlOrRd')
    axs[0].set_xticks(range(len(methods)))
    axs[0].set_yticks(range(len(methods)))
    axs[0].set_xticklabels(order, rotation=45, ha='left')
    axs[0].set_yticklabels(order)
    plt.colorbar(cax, label='Intersection / N Row')
    
    # plot barplot with set size on the right
    sns.barplot(data = n_assignments_per_method, y = "method", x = "count", ax = axs[1],
               palette = colors)
    axs[1].set_xlabel('')
    axs[1].set_ylabel('')
    axs[1].set_title('Number of cells\nwith single assignments')
    
    
def get_params_per_batch(batch_list, param, out_dir):
    res = pd.DataFrame({'gRNA': pd.read_csv(out_dir + 'PGMM_high-prior/batch1/gRNA_estimates.csv')['gRNA']})
    for batch in batch_list:
        batch_param = pd.read_csv(out_dir + 'PGMM_high-prior/batch' + str(batch) + '/gRNA_estimates.csv')
        batch_param = batch_param[['gRNA', param]]
        batch_param.columns = ['gRNA', 'batch'+str(batch)]
        res = res.merge(batch_param, on = 'gRNA')
    res = res.set_index('gRNA')
    return res