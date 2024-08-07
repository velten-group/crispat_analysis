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
    

def calculate_jaccard(pert_dict):
    '''
    Calculates pairwise Jaccard index
    
    Args:
        pert_dict (dictionary): dictionary which contains assigned cell_gRNAs for each method
        
    Returns:
        A pd DataFrame with the pairwise similarity scores
    '''
    matrix = np.zeros((len(pert_dict), len(pert_dict)))
    results = pd.DataFrame(matrix, index=pert_dict.keys(), columns=pert_dict.keys())

    for i in list(pert_dict.keys()):
        for j in list(pert_dict.keys()):
            set1 = set(pert_dict[i])
            set2 = set(pert_dict[j])    
            intersection_size = len(set1.intersection(set2))
            union_size = len(set1.union(set2))
            similarity = intersection_size / union_size
            results.loc[i, j] = similarity
    return results


def plot_intersection_heatmap(perturbations, method_order = None):
    '''
    Plot a heatmap with Jaccard index showing the intersecting assignments
    
    Args:
        perturbations (pd DataFrame): df with the assigned perturbations (needed columns: method, cell, gRNA)
        method_order (list): list defining the order of the rows and columns (default: alphabetic order)
        
    Returns:
        A matplotlib plot 
    '''
    # create dictionary with the cell-perturbation pairs per method
    pert_dict = {}
    
    if method_order != None:
        methods = method_order
    else:
        methods = perturbations['method'].unique()
        
    for method in methods:
        subset = perturbations[perturbations['method'] == method]
        cell_gRNA_list = subset['cell'] + '_' + subset['gRNA'] 
        pert_dict[method] = cell_gRNA_list.unique()
    
    # calculate matrix with intersection / total in row
    matrix = calculate_jaccard(pert_dict)
    intersections = pd.DataFrame(matrix, columns = methods, index = methods)
    
    # plot heatmap
    fig, axs = plt.subplots(figsize = (6, 3), nrows = 1, ncols = 1)
    cax = axs.matshow(intersections, cmap='YlOrRd')
    axs.set_xticks(range(len(methods)))
    axs.set_yticks(range(len(methods)))
    axs.set_xticklabels(methods, rotation=45, ha='left')
    axs.set_yticklabels(methods)
    plt.colorbar(cax, label='Jaccard index')
      
    
def plot_n_assigned_cells(perturbations, colors = None):
    '''
    Plots a barplot with the number of assigned cells and uniquely assigned cells per method
    
    Args:
        perturbations (pd DataFrame): df with the assigned perturbations (needed columns: method, cell, gRNA)
        colors (dictionary, optional): specifies the colors to use for each method in the barplot 
        
    Returns:
        A matplotlib plot 
    '''
    perturbations = perturbations[['cell', 'method', 'gRNA']].drop_duplicates()
    
    # calculate number of assignments per method
    n_total_assignments = perturbations[['method','cell']].drop_duplicates().groupby(['method']).size().reset_index(name='count')
    n_total_assignments['group'] = 'All assigned cells'

    # calculate number of cells with single gRNA assigned
    grouped_df = perturbations.groupby(['method', 'cell']).size().reset_index(name='grna_count')
    filtered_df = perturbations.merge(grouped_df, on=['method', 'cell'])
    filtered_df = filtered_df[filtered_df['grna_count'] == 1]
    n_single_assignments = filtered_df.groupby(['method']).size().reset_index(name='count')
    n_single_assignments['group'] = 'Uniquely assigned cells'

    combined_df = pd.concat([n_total_assignments, n_single_assignments])
    
    if colors == None:
        colors = sns.color_palette("husl", n_colors = n_total_assignments.shape[0])
    else:
        # order by guide assignment groups
        method_order = list(colors.keys())
        combined_df['method'] = pd.Categorical(combined_df['method'], categories=method_order, ordered=True)
        combined_df.sort_values(by='method', inplace=True)
    
    plt.figure(figsize = (3,3))
    scatter = sns.scatterplot(
        data = combined_df, x = 'count', y='method', 
        palette=colors, hue='method', s=50, style='group', 
        markers={'All assigned cells': 'X', 'Uniquely assigned cells': 'o'}
    )
    plt.grid(True, which='both', axis='y', color='grey', linestyle='-', linewidth=0.3)
    plt.grid(False, which='both', axis='x')

    # manually create legend
    handles, labels = scatter.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, label in enumerate(labels) if label in ['All assigned cells', 'Uniquely assigned cells']]
    new_labels = [label for label in labels if label in ['All assigned cells', 'Uniquely assigned cells']]
    scatter.legend(new_handles, new_labels, title='', loc = "upper center", bbox_to_anchor=(0.5, 1.25))

    plt.xlabel('Number of assigned cells')
    plt.ylabel('')
    plt.xlim(0, None)
    
    
def get_params_per_batch(batch_list, param, out_dir):
    res = pd.DataFrame({'gRNA': pd.read_csv(out_dir + 'PGMM_high-prior/batch1/gRNA_estimates.csv')['gRNA']})
    for batch in batch_list:
        batch_param = pd.read_csv(out_dir + 'PGMM_high-prior/batch' + str(batch) + '/gRNA_estimates.csv')
        batch_param = batch_param[['gRNA', param]]
        batch_param.columns = ['gRNA', 'batch'+str(batch)]
        res = res.merge(batch_param, on = 'gRNA')
    res = res.set_index('gRNA')
    return res