import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from torch.distributions import constraints
from scipy import sparse, stats, special
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time
from joblib import Parallel, delayed

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate
from pyro.distributions.mixture import MaskedMixture
from pyro.distributions.torch_distribution import TorchDistributionMixin 
    
    
@config_enumerate
def model(data, n_batches):
    '''
    Negative Binomial Model 
    '''
    # Global variables for each gRNA
    beta0 = pyro.sample("beta0", dist.Normal(0.0, 2.0))
    beta1 = pyro.sample("beta1", dist.LogNormal(1.0, 2.0))
    gamma = pyro.sample("gamma", dist.MultivariateNormal(torch.full((n_batches,), 0.0),
                                                         torch.eye(n_batches) * 1.0))
    pi = pyro.sample("pi", dist.Beta(1.0, 10.0))
    theta = pyro.sample("theta", dist.LogNormal(1.0,1.0))
    
    # Get data and confounders
    batch = torch.tensor(data.obs['batch'])
    seq_depth = torch.tensor(data.obs['total_counts'])
    obs = torch.tensor(data.X.toarray()).reshape(-1)

    with pyro.plate("cells", len(data)):
        # Local variables for each cell
        perturbation = pyro.sample("perturbation", dist.Bernoulli(pi)) 
        mu = torch.exp(beta0 + beta1 * perturbation + gamma[batch - 1] + torch.log(seq_depth + 1))
        sigma2 = mu + (mu ** 2) / theta
        r = (mu ** 2) / (sigma2 - mu)
        p = mu / sigma2
        pyro.sample("obs", dist.NegativeBinomial(r, 1-p), obs=obs)
        
        
def init_loc_fn(site, n_batches):
    '''
    Define initial parameter values
    '''
    if site["name"] == "beta0":
        return torch.tensor([0.0])
    if site["name"] == "beta1":
        return torch.tensor([3.0])
    if site["name"] == "gamma":
        return torch.zeros(n_batches)
    if site["name"] == "pi":
        return torch.tensor([0.01])
    if site["name"] == "theta":
        return torch.tensor([5.0])
    raise ValueError(site["name"])

    
def initialize(seed, optim, elbo, data, n_batches):
    '''
    Initialization for SVI 
    Args:
        seed: (str) seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data: (tensor) observed transformed gRNA counts
        conditioned_model: (pyro model) conditioned pyro mixture model
        n_batches: (int) number of batches
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["beta0", "beta1", "gamma", "pi", "theta"]),
        init_loc_fn=lambda site: init_loc_fn(site, n_batches)
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data, n_batches)


def plot_loss(losses, gRNA, output_dir):
    '''
    Saves a plot of the loss over the SVI steps
    Args:
        losses: (list) loss over the SVI steps
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    plt.figure(figsize=(8, 3), dpi=300).set_facecolor("white")
    plt.plot(losses)
    plt.xlabel("iters")
    plt.ylabel("loss")
    plt.yscale("log")
    plt.title("Convergence of SVI for " + gRNA)
    plt.savefig(output_dir+"loss_plots/loss_"+gRNA+".png", bbox_inches="tight")
    plt.close()

    

def fit_NegBinomial(gRNA, adata_crispr, batch_list, output_dir, seed):
    '''
    Fits Negative Binomial model for UMI counts of one gRNA 
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
        batch_list: (list) list of all batches
        output_dir: (str) directory in which the resulting plots will be saved
        seed: (int) seed used for pyro
    Returns:
        List of cells perturbed with the specified gRNA, as well as the inferred threshold
    '''
    # Set optimizer and elbo parameters
    optim = pyro.optim.Adam({"lr": 0.01, "betas": [0.8, 0.99]})
    elbo = TraceEnum_ELBO(num_particles = 1, max_plate_nesting=1)
    n_steps = 2500
    
    # Data used to fit the model
    selected_gRNA = adata_crispr[:,[gRNA]]

    # Only fit model for gRNAs with non-zero counts in at least 3 cells and with a maximum count of at least 2
    data = selected_gRNA.X.toarray() 
    data = torch.tensor(data[data != 0])
    if len(data) < 2:   
        print(gRNA + " has only " + str(len(data)) + " cells with non-zero counts, so no model is fitted for that gRNA")
        return(0, 0)
    if max(data) < 2:
        print("Max UMI count for " + gRNA + " is " + str(max(data).item()) + ", so no model is fitted for that gRNA")
        return(0, 0)

    # Initialization of SVI
    n_batches = len(batch_list)
    initialize(seed, optim, elbo, selected_gRNA, n_batches)

    # Train the model n_steps steps with early stopping when loss doesn't change at least 0.001 for 50 steps
    losses = []
    min_loss = 1.e10
    last_step = 0
    subsample_size = 15000

    for step in range(n_steps):
        subsample = sc.pp.subsample(selected_gRNA, n_obs=subsample_size, random_state=2023+step, copy = True)
        loss = svi.step(subsample, n_batches)
        if loss < min_loss - 0.001:
            min_loss = loss
            last_step = step
        losses.append(loss)
        #if (step - last_step) > 50:
        #    break
        
    # MAP estimates of the model
    map_estimates = global_guide(data)
    batch_effects = pd.DataFrame({'gRNA': gRNA,
                                  'param': ['gamma_' + str(batch) for batch in batch_list],
                                  'value': map_estimates['gamma'].detach().numpy()})
    beta0 = map_estimates["beta0"].item()
    beta1 = map_estimates["beta1"].item()
    pi = map_estimates["pi"].item()
    theta = map_estimates["theta"].item()
    estimates = pd.concat([pd.DataFrame({'gRNA': gRNA, 
                                         'param': ['beta0', 'beta1', 'pi', 'theta'], 
                                         'value': [map_estimates["beta0"].item(), map_estimates["beta1"].item(),
                                                   map_estimates["pi"].item(), map_estimates["theta"].item()]}), 
                           batch_effects])

    # create plot of the loss
    plot_loss(losses, gRNA, output_dir)
    
    return(losses[-1], estimates)


def get_perturbed_cells(adata_crispr, estimates, gRNA):
    '''
    Gets perturbed cells for one gRNA 
    Args:
        adata_crispr: (AnnData object) UMI counts of CRISPR gRNAs in all cells
        estimates: (pd DataFrame) parameter estimates from the Negative Binomial model
        gRNA: (str) gRNA name
    Returns:
        Perturbed cells
    '''
    # Get data and confounders
    selected_gRNA = adata_crispr[:,[gRNA]]
    data = pd.DataFrame({'gRNA_counts': selected_gRNA.X.toarray().reshape(-1), 
                         'batch': selected_gRNA.obs['batch'], 
                         'seq_depth': selected_gRNA.obs['total_counts']})
    data = data[data['gRNA_counts'] != 0]
    
    # get inferred parameters
    beta0 = estimates.loc[estimates['param'] == 'beta0', 'value'].item()
    beta1 = estimates.loc[estimates['param'] == 'beta1', 'value'].item()
    gamma = list(estimates.loc[estimates['param'].str.startswith('gamma_'), 'value'])
    gamma_cells = np.array([gamma[batch - 1] for batch in data['batch']])
    pi = estimates.loc[estimates['param'] == 'pi', 'value'].item()
    theta = estimates.loc[estimates['param'] == 'theta', 'value'].item()
    
    # calculate mean per cell
    data['mu0'] = np.exp(beta0 + gamma_cells + np.log(data['seq_depth']))
    data['mu1'] = np.exp(beta0 + beta1 + gamma_cells + np.log(data['seq_depth']))
    
    # calculate sigma2 per cell
    data['sig0'] = data['mu0'] + (data['mu0']**2) / theta
    data['sig1'] = data['mu1'] + (data['mu1']**2) / theta
    
    # get probabilities for the two mixture components
    data['p0'] = nbinom.pmf(data['gRNA_counts'], (data['mu0']**2) / (data['sig0'] - data['mu0']), 
                            data['mu0']/data['sig0']) * (1 - pi)
    data['p1'] = nbinom.pmf(data['gRNA_counts'], (data['mu1']**2) / (data['sig1'] - data['mu1']), 
                            data['mu1']/data['sig1']) * pi
    data['pert_probability'] = data['p1'] / (data['p0'] + data['p1'])
    data['perturbation'] =  np.where(data['pert_probability'] >= 0.8, 1, 0)
    
    # filter for perturbed cells
    perturbed_cells = data[data['perturbation'] == 1]
    perturbed_cells['gRNA'] = gRNA
    perturbed_cells.index.name = 'cell'
    perturbed_cells = perturbed_cells.reset_index()
    
    return perturbed_cells


def main():
    t_start = time.time()
    
    # Parse command line arguments
    argv = sys.argv[1:]
    try:
        options, args = getopt.getopt(argv, "c:", ["config=", "start_gRNA="])
    except:
        print('Incorrect arguments!')
        
    for name, value in options:
        if name in ('-c', '--config'):
            config_filename = value
        if name in ('--start_gRNA'):
            start_gRNA = int(value)
    
    # Parse config file
    config = json.load(open(config_filename))
    end_gRNA = start_gRNA + config['step'] - 1
    
    # If output_dir doesn't exist, the output folders are created
    output_dir = config['out_dir'] 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        os.makedirs(output_dir + "loss_plots/")
        print("The output directory " + output_dir +  " was created")

    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    # subset to specified batches
    adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(config['batch_list'])]

    # Add total number of gRNA counts per cell 
    adata_crispr.obs['total_counts'] = np.array(adata_crispr.X.sum(axis = 1)).flatten()

    gRNA_list = adata_crispr.var_names.tolist()
    
    if config['gRNA_list'] != 'all':
        if end_gRNA >= len(gRNA_list):
            end_gRNA = len(gRNA_list) - 1
        gRNA_list = gRNA_list[start_gRNA:(end_gRNA + 1)]
        
    # Fit Negative Binomial Model for each gRNA
    perturbations = pd.DataFrame({})
    losses = pd.DataFrame({'gRNA': [], 'loss': []})
    estimates = pd.DataFrame()
    print('Fit Negative Binomial Model for each gRNA')
    for gRNA in gRNA_list:
        #time.sleep(0.01)
        # fit the model to infer the parameters
        loss, map_estimates = fit_NegBinomial(gRNA, adata_crispr, config['batch_list'], output_dir, config['seed'])
        losses = pd.concat([losses, pd.DataFrame({'gRNA': [gRNA], 'loss': [loss]})])
        estimates = pd.concat([estimates, map_estimates])
        
        # get the perturbed cells
        perturbed_cells = get_perturbed_cells(adata_crispr, map_estimates, gRNA)
        perturbations = pd.concat([perturbations, perturbed_cells])
        
    # Save data frames with the results
    if config['gRNA_list'] == 'all':
        perturbations.to_csv(output_dir + 'perturbations.csv', index = False)
        losses.to_csv(output_dir + 'gRNA_losses.csv', index = False)
        estimates.to_csv(output_dir + 'gRNA_estimates.csv', index = False)
    else:
        perturbations.to_csv(output_dir + 'perturbations_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
        losses.to_csv(output_dir + 'gRNA_losses_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
        estimates.to_csv(output_dir + 'gRNA_estimates_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
        
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['negative_binomial'], 'time': [t_end - t_start]})
    run_time.to_csv(output_dir + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)

if __name__ == "__main__":
    main()
