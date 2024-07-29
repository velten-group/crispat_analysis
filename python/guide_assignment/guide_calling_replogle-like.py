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
from scipy.sparse import csr_matrix
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
from pyro.distributions.mixture import MaskedMixture
from pyro.distributions.torch_distribution import TorchDistributionMixin 
    
    
class cont_Poisson(dist.TorchDistribution):
    arg_constraints = {'rate': constraints.nonnegative}
    def __init__(self, rate, validate_args=None):
        self.rate = rate
        super().__init__(batch_shape=rate.shape, event_shape=torch.Size([]), 
                         validate_args=validate_args)

    def expand(self, batch_shape, _instance=None):
        new = self._get_checked_instance(cont_Poisson, _instance)
        new.rate = self.rate.expand(batch_shape)
        super(cont_Poisson, new).__init__(batch_shape=batch_shape, event_shape=torch.Size([]))
        return super(cont_Poisson, new).expand(batch_shape, _instance=new)

    def sample(self, sample_shape=torch.Size()):
        poisson = dist.Poisson(self.rate)
        return poisson.sample(sample_shape)
    
    def log_prob(self, value):
        prob = poisson_prob(value, self.rate)
        return torch.log(prob)
    

@config_enumerate
def model(data):
    '''
    Poisson-Gaussian Mixture Model 
    '''
    # Global variables
    weights = pyro.sample("weights", dist.Dirichlet(torch.tensor([0.9, 0.1])))
    loc = pyro.sample("mu", dist.Normal(3.0, 2.0))
    scale = pyro.sample("scale", dist.LogNormal(2.0, 1.0))
    lam = pyro.sample("lam", dist.LogNormal(0.0, 1.0))

    with pyro.plate("data", len(data)):
        # Local variables
        assignment = pyro.sample("assignment", dist.Categorical(weights)) 
        assignment = assignment == 1 #boolean mask needed for MaskedMixture

        poisson = cont_Poisson(lam)
        normal = dist.Normal(loc, scale)
        pyro.sample("obs", MaskedMixture(assignment, poisson, normal, validate_args = False), obs=data)
        
        
def init_loc_fn(site):
    '''
    Define initial parameter values
    '''
    if site["name"] == "weights":
        return torch.tensor([0.9, 0.1])
    if site["name"] == "scale":
        return torch.tensor([2.0])
    if site["name"] == "mu":
        return torch.tensor([3.0])
    if site["name"] == "lam":
        return torch.tensor([0.1])
    raise ValueError(site["name"])

    
def initialize(seed, optim, elbo, data):
    '''
    Initialization for SVI 
    Args:
        seed: (str) seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data: (tensor) observed transformed gRNA counts
        conditioned_model: (pyro model) conditioned pyro mixture model
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["weights", "mu", "scale", "lam"]),
        #init_loc_fn=init_loc_fn,
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data)


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
    
    
def plot_fitted_model(data, weights, mu, scale, lam, threshold, gRNA, output_dir):
    '''
    Saves a plot of the data histogram and the fitted mixture model
    Args:
        data: (tensor) observed transformed gRNA counts
        weights: (np array) estimated proportion for Poisson and Normal component
        mu: (float) MAP estimate for the mean of the Normal distribution
        scale: (float) MAP estimate for the scale of the Normal distribution
        lam: (float) MAP estimate for the lambda of the Poisson distribution
        threshold: (int) threshold for assigning it to the Normal component
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    X = np.arange(0, max(data)+2, 0.1)
    Y1 = weights[0] * torch.tensor([poisson_prob(k, torch.tensor(lam)) for k in X])
    Y2 = weights[1] * stats.norm.pdf(X, mu, scale)

    fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
    sns.histplot(data, binwidth=0.5, color='grey', stat = "proportion")
    ax.plot(X, Y1, "r-", label = "Poisson")
    ax.plot(X, Y2, "b-", label = "Normal")
    ax.plot(X, Y1 + Y2, "k--", label = "Mixture model")
    ax.axvline(threshold, c = "green", label = "Threshold")
    plt.legend()
    plt.title("Poisson-Gaussian mixture model for " + gRNA)
    plt.ylabel("Probability Density")
    plt.xlabel("Log2 " + gRNA + " UMI counts")
    plt.savefig(output_dir + "fitted_model_plots/fitted_model_" + gRNA + ".png", bbox_inches="tight")
    plt.close()
    
    
def poisson_prob(k, lam):
    '''
    Calculates the probability P(X = k) of a Poisson distribution with parameter lam
    Args:
        k: (float) value for which to calculate the probability
        lam: (tensor) lambda parameter of the Poisson distribution
    Returns:
        List of probabilities
    '''
    # gamma(k + 1) = k! is used to get the probability also for non-integer values
    prob = (lam ** k) * torch.exp(-lam) / special.gamma(k + 1)
    return prob
    
    
def prob_normal_component(X, weights, mu, scale, lam):
    '''
    Calculates the probability for belonging to the Gaussian component given observations
    Args:
        X: (list) list of values for which the probability is calculated
        weights: (np array) estimated proportion for Poisson and Normal component
        mu: (float) MAP estimate for the mean of the Normal distribution
        scale: (float) MAP estimate for the scale of the Normal distribution
        lam: (float) MAP estimate for the lambda of the Poisson distribution
    Returns:
        List of probabilities
    '''
    nominator = torch.tensor(stats.norm.pdf(X, mu, scale) * weights[1])
    denominator = nominator + torch.tensor([poisson_prob(k, torch.tensor(lam)) for k in X]) * weights[0]
    prob = nominator / denominator
    return prob
    

def fit_PGMM(gRNA, adata_crispr, output_dir, seed):
    '''
    Fits Poisson-Gaussian mixture model for log2 of non-zero UMI counts of one gRNA 
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
        output_dir: (str) directory in which the resulting plots will be saved
        seed: (int) seed used for pyro
    Returns:
        List of cells perturbed with the specified gRNA, as well as the inferred threshold
    '''
    # Set optimizer and elbo parameters
    optim = pyro.optim.Adam({"lr": 0.01, "betas": [0.8, 0.99]})
    elbo = TraceEnum_ELBO(num_particles = 1, max_plate_nesting=1)
    n_steps = 500
    
    # Data used to fit the model: log2 transformation of non-zero UMI counts for a given gRNA 
    # (rounded to integers due to discrete Poisson model)
    selected_guide = adata_crispr[:,[gRNA]].X
    data = selected_guide.toarray() 
    data = torch.tensor(np.log2(data[data != 0])).float()

    # Only fit model for gRNAs with non-zero counts in at least 2 cells and with a maximum count of at least 2
    if len(data) < 2:   
        print(gRNA + " has only " + str(len(data)) + " cells with non-zero counts, so no model is fitted for that gRNA")
        return([], 0, 0, 0)
    if max(data) < np.log2(2):
        print("Max UMI count for " + gRNA + " is " + str(max(data).item()) + ", so no model is fitted for that gRNA")
        return([], 0, 0, 0)
    
    # Choose the best among 10 random initializations.
    loss, seed = min((initialize(seed, optim, elbo, data), seed) for seed in range(10))

    # Initialization of SVI
    initialize(seed, optim, elbo, data)

    # Train the model n_steps steps with early stopping when loss doesn't change at least 0.001 for 50 steps
    losses = []
    min_loss = 1.e8
    last_step = 0

    for step in range(n_steps):
        loss = svi.step(data)
        if loss < min_loss - 0.001:
            min_loss = loss
            last_step = step
        losses.append(loss)
        if (step - last_step) > 50:
            break
        
    # MAP estimates of the model
    map_estimates = global_guide(data)
    weights = map_estimates["weights"].data.numpy()
    mu = map_estimates["mu"].item()
    scale = map_estimates["scale"].item()
    lam = map_estimates["lam"].item()
    estimates = pd.DataFrame({'gRNA': [gRNA], 
                              'weight_Poisson': [weights[0]], 
                              'weight_Normal': [weights[1]], 
                              'lambda': [lam], 
                              'mu': [mu], 
                              'scale': [scale]})

    # create plot of the loss
    plot_loss(losses, gRNA, output_dir)
    
    # threshold for which probability is higher to belong to the normal component
    X = np.arange(1, max(selected_guide.toarray())+1, 1)
    log_X = np.log2(X)
    df = pd.DataFrame({'t': X, 'prob_normal_component': prob_normal_component(log_X, weights, mu, scale, lam)})
    threshold = df.loc[(df.prob_normal_component > 0.5), 't'].min()
    
    # create plot of the mixture distribution
    plot_fitted_model(data, weights, mu, scale, lam, np.log2(threshold), gRNA, output_dir)
    
    # get cells with gRNA counts above the threshold
    perturbed_cells = adata_crispr.obs_names[selected_guide.toarray().reshape(-1) >= threshold].tolist()
    return(perturbed_cells, threshold, losses[-1], estimates)


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
        os.makedirs(output_dir + "fitted_model_plots/")
        os.makedirs(output_dir + "loss_plots/")
        print("The output directory " + output_dir +  " was created") 

    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(config['input_file'])
    # subset to specified batches
    adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(config['batch_list'])]

    gRNA_list = adata_crispr.var_names.tolist()
    
    if config['gRNA_list'] != 'all':
        gRNA_list = gRNA_list[config['start_gRNA']:config['end_gRNA'] + 1]
    
    # Downsampling of gRNA counts per cell to a maximum of 800 counts per cell
    #sc.pp.downsample_counts(adata_crispr, counts_per_cell = 800)

    # Fit Poisson-Gaussian Mixture Model (PGMM) for each gRNA
    perturbations = pd.DataFrame({'cell': [], 'gRNA': [], 'target_gene': []})
    thresholds = pd.DataFrame({'gRNA': [], 'threshold': []})
    losses = pd.DataFrame({'gRNA': [], 'loss': []})
    estimates = pd.DataFrame({'gRNA': [], 'weight_Poisson': [], 'weight_Normal': [], 'lambda': [], 'mu': [], 'scale': []})
    print('Fit Poisson-Gaussian Mixture Model for each gRNA: ')
    for gRNA in gRNA_list:
        #time.sleep(0.01)
        perturbed_cells, threshold, loss, map_estimates = fit_PGMM(gRNA, adata_crispr, output_dir, config['seed'])
        if len(perturbed_cells) != 0:
            df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
            perturbations = pd.concat([perturbations, df], ignore_index = True)
            thresholds = pd.concat([thresholds, pd.DataFrame({'gRNA': [gRNA], 'threshold': [threshold]})])
            losses = pd.concat([losses, pd.DataFrame({'gRNA': [gRNA], 'loss': [loss]})])
            estimates = pd.concat([estimates, map_estimates])

    # Save data frames with the results
    perturbations.to_csv(output_dir + 'perturbations.csv', index = False)
    thresholds.to_csv(output_dir + 'gRNA_thresholds.csv', index = False)
    losses.to_csv(output_dir + 'gRNA_losses.csv', index = False)
    estimates.to_csv(output_dir + 'gRNA_estimates.csv', index = False)
    
    # Save run time
    t_end = time.time()
    run_time = pd.DataFrame({'method': ['Replogle'], 'time': [t_end - t_start]})
    run_time.to_csv(output_dir + 'run_time.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    main()
