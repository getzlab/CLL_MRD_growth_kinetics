import pandas as pd
from helper import *
from matplotlib import pyplot as plt
import dalmatian
import numpy as np
def load_wbc_file(file):
    """
    Parse wbc file

    Args:
        wbc file

    Returns:
        variables as a list

    """

    df = pd.read_csv(file, sep='\t')

    # Get all time points
    times = df.dfd.tolist()

    # Get times that are in the interpolation range
    t = pd.notna(df['sample_id'])
    sample_time_index = [i for i, x in enumerate(t) if x]
    first_sample_index = sample_time_index[0]
    last_sample_index = sample_time_index[-1]
    df_interpolate = df.iloc[first_sample_index:last_sample_index + 1, :]
    times_interpolate = df_interpolate.dfd.tolist()
    wbc_interpolate = df_interpolate.wbc.tolist()

    # Get all wbc points
    wbc = df.wbc.tolist()

    # Get the sample list
    sample_list = df.sample_id[pd.notna(df.sample_id)].tolist()

    # Get the times at sample points
    times_sample = df.dfd[pd.notna(df.sample_id)].tolist()

    # Get the wbc at sample points
    wbc_sample = df.wbc[pd.notna(df.sample_id)].tolist()

    # Get the times that are not at sample points
    times_others = sorted(list(set(times) - set(times_sample)))

    return times_sample, times, times_interpolate, wbc, wbc_interpolate, sample_list, wbc_sample, times_others


# Function to load data
def load_data(patient_id, wbc_file, treatment_file):
    # Load WBC and treatment data
    wbc_df = pd.read_csv(wbc_file)
    treatment_df = pd.read_csv(treatment_file, sep='\t')

    # Get input files from Terra through dalmatian
    workspace = 'broad-firecloud-ibmwatson/TAG_CLL_Clonal_Kinetic_UMI_PrAN'
    wm = dalmatian.WorkspaceManager(workspace)
    participants = wm.get_participants()

    # Load additional data for the patient
    cluster_CCF_df = pd.read_csv(participants.loc[patient_id]['cluster_ccfs'], sep='\t')
    abundance_df = pd.read_csv(participants.loc[patient_id]['abundances_tsv'], sep='\t')
    mcmc_df = pd.read_csv(participants.loc[patient_id]['cell_population_mcmc_trace'], sep='\t')
    tree_df = pd.read_csv(participants.loc[patient_id]['tree_tsv'], sep='\t')

    return wbc_df, treatment_df, cluster_CCF_df, abundance_df, mcmc_df, tree_df





def plot_ccf(df, ax, times_sample, treatment):
    # Keep the necessary columns
    cols = ['Sample_ID', 'Cluster_ID', 'postDP_ccf_mean', 'postDP_ccf_CI_low', 'postDP_ccf_CI_high']
    df = df[cols]
    cluster_list = df.Cluster_ID.unique().tolist()
    number_samples = len(df.Sample_ID.unique())

    tick_list = ['T' + str(i) for i in range(number_samples)]
    x_axis = [i / 365 for i in times_sample]
    ax.set_xticks(x_axis)

    secax = ax.secondary_xaxis('top')
    secax.set_xlabel("Time (years)")
    ax.grid(True)

    for i in cluster_list:
        x = df[df.Cluster_ID == i].Sample_ID
        y = df[df.Cluster_ID == i].postDP_ccf_mean
        ci_low = df[df.Cluster_ID == i].postDP_ccf_CI_low
        ci_high = df[df.Cluster_ID == i].postDP_ccf_CI_high

        #         if x_scale == 'sample':
        #             x_axis = np.arange(0,number_samples)

        #         else:
        #             x_axis = times_sample
        ax.plot(x_axis, y, c=ClusterColors.get_hex_string(i), marker='o', label=i)

        #         ax.plot(x_axis, y,c= ClusterColors.get_hex_string(i), marker ='o')
        ax.fill_between(x_axis, ci_low, ci_high, color=ClusterColors.get_hex_string(i), alpha=0.1)

        ax.set_xlabel('Samples')
        ax.set_xticks(x_axis)

        ax.set_xticklabels(tick_list, fontsize=8)

        ax.set_ylabel('CCF')
        ax.legend()

    cmap = plt.get_cmap("Pastel1")
    xlim = ax.get_xlim()[1]
    for i, row in treatment.iterrows():
        treatment_name = row.tx
        start = row.tx_start / 365
        end = row.tx_end / 365
        if np.isnan(end):
            end = xlim
        length = end - start
        center = (start + end) / 2

        ax.axvspan(xmin=start, xmax=end, facecolor=cmap(i), alpha=0.2)


def get_abundance(abundance, mcmc, sample_list):
    """
    Parse abundance file and get abundance information and calculate interpolated abundance

    Args:
        cell_population_abundance file and mcmc_trace file
        mode: "interpolate_only" or "extrapolate"

    Returns:
        lists

    """

    cluster_map = {}

    # Get the cluster list from the mcmc file
    cluster_list = mcmc.Cluster_ID.unique().tolist()

    # Map the abundacne to each cluster in a dictionary
    for cluster_id in cluster_list:
        abundance['cluster'] = abundance['Cell_population'].str.split('_').str[-1].str[2:].astype(int)

        cluster_map[cluster_id] = abundance[abundance.cluster == cluster_id]

    ## To get abundance information from the cell_population_abundances.tsv
    cluster_abundance = {}
    for cluster_id, abundances in cluster_map.items():
        cluster_abundances = []
        ## iterate through the samples in the wbc file to make sure the order is correct
        for sample_name in sample_list:
            sample_abundances = float(abundances[abundances.Sample_ID == sample_name].Cell_abundance.iloc[0])
            cluster_abundances.append(sample_abundances / 100)
        cluster_abundance[cluster_id] = cluster_abundances


    return cluster_list, cluster_abundance


def calc_subclone(wbc, abundance, cluster_list, input_type="default"):
    """
    Calculate subclone abundance from cluster abundance

    Args:
        wbc and abundance

        input type: "default" or "mcmc (contains values of all iterations)"
    Returns:
        dictionaries

    """

    subclone_population = {}
    log_subclone = {}

    if input_type == "default":
        for cluster_id in cluster_list:
            subclone_population[cluster_id] = [_wbc * (_abundance + 1e-4) for _wbc, _abundance in
                                               zip(wbc, abundance[cluster_id])]

            log_subclone[cluster_id] = [np.log(item) for item in subclone_population[cluster_id]]

    if input_type == "mcmc":
        for cluster_id in cluster_list:

            subclone_population_per_iter = {}
            log_subclone_per_iter = {}

            for iter_idx in range(250):
                subclone_population_per_iter[iter_idx] = [_wbc * (_abundance + 1e-4) for _wbc, _abundance in
                                                          zip(wbc, abundance[cluster_id][iter_idx])]

                log_subclone_per_iter[iter_idx] = [np.log(item) for item in subclone_population_per_iter[iter_idx]]

            subclone_population[cluster_id] = subclone_population_per_iter
            log_subclone[cluster_id] = log_subclone_per_iter

    return subclone_population, log_subclone


def get_all_abundance(cluster_list, mcmc_df, sample_list, times_sample):
    """
    Get cell abundance and interpolated cell abundance from all mcmc iterations

    Args:
        cluster_list, mcmc_df, sample_list
    Returns:
        a dictionary with key: cluster and value: 250 iterations of cell abundance

    """
    all_abundance = {}

    for cluster in cluster_list:

        abundance_per_iter = {}

        for iter_idx in range(250):

            df_mcmc_iter_clust = mcmc_df[(mcmc_df.Iteration == iter_idx) & (mcmc_df.Cluster_ID == cluster)]

            cluster_abundances = []
            ## iterate through the samples in the wbc file to make sure the order is correct
            for sample_name in sample_list:
                sample_abundance = float(df_mcmc_iter_clust[df_mcmc_iter_clust.Sample_ID == sample_name].Abundance.iloc[0])
                cluster_abundances.append(sample_abundance / 100)

            abundance_per_iter[iter_idx] = cluster_abundances

        all_abundance[cluster] = abundance_per_iter

    return all_abundance

# New model

import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MultiClusterLinearRegression:
    def __init__(self, n_clusters, X, y):
        """
        Initialize the multi-cluster linear regression model.

        Parameters:
        - n_clusters: Number of clusters.
        - X: Input features (1D array of shape (n_samples,)).
        - y: Target values (1D array of shape (n_samples,)). None values indicate missing data.
        """
        self.n_clusters = n_clusters
        self.params = None
        self.X = X
        self.y = y

    def calculate_weight_ratio(self):
        """
        Calculate the weight ratio for balancing exome and wbc contributions.

        Returns:
        - weight_ratio: Ratio of WBC counts to exome counts.
        """
        num_wbc = len(self.y)
        num_exomes = sum(1 for item in self.y if item is not None)
        weight_ratio = num_wbc / ((num_exomes) * (self.n_clusters - 1))
        return weight_ratio

    def objective(self, params, logsumexp_points):
        """
        Objective function to minimize for fitting the model.

        Parameters:
        - params: Model parameters (intercepts and slopes).
        - logsumexp_points: Logsumexp target values (1D array of shape (n_samples,)). None values indicate missing data.

        Returns:
        - Negative log-likelihood to minimize.
        """
        intercepts = params[:self.n_clusters]
        slopes = params[self.n_clusters:]

        # Calculate predicted y values for all clusters
        y_pred = np.outer(self.X, slopes) + intercepts

        # Calculate negative log-likelihood
        likelihood = 0

        # For exome points
        for i, yi in enumerate(self.y):
            if yi is not None:

                likelihood += self.calculate_weight_ratio() * np.sum((yi - y_pred[i]) ** 2)

        # For logsumexp points
        for i, lse in enumerate(logsumexp_points):
            if lse is not None:
                likelihood += (lse - logsumexp(y_pred[i])) ** 2

        return likelihood  # Minimize negative log-likelihood

    def fit(self, logsumexp_points, initial_guess=None):
        """
        Fit the model to the data.

        Parameters:
        - logsumexp_points: Logsumexp target values (1D array of shape (n_samples,)). None values indicate missing data.
        - initial_guess: Initial guess for model parameters. If None, random values are used.
        """
        if initial_guess is None:
            initial_guess = np.random.randn(2 * self.n_clusters)

        result = minimize(
            self.objective,
            initial_guess,
            args=(logsumexp_points),
            method='SLSQP'
        )

        self.params = result.x

    def calculate_likelihood(self, params, logsumexp_points):
        """
        Calculate the likelihood components for debugging or analysis.

        Parameters:
        - params: Model parameters (intercepts and slopes).
        - logsumexp_points: Logsumexp target values (1D array of shape (n_samples,)). None values indicate missing data.

        Returns:
        - likelihood_y: Likelihood contribution from exome points.
        - likelihood_logsumexp: Likelihood contribution from logsumexp points.
        """
        intercepts = params[:self.n_clusters]
        slopes = params[self.n_clusters:]
        y_pred = np.outer(self.X, slopes) + intercepts

        # Calculate likelihood for exome points
        likelihood_y = 0
        for i, yi in enumerate(self.y):
            if yi is not None:
                likelihood_y += np.sum((yi - y_pred[i]) ** 2)
                logger.info(f"Exome logsubclone: {yi}, Exome predicted: {y_pred[i]}, Likelihood_y: {likelihood_y}")

        logger.info(f"Squared sum_y: {likelihood_y}")

        # Calculate likelihood for logsumexp points
        likelihood_logsumexp = 0
        for i, lse in enumerate(logsumexp_points):
            if lse is not None:
                likelihood_logsumexp += (lse - logsumexp(y_pred[i])) ** 2
                logger.info(f"Log WBC: {lse}, WBC predicted: {logsumexp(y_pred[i])}")

        logger.info(f"Squared sum_logsumexp: {likelihood_logsumexp}")

        total_squared_sum = likelihood_y + likelihood_logsumexp
        logger.info(f"Total squared sum: {total_squared_sum}")

        return likelihood_y, likelihood_logsumexp

    def predict(self, x_to_predict):
        """
        Predict target values using the fitted model.

        Returns:
        - Predicted values for all clusters (2D array of shape (n_samples, n_clusters)).
        """
        if self.params is None:
            raise ValueError("Model has not been fitted yet. Call `fit` first.")

        intercepts = self.params[:self.n_clusters]
        slopes = self.params[self.n_clusters:]
        return np.outer(x_to_predict, slopes) + intercepts


def create_inputs(times_sliced, log_subclone_sample_mcmc, iter, index_samples_model, times_sample):
    """
    Create input features (X) and target values (y) for the model.

    Parameters:
    - times_sliced: List of time points (in days).
    - log_subclone_sample_mcmc: Dictionary of subclone data in the format {cluster: {iteration: [values]}}.
    - iter: The iteration number to use for extracting values.
    - index_samples_model: Indices of the samples to include in the target values.

    Returns:
    - X: Input features (time in years).
    - y: Target values (list with None for missing data).
    """
    # Convert times to years
    X = [t / 365 for t in times_sliced]

    # Reconstruct the dictionary to {iteration: {cluster: [values]}}
    log_subclone_cluster_iter = reconstruct_log_subclone_dict(log_subclone_sample_mcmc)

    # Extract values for the specified iteration
    values = list(log_subclone_cluster_iter[iter].values())
    array = np.array(values).T  # Transpose to match the desired shape

    # Select the samples corresponding to the specified indices
    selected_sample = array[index_samples_model]

    # Map the selected samples to the correct positions in y
    y = map_samples_to_y(times_sliced, times_sample[index_samples_model], selected_sample)

    return X, y

def reconstruct_log_subclone_dict(log_subclone_sample_mcmc):
    """
    Reconstruct the log_subclone_sample_mcmc dictionary to {iteration: {cluster: [values]}}.

    Parameters:
    - log_subclone_sample_mcmc: Dictionary of subclone data in the format {cluster: {iteration: [values]}}.

    Returns:
    - log_subclone_cluster_iter: Reconstructed dictionary.
    """
    log_subclone_cluster_iter = {}
    for cluster, iterations in log_subclone_sample_mcmc.items():
        for iteration, value in iterations.items():
            if iteration not in log_subclone_cluster_iter:
                log_subclone_cluster_iter[iteration] = {}
            log_subclone_cluster_iter[iteration][cluster] = value
    return log_subclone_cluster_iter


def map_samples_to_y(times_sliced, times_sample_selected, selected_sample):
    """
    Map the selected samples to the correct positions in y.

    Parameters:
    - times_sliced: List of time points (in days).
    - times_sample_selected: List of selected time points (in days).
    - selected_sample: List of values corresponding to the selected time points.

    Returns:
    - y: List with values at the correct positions and None elsewhere.
    """
    y = [None] * len(times_sliced)
    indices = [i for i, t in enumerate(times_sliced) if t in times_sample_selected]

    for idx, sample_idx in enumerate(indices):
        y[sample_idx] = selected_sample[idx]

    return y


