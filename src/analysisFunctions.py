# Functions generally used in the analysis Jupyter notebooks.
import numpy as np
from numba import njit, prange, set_num_threads
import scipy
from itertools import zip_longest
from scipy.ndimage import uniform_filter1d

# Leading eigenvalue -- death rate * (speciation rate / population size)
# https://www.sciencedirect.com/science/article/pii/S0025556408000047?casa_token=uJsOivxpXqYAAAAA:qKWLH0ko0rIRh_BE-E4lwA9ad8jgqxj0q4KTACKdnSaofKQyJKl4-MqXtyMCa7SuVzCviEgQ
def firstEigen(u, m, x):
    # u = death rate, m = mutation, x = community size
    return (u * (m / x)) 

def secondEigen(u, m, x):
    return (2 * u *((m/x) + ((1 - m)/(x * (x + 1)))))

# Welford's online algorithm
# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
# Also here: https://nbviewer.org/github/changyaochen/changyaochen.github.io/blob/master/assets/notebooks/welford.ipynb#correctness
def welford_mean_var(i_step, untb_richness, num_simulations, richness_mean, richness_variance):
    # Update mean and variance incrementally
    if richness_mean is None:
        richness_mean = untb_richness
        richness_variance = np.zeros_like(untb_richness)
    else:
        # Update the mean (incrementally)
        new_mean = richness_mean + (untb_richness - richness_mean) / (i_step + 1)
        # Update the variance (incrementally)
        richness_variance += (untb_richness - richness_mean) * (untb_richness - new_mean)
        # Update the mean for the next iteration
        richness_mean = new_mean
    return richness_mean, richness_variance

# Additional functions for used in Fig2_Analysis.ipynb
# Slope estimation
@njit
def linregress_wrapper(data_in):
    # Remember scaling relationship has form: y = Ax**b such that log(y) = b*log(x) + log(A)
    x, y = data_in
    
    G = np.empty((len(x), 2), dtype=np.float64)  # Preallocate matrix
    G[:, 0] = x  # x
    G[:, 1] = 1  # g1
    
    coeffs, _, _, _ = np.linalg.lstsq(G, y)
    #result = scipy.stats.linregress(x, y)
    return coeffs

@njit
def bin_data_max_log_bin(data_in, n_bins, h_cut, l_cut):
    a, b = data_in
    # Cut bins based on largest log interval in log space between frequencies
    # Get Spacing
    cut_ind_l = np.where(a <= l_cut)
    ind_l = cut_ind_l[0][-1]
    interval = a[ind_l] - a[ind_l + 1] 
    # Get bins based on cut intervals within frequency bounds
    cut_ind_h = np.where(a >= h_cut)
    ind_h = cut_ind_h[0][0]
    x_bins = np.arange((a[ind_l]), a[ind_h], -interval)
    # Digitize the data into bins
    indices = np.digitize(a, x_bins, right=True)
    # Calculate the median for each bin
    bin_avg = np.zeros_like(x_bins, dtype = np.float64)
    for i in range(0, len(x_bins)):
        bin_data = b[(indices == i)]
        if bin_data.size > 0:
            bin_avg[i] = np.mean(bin_data)
        else:
            bin_avg[i] = -999  # No data in this bin
    # filter not -1
    inc = np.where(bin_avg != -999)
    o_bin = x_bins[inc]
    o_bin_avg = bin_avg[inc]
    x_bins = x_bins[inc]
    bin_avg = bin_avg[inc]
    return x_bins, bin_avg, o_bin, o_bin_avg

# New functions for Analysis in Fig3_Analysis.ipynb
# New functions to filter short-lived taxa
def add_time_column(dyn_FadLad, time_values):
    # Lets add new column corresponding to boundaries
    t = dyn_FadLad.copy()
    dim_iter = t.shape

    # Add new column
    time_col = np.full(shape = dim_iter[0], fill_value=np.nan)
    new_dynMat = np.column_stack((t, time_col))

    for i in range(dim_iter[0]):
        t_index = t[i,1]
        new_dynMat[i,2] = time_values[t_index].astype(np.float32)

    return new_dynMat

@njit
def filter_short_species(dynMatTime, cutoff):
    tmp = dynMatTime.copy()
    iter_species = np.unique(dynMatTime[:,0])
    retained = 0
    for i in range(len(iter_species)):
        mask = (tmp[:,0] == iter_species[i])
        species_range = np.max(tmp[mask,2]) - np.min(tmp[mask,2])
        if species_range < cutoff:
            tmp = tmp[~mask]
        else:
            retained += 1
    print(f"Number of species removed: {len(iter_species) - retained}")
    return tmp

def shift_species_ID(filtered_divDynMatrix):
    replace_divDyn = filtered_divDynMatrix.copy()
    old_species_ID = np.unique(filtered_divDynMatrix[:,0])
    new_species_ID = np.arange(0,len(np.unique(filtered_divDynMatrix[:,0])))
    for i in range(len(old_species_ID)):
        old_mask = (replace_divDyn[:,0] == old_species_ID[i])
        replace_divDyn[old_mask,0] = new_species_ID[i]
    return replace_divDyn

def persistenceCurve(presenceMatrix):
    # Input is presence-absence matrix
    pres_shape = presenceMatrix.shape
    save_out = np.full(shape = (pres_shape[1],1), fill_value=0)
    for species in range(pres_shape[1]):
        save_out[species,0] = np.sum(presenceMatrix[:,species] > 0)
    srt = np.sort(save_out,axis = 0)[::-1]
    return srt

def persistenceCurve2(presenceFadLad):
    # Input is fad-lad array
    pres_species = np.unique(presenceFadLad[:,0])
    save_out = np.full(shape = (len(pres_species),1), fill_value=0)
    for i in range(len(pres_species)):
        mask = (presenceFadLad[:,0] == pres_species[i])
        save_out[i,0] = np.max(presenceFadLad[mask,1]) - np.min(presenceFadLad[mask,1])
    srt = np.sort(save_out,axis = 0)[::-1]
    return srt

# Additonal functions used for Fig4_Analysis.ipynb
def findDeathRate(tc, N, m):
    return (N / tc)/m

def findTC(N, M, D):
    return 1 / (D * (M / N))

# Smoothing Function with Convolution 
def smooth_log_convolve(y, kernel):
    return np.convolve(np.log10(y), kernel, mode='same')

# Mean and STD for data
def mean_std_uneven(x_array, t_array,  bin_edges, n_bins = 25):
    # To get correct means and confidence intervals you will need to bin it then average across bins
    all_times = np.concatenate(t_array)
    all_x = np.concatenate(x_array)

    if bin_edges is None:
        num_bins = n_bins
        bin_edges = np.linspace(np.min(all_times), np.max(all_times), num_bins + 1)
        bin_indices = np.digitize(all_times, bin_edges)
    else:
        bin_indices = np.digitize(all_times, bin_edges)

    # Compute mean divRT per bin
    binned_means = []
    binned_std = []
    bin_centers = []

    for i in range(1, len(bin_edges)):
        in_bin = bin_indices == i
        if np.any(in_bin):  # avoid empty bins
            mean_x = np.nanmean(all_x[in_bin])
            std_x = np.nanstd(all_x[in_bin])
            binned_means.append(mean_x)
            binned_std.append(std_x)
            # Bin center for reference
            bin_center = (bin_edges[i - 1] + bin_edges[i]) / 2
            bin_centers.append(bin_center)

    # bin_centers → x-axis (time bins)
    # binned_means → y-axis (mean divRT values)
    return np.asarray(binned_means), np.asarray(binned_std), np.asarray(bin_centers)

def moment_curve(data):
    # Get moments
    clean_data = data[~np.isnan(data)]
    mean = np.nanmean(clean_data)
    std = np.nanstd(clean_data)
    # Fit normal distribution
    mu, sigma = scipy.stats.norm.fit(clean_data)
    x = np.linspace(min(clean_data), max(clean_data), 100)
    pdf = scipy.stats.norm.pdf(x, mu, sigma)
    return pdf, x, mean, std

def get_list_moments(value_list):
    zipped = zip_longest(*value_list, fillvalue=None)
    means = [np.mean([x for x in group if x is not None]) for group in zipped]
    zipped = zip_longest(*value_list, fillvalue=None)
    std = [np.std([x for x in group if x is not None]) for group in zipped]
    return np.asarray(means), np.asarray(std)

def smooth_array(arr, window=5):
    return uniform_filter1d(arr, size=window)