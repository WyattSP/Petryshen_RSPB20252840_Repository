# Functions related to sampling system dynamics
import numpy as np
from numba import njit, prange, set_num_threads, guvectorize, vectorize, int64, float64
import scipy.sparse

# Want quantiles of species persistence
def get_persistence(untb_results):
    # Allocate array to save output
    # Size is number of unqiue species and empty row of persistence 
    tmp = np.full(shape = (int(np.max(untb_results)), 2), fill_value=np.nan)
    for i in range(tmp.shape[0]):
        tmp[i, 0] = i + 1
        tmp[i, 1] = np.any(untb_results == i + 1, axis = 1).astype(np.int_).sum()
    return tmp

def filter_sample(sample_array, n):
    mask = np.random.choice(sample_array.shape[1], size=n, replace=False)
    mask = np.sort(mask)
    new_arr = sample_array[:,mask]
    return new_arr

@njit
def regular_time_bins(species_dynamics, nbins):
    tbins = np.linspace(int(species_dynamics[0,0]), int(species_dynamics[0,-1]), nbins)
    tmp = np.full(shape = (species_dynamics.shape[0], len(tbins)), fill_value = 0)
    tmp[0,:] = tbins
    for k in np.arange(1, species_dynamics.shape[0]):
        for i in range(len(tbins) - 1):
            tmp[k,i] = np.nansum(species_dynamics[k,tbins[i]:tbins[i + 1]])
    return tmp

@njit
def numba_any(arr, value):    
    dims = arr.shape
    boolres = np.zeros(dims[0], dtype=np.bool_)
    l_mask = arr == value
    for i in range(dims[0]):
        for j in range(dims[1]):
            if l_mask[i, j]:
                boolres[i] = True
                break
    return np.sum(boolres) 

@njit(parallel = True)
def retrieve_species_dynamics_numba(untb_results, species_id):
    # If species_id is provided it will give full list of species dynamics from that list
    # Get value for all species with abundance > init_pop_cutoff
    n_species = species_id.shape[0]
    s_id = species_id
    tmp = np.full(shape = (untb_results.shape[0], n_species), fill_value = np.nan, dtype = np.int32)
    for i in prange(n_species):
        tmp[:,i] = np.sum(untb_results == s_id[i], axis=1) # Dynamics of those species along time vector
    return tmp

# Apply sampling transfer function
@vectorize([float64(int64, float64)], nopython=True, target="parallel")
def sample_prob(k, s):
    return 1 - ((1 - s) ** k)

@guvectorize([(float64[:], float64[:], float64[:])], '(n),(n)->(n)', nopython=True, target="parallel")
def observer_filter_inequality(spec_occ_arr, random_float_arr, output_arr):
    for i in range(spec_occ_arr.shape[0]):
        output_arr[i] = spec_occ_arr[i] > random_float_arr[i]

@njit
def PA_richness(presence_matrix):
    dims = presence_matrix.shape
    tmp = np.full(shape = dims[0], fill_value = np.nan)
    for i in range(dims[0]):
        tmp[i] = np.sum((presence_matrix[i,:] > 0))
    return tmp

@njit
def threshold_row_filter(array, threshold):
    c_array = array.copy()
    dims = array.shape

    col_iter = range(dims[1])
    row_iter = range(dims[0])

    for c in col_iter:

        result = array[:, c].copy()
        start = None  # To track the start of a sequence of 1s
        for i in row_iter:
            if result[i] == 1 and start is None:  # Start of a sequence
                start = i
            elif result[i] == 0 and start is not None:  # End of a sequence
                if i - start < threshold:  # Sequence length below threshold
                    result[start:i] = [0] * (i - start)  # Replace with 0s
                start = None  # Reset the start

        # Handle case where the row ends with a sequence of 1s
        if start is not None and dims[0] - start < threshold:
            result[start:] = [0] * (dims[0] - start)

        # Put back into array
        c_array[:,c] = result

    return c_array

# Functions below are from old analysis or don't work...
@njit(parallel = True)
def compute_richness(fad, lad, num_elements):
    richness = np.full(shape = num_elements,fill_value = np.nan)
    for i in range(num_elements):
        tf = i >= fad  # True if i is greater than or equal to fad
        tl = i <= lad  # True if i is less than or equal to lad
        richness[i] = np.logical_and(tf, tl).sum()
    return richness

@njit
def fads_lads(sparse_non_zero_element):
    tmpA = sparse_non_zero_element[0]
    tmpB = sparse_non_zero_element[1]
    l_tmpA = len(tmpA)
    # Initialize fad and lad arrays to store results
    fad = np.full(shape = l_tmpA // 2, fill_value=np.nan)
    lad = np.full(shape = l_tmpA // 2, fill_value=np.nan)
    for i in range(0, l_tmpA, 2):  
        fad[i // 2] = tmpB[i]
        lad[i // 2] = tmpB[i + 1]
    return fad, lad


# These are not functioning properly.
@njit(parallel = True)
def prescence_fill(presence_matrix):
    dims = presence_matrix.shape
    for r in range(dims[0]):
        pres = np.where(presence_matrix[r,:] == 1)[0]
        if len(pres) == 2:
            presence_matrix[r, pres[0]:pres[-1]] = 1
        else:
            next
    return presence_matrix

@njit(parallel = True)
def compute_richness_from_presence(prescence_fill_matrix):
    dims = prescence_fill_matrix.shape
    richness = np.full(shape = dims[0], fill_value = np.nan)
    for r in range(dims[0]):
        richness[r] = (prescence_fill_matrix[r,:] > 0).sum()
    return richness

# New Richness function
# These two functions are not working properly
#full_presence = prescence_fill(presence.T)
#t2_richness = compute_richness_from_presence(full_presence.T)

# Apply function to convert all species with single prescence element at the end of the matrix to have last value equal to 1
@njit
def fill_late_species(presence_matrix):
    dims = presence_matrix.shape
    tmp = presence_matrix.copy()
    for c in range(dims[1]):
        occ = np.where(presence_matrix[:,c] == 1)[0]
        if len(occ) == 1:
            tmp[-1, c] = 1
    return tmp

# Filter array for only fad and lad... 
@njit
def only_fad_lad(presence_matrix):
    dims = presence_matrix.shape
    tmp = np.full(shape = dims, fill_value = 0)
    for c in range(dims[1]):
        occ = np.where(presence_matrix[:,c] == 1)[0]
        if len(occ) >= 2:
            tmp[occ[0], c] = 1
            tmp[occ[-1], c] = 1
        elif len(occ) == 1:
            tmp[occ[0], c] = 1
    return tmp 

# This function used below
@njit
def regular_time_bins_presence_matrix_numba(species_dynamics, nbins):
    dims = species_dynamics.shape
    tbins = np.linspace(0, dims[0], nbins)
    tmp = np.full(shape = (nbins - 1, dims[1]), fill_value = 0, dtype=np.float32)
    for t in range(nbins - 1): # Iterate through new time bins
        s = int(np.floor(tbins[t])) # Start time average bin index (floor)
        e = int(np.floor(tbins[t + 1])) # End time average bin index (floor)
        for k in range(dims[1]): # Iterate over all species
            tmp[t, k] = np.nansum(species_dynamics[s:e, k])
    return tmp, np.floor(tbins)