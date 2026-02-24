# Function to calculate second-for-third extinction from Alroy (2015)
# Other functions related to diversity metrics
# Functions to get first and last appearance datums
import numpy as np
from numba import njit

def expand_species_time(matrix):
    """
    Given a matrix with:
    - Column 0: Species ID
    - Column 1: First appearance (FAD)
    - Column 2: Last appearance (LAD)
    - Column 3: Duration (LAD - FAD)

    Return a 2-column array:
    - Column 0: Species ID (repeated)
    - Column 1: Each time step the species occurs (from FAD to LAD, inclusive)
    """
    output_rows = []

    for row in matrix:
        species_id = int(row[0])
        fad = int(row[1])
        lad = int(row[2])

        for t in range(fad, lad + 1):  # Include LAD
            output_rows.append([species_id, t])

    return np.array(output_rows, dtype=int)

def FadLabTable(array):
    species = np.unique(array[:,0]) # Unique Species ID
    tmp_arr = np.full((len(species), 3), fill_value = 0, dtype=int)
    for i, sp in enumerate(species):
        fadlad = array[array[:, 0] == sp, 1]  # All time values for species
        fad = np.min(fadlad)
        lad = np.max(fadlad) if len(fadlad) > 1 else 0
        tmp_arr[i, :] = (sp, fad, lad)
    return tmp_arr

def naive_rates(richness_list):
    # Get Speciation and Extinction Rates (raw unbiased rates)
    # Negative Values are Originations
    turno = -np.diff(richness_list[:,1])  

    ext_man = []
    ori_man = []

    for val in turno:
        if val > 0:
            ori_man.append(val)
            ext_man.append(0)
        else:
            ori_man.append(0)
            ext_man.append(-val)

    # Add zero at the start to align with unRich time steps
    unOri = np.array([0] + ori_man)
    unExt = np.array([0] + ext_man)

    return unOri, unExt

def bin_times(timeVector, bins):
    # Simple scaling of original time axis to new time axis based on stage definitions from divDyn
    # A / cf = B
    # scale factor
    cf = (max(timeVector) - min(timeVector)) / np.max(bins)
    # Place into new divisions and bin
    new_taxis = (range(len(timeVector)) / cf)[::-1]
    # Map bin intervals
    stage_bins = np.digitize(new_taxis, bins)
    return stage_bins, new_taxis

def bin_occupancy(occupanyList, binMap):
    # Initalize new list
    new_list = []
    # Iterate over stages
    bin_set = np.unique(binMap)
    for s in bin_set:
        index = np.where(binMap == s)[0]
        tmp = None
        # Iterate over index mapped to each stage
        for i in index:
            if tmp is None:
                tmp = occupanyList[i]
            else:
                tmp = np.hstack((tmp, occupanyList[i]))
        new_list.append(tmp)
    return new_list

def convert_to_dictionary(inputMetrics):
    in_shape = inputMetrics.shape
    occupancy_list = []
    for i in range(in_shape[0]):
        occupancy_list.append(np.array((np.where(inputMetrics[i,:] > 0)[0]), dtype=int))
    return occupancy_list

# Extinction 
def extRateProp(input_dict):
    ext = np.full(len(input_dict) - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, len(input_dict)):
        a = input_dict[i - 1]
        b = input_dict[i]
        ab_l = a.shape[0] + b.shape[0]
        ext[i - 1] = np.size(np.setdiff1d(a, b)) / ab_l
    return ext

# Origination
def oriRateProp(input_dict):
    ori = np.full(len(input_dict) - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, len(input_dict)):
        a = input_dict[i]
        b = input_dict[i - 1]
        ab_l = a.shape[0] + b.shape[0]
        ori[i - 1] = np.size(np.setdiff1d(a,b)) / ab_l
    return ori

# Gapfiller
def extGf(input_dict):
    l_iter = len(input_dict)
    tmp_out = np.full(shape=l_iter - 3, fill_value=0, dtype=np.float64)
    for i in range(l_iter - 3):
        a_index = i
        b_index = i + 1
        c_index = i + 2
        d_index = i + 3
        # Gap-fillers
        g = np.setdiff1d(np.intersect1d(input_dict[a_index], input_dict[d_index]), input_dict[c_index])
        # Three-timers
        t3 = np.intersect1d(np.intersect1d(input_dict[a_index], input_dict[b_index]), input_dict[c_index])
        # Two-timers
        t2 = np.intersect1d(input_dict[a_index], input_dict[b_index])
        # Part-timers
        p = np.intersect1d(input_dict[a_index], input_dict[c_index])
        denominator = np.size(t3) + np.size(p) + np.size(g)
        # Avoid division by zero and log of non-positive numbers
        if denominator == 0:
            tmp_out[i] = 0  # Assign NaN if division is not possible
        else:
            out_e = (np.size(t2) + np.size(p)) / denominator
            tmp_out[i] = np.log(out_e)
    return tmp_out


# Adopted from Alroy 2015
def ext2n3(input_dict):
    l_iter = len(input_dict)
    tmp_out = np.full(shape = l_iter - 3, fill_value = 0, dtype=np.float64)
    for i in range(l_iter - 3):
        a_index = i
        b_index = i + 1
        c_index = i + 2
        d_index = i + 3
        # first-interval taxa
        s1 = np.setdiff1d(np.intersect1d(input_dict[a_index], 
                                            input_dict[b_index]), 
                                            np.concatenate([np.array(sublist) for sublist in input_dict[c_index:]]))
        # third-interval taxa
        s3 = np.intersect1d(input_dict[a_index], input_dict[d_index])
        # two-timers
        t2 = np.intersect1d(input_dict[a_index], input_dict[b_index])
        # part-timers
        p = np.intersect1d(input_dict[a_index], input_dict[c_index])

        denominator = (np.size(t2) + np.size(p))
        if denominator == 0:
            tmp_out[i] = 0  # Assign NaN if division is not possible
        else:
            out_e = (np.size(s1) - np.size(s3)) / denominator
            if out_e == 1:
                tmp_out[i] = 0
            else:
                tmp_out[i] = np.log(1 / (1 - out_e))
    return tmp_out

# Needs to be switch from above
def ori2n3(input_dict):
    l_iter = len(input_dict)
    tmp_out = np.full(shape = l_iter - 3, fill_value = 0, dtype=np.float64)
    for i in range(l_iter - 3):
        a_index = i
        b_index = i + 1
        c_index = i + 2
        d_index = i + 3
        # first-interval taxa
        s1 = np.setdiff1d(np.intersect1d(input_dict[a_index], 
                                            input_dict[b_index]), 
                                            np.concatenate([np.array(sublist) for sublist in input_dict[c_index:]]))
        # third-interval taxa
        s3 = np.intersect1d(input_dict[a_index], input_dict[d_index])
        # two-timers
        t2 = np.intersect1d(input_dict[a_index], input_dict[b_index])
        # part-timers
        p = np.intersect1d(input_dict[a_index], input_dict[c_index])
        out_e = (np.size(s1) - np.size(s3)) / (np.size(t2) + np.size(p))
        tmp_out[i] = np.log(1 / (1 - out_e))
    return 

# FAD LAD metrics based on bins
@njit
def getFadLad(occurence_matrix):
    in_shape = occurence_matrix.shape
    tmp = np.full(shape = (in_shape[1], 3), fill_value = np.nan)
    for col in range(in_shape[1]):
        # Species ID
        tmp[col, 0] = col
        # Occurence 
        occ_vector = np.where(occurence_matrix[:,col] > 0)[0]
        if len(occ_vector) > 1:
            tmp[col, 1] = occ_vector[0] # FAD
            tmp[col, 2] = occ_vector[-1] # LAD
    return tmp

def ExtinctionFadLad(richness_input):
    ext = np.full(shape = len(richness_input), fill_value = np.nan)
    for i in np.arange(1, len(richness_input)):
        a = richness_input[i - 1]
        b = richness_input[i]
        ext[i - 1] = (a - b) / (a + b)
    return ext

def OriginationFadLad(richness_input):
    ori = np.full(shape = len(richness_input), fill_value = np.nan)
    for i in np.arange(1, len(richness_input)):
        a = richness_input[i]
        b = richness_input[i - 1]
        ori[i - 1] = (a - b) / (a + b)
    return ori

# Definitions taken from Foote 2000
# https://geosci.uchicago.edu/~foote/REPRINTS/DeepTime.pdf

def orgRateMatrixProp(full_untb_matrix):
    # Nfl + Nft / Ntot
    m_shape = full_untb_matrix.shape
    ori = np.full(shape = m_shape[0] - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, m_shape[0]):
        a = np.unique(full_untb_matrix[i - 1,:])
        b = np.unique(full_untb_matrix[i,:])

        # Species in a, b 
        o = len(np.setdiff1d(b, a)) # In b not in a
        nt = len(np.union1d(a, b)) # Total in a and b

        ori[i] = o / nt
    return ori

def RateMatrixProp(full_untb_matrix):
    # Nfl + Nbl + Nft / Ntot
    m_shape = full_untb_matrix.shape
    ext = np.full(shape = m_shape[0] - 1, fill_value = 0, dtype=float)
    ori = np.full(shape = m_shape[0] - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, m_shape[0]):
        a = np.unique(full_untb_matrix[i - 1,:])
        b = np.unique(full_untb_matrix[i,:])
        # Species in a, b 
        o = len(np.setdiff1d(b, a, assume_unique=True)) # In b not in a
        e = len(np.setdiff1d(a, b, assume_unique=True)) # In a not in b
        nt = len(np.union1d(a, b)) # Total in a and b
        ori[i-1] = o / nt
        ext[i-1] = e / nt
    return ori, ext

def orgRateMatrixProp2(full_untb_matrix):
    # Nfl + Nft / Ntot
    m_shape = full_untb_matrix.shape
    ori = np.full(shape = m_shape[0] - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, m_shape[0] - 1):
        a = np.unique(full_untb_matrix[i - 1,:])
        b = np.unique(full_untb_matrix[i,:])
        c = np.unique(full_untb_matrix[i + 1,:])

        FL = len(np.setdiff1d(b, np.union1d(a, c))) # In b not in a, c
        bL = len(np.setdiff1d(np.union1d(a, b), c)) # In a, b not in c
        Ft = len(np.setdiff1d(np.union1d(b, c), a)) # In b, c not in a
        bt = len(np.union1d(np.union1d(a, b), c)) # In a, b, c

        #num = len(np.setdiff1d(a, b))
        #den = num + len(np.union1d(a, b))
        ori[i-1] = (FL + Ft) / (FL + bL + Ft + bt)
    return ori

def extRateMatrixProp2(full_untb_matrix):
    # Nfl + Nbl + Nft / Ntot
    m_shape = full_untb_matrix.shape
    ext = np.full(shape = m_shape[0] - 1, fill_value = 0, dtype=float)
    for i in np.arange(1, m_shape[0] - 1):
        #a = np.unique(full_untb_matrix[i - 1,:])
        #b = np.unique(full_untb_matrix[i,:])
        #num = len(np.setdiff1d(a, b))
        #den = num + len(np.union1d(a, b))

        a = np.unique(full_untb_matrix[i - 1,:])
        b = np.unique(full_untb_matrix[i,:])
        c = np.unique(full_untb_matrix[i + 1,:])

        FL = len(np.setdiff1d(b, np.union1d(a, c))) # In b not in a, c
        bL = len(np.setdiff1d(np.union1d(a, b), c)) # In a, b not in c
        Ft = len(np.setdiff1d(np.union1d(b, c), a)) # In b, c not in a
        bt = len(np.union1d(np.union1d(a, b), c)) # In a, b, c

        ext[i-1] = (FL + bL) / (FL + bL + Ft + bt)
    return ext

def RichnessRT(occ_table, bins):
    tmp = np.full(shape = len(bins) - 1, fill_value = np.nan)
    for i in np.arange(0, len(bins) - 1):
        tmp[i] = len(np.intersect1d(np.where(occ_table[:,1] >= bins[i])[0], np.where(occ_table[:,1] <= bins[i + 1])[0]))
    return tmp