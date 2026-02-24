# A number of different algorithms for Haar Fluctuation Analysis
# See Lovejoy 2016 or Spiridonov and Lovejoy 2022 for discussion on Haar Fluctuation Analysis
# Code was adapted from R code from Hébert et al. (2020)
# This was all adapted from mathmatica code from http://www.physics.mcgill.ca/~gang/software/index.html
from numba import jit
import numpy as np
import math

# HaarFluctuationAnalysis
@jit(nopython=True)
def Haar(trace, t, type="mean"):
    # Set sampling scale intervals
    mean_tdiff = np.mean(np.diff(t)) # Get max step size
    scales = [2 * mean_tdiff, max(t)//2]
    scales = 2 ** np.arange(np.log2(scales[0]), np.log2(scales[1]), 0.1)
    # Store tendency fluctuations
    out_h_fluc = np.full((len(scales), 3), np.nan)
    # Find sequence of non-overlapping lags with min displacement 1
    max_w = scales[-1]
    #lag_seq = np.arange(1, max_w + 1, 1)  # Keep all values if filtering based on e
    lag_seq = scales
    for i in range(len(lag_seq)-1):
        k = lag_seq[i]
        # Start/end of window
        j_start = np.arange(0, max_w * 2, k)
        j_end = np.arange(k, max_w * 2, k)
        # Store tendency fluctuations
        T_fluc_vec = np.full((len(j_end) - 1, 1), np.nan)
        for j in range(len(j_end) - 1):
            # Need to mask x with t for intervals
            E_lh = np.mean(trace[(t >= j_start[j]) & (t < j_end[j])])
            E_rh = np.mean(trace[(t >= j_start[j + 1]) & (t < j_end[j + 1])])
            T_fluc_vec[j] = abs(E_lh - E_rh) * 2  
        # Root-mean-square of tendency fluctuations
        if type == "rms":
            h_fluc = np.sqrt(np.mean(T_fluc_vec ** 2))
        elif type == "mean":
            h_fluc = np.mean(T_fluc_vec)
        elif type == "hG":
            h_fluc = np.sqrt(2 / np.pi) * np.std(np.vstack((T_fluc_vec, -T_fluc_vec)))
        # Calculate e
        e = np.floor((j_start[0] + j_end[0]) / 2) / j_end[0]
        # Save lag, h_fluc, e
        out_h_fluc[i,:] = [k, h_fluc, e]
    return out_h_fluc

@jit(nopython=True)
def Haar_alloc(trace, t, type="rms", returnFull = False):
    # Can remove boundscheck if desired in @jit
    # Set sampling scale intervals
    mean_tdiff = np.mean(np.diff(t)) # Get max step size
    scales = [2 * mean_tdiff, max(t)//2]
    scales = 2 ** np.arange(np.log2(scales[0]), np.log2(scales[1]), 0.1)

    # Store difference and tendency fluctuations
    out_h_fluc = np.full((len(scales), 3), np.nan)
    
    # Find sequence of non-overlapping lags with min displacement 1
    max_w = scales[-1]
    #scales = np.arange(1, max_w + 1, 1) # Keep all values if filtering based on e

    # Calculate difference and tendency fluctuations
    results_out = np.full(shape = (len(scales), len(scales)-1),fill_value = np.nan).T # allocate memory 
    for i in range(len(scales)-1):
        k = scales[i]
        j_start = np.arange(0, max_w * 2, k)
        j_end = np.arange(k, max_w * 2, k)
        for j in range(len(j_end) - 1):
            # Need to mask x with t for intervals
            E_lh = np.mean(trace[(t >= j_start[j]) & (t < j_end[j])])
            E_rh = np.mean(trace[(t >= j_start[j + 1]) & (t < j_end[j + 1])])
            results_out[j,i] = (abs(E_lh - E_rh) * 2)

    # Average difference/tendency fluctuations
    fluc_out = np.full(shape = results_out.shape[1],fill_value = np.nan)
    if type == "rms":
        for c in range(results_out.shape[1]):
            arr = results_out[~np.isnan(results_out[:,c]),c]
            if arr.size == 0:
                pass
            else:
                fluc_out[c] = np.sqrt(np.mean(arr ** 2))
    elif type == "mean":
        for c in range(results_out.shape[1]):
            arr = results_out[~np.isnan(results_out[:,c]),c]
            if arr.size == 0:
                pass
            else:
                fluc_out[c] = np.mean(arr)
    elif type == "hG":
        for c in range(results_out.shape[1]):
            arr = results_out[~np.isnan(results_out[:,c]),c]
            if arr.size == 0:
                pass
            else:
                fluc_out[c] = np.sqrt(2 / np.pi) * np.std(np.vstack((arr, -arr)))

    # Calculate e
    #e = np.floor((j_start[0] + j_end[0]) / 2) / j_end[0]

    # Save lag, h_fluc, e
    out_h_fluc[:,0] = scales
    out_h_fluc[:,1] = fluc_out
    #out_h_fluc[2,:] = [e]
    
    if returnFull:
        return out_h_fluc, results_out
    else:
        return out_h_fluc, None
    
# Function from Hebert et al. 2021
@jit(nopython=True)
def sdsmpl(array_in):
    s = np.std(array_in) 
    n = len(array_in)
    if n < 101:
        correction_factor = np.sqrt((n - 1) / 2) * np.exp(math.lgamma((n - 1) / 2) - math.lgamma(n / 2))
        return s * correction_factor
    else:
        return s
    
@jit(nopython=True)
def Haar_alloc_v2(trace, t, Scales = None, q = 1, returnFull = False):
    # Can remove boundscheck if desired in @jit
    vec_len = len(t)
    mean_tdiff = np.mean(np.diff(t))
    # Set sampling scale intervals
    if Scales is None:
        # Hashed out lines are translated from Hebert's RScales Haar function. I don't like them so added new below.
        #scales = 2 ** np.arange(np.log2(mean_tdiff),np.log2(np.max(t) + 1 - np.min(t)), 0.1)
        #scales = scales[scales <= (t[vec_len - 1] - t[0] + 1)]
        #scales = scales[scales > np.round(min(np.diff(t)),4)]
        scales = [2 * mean_tdiff, max(t)//2] # Set 2 * min and max fluc
        scales = 2 ** np.arange(np.log2(scales[0]), np.log2(scales[1]), 0.1) # Set all scale jumps
    else:
        scales = [Scales[0] * mean_tdiff, np.max(t) * Scales[1] - np.min(t) * Scales[1]]
        scales = 2 ** np.arange(np.log2(scales[0]), np.log2(scales[1]), np.float64(0.1))

    # Store difference and tendency fluctuations
    out_h_fluc = np.full((len(scales), 3), np.nan)
    
    # Find sequence of non-overlapping lags with min displacement 1
    max_w = scales[-1]
    #scales = np.arange(1, max_w + 1, 1) # Keep all values if filtering based on e

    # Calculate difference and tendency fluctuations
    results_out = np.full(shape = (len(scales), len(scales)-1),fill_value = np.nan).T # allocate memory 
    for i in range(len(scales)-1):
        k = scales[i]
        j_start = np.arange(0, max_w * 2, k)
        j_end = np.arange(k, max_w * 2, k)
        for j in range(len(j_end) - 1):
            # Need to mask x with t for intervals
            E_lh = np.mean(trace[(t >= j_start[j]) & (t < j_end[j])])
            E_rh = np.mean(trace[(t >= j_start[j + 1]) & (t < j_end[j + 1])])
            results_out[j,i] = ((E_lh - E_rh) * 2)

    # Average difference/tendency fluctuations
    fluc_out = np.full(shape = results_out.shape[1],fill_value = np.nan)
    if q == 1:
        for c in range(results_out.shape[1]):
            arr = results_out[~np.isnan(results_out[:,c]),c]
            if arr.size == 0:
                pass
            else:
                fluc_out[c] = np.sqrt(2 / np.pi) * sdsmpl(np.concatenate((arr, -arr)))
    else:
        for c in range(results_out.shape[1]):
            arr = results_out[~np.isnan(results_out[:,c]),c]
            if arr.size == 0:
                pass
            else:
                fluc_out[c] = np.mean(np.abs(arr)**q)**(1/q)

    # Calculate e
    #e = np.floor((j_start[0] + j_end[0]) / 2) / j_end[0]
    # Save lag, h_fluc, e
    out_h_fluc[:,0] = scales
    out_h_fluc[:,1] = fluc_out
    #out_h_fluc[2,:] = [e]
    
    if returnFull:
        return out_h_fluc, results_out
    else:
        return out_h_fluc, None


@jit(nopython=True)
def Haar_hebert(trace, t, Scales = None, q = 1, overlap = 1, returnFull = False, allq = False):
    # Can remove boundscheck if desired in @jit
    vec_len = len(t)
    mean_tdiff = np.mean(np.diff(t))
    # Set sampling scale intervals
    if Scales is None:
        # Hashed out lines are translated from Hebert's RScales Haar function. I don't like them so added new below.
        scales = 2 ** np.arange(np.log2(mean_tdiff),np.log2(np.max(t) + 1 - np.min(t)), 0.1) # Starts at mean of the t-step
        scales = scales[scales <= (t[vec_len - 1] - t[0] + 1)]
        scales = scales[scales > np.round(min(np.diff(t)),4)]
    else:
        scales = [Scales[0] * mean_tdiff, np.max(t) * Scales[1] - np.min(t) * Scales[1]]
        scales = 2 ** np.arange(np.log2(scales[0]), np.log2(scales[1]), 0.1)

    # Store difference and tendency fluctuations
    out_h_fluc = np.full((len(scales), 3), np.nan)

    if overlap >= 1:
        overlap = 0.99
    # Find sequence of non-overlapping lags with min displacement 1
    #max_w = scales[-1]
    #scales = np.arange(1, max_w + 1, 1) # Keep all values if filtering based on e
    e_vec = np.full(shape = len(scales), fill_value = np.nan)
    #e = np.floor((j_start[0] + j_end[0]) / 2) / j_end[0]
    # Calculate difference and tendency fluctuations
    results_out = np.full(shape = (len(scales), len(np.arange(t[0], t[-1], scales[0]))),fill_value = np.nan) # allocate memory 
    max_s = t[-1] # Max scale
    for i in range(len(scales)):
        current_scale = scales[i] # Current scale
        half_width = current_scale/2 # Current fluctuation
        max_fluc = int((max_s - t[0] + 1)//current_scale) # Max fluctuation
        # Set indexs
        step = np.float64((1 - overlap) * current_scale)
        j_start = np.arange(t[0], max_s, step)
        j_mid = np.arange(t[0] + half_width, max_s + current_scale, step)
        # Need to double check epsilon calculation
        e_vec[i] = ((j_start[1]// 2) - j_start[0]) / step
        # Iterate over intervals
        for j in range(max_fluc):
            # Need to mask x with t for intervals
            lh = trace[(t >= j_start[j]) & (t < (j_start[j] + half_width))]
            rh = trace[(t >= j_mid[j]) & (t < (j_mid[j] + half_width))]

            if len(lh) == 0 or len(rh) == 0:
                results_out[i,j] = np.nan
            else:
                results_out[i,j] = np.float64((np.mean(rh) - np.mean(lh)) * 2)

    # Average difference/tendency fluctuations
    fluc_out = np.full(shape = results_out.shape[0],fill_value = np.nan)
    if allq == False:
        if q == 1:
            for c in range(results_out.shape[0]):
                arr = results_out[c, ~np.isnan(results_out[c,:])]
                if arr.size == 0:
                    pass
                else:
                    fluc_out[c] = np.sqrt(2 / np.pi) * sdsmpl(np.concatenate((arr, -arr)))
        else:
            for c in range(results_out.shape[0]):
                arr = results_out[c, ~np.isnan(results_out[c,:])]
                if arr.size == 0:
                    pass
                else:
                    fluc_out[c] = np.mean(np.abs(arr)**q)**(1/q)
    else:
        q_vals = np.arange(0.1, 2.1, 0.1)
        q_fluc_matrix = np.full(shape = (len(q_vals), results_out.shape[0]), fill_value=np.nan)
        for iq in range(len(q_vals)):
            q = q_vals[iq]
            if q == 1:
                for c in range(results_out.shape[0]):
                    arr = results_out[c, ~np.isnan(results_out[c,:])]
                    if arr.size == 0:
                        pass
                    else:
                        q_fluc_matrix[iq,c] = np.sqrt(2 / np.pi) * sdsmpl(np.concatenate((arr, -arr)))
            else:
                for c in range(results_out.shape[0]):
                    arr = results_out[c, ~np.isnan(results_out[c,:])]
                    if arr.size == 0:
                        pass
                    else:
                        q_fluc_matrix[iq,c] = np.mean(np.abs(arr)**q)**(1/q)
        fluc_out = q_fluc_matrix[len(q_vals)//2, :]

    # Save lag, h_fluc, e
    out_h_fluc[:,0] = scales
    out_h_fluc[:,1] = fluc_out
    out_h_fluc[:,2] = e_vec
    
    if returnFull:
        return out_h_fluc, results_out, q_fluc_matrix
    else:
        return out_h_fluc, None, None

   
def arns_haar_fluct(x, t, a, eps, res, plot = False):
    """
    Calculates the haar fluctuation using the algorithm written by Arnscheidt and Rothman 2022 (?)
    
    Parameters:
        x (np.array): Time series of a process under investigation
        t (np.array): Vector of time corresponding to the time series under investigation
        a (float): ???
        eps (float): Cutoff threshold for values
        res (int): Resolution of the haar flucuation
        plot (boolean): If plot == True a scatter plot will be returned
    
    Returns:
        np.array: Returns the Haary tendency flucation along with a log-increasing time scale vector
    """
    nmax = len(x)
    x_cs = np.cumsum(x)
    
    Dx = []
    Dt = []
    width = np.arange(2,(nmax//2)*2,2)
    eps_out = []
    
    for iw in range(0,len(width)):
    
        dint = int(a*width[iw]) # spacing between successive haar fluctuations
        if dint==0:
            dint=1
    
        for i in range(0,(nmax-width[iw])//dint):
            i0 = i*dint
            tstep = t[i0+width[iw]]-t[i0]
            e = (t[i0+width[iw]//2]-t[i0])/tstep
            eps_out.append(e)
            if e>eps and e<(1-eps):
                Ds1 = x_cs[i0+width[iw]//2]-x_cs[i0]
                Ds2 = x_cs[i0+width[iw]]-x_cs[i0+width[iw]//2]
                Dx.append(np.abs((Ds1-Ds2)*2/width[iw])**2)
                Dt.append(tstep)
                
    Dx = np.array(Dx)
    Dt = np.array(Dt)
    
    # average over uniform log-spaced bins
    bin_edges = np.arange(np.floor(np.log10(np.min(Dt))),np.ceil(np.log10(np.max(Dt))),1/res)
    bin_center = bin_edges[1:]-np.diff(bin_edges)
    bin_edges_real = 10**bin_edges
    n = np.histogram(Dt,bin_edges_real)[0]
    vals = np.histogram(Dt,bin_edges_real,weights=Dx)[0]/n
    
    # quantiles
    vals_lower = np.empty_like(vals)
    vals_upper = np.empty_like(vals)
    
    for i in range(0,len(bin_center)):
        temp = Dx[(Dt>bin_edges_real[i])*(Dt<bin_edges_real[i+1])]
        if len(temp)>0:
            vals_lower[i] = np.quantile(temp,0.05)
            vals_upper[i] = np.quantile(temp,0.95)
        else:
            vals_lower[i] = np.nan
            vals_upper[i] = np.nan
    
    # truncate when n is a factor of 5 smaller than maximum
    inds = n>np.max(n)/5
    
    if plot == True:
       pass
        
    # Return v and t
    return(np.sqrt(vals[inds]), 10**bin_center[inds], eps_out)
