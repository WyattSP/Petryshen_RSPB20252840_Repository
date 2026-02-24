# Code adapted from Spiridonov and Lovejoy 2022
import numpy as np
import scipy

# Consider two time series and their Haar fluctuations with the mean removed.
# The normalized fluctuation correlation coefficients are given by:
# rho = (A**2 + B**2 - ((A-B)**2)) / ((2(A**2)**0.5) * ((B**2)**0.5))

## Equation from Spiridonov and Lovejoy (2022)
def CHFC(A, B):
    flucs = (A.shape[0])
    tmp = np.full(shape = (A.shape[0]), fill_value = np.nan)
    for i in range(flucs):
        a = A[i,~np.isnan(A[i,:])] - np.mean(A[i,~np.isnan(A[i,:])])
        b = B[i,~np.isnan(B[i,:])] - np.mean(B[i,~np.isnan(B[i,:])])
        N = len(a)

        if N <= 0 or a.shape != b.shape:
            tmp[i] = np.nan
        else:
            a_angle = (1 / N) * np.sum(a**2)
            b_angle = (1 / N) * np.sum(b**2)

            c = a - b
            c_angle = (1 / N) * np.sum(c**2)

            num = a_angle + b_angle - c_angle
            dem = 2*((a_angle**(1/2)) * (b_angle**(1/2)))
            tmp[i] = num / dem
    return tmp