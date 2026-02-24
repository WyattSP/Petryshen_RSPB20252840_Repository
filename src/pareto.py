# These functions are adapted from the R-Package PtProcess 
# Harte, D. (2010). PtProcess: An R Package for Modelling Marked Point Processes Indexed by Time. Journal of Statistical Software, 35(8), 1-32.
# https://doi.org/10.18637/jss.v035.i08.

import numpy as np 

# d is pdf; p is cdf; q is inverse cdf; r is random values

def dpareto(x, lam, a, log = False):
    """
    Calculates the probability density function (PDF) for a Pareto distribution.

    Parameters:
        x: vector of quantiles
        lam: shape parameter
        a: the random value takes values on the interval a <= x < inf. This is a scaler and is assumed to be constant for all values given in a function call
        log: logical; if TRUE, probabilities p are fiven as log(p)

    Returns:

    """
    x = np.asarray(x)
    lam = np.asarray(lam)
    if np.any(x < a):
        raise ValueError("Invalid data: x < a")
    if np.any(lam <= 0):
        raise ValueError("All values of lambda must be positive")
    if lam.size != 1 and lam.size != x.size:
        raise ValueError("Argument lambda is wrong length")
    if a <= 0:
        raise ValueError("Argument a must be positive")
    d = (lam / a) * (a / x) ** (lam + 1)
    if log:
        d = np.log(d)
    return d


def dtappareto(x, lam, theta, a, log = False):
    """
    Calculates the probability density function (PDF) for a tempered Pareto distribution.

    Parameters:
        x: vector of quantiles
        lam: shape parameter
        theta: tapering parameter
        a: the random value takes values on the interval a <= x < inf. This is a scaler and is assumed to be constant for all values given in a function call
        log: logical; if TRUE, probabilities p are fiven as log(p)

    Returns:

    """
    x = np.asarray(x)
    lam = np.asarray(lam)
    theta = np.asarray(theta)

    if np.any(x < a):
        raise ValueError("Invalid data: x < a")
    if np.any(lam <= 0):
        raise ValueError("All values of lambda must be positive")
    if lam.size != 1 and lam.size != x.size:
        raise ValueError("Argument lambda is wrong length")
    if np.any(theta <= 0):
        raise ValueError("All values of theta must be positive")
    if theta.size != 1 and theta.size != x.size:
        raise ValueError("Argument theta is wrong length")
    if a <= 0:
        raise ValueError("Argument a must be positive")

    d = (lam/x + 1/theta) * ((a/x)**lam * np.exp((a-x)/theta))
    if log:
        d = np.log(d)
    return d

def rpareto(n, lam, a):
    """
    Returns random variables for a Pareto distribution.

    Parameters:
        n: number of observations to simulate
        lam: shape parameter
        a: the random value takes values on the interval a <= x < inf. This is a scaler and is assumed to be constant for all values given in a function call

    Returns:

    """
    lam = np.asarray(lam)

    if np.any(lam <= 0):
        raise ValueError("All values of lambda must be positive")
    if lam.size != 1 and lam.size != n:
        raise ValueError("Argument theta is wrong length")
    if a <= 0:
        raise ValueError("Argument a must be positive")
    #Ensures lam is broadcasted to length n
    lam = lam if lam.size != 1 else lam * np.ones(n)
    return(a * np.exp(np.random.exponential(scale = 1/lam, size = n)))

def rtappareto(n, lam, theta, a):
    """
    Calculates the probability density function (PDF) for a tempered Pareto distribution.

    Parameters:
        n: number of observations to simulate
        lam: shape parameter
        theta: tapering parameter
        a: the random value takes values on the interval a <= x < inf. This is a scaler and is assumed to be constant for all values given in a function call

    Returns:

    """
    lam = np.asarray(lam)
    theta = np.asarray(theta)

    if np.any(theta <= 0):
        raise ValueError("All values of theta must be positive")
    if theta.size != 1 and theta.size != n:
        raise ValueError("Argument theta is wrong length")
    if np.any(lam <= 0):
        raise ValueError("All values of lambda must be positive")
    if lam.size != 1 and lam.size != n:
        raise ValueError("Argument lambda is wrong length")
    if a <= 0:
        raise ValueError("Argument a must be positive")

    # Ensure both lam and theta are broadcasted to length n
    lam = lam if lam.size != 1 else lam * np.ones(n)
    theta = theta if theta.size != 1 else theta * np.ones(n)

    pareto_sample = rpareto(n, lam, a)
    tempered_sample = a + np.random.exponential(scale=theta, size=n)

    return np.minimum(pareto_sample, tempered_sample)