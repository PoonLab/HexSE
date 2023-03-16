import scipy
import numpy as np
import scipy.stats as ss

def discretize(alpha, ncat, dist, scale=None):
    """
    Divide a distribution into a number of intervals with equal probability and get the mid point of those intervals
    From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
    :param alpha: int, shape parameter (in gamma), or mean (in lognorm)
    :param ncat: int, Number of categories
    :param dist: string, distribution of probabilities, can be a string
    :param scale: int, scale parameter of the distributions
    :return: array with ncat number of values
    """

    # Handle if distribution was specified from input file (parsed as a string)
    if dist == 'ss.gamma' or dist == 'gamma':
        dist = ss.gamma
    elif dist == 'ss.lognorm' or dist == 'lognorm':
        dist = ss.lognorm

    if dist == ss.gamma:
        # In gamma, mean = shape*scale. 
        # If neutral evolution, omega = 1; therefore scale should be 1/shape
        if scale is None:
            scale=1/alpha
        # Raise an exception when scale is a non positive number
        if scale <= 0:
            raise ValueError("Scale cannot be zero or negative")
        dist = dist(alpha, scale=scale)

    elif dist == ss.lognorm:
        if scale is None:
            scale=1  # default to 1
        dist = dist(s=alpha, scale=scale)

    quantiles = dist.ppf(np.arange(0, ncat) / ncat)
    rates = np.zeros(ncat, dtype=np.double)

    for i in range(ncat-1):
        rates[i] = (ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[i], quantiles[i+1])[0])

    rates[ncat-1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[ncat-1], np.inf)[0]

    return rates