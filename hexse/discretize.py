import scipy
import numpy as np
import scipy.stats as ss

def discretize(alpha, ncat, dist):
    """
    Divide a distribution into a number of intervals with equal probability and get the mid point of those intervals
    From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
    :param alpha: shape parameter
    :param ncat: Number of categories
    :param dist: distribution of probabilities, can be a string
    :return: array with ncat number of values
    """

    # Handle if distribution was specified from input file (parsed as a string)
    if dist == 'ss.gamma' or dist == 'gamma':
        dist = ss.gamma
    elif dist == 'ss.lognorm' or dist == 'lognorm':
        dist = ss.lognorm

    if dist == ss.gamma:
        dist = dist(alpha, scale=0.4)
    elif dist == ss.lognorm:
        dist = dist(s=alpha, scale=0.5)  # scale=np.exp(0.05 * alpha**2)

    quantiles = dist.ppf(np.arange(0, ncat) / ncat)
    rates = np.zeros(ncat, dtype=np.double)

    for i in range(ncat-1):
        rates[i] = (ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[i], quantiles[i+1])[0])

    rates[ncat-1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[ncat-1], np.inf)[0]

    return rates
