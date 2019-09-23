# Create a tree with application of mutation rates according to the landscape of possibilities
import scipy
import scipy.stats as ss
import numpy as np

TRANSITIONS_DICT = {'A':'G', 'G':'A', 'T':'C', 'C':'T'}

class Event_tree:

    def __init__(self, nt_frequencies, omega_values = None):
        """

        :param nt_frequencies: <Dict> Frequency of nucleotides in a given sequence, with nucleotide as keys
        :param kappa: transition/transversion rate ratio
        :param omega_values: <Dict> Numeric values drawn from gamma distribution. Applied in case that a mutation is non_synonymous

        """
        self.nt_frequencies = nt_frequencies
        self.omega_values = omega_values

    def get_omega_values(self):
        return self.omega

    def get_nt_frequencies(self):
        return self.nt_frequencies

    def draw_omega_values(self, alpha, ncat):
        """
        Draw ncat number of omega values from a discretized gamma distribution
        :param alpha: shape parameter
        :param ncat: Number of categories (expected omegas)
        :return: dictionary with number of categories as keys (e.i. {0: 0.29, 1: 0.65, 2: 1.06})
        """
        values = self.discretize(alpha = alpha, ncat = ncat)
        omega_values = list(values)
        return omega_values

    def discretize(self, alpha, ncat, dist=ss.gamma):
        """
        Divide the gamma distribution into a number of intervals with equal probability and get the mid point of those intervals
        From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
        :param alpha: shape parameter
        :param ncat: Number of categories
        :param dist: function from scipy stats
        :return: array with ncat number of values
        """
        if dist == ss.gamma:
            dist = dist(alpha, scale=1 / alpha)
        elif dist == ss.lognorm:
            dist = dist(s=alpha, scale=np.exp(0.5 * alpha ** 2))
        quantiles = dist.ppf(np.arange(0, ncat) / ncat)
        rates = np.zeros(ncat, dtype=np.double) # return a new array of shape ncat and type double
        for i in range(ncat - 1):
            rates[i] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                   quantiles[i], quantiles[i + 1])[0]
        rates[ncat - 1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                      quantiles[ncat - 1], np.inf)[0]
        return rates

    def check_transversion(self, from_nt, to_nt):
        trv_dict = {'is_trv':True}
        if TRANSITIONS_DICT[from_nt] == to_nt:
            trv_dict['is_trv'] = False
        return trv_dict

    def create_event_tree(self):
        event_tree = {'to_nt': {'A':{},
                                'T':{},
                                'C':{},
                                'G':{}}
                        }

        for to_nt in self.nt_frequencies.keys():
            if to_nt in event_tree['to_nt'].keys():
                # Add stationary frequencies to every nucleotide in the event tree
                event_tree['to_nt'][to_nt]['stationary_frequency'] = self.nt_frequencies[to_nt]
                # Update nucleotides with possible mutations
                event_tree['to_nt'][to_nt].update([( 'from_nt',{'A':{},'T':{},'C':{},'G':{}})])
                # For possible mutations, check if they are a transition or a transversion
                for from_nt in event_tree['to_nt'][to_nt]['from_nt'].keys():
                    # Nucleotide cannot change to itself
                    if from_nt == to_nt:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = None
                    else:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = self.check_transversion(from_nt, to_nt)
                        # create key to a dictionary that will store information about nucleotides affected by nonsyn mutations
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_nonsyn', {})])

        return event_tree

s_frequencies = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
events = Event_tree(s_frequencies)
dict = events.create_event_tree()
print(dict['to_nt']['T']['from_nt']['C'])