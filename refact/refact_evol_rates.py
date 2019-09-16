from src.sequence_info import Sequence
from src.sequence_info import NUCLEOTIDES
from src.sequence_info import COMPLEMENT_DICT
from src.sequence_info import CODON_DICT
from scipy.stats import gamma


class Rates(list):
    """
    Create a dictionary with keys as possible substitutions for every nucleotide, and mutation rates for each substitution
    e.g., {'A':, 'T':, 'C':, 'G':}
    """

    def __init__(self, mu, kappa, nucleotide):
        """

        :param nucleotide: Nucleotide object
        :param mu: Global substitution rate ()
        :param kappa: transition transvertion rate ratio
        """
        self.mu =  mu
        self.kappa = kappa
        self.nucleotide = nucleotide


    def get_letter