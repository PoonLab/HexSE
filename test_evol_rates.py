import unittest
from src.evol_rates import Rates
from src.sequence_info import Sequence
from scipy.stats import gamma
import numpy as np


class TestDrawOmegaValues(unittest.TestCase):
    """
    Scipy.stats generates random numbers using numpy.random.
    Two (pseudo)random number generators will generate the same number if they have the same seed (starting value).
    """

    def testSixCodons(self):

        # Sets the seed for the random number generator
        gamma_random_value = gamma(1)
        gamma_random_value.random_state = np.random.seed(seed=1000)

        # Create Sequence and Rate objects
        s = Sequence('ATGAAAGGGTTTACATAG')
        mu = 0.00001
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {(0, 17): [0.04156146032551616, 0.5061606240598411, 0.2654408488563643,
                              1.8435206051609963, 0.23203590181428196, 1.3566172496923332]}
        result = r.draw_omega_values()
        self.assertEqual(expected, result)

    def testTwoOrfs(self):

        # Sets the seed for the random number generator
        gamma_random_value = gamma(1)
        gamma_random_value.random_state = np.random.seed(seed=1000)

        # Create Sequence and Rate objects
        s = Sequence('AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC',
                     [(0, 8), (8, 0), (1, 9), (9, 1), (2, 10), (10, 2)])
        mu = 0.00001
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {(0, 8): [2.6757297014547543, 0.5368803750979751, 0.02940992249316154],
                    (1, 9): [4.018911610141678, 0.4149666303984266, 1.2265156417930274],
                    (2, 10): [0.4492243301978381, 0.03573692247459582, 1.9314233637688747],
                    (8, 0): [0.0774824520810144, 0.2805465373886369, 0.14306786624383916],
                    (9, 1): [1.0707641958670537, 1.4510803545625586, 0.8076319470011387],
                    (10, 2): [2.1636316217836407, 2.3454677028470754, 0.01047638676086102]}
        result = r.draw_omega_values()
        self.assertEqual(expected, result)

    def testNoOrfs(self):

        # Sets the seed for the random number generator
        gamma_random_value = gamma(1)
        gamma_random_value.random_state = np.random.seed(seed=100)

        # Create Sequence and Rate objects
        s = Sequence('ATTGCGACGATGACGCGCAGA')
        mu = 0.00001
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {}
        result = r.draw_omega_values()
        self.assertEqual(expected, result)


class TestGetFrequencyRates(unittest.TestCase):
    """
    Tests get_frequency_rates
    """
    def testShortSeq(self):

        # Create Sequence and Rate objects
        s = Sequence('AAAAAAAAA')
        mu = 0.00001
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        result = r.get_frequency_rates()
        self.assertEqual(expected, result)

    def testSimpleUse(self):

        # Create Sequence and Rate objects
        s = Sequence("GTACGATCGATCGATGCTAGC")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        result = r.get_frequency_rates()
        self.assertEqual(expected, result)

    def testLongerSeq(self):

        # Create Sequence and Rate objects
        s = Sequence("TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACT" 
                     "ACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCA" 
                     "ACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACA" 
                     "GCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCG" 
                     "CTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGT" 
                     "ACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCT" 
                     "TGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCT" 
                     "AGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = {'A': 0.25, 'C': 0.25, 'T': 0.22, 'G': 0.28}
        result = r.get_frequency_rates()
        self.assertEqual(expected, result)


class TestGetSynSubs(unittest.TestCase):
    """
    Tests get_syn_subs
    """

    def testNoORFs(self):

        # Create Sequence and Rate objects
        s = Sequence("ATGTTTCCC")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = [{'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},

                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},

                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []}]
        result = r.get_syn_subs()
        self.assertEqual(expected, result)

    def testSmallOrf(self):

        # Create Sequence and Rate objects
        s = Sequence("ATGGGGTGA")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = [{'A': [True],  'C': [False], 'G': [False],  'T': [False]},  # A
                    {'A': [False], 'C': [False], 'G': [False],  'T': [True]},   # T
                    {'A': [False], 'C': [False], 'G': [True],   'T': [False]},  # G

                    {'A': [False], 'C': [False], 'G': [True],   'T': [False]},  # G
                    {'A': [False], 'C': [False], 'G': [True],   'T': [False]},  # G
                    {'A': [True],  'C': [True],  'G': [True],   'T': [True]},   # G

                    {'A': [False], 'C': [False], 'G': [False],  'T': [True]},   # T
                    {'A': [True],  'C': [False], 'G': [True],   'T': [False]},  # G
                    {'A': [True],  'C': [False], 'G': [False],  'T': [False]}]  # A

        result = r.get_syn_subs()
        self.assertEqual(expected, result)

    def testSimpleUse(self):

        # Create Sequence and Rate objects
        s = Sequence("ATGCTGAAGTAG")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)
        expected = [{'A': [True],  'C': [False], 'G': [False], 'T': [False]},   # A
                    {'A': [False], 'C': [False], 'G': [False], 'T': [True]},    # T
                    {'A': [False], 'C': [False], 'G': [True],  'T': [False]},   # G

                    {'A': [False], 'C': [True],  'G': [False], 'T': [True]},    # C
                    {'A': [False], 'C': [False], 'G': [False], 'T': [True]},    # T
                    {'A': [True],  'C': [True],  'G': [True],  'T': [True]},    # G

                    {'A': [True],  'C': [False], 'G': [False], 'T': [False]},   # A
                    {'A': [True],  'C': [False], 'G': [False], 'T': [False]},   # A
                    {'A': [True],  'C': [False], 'G': [True],  'T': [False]},   # G

                    {'A': [False], 'C': [False], 'G': [False], 'T': [True]},    # T
                    {'A': [True],  'C': [False], 'G': [False], 'T': [False]},   # A
                    {'A': [True],  'C': [False], 'G': [True],  'T': [False]}]   # G

        result = r.get_syn_subs()
        self.assertEqual(expected, result)

    def testMultipleORFs(self):

        # Create Sequence and Rate objects
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065,   'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064,   'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)

        expected = [                                                                                # (0, 23)  (5, 19)
            {'A': [True],         'C': [False],        'G': [False],        'T': [False]},          # A
            {'A': [False],        'C': [False],        'G': [False],        'T': [True]},           # T
            {'A': [False],        'C': [False],        'G': [True],         'T': [False]},          # G

            {'A': [True],         'C': [True],         'G': [False],        'T': [False]},          # A
            {'A': [False],        'C': [False],        'G': [True],         'T': [False]},          # G
            {'A': [True, True],   'C': [False, False], 'G': [True, False],  'T': [False, False]},   # A         A   5

            {'A': [False, False], 'C': [False, False], 'G': [False, False], 'T': [True, True]},     # T         T
            {'A': [False, False], 'C': [False, False], 'G': [True, True],   'T': [False, False]},   # G         G
            {'A': [False, False], 'C': [False, False], 'G': [True, True],   'T': [False, False]},   # G         G   8

            {'A': [False, False], 'C': [True, True],   'G': [False, False], 'T': [False, False]},   # C         C
            {'A': [True, True],   'C': [False, True],  'G': [False, True],  'T': [False, True]},    # A         A
            {'A': [False, False], 'C': [True, True],   'G': [False, False], 'T': [True, False]},    # C         C   11

            {'A': [True, True],   'C': [False, False], 'G': [False, False], 'T': [False, False]},   # A         A
            {'A': [True, True],   'C': [False, False], 'G': [False, True],  'T': [False, False]},   # A         A
            {'A': [True, False],  'C': [False, False], 'G': [True, True],   'T': [False, False]},   # G         G   14

            {'A': [False, False], 'C': [False, False], 'G': [False, False], 'T': [True, True]},     # T         T
            {'A': [False, True],  'C': [False, True],  'G': [True, True],   'T': [False, True]},    # G         G
            {'A': [False, False], 'C': [True, False],  'G': [False, False], 'T': [True, True]},     # T         T   17

            {'A': [True, True],   'C': [False, False], 'G': [False, True],  'T': [False, False]},   # A         A
            {'A': [True, True],   'C': [False, False], 'G': [False, True],  'T': [False, False]},   # A         A
            {'A': [False],        'C': [True],         'G': [False],        'T': [True]},           # C

            {'A': [False],        'C': [False],        'G': [False],        'T': [True]},           # T
            {'A': [True],         'C': [False],        'G': [False],        'T': [False]},          # A
            {'A': [True],         'C': [False],        'G': [True],         'T': [False]}]          # G

        result = r.get_syn_subs()
        self.assertEqual(expected, result)
