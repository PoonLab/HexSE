import unittest
from src.evol_rates import Rates
from src.sequence_info import Sequence


class TestDrawOmegaValues(unittest.TestCase):
    pass


class TestGetFrequencyRates(unittest.TestCase):
    """
    Tests get_frequency_rates
    """
    def testShortSeq(self):
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
    maxDiff = None

    def testNoORFs(self):
        s = Sequence("ATGTTTCCC")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065, 'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064, 'G': 0.00001}}
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
        s = Sequence("ATGGGGTGA")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065, 'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064, 'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)
        expected = [{'A': [True],  'C': [False],  'G': [False], 'T': [False]},  # A
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
        s = Sequence("ATGCTGAAGTAG")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065, 'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064, 'G': 0.00001}}
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
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG")
        mu = 0.00005
        bias = {'A': {'C': 0.001, 'G': 0.065, 'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064, 'G': 0.00001}}
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
