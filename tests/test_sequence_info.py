import random
import unittest

from src.sequence_info import Sequence
from src.sequence_info import Nucleotide

MAX_DIFF = None
KAPPA = 0.3
GLOBAL_RATE = 1

OMEGA_VALUES_4 = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]

OMEGA_VALUES_5 = [0.25496553479835254, 0.5492012532921686, 0.842647578837806, 1.2308533419685328, 2.1223322911031457]

MU_VALUES_4 = {'mu1': 0.325970770628452,
               'mu2': 0.7739880200789496,
               'mu3': 1.481877317174674,
               'mu4': 4.351175963585551}

MU_VALUES_3 = {'mu1': 0.3965034306394888,
               'mu2': 1.0832728336432287,
               'mu3': 3.7199827893060657}

EMPTY_EVENT_TREE_MU_4 = {'to_nt': {'A': {'from_nt': {'A': None,
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}}}},
                                   'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'C': None,
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}}}},
                                   'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'G': None,
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}}}},
                                   'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {},
                                                                        'mu4': {}}},
                                                     'T': None}}}}

EMPTY_EVENT_TREE_MU_3 = {'to_nt': {'A': {'from_nt': {'A': None,
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}}}},
                                   'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'C': None,
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}}}},
                                   'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'G': None,
                                                     'T': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}}}},
                                   'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'C': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'G': {'category': {'mu1': {},
                                                                        'mu2': {},
                                                                        'mu3': {}}},
                                                     'T': None}}}}


class TestSequence1(unittest.TestCase):
    """
    Sequence: GTACGATCGATCGATGCTAGC
    ORFs:
        (0, 21) (+)
    Notes:
        No START or STOP codons
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9001)     # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, KAPPA, GLOBAL_RATE, pi1, OMEGA_VALUES_4, MU_VALUES_4)

        self.seq1_codons = self.sequence1.find_codons('+0', [(0, 21)])

    def testReverseComplement(self):
        s = str(self.sequence1)
        expected = 'GCTAGCATCGATCGATCGTAC'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'CATGCTAGCTAGCTACGATCG'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence1 = self.sequence1.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence1, new_sequence1)
        self.assertEqual(self.sequence1.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence1.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence1.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
        self.assertEqual(self.sequence1.omega_values, new_sequence1.omega_values)
        self.assertEqual(self.sequence1.cat_values, new_sequence1.cat_values)
        self.assertEqual(self.sequence1.is_circular, new_sequence1.is_circular)
        self.assertEqual(self.sequence1.total_omegas, new_sequence1.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence1.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence1.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq1_nt in enumerate(self.sequence1.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(seq1_nt, new_nt)

            self.assertEqual(seq1_nt.state, new_nt.state)
            self.assertEqual(seq1_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq1_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq1_nt.rates, new_nt.rates)
            self.assertEqual(seq1_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq1_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq1_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq1_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq1_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq1_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq1_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.24, 'C': 0.24, 'G': 0.29, 'T': 0.24}
        result = Sequence.get_frequency_rates(str(self.sequence1))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_4
        result = self.sequence1.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator
        # Only 4 possible omegas because there is only 1 ORF
        exp_total_omegas = {(0, 0, 0, 1, 0): 1.9817219453273531,
                            (0, 0, 1, 0, 0): 1.0699896623909886,
                            (0, 1, 0, 0, 0): 0.6550136761581515,
                            (1, 0, 0, 0, 0): 0.29327471612351436}

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)

        exp_sub_rates = {'A': 0.12603317336204048, 'C': 0.030344325868823886, 'G': None, 'T': 0.030344325868823886}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omegas = {'A': (1, 0, 0, 0, 0), 'C': (0, 0, 1, 0, 0), 'G': None, 'T': (0, 0, 1, 0, 0)}
        self.assertEqual(exp_omegas, nt.omega_keys)
        exp_cat_keys = {'A': 'mu3', 'C': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.18672182509968824
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence1.total_omegas)

    def testIsTransv(self):
        nt = self.sequence1.nt_sequence[2]  # A
        result = self.sequence1.is_transv(nt.state, 'A')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence1.is_transv(nt.state, 'G')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence1.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

    def testCodonIterator(self):
        orf = self.sequence1.nt_sequence[0: 21]
        result = []
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                    ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]

        for codon in self.sequence1.codon_iterator(orf, 0, 21):
            result.append(codon)

        for pos, cdn in enumerate(result):
            for idx, nt in enumerate(cdn):
                self.assertEqual(expected[pos][idx], nt.state)

    def testFindCodons(self):
        # Test forward strand ORF
        expected = ['GTA', 'CGA', 'TCG', 'ATC', 'GAT', 'GCT', 'AGC']
        result = self.sequence1.find_codons('+0', [(0, 21)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStopStartCodon(self):
        nt = self.sequence1.nt_sequence[0]  # G in codon GTA
        expected = False
        result = self.sequence1.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {(0, 0, 0, 1, 0): 0.39634438906547004,
                                                            (0, 0, 0, 0, 1): 0.1999999999999997}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {(0, 0, 0, 0, 1): 0.1999999999999997}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                            (0, 1, 0, 0, 0): 0.1310027352316301,
                                                            (1, 0, 0, 0, 0): 0.058654943224702785}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {(0, 1, 0, 0, 0): 0.1310027352316301}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {(0, 1, 0, 0, 0): 0.1310027352316301}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {(1, 0, 0, 0, 0): 0.058654943224702785,
                                                            (0, 0, 0, 0, 1): 0.1999999999999997}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}}}}}},
                     'T': {'from_nt': {'A': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 0, 0, 1):  0.1999999999999997}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(1, 0, 0, 0, 0): 0.058654943224702785}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}}}},
                                       'T': None,
                                       'C': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(1, 0, 0, 0, 0): 0.058654943224702785}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(1, 0, 0, 0, 0): 0.058654943224702785,
                                                                       (0, 0, 0, 0, 1): 0.1999999999999997}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 0, 1):  0.1999999999999997}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {}}}},
                                       'G': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301,
                                                                       (1, 0, 0, 0, 0): 0.058654943224702785}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {}}}}}},
                     'C': {'from_nt': {'A': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (1, 0, 0, 0, 0): 0.058654943224702785}}}},
                                       'T': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 1, 0, 0): 0.2139979324781974,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(1, 0, 0, 0, 0): 0.058654943224702785,
                                                                       (0, 0, 0, 0, 1): 0.1999999999999997}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997}}}},
                                       'C': None,
                                       'G': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004,
                                                                       (0, 0, 1, 0, 0): 0.2139979324781974}}}}}},
                     'G': {'from_nt': {'A': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                       (0, 0, 1, 0, 0): 0.2139979324781974,
                                                                       (0, 0, 0, 1, 0): 0.39634438906547004,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}}}},
                                       'T': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004,
                                                                       (0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004,
                                                                       (0, 0, 0, 0, 1): 0.1999999999999997}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {}}}},
                                       'C': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.04701719357592987,
                                                             'omega': {(0, 1, 0, 0, 0): 0.1310027352316301}},
                                                     'mu2': {'prob': 0.11163806035536104,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004,
                                                                       (1, 0, 0, 0, 0): 0.058654943224702785}},
                                                     'mu3': {'prob': 0.21374220928782836,
                                                             'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                                     'mu4': {'prob': 0.6276025367808807,
                                                             'omega': {}}}},
                                       'G': None}}}}

        result = self.sequence1.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(9001)

        g0 = self.sequence1.nt_sequence[0]
        t1 = self.sequence1.nt_sequence[1]
        a2 = self.sequence1.nt_sequence[2]
        c3 = self.sequence1.nt_sequence[3]
        g4 = self.sequence1.nt_sequence[4]
        a5 = self.sequence1.nt_sequence[5]
        t6 = self.sequence1.nt_sequence[6]
        c7 = self.sequence1.nt_sequence[7]
        g8 = self.sequence1.nt_sequence[8]
        a9 = self.sequence1.nt_sequence[9]
        t10 = self.sequence1.nt_sequence[10]
        c11 = self.sequence1.nt_sequence[11]
        g12 = self.sequence1.nt_sequence[12]
        a13 = self.sequence1.nt_sequence[13]
        t14 = self.sequence1.nt_sequence[14]
        g15 = self.sequence1.nt_sequence[15]
        c16 = self.sequence1.nt_sequence[16]
        t17 = self.sequence1.nt_sequence[17]
        a18 = self.sequence1.nt_sequence[18]
        g19 = self.sequence1.nt_sequence[19]
        c20 = self.sequence1.nt_sequence[20]

        exp_event_tree = {'to_nt':
                          {'A': {'from_nt':
                                 {'A': None,
                                  'T': {'category': {'mu1': {(0, 0, 0, 1, 0): [t6]},
                                                     'mu2': {(0, 0, 0, 1, 0): [t14],
                                                             (0, 0, 0, 0, 1): [t17]},
                                                     'mu3': {(0, 0, 1, 0, 0): [t1]},
                                                     'mu4': {(0, 0, 0, 1, 0): [t10]}}},
                                  'C': {'category': {'mu1': {(0, 0, 0, 0, 1): [c11]},
                                                     'mu2': {(0, 0, 0, 0, 1): [c3],
                                                             (0, 1, 0, 0, 0): [c16],
                                                             (1, 0, 0, 0, 0): [c20]},
                                                     'mu3': {},
                                                     'mu4': {}}},
                                  'G': {'category': {'mu1': {(0, 1, 0, 0, 0): [g15]},
                                                     'mu2': {(0, 1, 0, 0, 0): [g4, g19]},
                                                     'mu3': {(1, 0, 0, 0, 0): [g0],
                                                             (0, 0, 0, 0, 1): [g8]},
                                                     'mu4': {(0, 0, 1, 0, 0): [g12]}}}}},
                           'T': {'from_nt':
                                 {'A': {'category': {'mu1': {(0, 0, 0, 0, 1): [a5]},
                                                     'mu2': {(0, 0, 0, 0, 1): [a2],
                                                             (0, 1, 0, 0, 0): [a13]},
                                                     'mu3': {(1, 0, 0, 0, 0): [a9]},
                                                     'mu4': {(0, 0, 0, 1, 0): [a18]}}},
                                  'T': None,
                                  'C': {'category': {'mu1': {(1, 0, 0, 0, 0): [c16]},
                                                     'mu2': {(1, 0, 0, 0, 0): [c7],
                                                             (0, 0, 0, 0, 1): [c20]},
                                                     'mu3': {(0, 0, 0, 0, 1): [c11]},
                                                     'mu4': {}}},
                                  'G': {'category': {'mu1': {(0, 0, 1, 0, 0): [g0]},
                                                     'mu2': {(0, 0, 1, 0, 0): [g4]},
                                                     'mu3': {(0, 0, 0, 0, 1): [g8],
                                                             (0, 1, 0, 0, 0): [g12, g15],
                                                             (1, 0, 0, 0, 0): [g19]},
                                                     'mu4': {}}}}},
                           'C': {'from_nt':
                                 {'A': {'category': {'mu1': {(0, 0, 1, 0, 0): [a9]},
                                                     'mu2': {},
                                                     'mu3': {(0, 0, 0, 0, 1): [a5],
                                                             (0, 1, 0, 0, 0): [a18]},
                                                     'mu4': {(0, 0, 0, 0, 1): [a2],
                                                             (1, 0, 0, 0, 0): [a13]}}},
                                  'T': {'category': {'mu1': {(0, 0, 1, 0, 0): [t1],
                                                             (0, 1, 0, 0, 0): [t6]},
                                                     'mu2': {(1, 0, 0, 0, 0): [t10],
                                                             (0, 0, 0, 0, 1): [t14]},
                                                     'mu3': {},
                                                     'mu4': {(0, 0, 0, 0, 1): [t17]}}},
                                  'C': None,
                                  'G': {'category': {'mu1': {(0, 0, 1, 0, 0): [g0]},
                                                     'mu2': {},
                                                     'mu3': {(0, 0, 0, 0, 1): [g8],
                                                             (0, 1, 0, 0, 0): [g12]},
                                                     'mu4': {(0, 0, 0, 1, 0): [g4, g19],
                                                             (0, 0, 1, 0, 0): [g15]}}}}},
                           'G': {'from_nt':
                                 {'A': {'category': {'mu1': {},
                                                     'mu2': {},
                                                     'mu3': {},
                                                     'mu4': {(0, 0, 0, 0, 1): [a2, a5],
                                                             (0, 0, 1, 0, 0): [a9],
                                                             (0, 0, 0, 1, 0): [a13],
                                                             (0, 1, 0, 0, 0): [a18]}}},
                                  'T': {'category': {'mu1': {(0, 0, 0, 1, 0): [t10]},
                                                     'mu2': {(0, 0, 0, 1, 0): [t1],
                                                             (0, 1, 0, 0, 0): [t14]},
                                                     'mu3': {(0, 0, 0, 1, 0): [t6],
                                                             (0, 0, 0, 0, 1): [t17]},
                                                     'mu4': {}}},
                                  'C': {'category': {'mu1': {(0, 1, 0, 0, 0): [c20]},
                                                     'mu2': {(0, 0, 0, 1, 0): [c3],
                                                             (1, 0, 0, 0, 0): [c11, c16]},
                                                     'mu3': {(0, 0, 0, 1, 0): [c7]},
                                                     'mu4': {}}},
                                  'G': None}}}}

        res_omega_key = self.sequence1.nt_in_event_tree(g0)
        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 61
        res_count = self.sequence1.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 10)
        nt.rates = {'A': None, 'C': 0.3132846693781597, 'G': 1.0442822312605322, 'T': 0.05572713744568437}
        nt.get_mutation_rate()
        exp_mutation_rate = 1.4132940380843764
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': 0.12603317336204048, 'C': 0.030344325868823886, 'G': None,                'T': 0.030344325868823886},
            {'A': 0.11416272554221871, 'C': 0.08370848515537624,  'G': 0.11043569122638641, 'T': None},
            {'A': None,                'C': 0.3132846693781597,   'G': 1.0442822312605322,  'T': 0.05572713744568437},
            {'A': 0.05572713744568437, 'C': None,                 'G': 0.11043569122638641, 'T': 0.0},
            {'A': 0.14702209411694175, 'C': 0.7501854178665989,   'G': None,                'T': 0.07204984868600821},
            {'A': None,                'C': 0.10669516683657652,  'G': 1.0442822312605322,  'T': 0.02346989548524854},
            {'A': 0.0465108069376564,  'C': 0.0512436750694675,   'G': 0.2114401535804069,  'T': None},
            {'A': 0.0,                 'C': None,                 'G': 0.2114401535804069,  'T': 0.054477868049197166},
            {'A': 0.42974442198065543, 'C': 0.1289233265941966,   'G': None,                'T': 0.1289233265941966},
            {'A': None,                'C': 0.025112545546612873, 'G': 1.1173711920673652,  'T': 0.03129099476574798},
            {'A': 0.6208431044413233,  'C': 0.054477868049197166, 'G': 0.0465108069376564,  'T': None},
            {'A': 0.02346989548524854, 'C': None,                 'G': 0.01634336041475915, 'T': 0.35565055612192176},
            {'A': 1.3501568570813995,  'C': 0.08444654209500271,  'G': None,                'T': 0.08444654209500271},
            {'A': None,                'C': 0.09187847247772883,  'G': 2.069477014804411,   'T': 0.0365020371600683},
            {'A': 0.11043569122638641, 'C': 0.1857571248189479,   'G': 0.0365020371600683,  'T': None},
            {'A': 0.0619194407089399,  'C': 0.40504705712441985,  'G': None,                'T': 0.08444654209500271},
            {'A': 0.0365020371600683,  'C': None,                 'G': 0.01634336041475915, 'T': 0.02294375645294939},
            {'A': 0.05572713744568437, 'C': 1.0442822312605322,   'G': 0.10669516683657652, 'T': None},
            {'A': None,                'C': 0.06988679345793328,  'G': 0.6840191432445981,  'T': 0.6208431044413233},
            {'A': 0.14702209411694175, 'C': 0.7501854178665989,   'G': None,                'T': 0.03780995200861215},
            {'A': 0.01634336041475915, 'C': None,                 'G': 0.01537310252084025, 'T': 0.1857571248189479}
        ]

        exp_total_rates = [0.18672182509968824,
                           0.30830690192398136,
                           1.4132940380843764,
                           0.16616282867207077,
                           0.9692573606695488,
                           1.1744472935823573,
                           0.3091946355875308,
                           0.2659180216296041,
                           0.6875910751690486,
                           1.173774732379726,
                           0.7218317794281769,
                           0.39546381202192943,
                           1.519049941271405,
                           2.1978575244422083,
                           0.3326948532054026,
                           0.5514130399283624,
                           0.07578915402777683,
                           1.2067045355427932,
                           1.3747490411438545,
                           0.9350174639921528,
                           0.2174735877545473]

        for pos, nt in enumerate(self.sequence1.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq1_codons[0]     # GTA
        nt = self.sequence1.nt_sequence[1]  # T
        expected = 1
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['G', 'T', 'A']
        exp_mutated = ['G', 'T', 'C']
        res_codon, res_mutated = self.seq1_codons[0].mutate_codon(2, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        exp_codon = ['G', 'A', 'T']
        exp_mutated = ['T', 'A', 'T']
        res_codon, res_mutated = self.seq1_codons[4].mutate_codon(0, 'T')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq1_codons[0]          # GTA = Valine

        # Mutation at wobble position
        expected = False                     # Synonymous
        result = codon.is_nonsyn(2, 'C')     # GTC = Valine
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(0, 'C')    # CTA = Leucine
        self.assertEqual(expected, result)

    def testIsStop(self):
        codon = self.seq1_codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.creates_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.creates_stop(1, 'A')
        self.assertEqual(expected, result)

    def testIsStart(self):
        codon = self.seq1_codons[0]               # GTA = Valine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq1_codons[2]     # TCG
        expected = True
        result = codon.creates_stop(1, 'A')     # TCG --> TAG
        self.assertEqual(expected, result)

        result = codon.creates_stop(1, 'G')     # TCG -> TGG
        expected = False
        self.assertEqual(expected, result)


class TestSequence2(unittest.TestCase):
    """
    Sequence: ATGAATAAACCCGTATGA
    ORFs (indexing relative to forward strand)
        (0, 18) (+) ATG AAT AAA CCC GTA TGA
        (3, 15) (-)
    Notes:
        2 overlapping ORFs
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4001)
        s2 = 'ATGAATAAACCCGTATGA'
        sorted_orfs = {'+0': [[(0, 18)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [[(3, 15)]]}
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, sorted_orfs, KAPPA, GLOBAL_RATE, pi2, OMEGA_VALUES_5, MU_VALUES_3)

        self.plus_0_codons = self.sequence2.find_codons('+0', [(0, 18)])
        self.minus_2_codons = self.sequence2.find_codons('-2', [(3, 15)])

    def testReverseComplement(self):
        s = str(self.sequence2)
        expected = 'TCATACGGGTTTATTCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTTATTTGGGCATACT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence2 = self.sequence2.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence2, new_sequence2)
        self.assertEqual(self.sequence2.orfs, new_sequence2.orfs)
        self.assertEqual(self.sequence2.kappa, new_sequence2.kappa)
        self.assertEqual(self.sequence2.global_rate, new_sequence2.global_rate)
        self.assertEqual(self.sequence2.pi, new_sequence2.pi)
        self.assertEqual(self.sequence2.omega_values, new_sequence2.omega_values)
        self.assertEqual(self.sequence2.cat_values, new_sequence2.cat_values)
        self.assertEqual(self.sequence2.is_circular, new_sequence2.is_circular)
        self.assertEqual(self.sequence2.total_omegas, new_sequence2.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence2.event_tree, new_sequence2.event_tree)
        self.assertCountEqual(self.sequence2.event_tree, new_sequence2.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq2_nt in enumerate(self.sequence2.nt_sequence):
            new_nt = new_sequence2.nt_sequence[pos]
            self.assertIsNot(seq2_nt, new_nt)

            self.assertEqual(seq2_nt.state, new_nt.state)
            self.assertEqual(seq2_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq2_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq2_nt.rates, new_nt.rates)
            self.assertEqual(seq2_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq2_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq2_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq2_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq2_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq2_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq2_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_3
        result = self.sequence2.create_event_tree()
        self.assertEqual(expected, result)

    def testGetFrequencyRates(self):
        expected = {'A': 0.44, 'C': 0.17, 'G': 0.17, 'T': 0.22}
        result = Sequence.get_frequency_rates(str(self.sequence2))
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(4001)   # Set seed value to initialize pseudo-random number generator

        exp_total_omegas = {(0, 0, 0, 1, 0, 1): 1.2308533419685328,
                            (0, 0, 0, 1, 1, 0): 2.61227979327204,
                            (0, 0, 1, 0, 0, 1): 0.842647578837806,
                            (0, 0, 1, 0, 1, 0): 1.7883781665873595,
                            (0, 1, 0, 0, 0, 1): 0.5492012532921686,
                            (0, 1, 0, 0, 1, 0): 1.1655875541762872,
                            (0, 1, 1, 0, 0, 0): 0.4627831063813346,
                            (0, 2, 0, 0, 0, 0): 0.3016220166176888,
                            (1, 0, 0, 0, 0, 1): 0.25496553479835254,
                            (1, 0, 0, 0, 1, 0): 0.5411215876209263,
                            (1, 0, 0, 1, 0, 0): 0.31382518059334646}

        # Synonymous mutation in (+) strand, non-synonymous mutation in (-) frame
        nt = self.sequence2.nt_sequence[11]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.2335164156624226, 'C': None, 'G': 0.030341674692605226, 'T': 0.10113891564201742}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 0, 0, 1, 0, 1), 'C': None, 'G': (0, 1, 0, 0, 0, 1), 'T': (0, 1, 0, 0, 0, 1)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu3', 'G': 'mu2', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.36499700599704527
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation would destroy a START codon
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (0, 0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0, 0), 'T': (0, 0, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation in STOP codon on reverse strand ORF
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (0, 0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0, 0), 'T': (0, 0, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Non-synonymous mutation in both (+) and (-) strands
        random.seed(4001)
        nt = self.sequence2.nt_sequence[10]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.037346151696315194, 'C': None, 'G': 0.037346151696315194, 'T': 0.17608224296983835}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 1, 0, 1, 0, 0), 'C': None, 'G': (0, 1, 0, 1, 0, 0), 'T': (0, 0, 0, 1, 1, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu2', 'G': 'mu2', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.25077454636246876
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {(0, 0, 0, 1, 0, 1): 1.2308533419685328,
                            (0, 0, 0, 1, 1, 0): 2.61227979327204,
                            (0, 0, 1, 0, 0, 1): 0.842647578837806,
                            (0, 0, 1, 0, 1, 0): 1.7883781665873595,
                            (0, 1, 0, 0, 0, 1): 0.5492012532921686,
                            (0, 1, 0, 0, 1, 0): 1.1655875541762872,
                            (0, 1, 0, 1, 0, 0): 0.6759861980279724,
                            (0, 1, 1, 0, 0, 0): 0.4627831063813346,
                            (0, 2, 0, 0, 0, 0): 0.3016220166176888,
                            (1, 0, 0, 0, 0, 1): 0.25496553479835254,
                            (1, 0, 0, 0, 1, 0): 0.5411215876209263,
                            (1, 0, 0, 1, 0, 0): 0.31382518059334646}
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

    def testIsTransv(self):
        nt = self.sequence2.nt_sequence[9]  # C
        result = self.sequence2.is_transv(nt.state, 'C')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence2.is_transv(nt.state, 'T')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence2.is_transv(nt.state, 'G')
        expected = True
        self.assertEqual(expected, result)

    def testFindCodons(self):

        # Test forward strand ORF
        expected = ['ATG', 'AAT', 'AAA', 'CCC', 'GTA', 'TGA']
        result = self.sequence2.find_codons('+0', [(0, 18)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Test reverse strand ORF
        expected = ['ATG', 'CCC', 'AAA', 'TAA']
        result = self.sequence2.find_codons('-2', [(3, 15)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-2')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence2.nt_sequence[5]  # T in codon AAT (+ strand) and TAA in (- strand)
        expected = True
        result = self.sequence2.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence2.nt_sequence[14]  # A in codon GTA (+ strand) and ATG (- strand)
        expected = True
        result = self.sequence2.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {(1, 0, 0, 0, 0, 1): 0.0026373110962595162}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {(0, 0, 0, 1, 1, 0): 0.26568047126179006,
                                                            (0, 1, 0, 0, 0, 1): 0.005680824902599456}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {(0, 1, 0, 0, 1, 0): 0.11854543739457091}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {(1, 0, 0, 1, 0, 0): 0.0319174163841851}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {(0, 1, 0, 0, 0, 1): 0.005680824902599456,
                                                            (1, 0, 0, 0, 0, 1): 0.0026373110962595162}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {(0, 1, 1, 0, 0, 0): 0.0470670998229415}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {(1, 0, 0, 0, 1, 0): 0.055034471720573526}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {(0, 0, 1, 0, 1, 0): 0.18188601210212918}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {(1, 0, 0, 0, 0, 1): 0.0026373110962595162,
                                                            (0, 0, 1, 0, 0, 1): 0.008716173390504544}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {(0, 0, 1, 0, 0, 1): 0.008716173390504544}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {(1, 0, 0, 0, 0, 1): 0.0026373110962595162,
                                                            (0, 2, 0, 0, 0, 0): 0.03067630034283003}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {}}}},
                            'G': None}}}}

        result = self.sequence2.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(4001)

        # All other nucleotides are part of START or STOP codons
        a6 = self.sequence2.nt_sequence[6]
        a7 = self.sequence2.nt_sequence[7]
        a8 = self.sequence2.nt_sequence[8]
        c9 = self.sequence2.nt_sequence[9]
        c10 = self.sequence2.nt_sequence[10]
        c11 = self.sequence2.nt_sequence[11]

        exp_event_tree = {'to_nt':
                          {'A': {'from_nt': {'A': None,
                                             'T': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}},
                                             'C': {'category': {'mu1': {},
                                                                'mu2': {(1, 0, 0, 0, 0, 1): [c9]},
                                                                'mu3': {(0, 0, 0, 1, 1, 0): [c10],
                                                                        (0, 1, 0, 0, 0, 1): [c11]}}},
                                             'G': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}}}},
                           'T': {'from_nt': {'A': {'category': {'mu1': {(0, 1, 0, 0, 1, 0): [a7]},
                                                                'mu2': {},
                                                                'mu3': {}}},
                                             'T': None,
                                             'C': {'category': {'mu1': {(1, 0, 0, 1, 0, 0): [c10]},
                                                                'mu2': {},
                                                                'mu3': {(0, 1, 0, 0, 0, 1): [c9],
                                                                        (1, 0, 0, 0, 0, 1): [c11]}}},
                                             'G': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}}}},
                           'C': {'from_nt': {'A': {'category': {'mu1': {(0, 1, 1, 0, 0, 0): [a6, a8]},
                                                                'mu2': {},
                                                                'mu3': {(1, 0, 0, 0, 1, 0): [a7]}}},
                                             'T': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}},
                                             'C': None,
                                             'G': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}}}},
                           'G': {'from_nt': {'A': {'category': {'mu1': {(0, 0, 1, 0, 1, 0): [a7]},
                                                                'mu2': {(1, 0, 0, 0, 0, 1): [a6],
                                                                        (0, 0, 1, 0, 0, 1): [a8]},
                                                                'mu3': {}}},
                                             'T': {'category': {'mu1': {},
                                                                'mu2': {},
                                                                'mu3': {}}},
                                             'C': {'category': {'mu1': {(0, 0, 1, 0, 0, 1): [c11]},
                                                                'mu2': {(1, 0, 0, 0, 0, 1): [c9],
                                                                        (0, 2, 0, 0, 0, 0): [c10]},
                                                                'mu3': {}}},
                                             'G': None}}}}

        res_omega_key = self.sequence2.nt_in_event_tree(c11)
        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 16
        res_count = self.sequence2.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('T', 2)
        nt.rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        nt.get_mutation_rate()
        exp_mutation_rate = 0
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                  'G': None,                 'T': 0.0},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': None,                 'C': 0.024221351790530223, 'G': 0.12152678443944394,  'T': 0.0},
            {'A': None,                 'C': 0.2657111150590803,   'G': 0.31200315446636473,  'T': 0.06100504924028973},
            {'A': None,                 'C': 0.02422135179053022,  'G': 0.4016395814157039,   'T': 0.0},
            {'A': 0.014086059105481004, 'C': None,                 'G': 0.014086059105481004, 'T': 0.34731326571907206},
            {'A': 0.4955994294630238,   'C': None,                 'G': 0.016663685768162066, 'T': 0.02115356932347423},
            {'A': 0.10419397971572163,  'C': None,                 'G': 0.01703974544729171,  'T': 0.161239458223735},
            {'A': 0.0,                  'C': 0.0,                  'G': None,                 'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                  'G': None,                 'T': 0.0},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           0,
                           0,
                           0,
                           0.14574813622997415,
                           0.6387193187657347,
                           0.4258609332062341,
                           0.37548538393003406,
                           0.5334166845546601,
                           0.2824731833867483,
                           0,
                           0,
                           0,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence2.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testSetTotalOmegas(self):
        random.seed(4001)
        # Last position in tuple is the synonymous flag
        chosen_omegas = [(0, 0, 0, 0, 0),   (1, 0, 0, 0, 0),    (0, 1, 0, 0, 0),
                         (0, 0, 1, 0, 0),   (0, 0, 0, 1, 0),    (0, 0, 0, 0, 1),
                         (1, 1, 0, 0, 0),   (0, 1, 1, 0, 0),    (0, 0, 1, 1, 0),
                         (0, 0, 0, 1, 1),   (2, 0, 0, 0, 0),    (0, 2, 0, 0, 0),
                         (0, 0, 2, 0, 0),   (0, 0, 0, 2, 0),    (0, 0, 0, 0, 2),
                         (2, 1, 0, 0, 0),   (0, 2, 1, 0, 0),    (0, 0, 2, 1, 0),
                         (0, 0, 0, 2, 1),   (1, 0, 0, 0, 2)]

        self.sequence2.total_omegas = {}

        exp_total_omegas = {(0, 0, 0, 1, 0): 1.2308533419685328,
                             (0, 0, 0, 1, 1): 1.2308533419685328,
                             (0, 0, 0, 2, 0): 1.5149999494351059,
                             (0, 0, 0, 2, 1): 1.5149999494351059,
                             (0, 0, 1, 0, 0): 0.842647578837806,
                             (0, 0, 1, 1, 0): 1.0371755885142062,
                             (0, 0, 2, 0, 0): 0.7100549421212166,
                             (0, 0, 2, 1, 0): 0.8739734984911726,
                             (0, 1, 0, 0, 0): 0.5492012532921686,
                             (0, 1, 1, 0, 0): 0.4627831063813346,
                             (0, 2, 0, 0, 0): 0.3016220166176888,
                             (0, 2, 1, 0, 0): 0.25416106202707195,
                             (1, 0, 0, 0, 0): 0.25496553479835254,
                             (1, 0, 0, 0, 2): 0.25496553479835254,
                             (1, 1, 0, 0, 0): 0.14002739125756325,
                             (2, 0, 0, 0, 0): 0.06500742393500993,
                             (2, 1, 0, 0, 0): 0.03570215869840277}

        for omega in chosen_omegas:
            self.sequence2.set_total_omegas(omega)

        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

    def testNtInPos(self):
        plus_0_codon = self.plus_0_codons[4]    # GTA
        nt = self.sequence2.nt_sequence[12]     # G
        expected = 0
        result = plus_0_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        minus_2_codon = self.minus_2_codons[0]  # ATG
        expected = 2
        result = minus_2_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['C', 'C', 'C']
        exp_mutated = ['C', 'C', 'G']
        res_codon, res_mutated = self.plus_0_codons[3].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        res_codon, res_mutated = self.minus_2_codons[1].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        plus_0_codon = self.plus_0_codons[2]     # AAA
        expected = False                         # Synonymous
        result = plus_0_codon.is_nonsyn(2, 'G')
        self.assertEqual(expected, result)

        expected = True
        result = plus_0_codon.is_nonsyn(2, 'C')  # Lys --> Asn
        self.assertEqual(expected, result)

    def testIsStop(self):
        codon = self.plus_0_codons[5]
        exp = True
        result = codon.is_stop()
        self.assertEqual(exp, result)

        codon = self.minus_2_codons[2]
        exp = False
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        plus_0_codon = self.plus_0_codons[0]    # ATG
        expected = True
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

        plus_0_codon = self.plus_0_codons[4]    # GTA
        expected = False
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

        minus_2_codon = self.minus_2_codons[0]  # ATG
        expected = True
        result = minus_2_codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.minus_2_codons[3]
        exp = True
        result = codon.creates_stop(2, 'G')
        self.assertEqual(exp, result)

        codon = self.plus_0_codons[1]
        exp = False
        result = codon.creates_stop(2, 'A')
        self.assertEqual(exp, result)


class TestSequence3(unittest.TestCase):
    """
    Sequence: ATGACGTGGTGA
    ORFs:
        (0, 12) (+)
    Notes:
        Simple ORF
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(555)

        s3 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, OMEGA_VALUES_4, MU_VALUES_3)

        self.seq3_codons = self.sequence3.find_codons('+0', [(0, 12)])

    def testReverseComplement(self):
        s = str(self.sequence3)
        expected = 'TCACCACGTCAT'   # Reverse and complement
        result = Sequence.complement(s, rev=True)
        self.assertEqual(expected, result)

        expected = 'TACTGCACCACT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence3 = self.sequence3.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence3, new_sequence3)
        self.assertEqual(self.sequence3.orfs, new_sequence3.orfs)
        self.assertEqual(self.sequence3.kappa, new_sequence3.kappa)
        self.assertEqual(self.sequence3.global_rate, new_sequence3.global_rate)
        self.assertEqual(self.sequence3.pi, new_sequence3.pi)
        self.assertEqual(self.sequence3.omega_values, new_sequence3.omega_values)
        self.assertEqual(self.sequence3.cat_values, new_sequence3.cat_values)
        self.assertEqual(self.sequence3.is_circular, new_sequence3.is_circular)
        self.assertEqual(self.sequence3.total_omegas, new_sequence3.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence3.event_tree, new_sequence3.event_tree)
        self.assertCountEqual(self.sequence3.event_tree, new_sequence3.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq3_nt in enumerate(self.sequence3.nt_sequence):
            new_nt = new_sequence3.nt_sequence[pos]
            self.assertIsNot(seq3_nt, new_nt)

            self.assertEqual(seq3_nt.state, new_nt.state)
            self.assertEqual(seq3_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq3_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq3_nt.rates, new_nt.rates)
            self.assertEqual(seq3_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq3_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq3_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq3_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq3_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq3_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq3_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.25, 'C': 0.08, 'G': 0.42, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence3))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_3
        result = self.sequence3.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):

        # Tests a non-synonymous mutation
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[8]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.053456076057738736, 'G': None, 'T': 0.307016589860621}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 0, 0, 0, 0), 'C': (0, 0, 1, 0, 0), 'G': None, 'T': (0, 1, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'T': 'mu3'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.36047266591835975
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        # Only 4 possible omegas because there is only 1 ORF
        exp_total_omegas = {(0, 0, 0, 1, 0): 1.9817219453273531,
                            (0, 0, 1, 0, 0): 1.0699896623909886,
                            (0, 1, 0, 0, 0): 0.6550136761581515,
                            (1, 0, 0, 0, 0): 0.29327471612351436}
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.1665314408685853, 'C': 0.13649237703904682, 'G': None, 'T': 0.04995943226057559}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 0, 0, 0, 1), 'C': (0, 0, 0, 0, 1), 'G': None, 'T': (0, 0, 0, 0, 1)}  # Synonymous
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.3529832501682077
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a mutation that would destroy a STOP codon
        random.seed(555)
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 0, 0, 0, 0), 'C': (0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0), 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'G': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

    def testIsTransv(self):
        nt = self.sequence3.nt_sequence[1]  # T
        result = self.sequence3.is_transv(nt.state, 'T')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence3.is_transv(nt.state, 'C')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence3.is_transv(nt.state, 'G')
        expected = True
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ACG', 'TGG', 'TGA']
        result = self.sequence3.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence3.nt_sequence[0]      # A in START codon
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[11]      # A in STOP codon (TGA)
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[4]
        expected = False
        result = self.sequence3.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[8]  # Last G in TGG codon
        expected = True     # Introduces a STOP codon
        result = self.sequence3.is_start_stop_codon(nt, 'A')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):

        expected = {'to_nt': {'A':
                              {'from_nt': {'A': None,
                                           'T': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {(1, 0, 0, 0, 0): 0.058654943224702785}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}},
                                           'C': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {(0, 1, 0, 0, 0): 0.1310027352316301}}}},
                                           'G': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {(0, 0, 0, 0, 1): 0.1999999999999997}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}}}},
                              'T':
                              {'from_nt': {'A': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}},
                                           'T': None,
                                           'C': {'prob': 0.625,
                                                 'cat':
                                                     {'mu1': {'prob': 0.07625419304108505,
                                                              'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                                      'mu2': {'prob': 0.20833135198746808,
                                                              'omega': {}},
                                                      'mu3': {'prob': 0.7154144549714468,
                                                              'omega': {}}}},
                                           'G': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {(0, 0, 0, 0, 1): 0.1999999999999997}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {(1, 0, 0, 0, 0): 0.058654943224702785}}}}}},
                              'C':
                              {'from_nt': {'A': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}},
                                           'T': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {(0, 0, 0, 1, 0): 0.39634438906547004}}}},
                                           'C': None,
                                           'G': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {(0, 1, 0, 0, 0): 0.1310027352316301}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {(0, 0, 0, 0, 1): 0.1999999999999997,
                                                                           (1, 0, 0, 0, 0): 0.058654943224702785}}}}}},
                              'G':
                              {'from_nt': {'A': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}}}},
                                           'T': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {(0, 0, 1, 0, 0): 0.2139979324781974}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}},
                                           'C': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.20833135198746808,
                                                                 'omega': {(1, 0, 0, 0, 0): 0.058654943224702785}},
                                                         'mu3': {'prob': 0.7154144549714468,
                                                                 'omega': {}}}},
                                           'G': None}}}}
        result = self.sequence3.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(555)

        a3 = self.sequence3.nt_sequence[3]
        c4 = self.sequence3.nt_sequence[4]
        g5 = self.sequence3.nt_sequence[5]
        t6 = self.sequence3.nt_sequence[6]
        g7 = self.sequence3.nt_sequence[7]
        g8 = self.sequence3.nt_sequence[8]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {(1, 0, 0, 0, 0): [t6]},
                                                                         'mu3': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 1, 0, 0, 0): [c4]}}},
                                                      'G': {'category': {'mu1': {(0, 0, 0, 0, 1): [g5]},
                                                                         'mu2': {}, 'mu3': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {(0, 0, 0, 1, 0): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {}}},
                                                      'T': None,
                                                      'C': {'category': {'mu1': {(0, 0, 0, 1, 0): [c4]},
                                                                         'mu2': {},
                                                                         'mu3': {}}},
                                                      'G': {'category': {'mu1': {(0, 0, 0, 0, 1): [g5]},
                                                                         'mu2': {(0, 0, 1, 0, 0): [g8]},
                                                                         'mu3': {(1, 0, 0, 0, 0): [g7]}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {(0, 0, 0, 1, 0): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 0, 1, 0): [t6]}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {(0, 1, 0, 0, 0): [g8]},
                                                                         'mu3': {(0, 0, 0, 0, 1): [g5],
                                                                                 (1, 0, 0, 0, 0): [g7]}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 1, 0, 0): [a3]}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {(0, 0, 1, 0, 0): [t6]},
                                                                         'mu3': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {(1, 0, 0, 0, 0): [c4]},
                                                                         'mu3': {}}},
                                                      'G': None}}}}

        res_omega_key = self.sequence3.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 16
        res_count = self.sequence3.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('C', 4)
        nt.rates = {'A': 0.0584793504496421, 'C': None, 'G': 0.0076247167865047885, 'T': 0.06286076399166855}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.12896483122781544
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None,                 'C': 0.0,                 'G': 0.0,                   'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                   'T': None},
            {'A': 0.0,                  'C': 0.0,                 'G': None,                  'T': 0.0},
            {'A': None,                 'C': 0.05893196624218927, 'G': 0.9950857822074713,    'T': 0.05893196624218927},
            {'A': 0.0584793504496421,   'C': None,                'G': 0.0076247167865047885, 'T': 0.06286076399166855},
            {'A': 0.1665314408685853,   'C': 0.4687178314525643,  'G': None,                  'T': 0.04995943226057559},
            {'A': 0.023827239957827467, 'C': 1.8429928824519723,  'G': 0.08693180501604358,   'T': None},
            {'A': 0.0,                  'C': 0.13746308896128004, 'G': None,                  'T': 0.13746308896128004},
            {'A': 0.0,                  'C': 0.08940437365191053, 'G': None,                  'T': 0.14604543242695323},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                   'T': None},
            {'A': 0.0,                  'C': 0.0,                 'G': None,                  'T': 0.0},
            {'A': None,                 'C': 0.0,                 'G': 0.0,                   'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           1.11294971469185,
                           0.12896483122781544,
                           0.6852087045817252,
                           1.9537519274258432,
                           0.2749261779225601,
                           0.23544980607886376,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence3.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq3_codons[2]         # TGG
        nt = self.sequence3.nt_sequence[6]  # T
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['T', 'G', 'A']
        exp_mutated = ['T', 'C', 'A']
        res_codon, res_mutated = self.seq3_codons[3].mutate_codon(1, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq3_codons[0]         # ATG
        expected = True
        result = codon.is_nonsyn(2, 'T')    # ATT
        self.assertEqual(expected, result)

        codon = self.seq3_codons[3]     # TGA
        expected = False
        result = codon.is_nonsyn(1, 'A')
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = False
        codon = self.seq3_codons[0]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = True
        codon = self.seq3_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq3_codons[2]     # TGG
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq3_codons[2]  # TGG
        exp = True
        result = codon.creates_stop(1, 'A')  # Creates TAG
        self.assertEqual(exp, result)

        exp = False
        result = codon.creates_stop(1, 'C')  # Creates TCG
        self.assertEqual(exp, result)


class TestSequence4(unittest.TestCase):
    """
    Sequence: ATGATGCCCTAA
    ORFs:
        (0, 12) (+)
    Notes:
        Internal methionine
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(5001)

        s4 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(4000)
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, OMEGA_VALUES_5, MU_VALUES_4)

        self.seq4_codons = self.sequence4.find_codons('+0', [(0, 12)])

    def testReverseComplement(self):
        s = str(self.sequence4)
        expected = 'TTAGGGCATCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTACGGGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence4 = self.sequence4.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence4, new_sequence4)
        self.assertEqual(self.sequence4.orfs, new_sequence4.orfs)
        self.assertEqual(self.sequence4.kappa, new_sequence4.kappa)
        self.assertEqual(self.sequence4.global_rate, new_sequence4.global_rate)
        self.assertEqual(self.sequence4.pi, new_sequence4.pi)
        self.assertEqual(self.sequence4.omega_values, new_sequence4.omega_values)
        self.assertEqual(self.sequence4.cat_values, new_sequence4.cat_values)
        self.assertEqual(self.sequence4.is_circular, new_sequence4.is_circular)
        self.assertEqual(self.sequence4.total_omegas, new_sequence4.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence4.event_tree, new_sequence4.event_tree)
        self.assertCountEqual(self.sequence4.event_tree, new_sequence4.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq4_nt in enumerate(self.sequence4.nt_sequence):
            new_nt = new_sequence4.nt_sequence[pos]
            self.assertIsNot(seq4_nt, new_nt)

            self.assertEqual(seq4_nt.state, new_nt.state)
            self.assertEqual(seq4_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq4_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq4_nt.rates, new_nt.rates)
            self.assertEqual(seq4_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq4_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq4_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq4_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq4_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq4_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq4_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.33, 'C': 0.25, 'G': 0.17, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence4))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_4
        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)
        # Only 5 possible omegas because there is only 1 ORF
        exp_total_omegas = {(0, 0, 0, 0, 1, 0): 2.1223322911031457,
                            (0, 0, 0, 1, 0, 0): 1.2308533419685328,
                            (0, 0, 1, 0, 0, 0): 0.842647578837806,
                            (0, 1, 0, 0, 0, 0): 0.5492012532921686,
                            (1, 0, 0, 0, 0, 0): 0.25496553479835254}

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (0, 0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0, 0), 'T': (0, 0, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu3', 'G': 'mu4', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence4.total_omegas)

        # Tests mutation in internal methionine
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.019536686677227796, 'G': 0.5420777234203797, 'T': 0.9142294839271567}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (1, 0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 1, 0), 'T': (0, 0, 0, 0, 1, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu2', 'G': 'mu2', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 1.4758438940247642
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence4.total_omegas)

    def testIsTransv(self):
        nt = self.sequence4.nt_sequence[2]  # G
        result = self.sequence4.is_transv(nt.state, 'G')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence4.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence4.is_transv(nt.state, 'A')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ATG', 'CCC', 'TAA']
        result = self.sequence4.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence4.nt_sequence[3]  # A in internal methione codon
        expected = False
        result = self.sequence4.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence4.nt_sequence[0]  # A in start codon
        expected = True
        result = self.sequence4.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A':
                              {'from_nt': {'A': None,
                                           'T': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(0, 0, 0, 1, 0, 0): 0.20514222366142196}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'C': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423,
                                                                           (0, 0, 0, 0, 0, 1): 0.16666666666666652}}}},
                                           'G': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}}}}}},
                              'T':
                              {'from_nt': {'A': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'T': None,
                                           'C': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {(0, 1, 0, 0, 0, 0): 0.09153354221536135}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {(0, 0, 0, 0, 0, 1): 0.16666666666666652}}}},
                                           'G': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}}}},
                              'C':
                              {'from_nt': {'A': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'T': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {(0, 1, 0, 0, 0, 0): 0.09153354221536135}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'C': None,
                                           'G': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}}}},
                              'G':
                              {'from_nt': {'A': {'prob': 0.625,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'T': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {(0, 1, 0, 0, 0, 0): 0.09153354221536135}}}},
                                           'C': {'prob': 0.18749999999999997,
                                                 'cat': {'mu1': {'prob': 0.04701719357592987,
                                                                 'omega': {}},
                                                         'mu2': {'prob': 0.11163806035536104,
                                                                 'omega': {}},
                                                         'mu3': {'prob': 0.21374220928782836,
                                                                 'omega': {(0, 0, 0, 0, 1, 0): 0.35372204851719063,
                                                                           (0, 0, 0, 0, 0, 1): 0.16666666666666652}},
                                                         'mu4': {'prob': 0.6276025367808807,
                                                                 'omega': {}}}},
                                           'G': None}}}}

        result = self.sequence4.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(5001)

        a3 = self.sequence4.nt_sequence[3]
        t4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        c6 = self.sequence4.nt_sequence[6]
        c7 = self.sequence4.nt_sequence[7]
        c8 = self.sequence4.nt_sequence[8]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 0, 1, 0, 0): [t4]},
                                                                         'mu4': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 1, 0, 0, 0): [c6]},
                                                                         'mu4': {(0, 0, 1, 0, 0, 0): [c7],
                                                                                 (0, 0, 0, 0, 0, 1): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {(0, 0, 1, 0, 0, 0): [g5]}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {(1, 0, 0, 0, 0, 0): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'T': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {(0, 1, 0, 0, 0, 0): [c6]},
                                                                         'mu3': {(0, 0, 1, 0, 0, 0): [c7]},
                                                                         'mu4': {(0, 0, 0, 0, 0, 1): [c8]}}},
                                                      'G': {'category': {'mu1': {(0, 0, 1, 0, 0, 0): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 1, 0, 0, 0): [a3]},
                                                                         'mu4': {}}},
                                                      'T': {'category': {'mu1': {(0, 1, 0, 0, 0, 0): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {(0, 0, 1, 0, 0, 0): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(1, 0, 0, 0, 0, 0): [a3]},
                                                                         'mu4': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {(0, 1, 0, 0, 0, 0): [t4]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 0, 0, 1, 0): [c6, c7],
                                                                                 (0, 0, 0, 0, 0, 1): [c8]},
                                                                         'mu4': {}}},
                                                      'G': None}}}}

        res_omega_key = self.sequence4.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 18
        res_count = self.sequence4.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 10)
        nt.rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        nt.get_mutation_rate()
        exp_mutation_rate = 0
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None,                'C': 0.0,                  'G': 0.0,                 'T': 0.0},
            {'A': 0.0,                 'C': 0.0,                  'G': 0.0,                 'T': None},
            {'A': 0.0,                 'C': 0.0,                  'G': None,                'T': 0.0},
            {'A': None,                'C': 0.12362133301173836,  'G': 0.12468312208406625, 'T': 0.008228019874329523},
            {'A': 0.13679802361738583, 'C': 0.04475588894143996,  'G': 0.17922534693719583, 'T': None},
            {'A': 0.6233063414381459,  'C': 0.014008602512739905, 'G': None,                'T': 0.014008602512739905},
            {'A': 0.0936525250088927,  'C': None,                 'G': 0.23587770612698317, 'T': 0.10626879766512083},
            {'A': 0.2749880918109467,  'C': None,                 'G': 0.23587770612698317, 'T': 0.31217508336297567},
            {'A': 0.32633819726891633, 'C': None,                 'G': 0.11114079878810054, 'T': 1.0877939908963878},
            {'A': 0.0,                 'C': 0.0,                  'G': 0.0,                 'T': None},
            {'A': None,                'C': 0.0,                  'G': 0.0,                 'T': 0.0},
            {'A': None,                'C': 0.0,                  'G': 0.0,                 'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           0.25653247497013415,
                           0.36077925949602163,
                           0.6513235464636258,
                           0.4357990288009967,
                           0.8230408813009056,
                           1.5252729869534047,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence4.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq4_codons[3]         # TAA
        nt = self.sequence4.nt_sequence[9]  # T
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['A', 'T', 'G']
        exp_mutated = ['T', 'T', 'G']
        res_codon, res_mutated = self.seq4_codons[1].mutate_codon(0, 'T')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq4_codons[1]         # ATG (Met)
        expected = True
        result = codon.is_nonsyn(2, 'A')    # ATA (Ile)
        self.assertEqual(expected, result)

        codon = self.seq4_codons[2]         # CCC (Pro)
        expected = False
        result = codon.is_nonsyn(2, 'A')    # CCA (Pro)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.seq4_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq4_codons[2]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq4_codons[1]     # internal methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

        codon = self.seq4_codons[0]     # First methionine in ORF
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        exp = True
        codon = self.seq4_codons[3]
        result = codon.creates_stop(2, 'G')  # TAA to TAG (STOP to STOP)
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq4_codons[1]
        result = codon.creates_stop(2, 'A')  # Creates ATA
        self.assertEqual(exp, result)


class TestSequence5(unittest.TestCase):
    """
    Sequence: ATGAATGCCTGACTAA
    ORFs:
        (0, 12) (+)
        (4, 16) (+)
    Notes:
        2 overlapping ORFs in the forward strand
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9991)

        s5 = 'ATGAATGCCTGACTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [[(4, 16)]], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, OMEGA_VALUES_4, MU_VALUES_4)

        self.plus_0_codons = self.sequence5.find_codons('+0', [(0, 12)])
        self.plus_1_codons = self.sequence5.find_codons('+1', [(4, 16)])

    def testReverseComplement(self):
        s = str(self.sequence5)
        expected = 'TTAGTCAGGCATTCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTTACGGACTGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence5 = self.sequence5.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence5, new_sequence5)
        self.assertEqual(self.sequence5.orfs, new_sequence5.orfs)
        self.assertEqual(self.sequence5.kappa, new_sequence5.kappa)
        self.assertEqual(self.sequence5.global_rate, new_sequence5.global_rate)
        self.assertEqual(self.sequence5.pi, new_sequence5.pi)
        self.assertEqual(self.sequence5.omega_values, new_sequence5.omega_values)
        self.assertEqual(self.sequence5.cat_values, new_sequence5.cat_values)
        self.assertEqual(self.sequence5.is_circular, new_sequence5.is_circular)
        self.assertEqual(self.sequence5.total_omegas, new_sequence5.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence5.event_tree, new_sequence5.event_tree)
        self.assertCountEqual(self.sequence5.event_tree, new_sequence5.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq5_nt in enumerate(self.sequence5.nt_sequence):
            new_nt = new_sequence5.nt_sequence[pos]
            self.assertIsNot(seq5_nt, new_nt)

            self.assertEqual(seq5_nt.state, new_nt.state)
            self.assertEqual(seq5_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq5_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq5_nt.rates, new_nt.rates)
            self.assertEqual(seq5_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq5_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq5_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq5_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq5_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq5_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq5_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.38, 'C': 0.19, 'G': 0.19, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence5))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_4
        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9991)

        exp_all_omegas = {(0, 0, 1, 0, 0): 1.0699896623909886,
                          (0, 0, 1, 0, 1): 1.0699896623909886,
                          (0, 1, 0, 0, 0): 0.6550136761581515,
                          (0, 1, 0, 1, 0): 1.298054976532153,
                          (1, 0, 0, 0, 0): 0.29327471612351436,
                          (1, 0, 0, 0, 1): 0.29327471612351436,
                          (1, 0, 1, 0, 0): 0.31380091449281217}

        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0), 'T': (0, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0, 0, 0, 0, 0), 'C': (0, 0, 0, 0, 0), 'G': (0, 0, 0, 0, 0), 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'G': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests a mutation that is synonymous in one frame and non-synonymous in the other
        random.seed(9991)
        nt = self.sequence5.nt_sequence[8]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.00544914215757548, 'C': None, 'G': 0.07273712404486865, 'T': 0.12273685164103774}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (1, 0, 0, 0, 1), 'C': None, 'G': (1, 0, 0, 0, 1), 'T': (0, 0, 0, 1, 1)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'G': 'mu4', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.20092311784348188
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_all_omegas = {(0, 0, 0, 1, 1): 1.9817219453273531,
                          (0, 0, 1, 0, 0): 1.0699896623909886,
                          (0, 0, 1, 0, 1): 1.0699896623909886,
                          (0, 1, 0, 0, 0): 0.6550136761581515,
                          (0, 1, 0, 1, 0): 1.298054976532153,
                          (1, 0, 0, 0, 0): 0.29327471612351436,
                          (1, 0, 0, 0, 1): 0.29327471612351436,
                          (1, 0, 1, 0, 0): 0.31380091449281217}
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests a nucleotide involved in multiple non-synonymous codons
        random.seed(9991)
        nt = self.sequence5.nt_sequence[7]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0015980956193796234, 'C': None, 'G': 0.0729690937192296, 'T': 0.10306656193502935}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (2, 0, 0, 0, 0), 'C': None, 'G': (0, 0, 0, 2, 0), 'T': (0, 1, 1, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.17763375127363856
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_all_omegas = {(0, 0, 0, 1, 1): 1.9817219453273531,
                          (0, 0, 0, 2, 0): 3.927221868592029,
                          (0, 0, 1, 0, 0): 1.0699896623909886,
                          (0, 0, 1, 0, 1): 1.0699896623909886,
                          (0, 1, 0, 0, 0): 0.6550136761581515,
                          (0, 1, 0, 1, 0): 1.298054976532153,
                          (0, 1, 1, 0, 0): 0.7008578622139409,
                          (1, 0, 0, 0, 0): 0.29327471612351436,
                          (1, 0, 0, 0, 1): 0.29327471612351436,
                          (1, 0, 1, 0, 0): 0.31380091449281217,
                          (2, 0, 0, 0, 0): 0.08601005911732792}
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

    def testIsTransv(self):
        nt = self.sequence5.nt_sequence[4]  # A
        result = self.sequence5.is_transv(nt.state, 'A')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence5.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence5.is_transv(nt.state, 'G')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        # Check codons in +0 frame
        expected = ['ATG', 'AAT', 'GCC', 'TGA']
        result = self.sequence5.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Check Codons in +1 frame
        expected = ['ATG', 'CCT', 'GAC', 'TAA']
        result = self.sequence5.find_codons('+1', [(4, 16)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+1')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence5.nt_sequence[4]  # A in START codon in ORF (4, 16)
        expected = True
        result = self.sequence5.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

        nt = self.sequence5.nt_sequence[15]     # last A in STOP codon in ORF (4, 16)
        expected = True
        result = self.sequence5.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence5.nt_sequence[3]  # First codon in AAT in ORF (0, 12)
        expected = False
        result = self.sequence5.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'from_nt':
                                    {'A': None,
                                     'T': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'C': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {(0, 0, 1, 0, 1): 0.029787448030113194,
                                                                     (0, 1, 0, 0, 0): 0.10928919466474106}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {(0, 1, 0, 1, 0): 0.21658079545427042}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'G': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}}}},
                              'T': {'from_nt':
                                    {'A': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {(0, 1, 0, 0, 0): 0.10928919466474106}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'T': None,
                                     'C': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {(1, 0, 1, 0, 0): 0.052357760575518504}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {(0, 0, 1, 0, 1): 0.029787448030113194}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {(0, 0, 0, 0, 1): 0.16685024854099906}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'G': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}}}},
                              'C': {'from_nt':
                                    {'A': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {(1, 0, 0, 0, 0): 0.04893295927599932}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'T': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'C': None,
                                     'G': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}}}},
                              'G': {'from_nt':
                                    {'A': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {(0, 1, 0, 0, 0): 0.10928919466474106}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'T': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'C': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.04701719357592987,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.11163806035536104,
                                                           'omega': {(1, 0, 0, 0, 1): 0.00816447641704707,
                                                                     (0, 0, 1, 0, 0): 0.17852804110623613}},
                                                   'mu3': {'prob': 0.21374220928782836,
                                                           'omega': {(1, 0, 1, 0, 0): 0.052357760575518504}},
                                                   'mu4': {'prob': 0.6276025367808807,
                                                           'omega': {}}}},
                                     'G': None}}}}

        result = self.sequence5.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(9991)

        # Only 4 possible substitution sites because the
        a3 = self.sequence5.nt_sequence[3]
        c7 = self.sequence5.nt_sequence[7]
        c8 = self.sequence5.nt_sequence[8]
        c12 = self.sequence5.nt_sequence[12]
        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'C': {'category': {'mu1': {(0, 0, 1, 0, 1): [c8],
                                                                                 (0, 1, 0, 0, 0): [c12]},
                                                                         'mu2': {(0, 1, 0, 1, 0): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 1, 0, 0, 0): [a3]},
                                                                         'mu4': {}}},
                                                      'T': None,
                                                      'C': {'category': {'mu1': {(1, 0, 1, 0, 0): [c7]},
                                                                         'mu2': {(0, 0, 1, 0, 1): [c8]},
                                                                         'mu3': {(0, 0, 0, 0, 1): [c12]},
                                                                         'mu4': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(1, 0, 0, 0, 0): [a3]},
                                                                         'mu4': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 1, 0, 0, 0): [a3]},
                                                                         'mu4': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {(1, 0, 0, 0, 1): [c8],
                                                                                 (0, 0, 1, 0, 0): [c12]},
                                                                         'mu3': {(1, 0, 1, 0, 0): [c7]},
                                                                         'mu4': {}}},
                                                      'G': None}}}}

        res_omega_key = self.sequence5.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 12
        res_count = self.sequence5.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 3)
        nt.rates = {'A': None, 'C': 0.04954407504576764, 'G': 0.3688469654724257, 'T': 0.1106540896417277}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.529045130159921
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None,                 'C': 0.0,                 'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                 'G': None,                 'T': 0.0},
            {'A': None,                 'C': 0.04954407504576764, 'G': 0.3688469654724257,   'T': 0.1106540896417277},
            {'A': None,                 'C': 0.0,                 'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                 'G': None,                 'T': 0.0},
            {'A': 0.05726670307066566,  'C': None,                'G': 0.026505824065847356, 'T': 0.019435085925015638},
            {'A': 0.01988076522440186,  'C': None,                'G': 0.012938493661684328, 'T': 0.15735024425679955},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                 'G': None,                 'T': 0.0},
            {'A': None,                 'C': 0.0,                 'G': 0.0,                  'T': 0.0},
            {'A': 0.012170372828998532, 'C': None,                'G': 0.04720507327703986,  'T': 0.2815566902631881},
            {'A': 0.0,                  'C': 0.0,                 'G': 0.0,                  'T': None},
            {'A': None,                 'C': 0.0,                 'G': 0.0,                  'T': 0.0},
            {'A': None,                 'C': 0.0,                 'G': 0.0,                  'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           0.529045130159921,
                           0,
                           0,
                           0,
                           0.10320761306152866,
                           0.19016950314288575,
                           0,
                           0,
                           0,
                           0.3409321363692265,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence5.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        plus_0_codon = self.plus_0_codons[2]    # GCC
        nt = self.sequence5.nt_sequence[7]      # First C
        expected = 1
        result = plus_0_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        plus_1_codon = self.plus_1_codons[1]    # CCT
        expected = 0
        result = plus_1_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['G', 'C', 'C']
        exp_mutated = ['G', 'C', 'G']
        res_codon, res_mutated = self.plus_0_codons[2].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        exp_codon = ['G', 'A', 'C']
        exp_mutated = ['G', 'A', 'G']
        res_codon, res_mutated = self.plus_1_codons[2].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        plus_0_codon = self.plus_0_codons[1]     # AAT (Asn)
        expected = True
        result = plus_0_codon.is_nonsyn(0, 'C')  # CAT (His)
        self.assertEqual(expected, result)

        plus_1_codon = self.plus_1_codons[1]     # CCT (Pro)
        expected = False
        result = plus_1_codon.is_nonsyn(2, 'C')  # CCT (Pro)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.plus_1_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.plus_0_codons[0]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        plus_1_codon = self.plus_1_codons[0]    # ATG
        expected = True
        result = plus_1_codon.is_start()
        self.assertEqual(expected, result)

        plus_0_codon = self.plus_0_codons[3]    # TGA
        expected = False
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.plus_0_codons[3]
        exp = True
        result = codon.creates_stop(1, 'A')
        self.assertEqual(exp, result)

        codon = self.plus_1_codons[2]
        exp = False
        result = codon.creates_stop(0, 'T')
        self.assertEqual(exp, result)


class TestSequence6(unittest.TestCase):
    """
    Sequence: ATGATGGCCCTAA
    ORFs:
        [(0, 5), (6, 13)] (+)
    Notes:
        Spliced ORF (5th nucleotide is excluded)
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4000)

        s6 = 'ATGATGGCCCTAA'
        sorted_orfs = {'+0': [[(0, 5), (6, 13)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, OMEGA_VALUES_5, MU_VALUES_3)
        self.seq6_codons = self.sequence6.find_codons('+0', [(0, 5), (6, 13)])

    def testReverseComplement(self):
        s = str(self.sequence6)
        expected = 'TTAGGGCCATCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTACCGGGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence6 = self.sequence6.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence6, new_sequence6)
        self.assertEqual(self.sequence6.orfs, new_sequence6.orfs)
        self.assertEqual(self.sequence6.kappa, new_sequence6.kappa)
        self.assertEqual(self.sequence6.global_rate, new_sequence6.global_rate)
        self.assertEqual(self.sequence6.pi, new_sequence6.pi)
        self.assertEqual(self.sequence6.omega_values, new_sequence6.omega_values)
        self.assertEqual(self.sequence6.cat_values, new_sequence6.cat_values)
        self.assertEqual(self.sequence6.is_circular, new_sequence6.is_circular)
        self.assertEqual(self.sequence6.total_omegas, new_sequence6.total_omegas)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence6.event_tree, new_sequence6.event_tree)
        self.assertCountEqual(self.sequence6.event_tree, new_sequence6.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq6_nt in enumerate(self.sequence6.nt_sequence):
            new_nt = new_sequence6.nt_sequence[pos]
            self.assertIsNot(seq6_nt, new_nt)
            self.assertEqual(seq6_nt.state, new_nt.state)
            self.assertEqual(seq6_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq6_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq6_nt.rates, new_nt.rates)
            self.assertEqual(seq6_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq6_nt.cat_keys, new_nt.cat_keys)
            self.assertEqual(seq6_nt.omega_in_event_tree, new_nt.omega_in_event_tree)
            self.assertEqual(seq6_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(str(seq6_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq6_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq6_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.31, 'C': 0.23, 'G': 0.23, 'T': 0.23}
        result = Sequence.get_frequency_rates(str(self.sequence6))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE_MU_3
        result = self.sequence6.create_event_tree()
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        # Only 5 possible omegas because there is only 1 ORF
        exp_total_omegas = {(0, 0, 0, 0, 1, 0): 2.1223322911031457,
                            (0, 0, 0, 1, 0, 0): 1.2308533419685328,
                            (0, 0, 1, 0, 0, 0): 0.842647578837806,
                            (0, 1, 0, 0, 0, 0): 0.5492012532921686,
                            (1, 0, 0, 0, 0, 0): 0.25496553479835254}

        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': 0.2491527517379426, 'C': 0.07474582552138279, 'G': None, 'T': 0.07474582552138279}
        self.assertEqual(exp_sub_rates, nt.rates)
        # No omega keys because nt treated as synonymous mutation
        exp_omega_keys = {'A': (0, 0, 0, 0, 0, 0), 'C': (0, 0, 0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu2', 'C': 'mu2', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.3986444027807081
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

        # Tests internal methionine
        random.seed(4000)
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 0.12400154884247461, 'G': 1.4194135069676665, 'T': 0.1900007865404335}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (0, 0, 0, 1, 0, 0), 'G': (0, 0, 0, 1, 0, 0), 'T': (0, 1, 0, 0, 0, 0)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu3'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 1.7334158423505746
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

    def testIsTransv(self):
        nt = self.sequence6.nt_sequence[6]  # G
        result = self.sequence6.is_transv(nt.state, 'G')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence6.is_transv(nt.state, 'T')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence6.is_transv(nt.state, 'A')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ATG', 'CCC', 'TAA']
        result = self.sequence6.find_codons('+0', [(0, 5), (6, 13)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStopStartCodon(self):
        nt = self.sequence6.nt_sequence[5]
        expected = False
        result = self.sequence6.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence6.nt_sequence[2]
        expected = True
        result = self.sequence6.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'from_nt':
                                    {'A': None,
                                     'T': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'C': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423,
                                                                     (0, 0, 0, 0, 0, 1): 0.16666666666666652}}}},
                                     'G': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {(0, 1, 0, 0, 0, 0): 0.09153354221536135}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {(0, 0, 0, 0, 0, 0): 0.16666666666666652}}}}}},
                              'T': {'from_nt':
                                    {'A': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'T': None,
                                     'C': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385,
                                                                     (0, 0, 0, 0, 0, 1): 0.16666666666666652}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'G': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(0, 0, 0, 0, 0, 0): 0.16666666666666652,
                                                                     (0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}}}},
                              'C': {'from_nt':
                                    {'A': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(0, 0, 1, 0, 0, 0): 0.14044126313963423}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'T': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'C': None,
                                     'G': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(0, 1, 0, 0, 0, 0): 0.09153354221536135}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {(0, 0, 0, 0, 0, 0): 0.16666666666666652}}}}}},
                              'G': {'from_nt':
                                    {'A': {'prob': 0.625,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {}}}},
                                     'T': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {(0, 0, 0, 1, 0, 0): 0.20514222366142196}}}},
                                     'C': {'prob': 0.18749999999999997,
                                           'cat': {'mu1': {'prob': 0.07625419304108505,
                                                           'omega': {(0, 0, 0, 0, 0, 1): 0.16666666666666652}},
                                                   'mu2': {'prob': 0.20833135198746808,
                                                           'omega': {(1, 0, 0, 0, 0, 0): 0.042494255799725385}},
                                                   'mu3': {'prob': 0.7154144549714468,
                                                           'omega': {(0, 0, 0, 0, 1, 0): 0.35372204851719063}}}},
                                     'G': None}}}}
        result = self.sequence6.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(4000)

        a3 = self.sequence6.nt_sequence[3]
        t4 = self.sequence6.nt_sequence[4]
        g5 = self.sequence6.nt_sequence[5]
        g6 = self.sequence6.nt_sequence[6]
        c7 = self.sequence6.nt_sequence[7]
        c8 = self.sequence6.nt_sequence[8]
        c9 = self.sequence6.nt_sequence[9]
        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {(1, 0, 0, 0, 0, 0): [t4]},
                                                                         'mu3': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 1, 0, 0, 0): [c7, c8],
                                                                                 (0, 0, 0, 0, 0, 1): [c9]}}},
                                                      'G': {'category': {'mu1': {(0, 1, 0, 0, 0, 0): [g6]},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 0, 0, 0, 0): [g5]}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {(0, 0, 1, 0, 0, 0): [a3]},
                                                                         'mu3': {}}},
                                                      'T': None,
                                                      'C': {'category': {'mu1': {(0, 0, 1, 0, 0, 0): [c8]},
                                                                         'mu2': {(1, 0, 0, 0, 0, 0): [c7],
                                                                                 (0, 0, 0, 0, 0, 1): [c9]},
                                                                         'mu3': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {(0, 0, 0, 0, 0, 0): [g5],
                                                                                 (0, 0, 1, 0, 0, 0): [g6]},
                                                                         'mu3': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {(0, 0, 1, 0, 0, 0): [a3]},
                                                                         'mu3': {}}},
                                                      'T': {'category': {'mu1': {(1, 0, 0, 0, 0, 0): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {(0, 1, 0, 0, 0, 0): [g6]},
                                                                         'mu3': {(0, 0, 0, 0, 0, 0): [g5]}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {(1, 0, 0, 0, 0, 0): [a3]},
                                                                         'mu3': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {(0, 0, 0, 1, 0, 0): [t4]}}},
                                                      'C': {'category': {'mu1': {(0, 0, 0, 0, 0, 1): [c9]},
                                                                         'mu2': {(1, 0, 0, 0, 0, 0): [c7]},
                                                                         'mu3': {(0, 0, 0, 0, 1, 0): [c8]}}},
                                                      'G': None}}}}

        res_omega_key = self.sequence6.nt_in_event_tree(g5)
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

    def testCountNtsOnEventTree(self):
        exp_count = 21
        res_count = self.sequence6.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('C', 10)
        nt.rates = {'A': 0.25667881246211854, 'C': None, 'G': 0.02735873671412473, 'T': 0.2491527517379426}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.5331903009141858
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': 0.0,                  'C': 0.0,                  'G': None,                 'T': 0.0},
            {'A': None,                 'C': 0.08489200243559196,  'G': 0.08562114358233551,  'T': 0.08489200243559196},
            {'A': 0.019057609378003708, 'C': 0.023251783125747112, 'G': 0.3159339741315129,   'T': None},
            {'A': 0.8555960415403951,   'C': 0.25667881246211854,  'G': None,                 'T': 0.07474582552138279},
            {'A': 0.05008484163962589,  'C': 0.04105050105470119,  'G': None,                 'T': 0.0629843889038263},
            {'A': 0.21628977986016748,  'C': None,                 'G': 0.019057609378003708, 'T': 0.06352536459334569},
            {'A': 0.21628977986016748,  'C': None,                 'G': 0.5447577321303627,   'T': 0.07684591084072731},
            {'A': 0.25667881246211854,  'C': None,                 'G': 0.02735873671412473,  'T': 0.2491527517379426},
            {'A': 0.0,                  'C': 0.0,                  'G': 0.0,                  'T': None},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0},
            {'A': None,                 'C': 0.0,                  'G': 0.0,                  'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           0.2554051484535194,
                           0.35824336663526374,
                           1.1870206795238967,
                           0.15411973159815337,
                           0.2988727538315169,
                           0.8378934228312576,
                           0.5331903009141858,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence6.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq6_codons[0]         # ATG
        nt = self.sequence6.nt_sequence[0]  # A
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codon = self.seq6_codons[3]          # TAA
        nt = self.sequence6.nt_sequence[12]  # A
        expected = 2
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['A', 'T', 'G']
        exp_mutated = ['C', 'T', 'G']
        res_codon, res_mutated = self.seq6_codons[1].mutate_codon(0, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq6_codons[3]             # TAA (STOP)
        expected = True
        result = codon.is_nonsyn(1, 'C')        # TAC (Tyr)
        self.assertEqual(expected, result)

        expected = False
        result = codon.is_nonsyn(1, 'G')        # TGA (STOP)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.seq6_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq6_codons[2]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq6_codons[0]     # First methionine
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

        codon = self.seq6_codons[1]     # Internal methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq6_codons[3]
        exp = True
        result = codon.creates_stop(2, 'G')
        self.assertEqual(exp, result)

        exp = False
        result = codon.creates_stop(0, 'G')
        self.assertEqual(exp, result)


if __name__ == '__main__':
    unittest.main()
