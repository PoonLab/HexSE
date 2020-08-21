import random
import unittest

from src.sequence_info import Sequence

MAX_DIFF = None
KAPPA = 0.3
GLOBAL_RATE = 1

OMEGA_VALUES_4 = {'omega1':  0.29327471612351436,
                  'omega2': 0.6550136761581515,
                  'omega3': 1.0699896623909886,
                  'omega4': 1.9817219453273531}

OMEGA_VALUES_5 = {'omega1': 0.25496553479835254,
                  'omega2': 0.5492012532921686,
                  'omega3': 0.842647578837806,
                  'omega4': 1.2308533419685328,
                  'omega5': 2.1223322911031457}

MU_VALUES_4 = {'mu1': 0.325970770628452,
               'mu2': 0.7739880200789496,
               'mu3': 1.481877317174674,
               'mu4': 4.351175963585551}

MU_VALUES_3 = {'mu1': 0.3965034306394888,
               'mu2': 1.0832728336432287,
               'mu3': 3.7199827893060657}

EMPTY_EVENT_TREE_MU_4 = {'to_nt': {'A': {'from_nt': {'A': None,
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}}}},
                                   'C': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'C': None,
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}}}},
                                   'G': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'G': None,
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}}}},
                                   'T': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []},
                                                                        'mu4': {'syn_mutations': []}}},
                                                     'T': None}}}}

EMPTY_EVENT_TREE_MU_3 = {'to_nt': {'A': {'from_nt': {'A': None,
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}}}},
                                   'C': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'C': None,
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}}}},
                                   'G': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'G': None,
                                                     'T': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}}}},
                                   'T': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'C': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
                                                     'G': {'category': {'mu1': {'syn_mutations': []},
                                                                        'mu2': {'syn_mutations': []},
                                                                        'mu3': {'syn_mutations': []}}},
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

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.12603317336204048, 'C': 0.030344325868823886,
                         'G': None,                'T': 0.030344325868823886}
        exp_omegas = {'A': ['omega1'], 'C': ['omega3'], 'G': None, 'T': ['omega3']}
        exp_cat_keys = {'A': 'mu3', 'C': 'mu1', 'T': 'mu1'}
        exp_total_rate = 0.18672182509968824
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301,
                                                            'omega1': 0.058654943224702785}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974,
                                                            'omega2': 0.1310027352316301}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004,
                                                            'omega3': 0.2139979324781974}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974,
                                                            'omega4': 0.39634438906547004,
                                                            'omega2': 0.1310027352316301}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004,
                                                            'omega2': 0.1310027352316301}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
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
                                  'T': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega4': [t6]},
                                                     'mu2': {'syn_mutations': [t17],
                                                             'omega4': [t14]},
                                                     'mu3': {'syn_mutations': [],
                                                             'omega3': [t1]},
                                                     'mu4': {'syn_mutations': [],
                                                             'omega4': [t10]}}},
                                  'C': {'category': {'mu1': {'syn_mutations': [c11]},
                                                     'mu2': {'syn_mutations': [c3],
                                                             'omega2': [c16],
                                                             'omega1': [c20]},
                                                     'mu3': {'syn_mutations': []},
                                                     'mu4': {'syn_mutations': []}}},
                                  'G': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega2': [g15]},
                                                     'mu2': {'syn_mutations': [],
                                                             'omega2': [g4, g19]},
                                                     'mu3': {'syn_mutations': [g8],
                                                             'omega1': [g0]},
                                                     'mu4': {'syn_mutations': [],
                                                             'omega3': [g12]}}}}},
                           'T': {'from_nt':
                                 {'A': {'category': {'mu1': {'syn_mutations': [a5]},
                                                     'mu2': {'syn_mutations': [a2],
                                                             'omega2': [a13]},
                                                     'mu3': {'syn_mutations': [],
                                                             'omega1': [a9]},
                                                     'mu4': {'syn_mutations': [],
                                                             'omega4': [a18]}}},
                                  'T': None,
                                  'C': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega1': [c16]},
                                                     'mu2': {'syn_mutations': [c20],
                                                             'omega1': [c7]},
                                                     'mu3': {'syn_mutations': [c11]},
                                                     'mu4': {'syn_mutations': []}}},
                                  'G': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega3': [g0]},
                                                     'mu2': {'syn_mutations': [],
                                                             'omega3': [g4]},
                                                     'mu3': {'syn_mutations': [g8],
                                                             'omega2': [g12, g15],
                                                             'omega1': [g19]},
                                                     'mu4': {'syn_mutations': []}}}}},
                           'C': {'from_nt':
                                 {'A': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega3': [a9]},
                                                     'mu2': {'syn_mutations': []},
                                                     'mu3': {'syn_mutations': [a5],
                                                             'omega2': [a18]},
                                                     'mu4': {'syn_mutations': [a2],
                                                             'omega1': [a13]}}},
                                  'T': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega3': [t1],
                                                             'omega2': [t6]},
                                                     'mu2': {'syn_mutations': [t14],
                                                             'omega1': [t10]},
                                                     'mu3': {'syn_mutations': []},
                                                     'mu4': {'syn_mutations': [t17]}}},
                                  'C': None,
                                  'G': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega3': [g0]},
                                                     'mu2': {'syn_mutations': []},
                                                     'mu3': {'syn_mutations': [g8],
                                                             'omega2': [g12]},
                                                     'mu4': {'syn_mutations': [],
                                                             'omega4': [g4, g19],
                                                             'omega3': [g15]}}}}},
                           'G': {'from_nt':
                                 {'A': {'category': {'mu1': {'syn_mutations': []},
                                                     'mu2': {'syn_mutations': []},
                                                     'mu3': {'syn_mutations': []},
                                                     'mu4': {'syn_mutations': [a2, a5],
                                                             'omega3': [a9],
                                                             'omega4': [a13],
                                                             'omega2': [a18]}}},
                                  'T': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega4': [t10]},
                                                     'mu2': {'syn_mutations': [],
                                                             'omega4': [t1],
                                                             'omega2': [t14]},
                                                     'mu3': {'syn_mutations': [t17],
                                                             'omega4': [t6]},
                                                     'mu4': {'syn_mutations': []}}},
                                  'C': {'category': {'mu1': {'syn_mutations': [],
                                                             'omega2': [c20]},
                                                     'mu2': {'syn_mutations': [],
                                                             'omega4': [c3],
                                                             'omega1': [c11, c16]},
                                                     'mu3': {'syn_mutations': [],
                                                             'omega4': [c7]},
                                                     'mu4': {'syn_mutations': []}}},
                                  'G': None}}}}

        res_omega_key = self.sequence1.nt_in_event_tree(g0)
        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        # Only 4 possible omegas because there is only 1 ORF
        self.assertEqual(OMEGA_VALUES_4, self.sequence1.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 61
        res_count = self.sequence1.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


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

        # Synonymous mutation in (+) strand, non-synonymous mutation in (-) frame
        nt = self.sequence2.nt_sequence[11]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.2335164156624226, 'C': None, 'G': 0.030341674692605226, 'T': 0.10113891564201742}
        exp_omega_keys = {'A': ['omega4'], 'C': None, 'G': ['omega2'], 'T': ['omega2']}
        exp_cat_keys = {'A': 'mu3', 'G': 'mu2', 'T': 'mu2'}
        exp_total_rate = 0.36499700599704527
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Mutation would destroy a START codon
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omega_keys = {'A': None, 'C': [], 'G': [], 'T': []}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu1'}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Mutation in STOP codon on reverse strand ORF
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omega_keys = {'A': None, 'C': [], 'G': [], 'T': []}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu1'}
        exp_total_rate = 0.0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Non-synonymous mutation in both (+) and (-) strands
        random.seed(4001)
        nt = self.sequence2.nt_sequence[10]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.037346151696315194, 'C': None,
                         'G': 0.037346151696315194, 'T': 0.17608224296983835}
        exp_omega_keys = {'A': ['omega4', 'omega2'],  'C': None,
                          'G': ['omega2', 'omega4'],  'T': ['omega4', 'omega5']}
        exp_cat_keys = {'A': 'mu2', 'G': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.25077454636246876
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
                            'T':
                                {'prob': 0.18749999999999997,
                                 'cat': {'mu1': {'prob': 0.07625419304108505,
                                                 'omega': {'syn_mutations': 0.08764083251715918}},
                                         'mu2': {'prob': 0.20833135198746808,
                                                 'omega': {'syn_mutations': 0.08764083251715918}},
                                         'mu3': {'prob': 0.7154144549714468,
                                                 'omega': {'syn_mutations': 0.08764083251715918}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.08764083251715918}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.08764083251715918,
                                                            'omega1': 0.022345391732910336}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.08764083251715918,
                                                            'omega2_omega5': 0.10215306361964918,
                                                            'omega2': 0.048132455057992865}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.08764083251715918}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.08764083251715918}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.08764083251715918}}}}}},
                     'T': {'from_nt': {'A': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega2_omega5': 0.10215306361964918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega1_omega1': 0.0056973047534601695}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}},
                                       'T': None,
                                       'C': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega5': 0.1860029688703295}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega1_omega4': 0.02750390009204871,
                                                                       'omega1': 0.022345391732910336}}}},
                                       'G': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}}}},
                     'C': {'from_nt': {'A': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega2_omega3': 0.0405586967181372}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega1_omega5': 0.04742434643210488}}}},
                                       'T': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}},
                                       'C': None,
                                       'G': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}}}},
                     'G': {'from_nt': {'A': {'prob': 0.625,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega3_omega5': 0.15673495137522697}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega1': 0.022345391732910336,
                                                                       'omega3': 0.07385033532791384}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}},
                                       'T': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918}}}},
                                       'C': {'prob': 0.18749999999999997,
                                             'cat': {'mu1': {'prob': 0.07625419304108505,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega2_omega4': 0.05924399316528073,
                                                                       'omega3': 0.07385033532791384}},
                                                     'mu2': {'prob': 0.20833135198746808,
                                                             'omega': {'syn_mutations': 0.08764083251715918}},
                                                     'mu3': {'prob': 0.7154144549714468,
                                                             'omega': {'syn_mutations': 0.08764083251715918,
                                                                       'omega2': 0.048132455057992865}}}},
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
                          {'A':
                           {'from_nt': {'A': None,
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega1': [c9],
                                                                   'syn_mutations': []},
                                                           'mu3': {'omega2': [c11],
                                                                   'omega2_omega5': [c10],
                                                                   'syn_mutations': []}}},
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}}}},
                           'C':
                           {'from_nt': {'A': {'category': {'mu1': {'omega2_omega3': [a6, a8],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega1_omega5': [a7],
                                                                   'syn_mutations': []}}},
                                        'C': None,
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}}}},
                           'G': {'from_nt': {'A': {'category': {'mu1': {'omega5_omega3': [a7],
                                                                        'syn_mutations': []},
                                                                'mu2': {'omega1': [a6],
                                                                        'omega3': [a8],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'omega2_omega4': [c10],
                                                                        'omega3': [c11],
                                                                        'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'omega2': [c9],
                                                                        'syn_mutations': []}}},
                                             'G': None,
                                             'T': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}}}},
                           'T': {'from_nt': {'A': {'category': {'mu1': {'omega2_omega5': [a7],
                                                                        'syn_mutations': []},
                                                                'mu2': {'omega1_omega1': [a8],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega5': [c9],
                                                                        'syn_mutations': []},
                                                                'mu3': {'omega1': [c11],
                                                                        'omega4_omega1': [c10],
                                                                        'syn_mutations': []}}},
                                             'G': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'T': None}}}}

        res_omega_key = self.sequence2.nt_in_event_tree(c11)
        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        exp_total_omegas = {'omega1': 0.25496553479835254,
                            'omega1_omega1': 0.06500742393500993,
                            'omega1_omega5': 0.5411215876209263,
                            'omega2': 0.5492012532921686,
                            'omega2_omega3': 0.4627831063813346,
                            'omega2_omega4': 0.6759861980279724,
                            'omega2_omega5': 1.1655875541762872,
                            'omega3': 0.842647578837806,
                            'omega4_omega1': 0.31382518059334646,
                            'omega5': 2.1223322911031457,
                            'omega5_omega3': 1.7883781665873595}
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 17
        res_count = self.sequence2.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


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
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.053216889078518174, 'C': 0.06492879242812255, 'G': 0.008721332329709673, 'T': None}
        exp_omega_keys = {'A': ['omega2'], 'C': ['omega2'], 'G': ['omega1'], 'T': None}
        exp_cat_keys = {'A': 'mu2', 'C': 'mu1', 'G': 'mu1'}
        exp_total_rate = 0.1268670138363504

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation in TGG at position 2
        random.seed(555)
        nt = self.sequence3.nt_sequence[8]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.053456076057738736, 'G': None, 'T': 0.307016589860621}
        exp_omega_keys = {'A': [], 'C': ['omega3'], 'G': None, 'T': ['omega2']}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'T': 'mu3'}
        exp_total_rate = 0.36047266591835975

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.1665314408685853, 'C': 0.13649237703904682, 'G': None, 'T': 0.04995943226057559}
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.3529832501682077

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation that would destroy a STOP codon
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.1665314408685853, 'C': 0.13649237703904682, 'G': None, 'T': 0.04995943226057559}
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.3529832501682077

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega4': 0.39634438906547004}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega2': 0.1310027352316301}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega3': 0.2139979324781974}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.1999999999999997}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.1999999999999997,
                                                            'omega1': 0.058654943224702785}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.1999999999999997}}}},
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

        exp_event_tree = {'to_nt':
                          {'A':
                           {'from_nt': {'A': None,
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega2': [c4],
                                                                   'syn_mutations': []}}},
                                        'G': {'category': {'mu1': {'syn_mutations': [g5]},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega1': [t6],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}}}},
                           'C':
                           {'from_nt': {'A': {'category': {'mu1': {'omega4': [a3],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}},
                                        'C': None,
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega2': [g8],
                                                                   'syn_mutations': []},
                                                           'mu3': {'omega1': [g7],
                                                                   'syn_mutations': [g5]}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega4': [t6],
                                                                   'syn_mutations': []}}}}},
                           'G': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'omega3': [a3],
                                                                        'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega1': [c4],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'G': None,
                                             'T': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega3': [t6],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}}}},
                           'T': {'from_nt': {'A': {'category': {'mu1': {'omega4': [a3],
                                                                        'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'omega4': [c4],
                                                                        'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'G': {'category': {'mu1': {'syn_mutations': [g5]},
                                                                'mu2': {'omega3': [g8],
                                                                        'syn_mutations': []},
                                                                'mu3': {'omega1': [g7],
                                                                        'syn_mutations': []}}},
                                             'T': None}}}}

        res_omega_key = self.sequence3.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        # Only 4 possible omegas because there is only 1 ORF
        self.assertEqual(OMEGA_VALUES_4, self.sequence3.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 16
        res_count = self.sequence3.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


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

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.3113585720876178, 'G': 0.0, 'T': 0.06849001095359351}
        exp_omega_keys = {'A': None, 'C': ['omega5'], 'G': [], 'T': ['omega5']}
        exp_cat_keys = {'C': 'mu3', 'G': 'mu4', 'T': 'mu1'}
        exp_total_rate = 0.37984858304121133
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutation in internal methionine
        # ATG (TGT, GTG)
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.019536686677227796, 'G': 0.5420777234203797, 'T': 0.9142294839271567}
        exp_omega_keys = {'A': None, 'C': ['omega1'], 'G': ['omega5'], 'T': ['omega5']}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu2', 'T': 'mu4'}
        exp_total_rate = 1.4758438940247642
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega4': 0.20514222366142196}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega2': 0.09153354221536135}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega2': 0.09153354221536135}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega2': 0.09153354221536135}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega5': 0.35372204851719063}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
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

        exp_event_tree = {'to_nt':
                          {'A':
                           {'from_nt': {'A': None,
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega3': [c6],
                                                                   'syn_mutations': []},
                                                           'mu4': {'omega3': [c7],
                                                                   'syn_mutations': [c8]}}},
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'omega3': [g5],
                                                                   'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega4': [t4],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}}}},
                           'C':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega3': [a3],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': None,
                                        'G': {'category': {'mu1': {'omega3': [g5],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'omega2': [t4],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}}}},
                           'G':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega1': [a3],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega5': [c6, c7],
                                                                   'syn_mutations': [c8]},
                                                           'mu4': {'syn_mutations': []}}},
                                        'G': None,
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'omega2': [t4],
                                                                   'syn_mutations': []}}}}},
                           'T':
                           {'from_nt': {'A': {'category': {'mu1': {'omega1': [a3],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega2': [c6],
                                                                   'syn_mutations': []},
                                                           'mu3': {'omega3': [c7],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': [c8]}}},
                                        'G': {'category': {'mu1': {'omega3': [g5],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'T': None}}}}

        res_omega_key = self.sequence4.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        # Only 5 possible omegas because there is only 1 ORF
        self.assertEqual(OMEGA_VALUES_5, self.sequence4.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 18
        res_count = self.sequence4.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


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
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omega_keys = {'A': None, 'C': [], 'G': [], 'T': []}
        exp_cat_keys = {'C': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        exp_omega_keys = {'A': [], 'C': [], 'G': [], 'T': None}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'G': 'mu1'}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega3': 0.23109259363549015,
                                                            'omega2': 0.14146754366940048}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega2_omega4': 0.2803493358412066}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega2': 0.14146754366940048}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega1_omega3': 0.06777361479669504}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega3': 0.23109259363549015}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega1': 0.06334043886484732}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega2': 0.14146754366940048}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.04701719357592987,
                                                  'omega': {'syn_mutations': 0.21597647319236044}},
                                          'mu2': {'prob': 0.11163806035536104,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega1': 0.06334043886484732,
                                                            'omega3': 0.23109259363549015}},
                                          'mu3': {'prob': 0.21374220928782836,
                                                  'omega': {'syn_mutations': 0.21597647319236044,
                                                            'omega1_omega3': 0.06777361479669504}},
                                          'mu4': {'prob': 0.6276025367808807,
                                                  'omega': {'syn_mutations': 0.21597647319236044}}}},
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

        exp_event_tree = {'to_nt':
                          {'A':
                           {'from_nt': {'A': None,
                                        'C': {'category': {'mu1': {'omega2': [c12],
                                                                   'omega3': [c8],
                                                                   'syn_mutations': []},
                                                           'mu2': {'omega4_omega2': [c7],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}}}},
                           'C':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega1': [a3],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': None,
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}}}},
                           'G':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega2': [a3],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega1': [c8],
                                                                   'omega3': [c12],
                                                                   'syn_mutations': []},
                                                           'mu3': {'omega3_omega1': [c7],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'G': None,
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}}}},
                           'T':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega2': [a3],
                                                                   'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'C': {'category': {'mu1': {'omega3_omega1': [c7],
                                                                   'syn_mutations': []},
                                                           'mu2': {'omega3': [c8],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': [c12]},
                                                           'mu4': {'syn_mutations': []}}},
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []},
                                                           'mu4': {'syn_mutations': []}}},
                                        'T': None}}}}

        res_omega_key = self.sequence5.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        exp_all_omegas = {'omega1': 0.29327471612351436,
                          'omega2': 0.6550136761581515,
                          'omega3': 1.0699896623909886,
                          'omega3_omega1': 0.31380091449281217,
                          'omega4_omega2': 1.298054976532153}
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 12
        res_count = self.sequence5.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


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
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.2491527517379426, 'C': 0.07474582552138279, 'G': None, 'T': 0.07474582552138279}
        # No omega keys because nt treated as synonymous mutation
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu2', 'C': 'mu2', 'T': 'mu2'}
        exp_total_rate = 0.3986444027807081
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests internal methionine
        random.seed(4000)
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.12400154884247461, 'G': 1.4194135069676665, 'T': 0.1900007865404335}
        exp_omega_keys = {'A': None, 'C': ['omega4'], 'G': ['omega4'], 'T': ['omega2']}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu3'}
        exp_total_rate = 1.7334158423505746
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

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
        expected = {'to_nt':
                    {'A': {'from_nt':
                           {'A': None,
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}}}},
                            'G': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega2': 0.09153354221536135}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}}}},
                     'T': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': None,
                            'C': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}}}},
                     'C': {'from_nt':
                           {'A': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega3': 0.14044126313963423}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'C': None,
                            'G': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega2': 0.09153354221536135}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}}}},
                     'G': {'from_nt':
                           {'A': {'prob': 0.625,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652}}}},
                            'T': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega4': 0.20514222366142196}}}},
                            'C': {'prob': 0.18749999999999997,
                                  'cat': {'mu1': {'prob': 0.07625419304108505,
                                                  'omega': {'syn_mutations': 0.16666666666666652}},
                                          'mu2': {'prob': 0.20833135198746808,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega1': 0.042494255799725385}},
                                          'mu3': {'prob': 0.7154144549714468,
                                                  'omega': {'syn_mutations': 0.16666666666666652,
                                                            'omega5': 0.35372204851719063}}}},
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

        exp_event_tree = {'to_nt':
                          {'A':
                           {'from_nt': {'A': None,
                                        'C': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'omega3': [c7, c8],
                                                                   'syn_mutations': [c9]}}},
                                        'G': {'category': {'mu1': {'omega2': [g6],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': [g5]}}},
                                        'T': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega1': [t4],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}}}},
                           'C':
                           {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega3': [a3],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}},
                                        'C': None,
                                        'G': {'category': {'mu1': {'syn_mutations': []},
                                                           'mu2': {'omega2': [g6],
                                                                   'syn_mutations': []},
                                                           'mu3': {'syn_mutations': [g5]}}},
                                        'T': {'category': {'mu1': {'omega1': [t4],
                                                                   'syn_mutations': []},
                                                           'mu2': {'syn_mutations': []},
                                                           'mu3': {'syn_mutations': []}}}}},
                           'G': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega1': [a3],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'syn_mutations': [c9]},
                                                                'mu2': {'omega1': [c7],
                                                                        'syn_mutations': []},
                                                                'mu3': {'omega5': [c8],
                                                                        'syn_mutations': []}}},
                                             'G': None,
                                             'T': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'syn_mutations': []},
                                                                'mu3': {'omega4': [t4],
                                                                        'syn_mutations': []}}}}},
                           'T': {'from_nt': {'A': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega3': [a3],
                                                                        'syn_mutations': []},
                                                                'mu3': {'syn_mutations': []}}},
                                             'C': {'category': {'mu1': {'omega3': [c8],
                                                                        'syn_mutations': []},
                                                                'mu2': {'omega1': [c7],
                                                                        'syn_mutations': [c9]},
                                                                'mu3': {'syn_mutations': []}}},
                                             'G': {'category': {'mu1': {'syn_mutations': []},
                                                                'mu2': {'omega3': [g6],
                                                                        'syn_mutations': [g5]},
                                                                'mu3': {'syn_mutations': []}}},
                                             'T': None}}}}

        res_omega_key = self.sequence6.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

        # Only 5 possible omegas because there is only 1 ORF
        self.assertEqual(OMEGA_VALUES_5, self.sequence6.total_omegas)

    def testCountNtsOnEventTree(self):
        exp_count = 21
        res_count = self.sequence6.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)


# ==========================================
# Tests for Codon
# ==========================================
class TestCodon(unittest.TestCase):

    def setUp(self):
        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [[(0, 21)]]}, kappa, mu, pi1, dN_values, dS_values)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [[(0, 12)]], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values)

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.nt_seq5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values)

        s6 = 'ATGAATGCCTGACTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [[(4, 16)]], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        random.seed(4000)
        self.nt_seq6 = Sequence(s6, sorted_orfs, kappa, mu, pi6, dN_values, dS_values)

    def testNtInPos(self):
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]
        nt = self.nt_seq1.nt_sequence[0]
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[1]
        nt = self.nt_seq1.nt_sequence[3]
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence[4]
        expected = 1
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence[5]
        expected = 2
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        # Test Sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        # CAT GCT AGC TAG CTA CGA TCG
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        exp_codon = ['G', 'T', 'A']
        exp_mutated = ['G', 'T', 'C']
        res_codon, res_mutated = codons[0].mutate_codon(2, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        codons = self.nt_seq1.find_codons('-0', [(0, 12)])
        exp_codon = ['G', 'A', 'T']
        exp_mutated = ['T', 'A', 'T']
        res_codon, res_mutated = codons[0].mutate_codon(0, 'A')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        # Test Sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]                    # GTA = Valine

        # Mutation at wobble position
        expected = False                     # Synonymous
        result = codon.is_nonsyn(2, 'C')     # GTC = Valine
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(0, 'C')    # CTA = Leucine
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATG ACG TGG TGA
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[2]                   # TGG = Tryptophan

        # Mutation at second position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(1, 'A')    # TAG = STOP
        self.assertEqual(expected, result)

        # Mutation at wobble position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(2, 'A')    # TGA = STOP
        self.assertEqual(expected, result)

        # Testing mutation at position 2 in ACG
        codon = codons[1]                   # ACG = Threonine
        expected = False                    # Synonymous
        result = codon.is_nonsyn(2, 'T')    # ACT = Threonine
        self.assertEqual(expected, result)

        # Testing mutation at position 1 in TGA
        expected = False                    # Synonymous
        codon = codons[3]                   # TGA = STOP
        result = codon.is_nonsyn(1, 'A')    # TAA = STOP
        self.assertEqual(expected, result)

    def testIsStop(self):
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.introduces_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.introduces_stop(1, 'A')
        self.assertEqual(expected, result)

        # G to A mutation in middle position (TAA = STOP)
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[3]
        expected = True
        result = codon.introduces_stop(1, 'A')
        self.assertEqual(expected, result)

    def testIsStart(self):
        # Testing sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]               # GTA = Valine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATG ACG TGG TGA
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[0]               # first Methionine
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 5
        # ATG ATG CCC TAA
        codons = self.nt_seq5.find_codons('+0', [(0, 12)])
        codon = codons[1]               # second Methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
