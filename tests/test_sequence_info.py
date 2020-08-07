import random
import unittest

from src.sequence_info import Sequence

MAX_DIFF = None
KAPPA = 0.3
MU = 0.0005

OMEGA_VALUES_4 = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 0.0]
OMEGA_VALUES_5 = [0.3075971532102495, 0.5999467116400771, 0.8730832667988639, 1.2223993729945162, 0.0]  # shape = 2.5
NT_CATEGORIES_DICT_4 = {'cat1': {'value': 5.265401029875868e-07},
                        'cat2': {'value': 0.0010780892949207104},
                        'cat3': {'value': 0.09375337615809742},
                        'cat4': {'value': 0.0}}
NT_CATEGORIES_DICT_3 = {'cat1': {'value': 0.0876290510827283},  # shape = 0.6
                        'cat2': {'value': 0.556154048516039},
                        'cat3': {'value': 0.0}}


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
        self.sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, KAPPA, MU, pi1, OMEGA_VALUES_4, NT_CATEGORIES_DICT_4)

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
        self.assertEqual(self.sequence1.mu, new_sequence1.mu)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
        self.assertEqual(self.sequence1.omega_values, new_sequence1.omega_values)
        self.assertEqual(self.sequence1.nt_categories, new_sequence1.nt_categories)
        self.assertEqual(self.sequence1.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence1.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence1.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence1.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
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
        expected = {'to_nt':
                    {'A':
                        {'stationary_frequency': 0.24,
                         'from_nt':
                             {'A': None,
                              'T': {'is_trv': True,
                                    'kappa': 0.3,
                                    'is_nonsyn':
                                        {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                  'cat2': {'value': 0.0010780892949207104},
                                                  'cat3': {'value': 0.09375337615809742},
                                                  'cat4': {'value': 0.0}},
                                         'False': {'cat1': {'value': 5.265401029875868e-07},
                                                   'cat2': {'value': 0.0010780892949207104},
                                                   'cat3': {'value': 0.09375337615809742},
                                                   'cat4': {'value': 0.0}}}},
                              'C': {'is_trv': True,
                                    'kappa': 0.3,
                                    'is_nonsyn':
                                        {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                  'cat2': {'value': 0.0010780892949207104},
                                                  'cat3': {'value': 0.09375337615809742},
                                                  'cat4': {'value': 0.0}},
                                         'False': {'cat1': {'value': 5.265401029875868e-07},
                                                   'cat2': {'value': 0.0010780892949207104},
                                                   'cat3': {'value': 0.09375337615809742},
                                                   'cat4': {'value': 0.0}}}},
                              'G': {'is_trv': False,
                                    'kappa': 1,
                                    'is_nonsyn':
                                        {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                  'cat2': {'value': 0.0010780892949207104},
                                                  'cat3': {'value': 0.09375337615809742},
                                                  'cat4': {'value': 0.0}},
                                         'False': {'cat1': {'value': 5.265401029875868e-07},
                                                   'cat2': {'value': 0.0010780892949207104},
                                                   'cat3': {'value': 0.09375337615809742},
                                                   'cat4': {'value': 0.0}}}}}},
                     'T': {'stationary_frequency': 0.24,
                           'from_nt':
                               {'A': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'T': None,
                                'C': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'G': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}}}},
                     'C': {'stationary_frequency': 0.24,
                           'from_nt':
                               {'A': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'T': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'C': None,
                                'G': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}}}},
                     'G': {'stationary_frequency': 0.29,
                           'from_nt':
                               {'A': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'T': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'C': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn':
                                          {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                    'cat2': {'value': 0.0010780892949207104},
                                                    'cat3': {'value': 0.09375337615809742},
                                                    'cat4': {'value': 0.0}},
                                           'False': {'cat1': {'value': 5.265401029875868e-07},
                                                     'cat2': {'value': 0.0010780892949207104},
                                                     'cat3': {'value': 0.09375337615809742},
                                                     'cat4': {'value': 0.0}}}},
                                'G': None}}}}
        result = self.sequence1.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)
        exp_sub_rates = {'A': 4.252483383790958e-05, 'C': 4.654455031400801e-05,
                         'G': None,                  'T': 4.654455031400801e-05}
        exp_omegas = {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 1, 0)}
        exp_total_rate = 0.0001356139344659256
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
        self.sequence2 = Sequence(s2, sorted_orfs, KAPPA, MU, pi2, OMEGA_VALUES_5, NT_CATEGORIES_DICT_3)

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
        self.assertEqual(self.sequence2.mu, new_sequence2.mu)
        self.assertEqual(self.sequence2.pi, new_sequence2.pi)
        self.assertEqual(self.sequence2.omega_values, new_sequence2.omega_values)
        self.assertEqual(self.sequence2.nt_categories, new_sequence2.nt_categories)
        self.assertEqual(self.sequence2.is_circular, new_sequence2.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence2.event_tree, new_sequence2.event_tree)
        self.assertCountEqual(self.sequence2.event_tree, new_sequence2.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence2.nt_sequence):
            new_nt = new_sequence2.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testCreateEventTree(self):
        expected = \
            {'to_nt':
                 {'A': {'stationary_frequency': 0.44,
                        'from_nt':
                            {'A': None,
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'T': {'stationary_frequency': 0.22,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': None,
                             'C': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'C': {'stationary_frequency': 0.17,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': None,
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'G': {'stationary_frequency': 0.17,
                        'from_nt':
                            {'A': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': None}}}}

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
        exp_sub_rates = {'A': 3.117118401136016e-05,  'C': None,
                         'G': 1.5298641146821965e-05, 'T': 0.00010390394670453388}
        exp_omegas = {'A': (0, 0, 0, 1, 0), 'C': None, 'G': (0, 1, 0, 0, 0), 'T': (0, 0, 0, 1, 0)}
        exp_total_rate = 0.00015037377186271603
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Mutation would destroy a START codon
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Non-synonymous mutation in both (+) and (-) strands
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None,                   'C': 4.8402715953170835e-05,
                         'G': 0.00016134238651056946, 'T': 9.862117498842749e-05}
        exp_omegas = {'A': None, 'C': (0, 1, 0, 1, 0), 'G': (0, 1, 0, 1, 0), 'T': (0, 0, 0, 2, 0)}
        exp_total_rate = 0.0003083662774521678
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, MU, pi3, OMEGA_VALUES_4, NT_CATEGORIES_DICT_3)

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
        self.assertEqual(self.sequence3.mu, new_sequence3.mu)
        self.assertEqual(self.sequence3.pi, new_sequence3.pi)
        self.assertEqual(self.sequence3.omega_values, new_sequence3.omega_values)
        self.assertEqual(self.sequence3.nt_categories, new_sequence3.nt_categories)
        self.assertEqual(self.sequence3.is_circular, new_sequence3.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence3.event_tree, new_sequence3.event_tree)
        self.assertCountEqual(self.sequence3.event_tree, new_sequence3.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence3.nt_sequence):
            new_nt = new_sequence3.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
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
        expected = \
            {'to_nt':
                 {'A': {'stationary_frequency': 0.25,
                        'from_nt':
                            {'A': None,
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'T': {'stationary_frequency': 0.25,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': None,
                             'C': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'C': {'stationary_frequency': 0.08,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': None,
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'G': {'stationary_frequency': 0.42,
                        'from_nt':
                            {'A': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': None}}}}

        result = self.sequence3.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):

        # Tests a non-synonymous mutation
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 2.456301285593068e-05, 'C': 0.00013374870779887357,
                         'G': 2.456301285593068e-05, 'T': None}
        exp_omegas = {'A': (0, 1, 0, 0), 'C': (0, 0, 1, 0), 'G': (0, 1, 0, 0), 'T': None}
        exp_total_rate = 0.00018287473351073495
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation in TGG at position 2
        random.seed(555)
        nt = self.sequence3.nt_sequence[8]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 4.126586159796355e-05, 'G': None, 'T': 6.740934873063228e-05}
        exp_omegas = {'A': None, 'C': (0, 1, 0, 0), 'G': None, 'T': (0, 0, 1, 0)}
        exp_total_rate = 0.00010867521032859582
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.00021, 'C': 6.3e-05, 'G': None, 'T': 6.3e-05}
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0.00033600000000000004
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation that would destroy a STOP codon
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.00021, 'C': 6.3e-05, 'G': None, 'T': 6.3e-05}
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0.00033600000000000004
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, MU, pi4, OMEGA_VALUES_5, NT_CATEGORIES_DICT_4)

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
        self.assertEqual(self.sequence4.mu, new_sequence4.mu)
        self.assertEqual(self.sequence4.pi, new_sequence4.pi)
        self.assertEqual(self.sequence4.omega_values, new_sequence4.omega_values)
        self.assertEqual(self.sequence4.nt_categories, new_sequence4.nt_categories)
        self.assertEqual(self.sequence4.is_circular, new_sequence4.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence4.event_tree, new_sequence4.event_tree)
        self.assertCountEqual(self.sequence4.event_tree, new_sequence4.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence4.nt_sequence):
            new_nt = new_sequence4.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
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
        expected = {'to_nt':
                        {'A': {'stationary_frequency': 0.33,
                               'from_nt':
                                   {'A': None,
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'T': {'stationary_frequency': 0.25,
                               'from_nt':
                                   {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'T': None,
                                    'C': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': {'is_trv': True,
                                          'kappa': 0.3, 'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'C': {'stationary_frequency': 0.25,
                               'from_nt':
                                   {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'T': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': None,
                                    'G': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'G': {'stationary_frequency': 0.17,
                               'from_nt':
                                   {'A': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                                 'cat2': {'value': 0.0010780892949207104},
                                                                 'cat3': {'value': 0.09375337615809742},
                                                                 'cat4': {'value': 0.0}},
                                                        'False': {'cat1': {'value': 5.265401029875868e-07},
                                                                  'cat2': {'value': 0.0010780892949207104},
                                                                  'cat3': {'value': 0.09375337615809742},
                                                                  'cat4': {'value': 0.0}}}},
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': None}}}}

        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 4.321762170654376e-05}
        exp_omegas = {'A': None, 'C': (0, 0, 0, 0, 1), 'G': None, 'T': (0, 0, 1, 0, 0)}
        exp_total_rate = 4.321762170654376e-05
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutation in internal methionine
        # ATG (TGT, GTG)
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 6.0508768963228546e-05, 'G': 0.0, 'T': 0.0}
        exp_omegas = {'A': None, 'C': (0, 0, 0, 1, 0), 'G': (0, 0, 0, 0, 1), 'T': (0, 0, 0, 0, 1)}
        exp_total_rate = 6.0508768963228546e-05
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, MU, pi5, OMEGA_VALUES_4, NT_CATEGORIES_DICT_4)

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
        self.assertEqual(self.sequence5.mu, new_sequence5.mu)
        self.assertEqual(self.sequence5.pi, new_sequence5.pi)
        self.assertEqual(self.sequence5.omega_values, new_sequence5.omega_values)
        self.assertEqual(self.sequence5.nt_categories, new_sequence5.nt_categories)
        self.assertEqual(self.sequence5.is_circular, new_sequence5.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence5.event_tree, new_sequence5.event_tree)
        self.assertCountEqual(self.sequence5.event_tree, new_sequence5.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence5.nt_sequence):
            new_nt = new_sequence5.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.38, 'C': 0.19, 'G': 0.19, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence5))
        self.assertEqual(expected, result)

    def testEventTree(self):
        expected = {'to_nt':
                        {'A': {'stationary_frequency': 0.38,
                               'from_nt':
                                   {'A': None,
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'T': {'stationary_frequency': 0.25,
                               'from_nt':
                                   {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'T': None,
                                    'C': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': {'is_trv': True,
                                          'kappa': 0.3, 'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'C': {'stationary_frequency': 0.19,
                               'from_nt':
                                   {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'T': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': None,
                                    'G': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}}}},
                         'G': {'stationary_frequency': 0.19,
                               'from_nt':
                                   {'A': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                                 'cat2': {'value': 0.0010780892949207104},
                                                                 'cat3': {'value': 0.09375337615809742},
                                                                 'cat4': {'value': 0.0}},
                                                        'False': {'cat1': {'value': 5.265401029875868e-07},
                                                                  'cat2': {'value': 0.0010780892949207104},
                                                                  'cat3': {'value': 0.09375337615809742},
                                                                  'cat4': {'value': 0.0}}}},
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn':
                                              {'True': {'cat1': {'value': 5.265401029875868e-07},
                                                        'cat2': {'value': 0.0010780892949207104},
                                                        'cat3': {'value': 0.09375337615809742},
                                                        'cat4': {'value': 0.0}},
                                               'False': {'cat1': {'value': 5.265401029875868e-07},
                                                         'cat2': {'value': 0.0010780892949207104},
                                                         'cat3': {'value': 0.09375337615809742},
                                                         'cat4': {'value': 0.0}}}},
                                    'G': None}}}}
        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, MU, pi6, OMEGA_VALUES_5, NT_CATEGORIES_DICT_3)

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
        self.assertEqual(self.sequence6.mu, new_sequence6.mu)
        self.assertEqual(self.sequence6.pi, new_sequence6.pi)
        self.assertEqual(self.sequence6.omega_values, new_sequence6.omega_values)
        self.assertEqual(self.sequence6.nt_categories, new_sequence6.nt_categories)
        self.assertEqual(self.sequence6.is_circular, new_sequence6.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence6.event_tree, new_sequence6.event_tree)
        self.assertCountEqual(self.sequence6.event_tree, new_sequence6.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence6.nt_sequence):
            new_nt = new_sequence6.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
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
        expected = \
            {'to_nt':
                 {'A': {'stationary_frequency': 0.31,
                        'from_nt':
                            {'A': None,
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'T': {'stationary_frequency': 0.23,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': None,
                             'C': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'C': {'stationary_frequency': 0.23,
                        'from_nt':
                            {'A': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': None,
                             'G': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}}}},
                  'G': {'stationary_frequency': 0.23,
                        'from_nt':
                            {'A': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn':
                                       {'True': {'cat1': {'value': 0.0876290510827283},
                                                 'cat2': {'value': 0.556154048516039},
                                                 'cat3': {'value': 0.0}},
                                        'False': {'cat1': {'value': 0.0876290510827283},
                                                  'cat2': {'value': 0.556154048516039},
                                                  'cat3': {'value': 0.0}}}},
                             'G': None}}}}

        result = self.sequence6.create_event_tree()
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.000115, 'C': 3.45e-05, 'G': None, 'T': 3.45e-05}
        # No omega keys because nt treated as synonymous mutation
        exp_omegas = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0.000184
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests internal methionine
        random.seed(4000)
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 5.6841570844245e-05, 'G': 0.00018947190281415002, 'T': 5.6841570844245e-05}
        exp_omegas = {'A': None, 'C': (0, 0, 0, 1, 0), 'G': (0, 0, 0, 1, 0), 'T': (0, 0, 0, 1, 0)}
        exp_total_rate = 0.00030315504450264
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
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

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


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
