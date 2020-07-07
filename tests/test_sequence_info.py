import random
import unittest

from src.sequence_info import Sequence


# ==========================================
# Tests for Sequence
# ==========================================
class TestSequenceInfo(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None
        random.seed(9001)     # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]

        self.sequence1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, dN_values, dS_values)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, dN_values, dS_values)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, dN_values, dS_values)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values)

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.sequence5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values)

    def testDeepcopy(self):
        random.seed(9001)

        # Check that Sequences are different objects with the same attributes
        new_sequence1 = self.sequence1.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence1, new_sequence1)
        self.assertEqual(self.sequence1.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence1.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence1.mu, new_sequence1.mu)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
        self.assertEqual(self.sequence1.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence1.dS_values, new_sequence1.dS_values)
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
        random.seed(9001)
        result = Sequence.get_frequency_rates('AAAAAAAAA')
        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        self.assertEqual(expected, result)

        result = Sequence.get_frequency_rates('GTACGATCGATCGATGCTAGC')
        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        self.assertEqual(expected, result)

        result = Sequence.get_frequency_rates('ATGACGTGGTGA')
        expected = {'A': 0.25, 'C': 0.08, 'T': 0.25, 'G': 0.42}
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        result = self.sequence1.get_substitution_rates(nt)
        expected = ({'A': 3.974321933436628e-05, 'C': 0.00015870631784843496, 'G': None, 'T': 0.00015870631784843496},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 1, 0)},   # dN values
                    {'A': (0, 0, 1, 0), 'C': (1, 0, 0, 0), 'G': None, 'T': (1, 0, 0, 0)})   # dS values
        self.assertEqual(expected, result)

        random.seed(9001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence1.nt_sequence[1]  # Second nucleotide is T
        result = self.sequence1.get_substitution_rates(nt)
        expected = ({'A': 9.867282041635764e-06, 'C': 0.0004378105319956827, 'G': 0.0001313431595987048, 'T': None},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': (0, 0, 1, 0), 'T': None},   # dN values
                    {'A': (0, 0, 1, 0), 'C': (1, 0, 0, 0), 'G': (1, 0, 0, 0), 'T': None})   # dS values
        self.assertEqual(expected, result)

        # Tests a sequence composed only of pyrimidines
        random.seed(900)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[1]
        result = self.sequence2.get_substitution_rates(nt)
        expected = ({'A': 8.539746787822021e-05, 'C': 0.000465, 'G': 8.539746787822021e-05, 'T': None},
                    {'A': (0, 1, 0, 0), 'C': (1, 0, 0, 0), 'G': (0, 1, 0, 0), 'T': None},   # dN values
                    {'A': (0, 0, 1, 0), 'C': (1, 0, 0, 0), 'G': (0, 0, 1, 0), 'T': None})   # dS values
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (non-syn mutations)
        random.seed(7)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[11]
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 0.0006622905526493184, 'C': 6.56715797993524e-05, 'G': None, 'T': 0.0},
                    {'A': (0, 0, 1, 1), 'C': (1, 0, 1, 0), 'G': None, 'T': (0, 1, 0, 0)},  # dN values
                    {'A': (1, 1, 0, 0), 'C': (2, 0, 0, 0), 'G': None, 'T': (1, 0, 0, 0)})  # dS values
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (syn and non-syn mutations)
        # Fwd: AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC
        # G->T and G->C substitutions induce nonsense mutations in the reverse strand
        random.seed(1000)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[7]
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
                    {'A': (0, 0, 0, 1), 'C': (0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0)},   # dN values
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0)})   # dS values
        # Substitution rates are 0 because a STOP codon would be introduced
        self.assertEqual(expected, result)

        # Tests a nucleotide position that would result in a synonymous mutation
        random.seed(555)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence4.nt_sequence[9]
        result = self.sequence4.get_substitution_rates(nt)
        expected = ({'A': 2.2956308569414034e-05, 'C': 0.000125,     'G': 1.679018660974705e-05, 'T': None},
                    {'A': (0, 1, 0, 0),           'C': (0, 1, 0, 0), 'G': (1, 0, 0, 0),          'T': None},  # dN values
                    {'A': (0, 0, 1, 0),           'C': (0, 1, 0, 0), 'G': (0, 1, 0, 0),          'T': None})  # dS values
        self.assertEqual(expected, result)

        # Tests a nonsense mutation
        random.seed(555)
        nt = self.sequence4.nt_sequence[8]
        result = self.sequence4.get_substitution_rates(nt)
        expected = ({'A': 0.0, 'C': 3.856659839661558e-05, 'G': None, 'T': 6.3e-05},
                    {'A': (0, 0, 0, 0), 'C': (0, 1, 0, 0), 'G': None, 'T': (0, 1, 0, 0)},   # dN values
                    {'A': (0, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 1, 0, 0)})   # dS values
        self.assertEqual(expected, result)

    def testIsTransv(self):
        result = self.sequence1.is_transv('A', 'A')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence3.is_transv('A', 'G')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence2.is_transv('C', 'A')
        expected = True
        self.assertEqual(expected, result)

    def testCodonIterator(self):

        # Tests iterating over the forward strand
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        results = []
        for codon in self.sequence1.codon_iterator(orf, 0, 21):
            results.append(codon)
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                    ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]
        self.assertEqual(expected, results)

        # Tests iterating over the reverse strand
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        results = []
        for codon in self.sequence1.codon_iterator(orf, 21, 0):
            results.append(codon)
        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'],
                    ['C', 'T', 'A'], ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        self.assertEqual(expected, results)

    def testFindCodons(self):
        """
        Note: For reverse strand orfs, find_codons reverses the sequence, without complementing the nucleotides.
        is_nonsyn() handles finding the complement of the nucleotides.
        """

        # Finds codons in the forward strand
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'], ['A', 'T', 'C'],
                    ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]

        result = self.sequence1.find_codons('+0', (0, 21))
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.state)

        # Find codons in the reverse strand
        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'], ['C', 'T', 'A'],
                    ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        result = self.sequence1.find_codons('-0', (21, 0))
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.state)

    def testCreateEventTree(self):
        result = self.sequence4.create_event_tree()
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.25},
                              'C': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': None,
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1}},
                                    'stationary_frequency': 0.08},
                              'G': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': None,
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.42},
                              'T': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': None},
                                    'stationary_frequency': 0.25}}}
        self.assertEqual(expected, result)

    def testNtsOnTips(self):
        # Get the Nucleotide objects
        a0 = self.sequence4.nt_sequence[0]
        t1 = self.sequence4.nt_sequence[1]
        g2 = self.sequence4.nt_sequence[2]
        a3 = self.sequence4.nt_sequence[3]
        c4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        t6 = self.sequence4.nt_sequence[6]
        g7 = self.sequence4.nt_sequence[7]
        g8 = self.sequence4.nt_sequence[8]
        t9 = self.sequence4.nt_sequence[9]
        g10 = self.sequence4.nt_sequence[10]
        a11 = self.sequence4.nt_sequence[11]

        result = self.sequence4.get_nts_on_tips()

        expected = {'to_nt':
                    {'A':
                        {'stationary_frequency': 0.25,
                         'from_nt':
                            {'A': None,
                             'T': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6, t9]},
                                                 'dS': {(0, 0, 0, 1): [t6],
                                                        (1, 0, 0, 0): [t9]}},
                                   'is_syn': [],
                                   'nts_in_subs': {t9, t6},
                                   'number_of_events': 0},
                             'C': {'is_trv': True,
                                   'kappa': 0.3,
                                   'is_nonsyn': {'dN': {(0, 0, 1, 0): [c4]},
                                                 'dS': {(0, 0, 1, 0): [c4]}},
                                   'is_syn': [],
                                   'nts_in_subs': {c4},
                                   'number_of_events': 0},
                             'G': {'is_trv': False,
                                   'kappa': 1,
                                   'is_nonsyn': {'dN': {},
                                                 'dS': {}},
                                   'is_syn': [g5],
                                   'nts_in_subs': {g5},
                                   'number_of_events': 1}},
                         'events_for_nt': 1},
                     'T': {'stationary_frequency': 0.25,
                           'from_nt':
                               {'A': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3, a11]},
                                                    'dS': {(1, 0, 0, 0): [a3],
                                                           (0, 0, 1, 0): [a11]}},
                                      'is_syn': [],
                                      'nts_in_subs': {a3, a11},
                                      'number_of_events': 0},
                                'T': None,
                                'C': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn': {'dN': {(0, 0, 0, 1): [c4]},
                                                    'dS': {(0, 1, 0, 0): [c4]}},
                                      'is_syn': [],
                                      'nts_in_subs': {c4},
                                      'number_of_events': 0},
                                'G': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(0, 1, 0, 0): [g7, g8],
                                                           (1, 0, 0, 0): [g10]},
                                                    'dS': {(0, 0, 0, 1): [g7, g10],
                                                           (0, 0, 1, 0): [g8]}},
                                      'is_syn': [g5],
                                      'nts_in_subs': {g7, g5, g10, g8},
                                      'number_of_events': 1}},
                           'events_for_nt': 1},
                     'C': {'stationary_frequency': 0.08,
                           'from_nt':
                               {'A': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(1, 0, 0, 0): [a3],
                                                           (0, 0, 1, 0): [a11]},
                                                    'dS': {(0, 0, 1, 0): [a3],
                                                           (1, 0, 0, 0): [a11]}},
                                      'is_syn': [],
                                      'nts_in_subs': {a3, a11},
                                      'number_of_events': 0},
                                'T': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn': {'dN': {(0, 1, 0, 0): [t6],
                                                           (0, 0, 0, 1): [t9]},
                                                    'dS': {(0, 1, 0, 0): [t6],
                                                           (1, 0, 0, 0): [t9]}},
                                      'is_syn': [],
                                      'nts_in_subs': {t9, t6},
                                      'number_of_events': 0},
                                'C': None,
                                'G': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(0, 0, 1, 0): [g7],
                                                           (0, 0, 0, 1): [g8, g10]},
                                                    'dS': {(0, 1, 0, 0): [g7],
                                                           (0, 0, 1, 0): [g8, g10]}},
                                      'is_syn': [g5],
                                      'nts_in_subs': {g7, g5, g10, g8},
                                      'number_of_events': 1}},
                           'events_for_nt': 1},
                     'G': {'stationary_frequency': 0.42,
                           'from_nt':
                               {'A': {'is_trv': False,
                                      'kappa': 1,
                                      'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3],
                                                           (0, 1, 0, 0): [a11]},
                                                    'dS': {(1, 0, 0, 0): [a3],
                                                           (0, 0, 1, 0): [a11]}},
                                      'is_syn': [],
                                      'nts_in_subs': {a3, a11},
                                      'number_of_events': 0},
                                'T': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6],
                                                           (0, 1, 0, 0): [t9]},
                                                    'dS': {(0, 1, 0, 0): [t6],
                                                           (1, 0, 0, 0): [t9]}},
                                      'is_syn': [],
                                      'nts_in_subs': {t9, t6},
                                      'number_of_events': 0},
                                'C': {'is_trv': True,
                                      'kappa': 0.3,
                                      'is_nonsyn': {'dN': {(0, 0, 1, 0): [c4]},
                                                    'dS': {(1, 0, 0, 0): [c4]}},
                                      'is_syn': [],
                                      'nts_in_subs': {c4},
                                      'number_of_events': 0},
                                'G': None},
                           'events_for_nt': 0}},
                    'total_events': 3}
        self.assertEqual(expected, result)


# ==========================================
# Tests for Nucleotide
# ==========================================
class TestNucleotide(unittest.TestCase):

    def setUp(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, dN_values, dS_values).nt_sequence

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.nt_seq2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, dN_values, dS_values).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, dN_values, dS_values).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values).nt_sequence

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.sequence5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values)

    def testGetMutationRate(self):

        seq1_expected_rates = [0.0003571558550312362,               # g0
                               0.0005827273266040434,               # t1
                               0.000192,                            # a2
                               7.2e-05,                             # c3
                               0.00034766669196033067,              # g4
                               0.000192,                            # a5
                               0.00025618838497572276,              # t6
                               0.0010541264093863016,               # c7
                               0.000232,                            # g8
                               0.0003079830809958123,               # a9
                               0.0002408033470208297,               # t10
                               0.0002873431595987048,               # c11
                               0.0001337127084216265,               # g12
                               0.0008629850479040508,               # a13
                               0.00016147550956710228,              # t14
                               0.00012183126917068458,              # g15
                               0.0004350559820278692,               # c16
                               0.000192,                            # t17
                               0.00016586728204163575,              # a18
                               0.0006069065088774146,               # g19
                               0.00021652268848481484]              # c20
        for pos, exp_seq1_rate in enumerate(seq1_expected_rates):
            nt = self.nt_seq1[pos]
            self.assertEqual(exp_seq1_rate, nt.mutation_rate)

        seq2_expected_rates = [0.0005891339975452602,               # t0
                               0.000744,                            # t1
                               0.0005656952120995977,               # t2
                               0.0006931397654833255,               # t3
                               0.004506792415479317,                # t4
                               0.0005785558496056398,               # t5
                               0.00010811765387774733,              # t6
                               0.0016808290584385383,               # t7
                               0.000744,                            # t8
                               0.00028363540788326586,              # t9
                               0.0006008687472494076,               # t10
                               0.0008628671796935054,               # t11
                               0.000744,                            # t12
                               0.000744]                            # t13
        for pos, exp_seq2_rate in enumerate(seq2_expected_rates):
            nt = self.nt_seq2[pos]
            self.assertEqual(exp_seq2_rate, nt.mutation_rate)

        seq3_expected_rates = [0.00022400000000000002,              # a0
                               0.00022400000000000002,              # a1
                               0.00022400000000000002,              # t2
                               0.00010471857777732271,              # t3
                               0.00035664546468782717,              # c4
                               0,                                   # a5
                               0,                                   # t6
                               0,                                   # g7
                               0.00024827025218117723,              # a8
                               0.0004090398495395135,               # a9
                               0.00040643926925193344,              # c10
                               2.239496575855691e-05,               # g11
                               0.00016889260223040578,              # a12
                               0.0005341739039971021,               # a13
                               0.0003126988468504775,               # a14
                               0.0002927806593074979,               # a15
                               0.0002113948905533971,               # t16
                               0.0011580290060315436,               # c17
                               0.00019231134008645506,              # t18
                               4.039496575855691e-05,               # g19
                               5.103283697304041e-05,               # t20
                               0.0007493685764313902,               # t21
                               0.00018081008501704828,              # c22
                               1.759038285768209e-05,               # g23
                               0.001246531333312454,                # c24
                               0.00036990181856380496,              # t25
                               0.00015337806440379493,              # t26
                               0.0004992019752883201,               # c27
                               0.00021534344906991846,              # a28
                               0.00011255202127871845,              # t29
                               0.00019631214362373376,              # t30
                               0.0002702461894726866,               # c31
                               0.0001291019253726049,               # a32
                               0.00015098256672753418,              # t33
                               0.00022400000000000002,              # t34
                               9.6e-05,                             # g35
                               0.0003316284761430394,               # c36
                               0.000264,                            # c37
                               0.0007239704854711276,               # c38
                               0.0003950968444482191,               # c39
                               0.000264,                            # c40
                               0.00014158570860052833,              # a41
                               0.001559979463247919,                # c42
                               0.00022400000000000002,              # a43
                               0.0001239781908189825,               # a44
                               0.0001534146175902228,               # t45
                               0.0003061786766654374,               # c46
                               0.0005437043017094783,               # t47
                               0.00032839160089398977,              # a48
                               8.36715797993524e-05,                # g49
                               9.6e-05,                             # g50
                               0.000264,                            # c51
                               0.000264,                            # c52
                               0.00022400000000000002,              # t53
                               0.00022400000000000002,              # a54
                               0.000264,                            # c55
                               0.000264,                            # c56
                               0.000264]                            # c57
        for pos, exp_seq3_rate in enumerate(seq3_expected_rates):
            nt = self.nt_seq3[pos]
            self.assertEqual(exp_seq3_rate, nt.mutation_rate)

        seq4_expected_rates = [0,                                   # a0
                               0,                                   # t1
                               0,                                   # g2
                               0.0004359306618420145,               # a3
                               0.00010530422972550672,              # c4
                               0.00033600000000000004,              # g5
                               4.083853581556804e-05,               # t6
                               0.0002928505292977334,               # g7
                               7.258214174242901e-05,               # g8
                               0.00020985656772831745,              # t9
                               5.6415027008750096e-05,              # g10
                               0.00025848242360752887]              # a11
        for pos, exp_seq4_rate in enumerate(seq4_expected_rates):
            nt = self.nt_seq4[pos]
            self.assertEqual(exp_seq4_rate, nt.mutation_rate)


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
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, dN_values, dS_values)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values)

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.nt_seq5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values)

    def testNtInPos(self):
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[0]
        nt = self.nt_seq1.nt_sequence[0]
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codons = self.nt_seq1.find_codons('+0', (0, 12))
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
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        exp_codon = ['G', 'T', 'A']
        exp_mutated = ['G', 'T', 'C']
        res_codon, res_mutated = codons[0].mutate_codon(2, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        codons = self.nt_seq1.find_codons('-0', (12, 0))
        exp_codon = ['G', 'A', 'T']
        exp_mutated = ['T', 'A', 'T']
        res_codon, res_mutated = codons[0].mutate_codon(0, 'A')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        # Test Sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', (0, 12))
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
        codons = self.nt_seq4.find_codons('+0', (0, 12))
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
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.is_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.is_stop(1, 'A')
        self.assertEqual(expected, result)

        # G to A mutation in middle position (TAA = STOP)
        codons = self.nt_seq4.find_codons('+0', (0, 12))
        codon = codons[3]
        expected = True
        result = codon.is_stop(1, 'A')
        self.assertEqual(expected, result)

    def testIsStart(self):
        # Testing sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[0]               # GTA = Valine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATG ACG TGG TGA
        codons = self.nt_seq4.find_codons('+0', (0, 12))
        codon = codons[0]               # first Methionine
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 5
        # ATG ATG CCC TAA
        codons = self.nt_seq5.find_codons('+0', (0, 12))
        codon = codons[1]               # second Methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
