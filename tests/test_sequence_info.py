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
        self.assertEqual(str(self.sequence1.event_tree), str(new_sequence1.event_tree))

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence1.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.dN_values, new_nt.dN_values)
            self.assertEqual(nt.dS_values, new_nt.dS_values)
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
        nt = self.sequence3.nt_sequence[6]
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 6.860859163125469e-05, 'C': 0.000946010880218476, 'G': 4.2000000000000004e-05, 'T': None},
                    {'A': (0, 0, 1, 0), 'C': (0, 0, 0, 1), 'G': (1, 0, 0, 0), 'T': None},   # dN values
                    {'A': (0, 1, 0, 0), 'C': (1, 0, 0, 0), 'G': (1, 0, 0, 0), 'T': None})   # dS values
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (syn and non-syn mutations)
        # Fwd: AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC
        # G->T and G->C substitutions induce nonsense mutations in the reverse strand
        random.seed(1000)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[7]
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 0.0007509006546280129, 'C': 0.0, 'G': None, 'T': 0.0},
                    {'A': (0, 0, 0, 2), 'C': (1, 0, 0, 0), 'G': None, 'T': (0, 1, 0, 0)},   # dN values
                    {'A': (1, 0, 1, 0), 'C': (0, 0, 0, 1), 'G': None, 'T': (0, 0, 0, 1)})   # dS values
        self.assertEqual(expected, result)

        # Tests a nucleotide position that would result in a synonymous mutation
        random.seed(555)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence4.nt_sequence[5]
        result = self.sequence4.get_substitution_rates(nt)
        expected = ({'A': 0.00021, 'C': 6.3e-05, 'G': None, 'T': 6.3e-05},
                    {'A': (0, 0, 0, 0), 'C': (0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0)},   # dN values
                    {'A': (0, 0, 0, 0), 'C': (0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0)})   # dS values
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
        expected = {'to_nt': {'A':
                                  {'stationary_frequency': 0.25,
                                   'from_nt':
                                       {'A': None,
                                        'T': {'is_trv': True,
                                              'kappa': 0.3,
                                              'is_nonsyn': {'dN': {(0, 0, 1, 0): [t1, t9],
                                                                  (0, 1, 0, 0): [t6]},
                                                            'dS': {(0, 0, 1, 0): [t1],
                                                                  (1, 0, 0, 0): [t6],
                                                                  (0, 0, 0, 1): [t9]}},
                                              'is_syn': [],
                                              'nts_in_subs': [t1, t9, t6],
                                              'number_of_events': 0},
                                        'C': {'is_trv': True,
                                              'kappa': 0.3,
                                              'is_nonsyn': {'dN': {(0, 1, 0, 0): [c4]},
                                                            'dS': {(0, 0, 1, 0): [c4]}},
                                              'is_syn': [],
                                              'nts_in_subs': [c4],
                                              'number_of_events': 0},
                                        'G': {'is_trv': False,
                                              'kappa': 1,
                                              'is_nonsyn': {'dN': {(0, 0, 0, 1): [g2]},
                                                            'dS': {(0, 0, 0, 1): [g2]}},
                                              'is_syn': [g5],
                                              'nts_in_subs': [g5, g2],
                                              'number_of_events': 1}},
                                   'events_for_nt': 1},
                              'T': {'stationary_frequency': 0.25,
                                    'from_nt':
                                        {'A': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [a0, a11],
                                                                   (0, 0, 0, 1): [a3]},
                                                             'dS': {(1, 0, 0, 0): [a0],
                                                                   (0, 0, 1, 0): [a3],
                                                                   (0, 0, 0, 1): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': [a0, a11, a3],
                                               'number_of_events': 0},
                                         'T': None,
                                         'C': {'is_trv': False,
                                               'kappa': 1,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [c4]},
                                                             'dS': {(1, 0, 0, 0): [c4]}},
                                               'is_syn': [],
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [g2, g10],
                                                                   (0, 1, 0, 0): [g7],
                                                                   (0, 0, 1, 0): [g8]},
                                                             'dS': {(0, 1, 0, 0): [g2],
                                                                   (0, 0, 1, 0): [g7],
                                                                   (1, 0, 0, 0): [g8, g10]}},
                                               'is_syn': [g5],
                                               'nts_in_subs': [g5, g2, g10, g7, g8],
                                               'number_of_events': 1}},
                                    'events_for_nt': 1},
                              'C': {'stationary_frequency': 0.08,
                                    'from_nt':
                                        {'A': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(1, 0, 0, 0): [a0, a11],
                                                                   (0, 0, 1, 0): [a3]},
                                                             'dS': {(0, 0, 1, 0): [a0],
                                                                   (0, 1, 0, 0): [a3],
                                                                   (1, 0, 0, 0): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': [a0, a11, a3],
                                               'number_of_events': 0},
                                         'T': {'is_trv': False,
                                               'kappa': 1,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [t1],
                                                                    (0, 0, 0, 1): [t6],
                                                                    (1, 0, 0, 0): [t9]},
                                                             'dS': {(1, 0, 0, 0): [t1],
                                                                    (0, 0, 1, 0): [t6, t9]}},
                                               'is_syn': [],
                                               'nts_in_subs': [t1, t6, t9],
                                               'number_of_events': 0},
                                         'C': None,
                                         'G': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 1, 0, 0): [g2],
                                                                    (0, 0, 1, 0): [g7, g8],
                                                                    (1, 0, 0, 0): [g10]},
                                                             'dS': {(0, 1, 0, 0): [g2, g10],
                                                                    (1, 0, 0, 0): [g7],
                                                                    (0, 0, 1, 0): [g8]}},
                                               'is_syn': [g5],
                                               'nts_in_subs': [g5, g2, g7, g8, g10],
                                               'number_of_events': 1}},
                                    'events_for_nt': 1},
                              'G': {'stationary_frequency': 0.42,
                                    'from_nt':
                                        {'A': {'is_trv': False,
                                               'kappa': 1,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [a0],
                                                                    (0, 1, 0, 0): [a3, a11]},
                                                             'dS': {(1, 0, 0, 0): [a0],
                                                                    (0, 0, 0, 1): [a3],
                                                                    (0, 0, 1, 0): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': [a0, a3, a11],
                                               'number_of_events': 0},
                                         'T': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [t1, t9],
                                                                    (1, 0, 0, 0): [t6]},
                                                             'dS': {(0, 1, 0, 0): [t1],
                                                                    (0, 0, 0, 1): [t6, t9]}},
                                               'is_syn': [],
                                               'nts_in_subs': [t1, t9, t6],
                                               'number_of_events': 0},
                                         'C': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [c4]},
                                                             'dS': {(1, 0, 0, 0): [c4]}},
                                               'is_syn': [],
                                               'nts_in_subs': [c4],
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
                               0.00032332685084578215,              # a5
                               0.00037456990815018656,              # t6
                               3.2395755567441353e-05,              # g7
                               0.00046415078872206957,              # a8
                               0.00011459615422288484,              # a9
                               0.000789911823084979,                # c10
                               9.439553648884362e-05,               # g11
                               0.00039078461919204766,              # a12
                               0.00027328516704295104,              # a13
                               0.00036107074192617697,              # a14
                               0.00017655918550525822,              # a15
                               0.0005720601016146754,               # t16
                               8.831265417917062e-05,               # c17
                               0.00016889260223040578,              # t18
                               0.00047789163168352714,              # g19
                               0.0002231353803076036,               # t20
                               0.00017165793053820124,              # t21
                               0.001246531333312454,                # c22
                               5.1337700605613604e-05,              # g23
                               0.00017294221882591291,              # c24
                               0.0004717808856990442,               # t25
                               6.507886436308092e-05,               # t26
                               8.908832781046372e-05,               # c27
                               0.00019631214362373376,              # a28
                               6.084458179381724e-05,               # t29
                               0.00036619447647979927,              # t30
                               0.000264,                            # c31
                               0.00039776489914726987,              # a32
                               0.00018920258566998787,              # t33
                               0.000335233686198489,                # t34
                               0.00025490526599784136,              # g35
                               0.00013433937112896293,              # c36
                               0.000264,                            # c37
                               0.001559979463247919,                # c38
                               0.00015446075439169997,              # c39
                               0.000264,                            # c40
                               0.0001534146175902228,               # a41
                               0.0005518761878033891,               # c42
                               0.00022400000000000002,              # a43
                               0.0001854640930243867,               # a44
                               0.0005770369502640317,               # t45
                               0.00024122649834313912,              # c46
                               0.0002809517503178035,               # t47
                               0.0003734876564341409,               # a48
                               5.327634285597265e-06,               # g49
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

        seq4_expected_rates = [0.0002820631384783016,               # a0
                               6.170238595738872e-05,               # t1
                               0.0008225834580011948,               # g2
                               0.00025848242360752887,              # a3
                               0.00023031258488830883,              # c4
                               0.00033600000000000004,              # g5
                               5.6934597457030823e-05,              # t6
                               0.00010156659839661558,              # g7
                               0.00038120878112926263,              # g8
                               7.400789115706304e-05,               # t9
                               0.0004887048960983141,               # g10
                               0.0003355134375527915]               # a11
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

    def testIsNonSyn(self):
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[0]    # GTA = Valine
        nt = self.nt_seq1.nt_sequence[20]  # C

        # Mutation at wobble position
        expected = False    # Synonymous
        result = codon.is_nonsyn(2, nt.state)
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True     # Non-synonymous
        result = codon.is_nonsyn(0, nt.state)
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATGACGTGGTGA
        codons = self.nt_seq4.find_codons('+0', (0, 12))
        codon = codons[2]   # TGG = Trp

        # Mutation at second position TGG --> TGA
        expected = True     # Non-synonymous
        result = codon.is_nonsyn(1, 'A')
        self.assertEqual(expected, result)

        # Mutation at last position TGG --> TGA
        expected = True    # Non-synonymous
        result = codon.is_nonsyn(2, nt.state)
        self.assertEqual(expected, result)

        # Testing mutation at position 2 in ACG ACG --> ACT
        expected = False
        codon = codons[1]
        result = codon.is_nonsyn(2, 'T')
        self.assertEqual(expected, result)

        # Testing mutation at position 10 TGA --> TAA
        expected = False
        codon = codons[3]
        result = codon.is_nonsyn(1, 'A')
        self.assertEqual(expected, result)

    def testIsStop(self):
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.is_stop(0, 'G')
        self.assertEqual(expected, result)

        # A to C mutation in middle position (TAG = STOP)
        expected = True
        result = codon.is_stop(1, 'A')
        self.assertEqual(expected, result)

        # TGA --> TAA
        codons = self.nt_seq4.find_codons('+0', (0, 12))
        codon = codons[3]
        expected = True
        result = codon.is_stop(1, 'A')
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
