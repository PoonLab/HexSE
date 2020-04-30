import random
import unittest

from src.sequence_info import DoubleLinkedList
from src.sequence_info import Nucleotide
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
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.sequence1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, omegas)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)

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
        nt = self.sequence1.nt_sequence.nucleotide_at_pos(0)  # First nucleotide is G
        result = self.sequence1.get_substitution_rates(nt)
        expected = ({'A': 4.252483383790958e-05,
                     'C': 4.654455031400801e-05,
                     'G': None,
                     'T': 4.654455031400801e-05},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 1, 0)})
        self.assertEqual(expected, result)

        random.seed(9001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence1.nt_sequence.nucleotide_at_pos(1)  # Second nucleotide is T
        result = self.sequence1.get_substitution_rates(nt)
        expected = ({'A': 1.0557889780446518e-05,
                     'C': 0.00012839875948691864,
                     'G': 3.851962784607559e-05,
                     'T': None},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': (0, 0, 1, 0), 'T': None})
        self.assertEqual(expected, result)

        # Tests a sequence composed only of pyrimidines
        random.seed(900)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence.nucleotide_at_pos(1)
        result = self.sequence2.get_substitution_rates(nt)
        expected = ({'A': 9.137440782406214e-05,
                     'C': 0.0004975451930118097,
                     'G': 4.091182289923025e-05,
                     'T': None},
                    {'A': (0, 1, 0, 0), 'C': (0, 0, 1, 0), 'G': (1, 0, 0, 0), 'T': None})
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (non-syn mutations)
        random.seed(7)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence.nucleotide_at_pos(6)
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 4.493956582042153e-05,
                     'C': 9.170191466214123e-05,
                     'G': 8.323232170374884e-05,
                     'T': None},
                    {'A': (0, 0, 1, 0), 'C': (0, 1, 0, 0), 'G': (0, 0, 0, 1), 'T': None})
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (syn and non-syn mutations)
        # Fwd: AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC
        # G->T and G->C substitutions induce nonsense mutations in the reverse strand
        random.seed(1000)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence.nucleotide_at_pos(7)
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 3.4871336457097086e-05,
                     'C': 0.0,
                     'G': None,
                     'T': 0.0},
                    {'A': (1, 0, 0, 1), 'C': (0, 0, 0, 1), 'G': None, 'T': (0, 0, 1, 0)})
        self.assertEqual(expected, result)

        # Tests a nucleotide position that would result in a synonymous mutation
        random.seed(555)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence4.nt_sequence.nucleotide_at_pos(5)
        result = self.sequence4.get_substitution_rates(nt)
        expected = ({'A': 0.00021,
                     'C': 6.3e-05,
                     'G': None,
                     'T': 6.3e-05},
                    {'A': (0, 0, 0, 0), 'C': (0, 0, 0, 0), 'G': None, 'T': (0, 0, 0, 0)})
        self.assertEqual(expected, result)

        # Tests a nonsense mutation
        random.seed(555)
        nt = self.sequence4.nt_sequence.nucleotide_at_pos(8)
        result = self.sequence4.get_substitution_rates(nt)
        expected = ({'A': 0.0,
                     'C': 4.126586159796355e-05,
                     'G': None,
                     'T': 6.740934873063228e-05},
                    {'A': (0, 0, 0, 0), 'C': (0, 1, 0, 0), 'G': None, 'T': (0, 0, 1, 0)})
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

    def testGetStringSequence(self):
        expected = 'GTACGATCGATCGATGCTAGC'
        result = self.sequence1.get_string_sequence()
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        result = self.sequence4.create_event_tree()
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'T': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.25},
                              'C': {'from_nt': {'A': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': None,
                                                'G': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1}},
                                    'stationary_frequency': 0.08},
                              'G': {'from_nt': {'A': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'C': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': None,
                                                'T': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.42},
                              'T': {'from_nt': {'A': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'G': {'is_nonsyn': {},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': None},
                                    'stationary_frequency': 0.25}}}
        self.assertEqual(expected, result)

    def testNtsOnTips(self):
        # Get the Nucleotide objects
        a0 = self.sequence4.nt_sequence.nucleotide_at_pos(0)
        t1 = self.sequence4.nt_sequence.nucleotide_at_pos(1)
        g2 = self.sequence4.nt_sequence.nucleotide_at_pos(2)
        a3 = self.sequence4.nt_sequence.nucleotide_at_pos(3)
        c4 = self.sequence4.nt_sequence.nucleotide_at_pos(4)
        g5 = self.sequence4.nt_sequence.nucleotide_at_pos(5)
        t6 = self.sequence4.nt_sequence.nucleotide_at_pos(6)
        g7 = self.sequence4.nt_sequence.nucleotide_at_pos(7)
        g8 = self.sequence4.nt_sequence.nucleotide_at_pos(8)
        t9 = self.sequence4.nt_sequence.nucleotide_at_pos(9)
        g10 = self.sequence4.nt_sequence.nucleotide_at_pos(10)
        a11 = self.sequence4.nt_sequence.nucleotide_at_pos(11)

        result = self.sequence4.get_nts_on_tips()
        expected = {'to_nt': {'A': {'stationary_frequency': 0.25,
                                    'from_nt': {'A': None,
                                                'T': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(1, 0, 0, 0): [t1],
                                                                    (0, 1, 0, 0): [t6],
                                                                    (0, 0, 0, 1): [t9]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [t1, t6, t9],
                                                      'number_of_events': 0},
                                                'C': {'is_trv': True, 'kappa': 0.3,
                                                      'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                                      'is_syn': [], 'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': {'is_trv': False,
                                                      'kappa': 1,
                                                      'is_nonsyn': {(0, 0, 1, 0): [g2]},
                                                      'is_syn': [g5, g7, g8, g10],
                                                      'nts_in_subs': [g5, g7, g8, g10, g2],
                                                      'number_of_events': 4}},
                                    'events_for_nt': 4},
                              'T': {'stationary_frequency': 0.25,
                                    'from_nt': {'A': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(0, 0, 1, 0): [a0],
                                                                    (0, 1, 0, 0): [a3],
                                                                    (1, 0, 0, 0): [a11]},
                                                      'is_syn': [], 'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0}, 'T': None,
                                                'C': {'is_trv': False,
                                                      'kappa': 1,
                                                      'is_nonsyn': {(0, 1, 0, 0): [c4]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(0, 0, 1, 0): [g2],
                                                                    (0, 1, 0, 0): [g7],
                                                                    (0, 0, 0, 1): [g8, g10]},
                                                      'is_syn': [g5],
                                                      'nts_in_subs': [g5, g2, g7, g8, g10],
                                                      'number_of_events': 1}}, 'events_for_nt': 1},
                              'C': {'stationary_frequency': 0.08,
                                    'from_nt': {'A': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(1, 0, 0, 0): [a0, a3, a11]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0},
                                                'T': {'is_trv': False,
                                                      'kappa': 1,
                                                      'is_nonsyn': {(0, 0, 1, 0): [t1, t9],
                                                                    (0, 0, 0, 1): [t6]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [t1, t9, t6],
                                                      'number_of_events': 0},
                                                'C': None,
                                                'G': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(0, 0, 1, 0): [g2, g7, g10],
                                                                    (0, 1, 0, 0): [g8]},
                                                      'is_syn': [g5],
                                                      'nts_in_subs': [g5, g2, g7, g10, g8],
                                                      'number_of_events': 1}},
                                    'events_for_nt': 1},
                              'G': {'stationary_frequency': 0.42,
                                    'from_nt': {'A': {'is_trv': False,
                                                      'kappa': 1,
                                                      'is_nonsyn': {(0, 0, 1, 0): [a0],
                                                                    (0, 0, 0, 1): [a3, a11]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0},
                                                'T': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(1, 0, 0, 0): [t1],
                                                                    (0, 1, 0, 0): [t6, t9]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [t1, t6, t9],
                                                      'number_of_events': 0},
                                                'C': {'is_trv': True,
                                                      'kappa': 0.3,
                                                      'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                                      'is_syn': [],
                                                      'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': None},
                                    'events_for_nt': 0}},
                              'total_events': 6}
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
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas).nt_sequence

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.nt_seq2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, omegas).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas).nt_sequence

    def testGetMutationRate(self):

        seq1_expected_rates = [0.0001356139344659256,               # g0
                               0.00014951453904781168,              # t1
                               0.000192,                            # a2
                               7.451962784607558e-05,               # c3
                               0.00021445050151207423,              # g4
                               0.000192,                            # a5
                               0.0002212856212025476,               # t6
                               0.00014994363117076288,              # c7
                               0.000232,                            # g8
                               0.00028496761812266927,              # a9
                               0.00014070176132674722,              # t10
                               0.00022734199003178473,              # c11
                               0.0003623873272993538,               # g12
                               0.00028688415106580445,              # a13
                               0.00020189987981223125,              # t14
                               0.0001939393378160447,               # g15
                               0.0002868841510658045,               # c16
                               0.000192,                            # t17
                               9.729308612259077e-05,               # a18
                               0.00024823760167470935,              # g19
                               0.0001690775176265221]               # c20
        for pos, exp_seq1_rate in enumerate(seq1_expected_rates):
            nt = self.nt_seq1.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq1_rate, nt.mutation_rate)

        seq2_expected_rates = [0.0011116760853799926,               # t0
                               0.001347214473853928,                # t1
                               0.0005972862307232924,               # t2
                               0.00045373477726983017,              # t3
                               0.0004947567402163136,               # t4
                               0.0008907137692767086,               # t5
                               4.1038013720301285e-05,              # t6
                               0.0005620865122741429,               # t7
                               0.000744,                            # t8
                               0.0012893253237744474,               # t9
                               0.0011042495202253436,               # t10
                               0.0006477488156481243,               # t11
                               0.000744,                            # t12
                               0.000744]                            # t13
        for pos, exp_seq2_rate in enumerate(seq2_expected_rates):
            nt = self.nt_seq2.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq2_rate, nt.mutation_rate)

        seq3_expected_rates = [0.00022400000000000002,              # a0
                               0.00022400000000000002,              # a1
                               0.00022400000000000002,              # t2
                               0.0001315300271379712,               # t3
                               0.00015105972458205792,              # c4
                               9.757019434373241e-05,               # a5
                               8.088657273312198e-05,               # t6
                               0.00012722531971401767,              # g7
                               0.0006223491730873115,               # a8
                               0.00037299093212676593,              # a9
                               8.516109776750775e-05,               # c10
                               9.049874011178012e-05,               # g11
                               0.00011998670479161041,              # a12
                               0.00021917390785691676,              # a13
                               4.814164729637698e-05,               # a14
                               0.00028608126335599195,              # a15
                               0.0002950595025379696,               # t16
                               0.00041242442849848566,              # c17
                               0.0001518013563596832,               # t18
                               8.777987208515278e-05,               # g19
                               0.00013306374507634475,              # t20
                               0.00013660832003822845,              # t21
                               0.00025633516090024886,              # c22
                               3.5155235748425486e-05,              # g23
                               0.00014141290207582652,              # c24
                               0.00020481970153202317,              # t25
                               4.007349139922032e-05,               # t26
                               0.00034986962921354855,              # c27
                               0.00033246222114311424,              # a28
                               0.0005098872686737637,               # t29
                               0.0004439057157533272,               # t30
                               0.00023248158673646788,              # c31
                               0.00021987380218631159,              # a32
                               0.0003673202039866725,               # t33
                               0.00019431753807718763,              # t34
                               8.873813855672038e-05,               # g35
                               0.0001610026629021978,               # c36
                               0.000264,                            # c37
                               0.00011587191489684776,              # c38
                               0.0003070667075580456,               # c39
                               0.000264,                            # c40
                               0.0002198738021863116,               # a41
                               0.00019346492182427744,              # c42
                               0.00022400000000000002,              # a43
                               0.00020244481076453242,              # a44
                               9.831556415490115e-05,               # t45
                               0.00022901709844811395,              # c46
                               0.00018962666521056838,              # t47
                               0.00037299093212676593,              # a48
                               4.0949939906115615e-05,              # g49
                               9.6e-05,                             # g50
                               0.000264,                            # c51
                               0.000264,                            # c52
                               0.00022400000000000002,              # t53
                               0.00022400000000000002,              # a54
                               0.000264,                            # c55
                               0.000264,                            # c56
                               0.000264]                            # c57
        for pos, exp_seq3_rate in enumerate(seq3_expected_rates):
            nt = self.nt_seq3.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq3_rate, nt.mutation_rate)

        seq4_expected_rates = [0.00028327605787648166,              # a0
                               0.00019843633299446633,              # t1
                               0.0004743947942133541,               # g2
                               0.00019843633299446633,              # a3
                               9.562805035526817e-05,               # c4
                               0.00033600000000000004,              # g5
                               0.0001465643347153617,               # t6
                               0.00013481869746126456,              # g7
                               0.0002496969651112465,               # g8
                               0.00024818789308831134,              # t9
                               0.0001661143441535868,               # g10
                               0.0001852884854149908]              # a11
        for pos, exp_seq4_rate in enumerate(seq4_expected_rates):
            nt = self.nt_seq4.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq4_rate, nt.mutation_rate)


# ==========================================
# Tests for DoubleLinkedList
# ==========================================
class TestDoubleLinkedList(unittest.TestCase):

    def setUp(self):
        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas).nt_sequence

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.nt_seq2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, omegas).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas).nt_sequence

    def testInsertNt(self):
        seq = DoubleLinkedList()

        # Test adding first nucleotide
        seq.insert_nt('A', 0)
        exp_head = Nucleotide('A', 0)
        exp_tail = exp_head

        self.assertEqual(exp_head.state, seq.head.state)
        self.assertEqual(exp_head.pos_in_seq, seq.head.pos_in_seq)
        self.assertEqual(exp_tail.state, seq.tail.state)
        self.assertEqual(exp_tail.pos_in_seq, seq.tail.pos_in_seq)

        # Test adding second nucleotide
        seq.insert_nt('T', 1)
        exp_tail = Nucleotide('T', 1)

        self.assertEqual(exp_head.state, seq.head.state)
        self.assertEqual(exp_head.pos_in_seq, seq.head.pos_in_seq)
        self.assertEqual(exp_tail.state, seq.tail.state)
        self.assertEqual(exp_tail.pos_in_seq, seq.tail.pos_in_seq)

    def testNucleotideAtPos(self):
        # Tests first Nucleotide
        result_nt = self.nt_seq1.nucleotide_at_pos(0)
        expected_state = 'G'
        self.assertEqual(expected_state, result_nt.state)

        # Tests last Nucleotide
        result_nt = self.nt_seq2.nucleotide_at_pos(13)
        expected_state = 'T'
        self.assertEqual(expected_state, result_nt.state)

        # Tests middle Nucleotide
        result_nt = self.nt_seq3.nucleotide_at_pos(5)
        expected_state = 'A'
        self.assertEqual(expected_state, result_nt.state)

    def testSliceSequence(self):
        # Slicing in the middle
        nts = self.nt_seq1.slice_sequence(1, 4)
        expected = ['T', 'A', 'C', 'G']
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Slicing from the start
        nts = self.nt_seq1.slice_sequence(0, 3)
        expected = ['G', 'T', 'A', 'C']
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Slicing from the end
        nts = self.nt_seq1.slice_sequence(18, 20)
        expected = ['A', 'G', 'C']
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Slicing out of bounds (right)
        nts = self.nt_seq1.slice_sequence(21, 23)
        expected = []
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Slicing out of bounds (left)
        nts = self.nt_seq1.slice_sequence(-1, -5)
        expected = []
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Slicing in reverse order
        nts = self.nt_seq1.slice_sequence(5, 1)
        expected = ['T', 'A', 'C', 'G', 'A']
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)

        # Accessing a single nucleotide
        nts = self.nt_seq1.slice_sequence(0, 0)
        expected = ['G']
        result = [nt.state for nt in nts]
        self.assertEqual(expected, result)


# ==========================================
# Tests for Codon
# ==========================================
class TestCodon(unittest.TestCase):

    def setUp(self):
        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)

    def testNtInPos(self):
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[0]
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(0)
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[1]
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(3)
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(4)
        expected = 1
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(5)
        expected = 2
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testIsNonSyn(self):
        codons = self.nt_seq1.find_codons('+0', (0, 12))
        codon = codons[0]    # GTA = Valine
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(20)  # C

        # Mutation at wobble position
        expected = False    # Synonymous
        result = codon.is_nonsyn(2, nt.state)
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True     # Non-synonymous
        result = codon.is_nonsyn(0, nt.state)
        self.assertEqual(expected, result)

        # Testing sequence ATGACGTGGTGA
        codons = self.nt_seq4.find_codons('+0', (0, 12))
        codon = codons[2]   # TGG = Trp
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(0)  # A

        # Mutation at second position TGA --> TAG
        expected = True     # Non-synonymous
        result = codon.is_nonsyn(1, 'A')
        self.assertEqual(expected, result)

        # Mutation at last position
        expected = True     # Non-synonymous
        result = codon.is_nonsyn(2, 'A')
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


if __name__ == '__main__':
    unittest.main()
