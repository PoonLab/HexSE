import random
import unittest

from refact.new_sequence_info import Sequence
from refact.new_sequence_info import DoubleLinkedList
from refact.new_sequence_info import Nucleotide


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
        self.sequence1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, {'+0': [(0, 11)]}, kappa, mu, pi2, omegas)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
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
        random.seed(1000)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence.nucleotide_at_pos(7)
        result = self.sequence3.get_substitution_rates(nt)
        expected = ({'A': 3.4871336457097086e-05,
                     'C': 3.81675959142053e-05,
                     'G': None,
                     'T': 1.0461400937129127e-05},
                    {'A': (1, 0, 0, 1), 'C': (0, 0, 1, 1), 'G': None, 'T': (1, 0, 0, 1)})
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
        for codon in self.sequence1.codon_iterator(orf, 0, 20):
            results.append(codon)
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                    ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]
        self.assertEqual(expected, results)

        # Tests iterating over the reverse strand
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        results = []
        for codon in self.sequence1.codon_iterator(orf, 20, 0):
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
        result = self.sequence1.find_codons('+0', (0, 20))
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.get_state())

        # Find codons in the reverse strand
        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'], ['C', 'T', 'A'],
                    ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        result = self.sequence1.find_codons('-0', (20, 0))
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.get_state())

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

        expected = {'to_nt': {'A': {'events_for_nt': 2,
                                    'from_nt': {'A': None,
                                                'C': {'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': {'is_nonsyn': {(0, 0, 0, 1): [g8],
                                                                    (0, 0, 1, 0): [g2, g7]},
                                                      'is_syn': [g5, g10],
                                                      'is_trv': False,
                                                      'kappa': 1,
                                                      'nts_in_subs': [g5, g10, g2, g7, g8],
                                                      'number_of_events': 2},
                                                'T': {'is_nonsyn': {(0, 1, 0, 0): [t6, t9],
                                                                    (1, 0, 0, 0): [t1]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [t1, t6, t9],
                                                      'number_of_events': 0}},
                                    'stationary_frequency': 0.25},
                              'C': {'events_for_nt': 1,
                                    'from_nt': {'A': {'is_nonsyn': {(1, 0, 0, 0): [a0, a3, a11]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0},
                                                'C': None,
                                                'G': {'is_nonsyn': {(0, 0, 0, 1): [g8],
                                                                    (0, 0, 1, 0): [g2],
                                                                    (0, 1, 0, 0): [g7],
                                                                    (1, 0, 0, 0): [g10]},
                                                      'is_syn': [g5],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [g5, g2, g7, g8, g10],
                                                      'number_of_events': 1},
                                                'T': {'is_nonsyn': {(0, 0, 0, 1): [t6],
                                                                    (0, 0, 1, 0): [t1, t9]},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1,
                                                      'nts_in_subs': [t1, t9, t6],
                                                      'number_of_events': 0}},
                                    'stationary_frequency': 0.08},
                              'G': {'events_for_nt': 0,
                                    'from_nt': {'A': {'is_nonsyn': {(0, 0, 0, 1): [a3],
                                                                    (0, 0, 1, 0): [a0],
                                                                    (0, 1, 0, 0): [a11]},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1,
                                                      'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0},
                                                'C': {'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': None,
                                                'T': {'is_nonsyn': {(0, 0, 0, 1): [t9],
                                                                    (0, 1, 0, 0): [t6],
                                                                    (1, 0, 0, 0): [t1]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [t1, t6, t9],
                                                      'number_of_events': 0}},
                                    'stationary_frequency': 0.42},
                              'T': {'events_for_nt': 1,
                                    'from_nt': {'A': {'is_nonsyn': {(0, 0, 1, 0): [a0],
                                                                    (0, 1, 0, 0): [a3],
                                                                    (1, 0, 0, 0): [a11]},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [a0, a3, a11],
                                                      'number_of_events': 0},
                                                'C': {'is_nonsyn': {(0, 1, 0, 0): [c4]},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1,
                                                      'nts_in_subs': [c4],
                                                      'number_of_events': 0},
                                                'G': {'is_nonsyn': {(0, 0, 0, 1): [g10],
                                                                    (0, 0, 1, 0): [g2, g8],
                                                                    (0, 1, 0, 0): [g7]},
                                                      'is_syn': [g5],
                                                      'is_trv': True,
                                                      'kappa': 0.3,
                                                      'nts_in_subs': [g5, g2, g8, g7, g10],
                                                      'number_of_events': 1},
                                                'T': None},
                                    'stationary_frequency': 0.25}},
                    'total_events': 4}

        result = self.sequence4.get_nts_on_tips()
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
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas).nt_sequence

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.nt_seq2 = Sequence(s2, {'+0': [(0, 11)]}, kappa, mu, pi2, omegas).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas).nt_sequence

    def testGetMutationRate(self):

        seq1_expected_rates = [0.0001356139344659256,               # g0
                               0.00014951453904781168,              # t1
                               0.000192,                            # a2
                               0.00020291838733299422,              # c3
                               0.0002541108558198061,               # g4
                               0.000192,                            # a5
                               0.00033272911581276056,              # t6
                               0.00028496761812266927,              # c7
                               0.000232,                            # g8
                               0.00017555974417030556,              # a9
                               0.00033272911581276056,              # t10
                               0.00019451962784607557,              # c11
                               0.00022772643797867983,              # g12
                               0.00025892241300017537,              # a13
                               0.00015413838212213998,              # t14
                               0.0003466516825378471,               # g15
                               0.00014505458381268202,              # c16
                               0.000192,                            # t17
                               0.00019049887967468768,              # a18
                               0.00021445050151207423,              # g19
                               0.0002298616178778603]               # c20
        for pos, exp_seq1_rate in enumerate(seq1_expected_rates):
            nt = self.nt_seq1.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq1_rate, nt.mutation_rate)

        seq2_expected_rates = [0.0008149072272842057,               # t0
                               0.0005041973621946621,               # t1
                               0.000782362034272396,                # t2
                               0.0002686589737207266,               # t3
                               0.0009232589622885183,               # t4
                               0.000705637965727605,                # t5
                               2.8377150119088977e-05,              # t6
                               0.0014744011273235508,               # t7
                               0.000744,                            # t8
                               0.0006724059786107684,               # t9
                               0.0004873301750616648,               # t10
                               0.0006477488156481243,               # t11
                               0.000744,                            # t12
                               0.000744]                            # t13
        for pos, exp_seq2_rate in enumerate(seq2_expected_rates):
            nt = self.nt_seq2.nucleotide_at_pos(pos)
            self.assertEqual(exp_seq2_rate, nt.mutation_rate)

        seq3_expected_rates = [0.00022400000000000002,              # a0
                               0.00022400000000000002,              # a1
                               0.00022400000000000002,              # t2
                               0.00024534841251567483,              # t3
                               0.00025038766525818244,              # c4
                               4.852607794352704e-05,               # a5
                               0.0003498912125648934,               # t6
                               0.00013610645494507173,              # g7
                               0.0005810105024816592,               # a8
                               0.00013153002713797121,              # a9
                               0.0002225470433518085,               # c10
                               6.869772697180496e-05,               # g11
                               0.00021917390785691676,              # a12
                               6.230625456314679e-05,               # a13
                               0.00026568203594183745,              # a14
                               0.0002950595025379696,               # a15
                               0.00026326621205451876,              # t16
                               0.0002225470433518085,               # c17
                               0.0002222486929538023,               # t18
                               5.854642287352648e-05,               # g19
                               9.096965281074431e-05,               # t20
                               0.00020705565663234752,              # t21
                               9.839680507721865e-05,               # c22
                               0.00025088133938514987,              # g23
                               0.0001473365136783115,               # c24
                               0.00024534841251567483,              # t25
                               0.0003253037295286643,               # t26
                               0.00048826322106726097,              # c27
                               0.0004439057157533272,               # a28
                               0.00024299220598339,                 # t29
                               0.00024534841251567483,              # t30
                               0.0002709289765767079,               # c31
                               0.00030207614850020467,              # a32
                               0.0001744336288891136,               # t33
                               0.00019431753807718763,              # t34
                               0.0001499533768135257,               # g35
                               0.0004123717862371957,               # c36
                               0.000264,                            # c37
                               0.00019346492182427744,              # c38
                               0.0003918304749186703,               # c39
                               0.000264,                            # c40
                               9.831556415490115e-05,               # a41
                               0.00022348856971245558,              # c42
                               0.00022400000000000002,              # a43
                               0.00013660832003822845,              # a44
                               0.00037299093212676593,              # t45
                               0.00022901709844811395,              # c46
                               0.00021987380218631159,              # t47
                               0.00018158104630298428,              # a48
                               0.00013134199003178473,              # g49
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

        seq4_expected_rates = [0.00019843633299446633,              # a0
                               0.0002190610826032811,               # t1
                               0.00024622808232180766,              # g2
                               0.0002481878930883114,               # a3
                               9.036091318349602e-05,               # c4
                               0.00033600000000000004,              # g5
                               0.00019631589480920673,              # t6
                               0.0005594863981901488,               # g7
                               0.0003680226187735122,               # g8
                               7.222015422600176e-05,               # t9
                               0.0003186752103285959,               # g10
                               0.00018075429532547538]              # a11
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
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas).nt_sequence

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        self.nt_seq2 = Sequence(s2, {'+0': [(0, 11)]}, kappa, mu, pi2, omegas).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas).nt_sequence

    def testInsertNt(self):
        seq = DoubleLinkedList()

        # Test adding first nucleotide
        seq.insert_nt('A', 0)
        exp_head = Nucleotide('A', 0)
        exp_current = exp_head

        self.assertEqual(exp_head.state, seq.head.state)
        self.assertEqual(exp_head.pos_in_seq, seq.head.pos_in_seq)
        self.assertEqual(exp_current.state, seq.current_nt.state)
        self.assertEqual(exp_current.pos_in_seq, seq.current_nt.pos_in_seq)

        # Test adding second nucleotide
        seq.insert_nt('T', 1)
        exp_current = Nucleotide('T', 1)

        self.assertEqual(exp_head.state, seq.head.state)
        self.assertEqual(exp_head.pos_in_seq, seq.head.pos_in_seq)
        self.assertEqual(exp_current.state, seq.current_nt.state)
        self.assertEqual(exp_current.pos_in_seq, seq.current_nt.pos_in_seq)

    def testNucleotideAtPos(self):
        # Tests first Nucleotide
        result_nt = self.nt_seq1.nucleotide_at_pos(0)
        expected_state = 'G'
        self.assertEqual(expected_state, result_nt.get_state())

        # Tests last Nucleotide
        result_nt = self.nt_seq2.nucleotide_at_pos(13)
        expected_state = 'T'
        self.assertEqual(expected_state, result_nt.get_state())

        # Tests middle Nucleotide
        result_nt = self.nt_seq3.nucleotide_at_pos(5)
        expected_state = 'A'
        self.assertEqual(expected_state, result_nt.get_state())

    def testSliceSequence(self):
        # Slicing in the middle
        nts = self.nt_seq1.slice_sequence(1, 4)
        expected = ['T', 'A', 'C', 'G']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Slicing from the start
        nts = self.nt_seq1.slice_sequence(0, 3)
        expected = ['G', 'T', 'A', 'C']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Slicing from the end
        nts = self.nt_seq1.slice_sequence(18, 20)
        expected = ['A', 'G', 'C']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Slicing out of bounds (right)
        nts = self.nt_seq1.slice_sequence(21, 23)
        expected = []
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Slicing out of bounds (left)
        nts = self.nt_seq1.slice_sequence(-1, -5)
        expected = []
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Slicing in reverse order
        nts = self.nt_seq1.slice_sequence(5, 1)
        expected = ['T', 'A', 'C', 'G', 'A']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

        # Accessing a single nucleotide
        nts = self.nt_seq1.slice_sequence(0, 0)
        expected = ['G']
        result = [nt.get_state() for nt in nts]
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
        self.nt_seq1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)

    def testNtInPos(self):
        codons = self.nt_seq1.find_codons('+0', (0, 11))
        codon = codons[0]
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(0)
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codons = self.nt_seq1.find_codons('+0', (0, 11))
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
        codons = self.nt_seq1.find_codons('+0', (0, 11))
        codon = codons[0]    # GTA = Valine
        nt = self.nt_seq1.nt_sequence.nucleotide_at_pos(20)  # C

        # At wobble position
        expected = False
        result = codon.is_nonsyn(2, nt.state)
        self.assertEqual(expected, result)

        # At first position
        expected = True
        result = codon.is_nonsyn(0, nt.state)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
