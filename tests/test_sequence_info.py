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
        self.sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)

    def testGetFrequencyRates(self):
        result = Sequence.get_frequency_rates('AAAAAAAAA')
        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        self.assertEqual(expected, result)

        result = Sequence.get_frequency_rates('GTACGATCGATCGATGCTAGC')
        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence1.nt_sequence.slice_sequence(1, 1)  # First nucleotide is G
        result = self.sequence1.get_substitution_rates(nt[0])
        expected = ({'A': 1.0557889780446518e-05,
                     'C': 0.00012839875948691864,
                     'G': 3.851962784607559e-05,
                     'T': None},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': (0, 0, 1, 0), 'T': None})
        self.assertEqual(expected, result)

        # Tests a sequence composed only of pyrimidines
        random.seed(900)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence.slice_sequence(1, 1)
        result = self.sequence2.get_substitution_rates(nt[0])
        expected = ({'A': 9.137440782406214e-05,
                     'C': 0.0004975451930118097,
                     'G': 4.091182289923025e-05,
                     'T': None},
                    {'A': (0, 1, 0, 0), 'C': (0, 0, 1, 0), 'G': (1, 0, 0, 0), 'T': None})
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (non-syn mutations)
        random.seed(7)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence.slice_sequence(6, 6)
        result = self.sequence3.get_substitution_rates(nt[0])
        expected = ({'A': 4.493956582042153e-05,
                     'C': 9.170191466214123e-05,
                     'G': 8.323232170374884e-05,
                     'T': None},
                    {'A': (0, 0, 1, 0), 'C': (0, 1, 0, 0), 'G': (0, 0, 0, 1), 'T': None})
        self.assertEqual(expected, result)

        # Tests a nucleotide involved in multiple orfs (syn and non-syn mutations)
        random.seed(1000)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence.slice_sequence(7, 7)
        result = self.sequence3.get_substitution_rates(nt[0])
        expected = ({'A': 3.4871336457097086e-05,
                     'C': 3.81675959142053e-05,
                     'G': None,
                     'T': 1.0461400937129127e-05},
                    {'A': (1, 0, 0, 1), 'C': (0, 0, 1, 1), 'G': None, 'T': (1, 0, 0, 1)})
        self.assertEqual(expected, result)

        # Tests a nucleotide position that would result in a synonymous mutation
        random.seed(555)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence4.nt_sequence.slice_sequence(5, 5)
        result = self.sequence4.get_substitution_rates(nt[0])
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
