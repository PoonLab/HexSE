import unittest
import random

from refact.new_sequence_info import Sequence
from refact.new_sequence_info import DoubleLinkedList
from refact.new_sequence_info import Nucleotide
from refact.new_sequence_info import Codon


def make_dll(original_sequence):
    s = Sequence(original_sequence, rcseq=None, sorted_orfs=None, mu=0.5, pi=None, kappa=0.3)
    s.nt_sequence = DoubleLinkedList()
    for pos, nt in enumerate(s.original_seq):
        s.nt_sequence.insert_nt(nt, pos)
    return s


class TestSequences(unittest.TestCase):
    def setUp(self):
        s1 = make_dll('ATTTTTTTTTTTTT')
        self.seq1 = s1.nt_sequence
        self.sequence1 = s1

        s2 = make_dll('TTTTTTTTTTTTTG')
        self.seq2 = s2.nt_sequence
        self.sequence2 = s2

        s3 = make_dll('TTTTTTCTTTTTTT')
        self.seq3 = s3.nt_sequence
        self.sequence3 = s3

        s4 = make_dll('ATGCCGTATGC')
        self.seq4 = s4.nt_sequence
        self.sequence4 = s4

        s5 = make_dll('ATGCCCTGA')
        self.seq5 = s5.nt_sequence
        self.sequence5 = s5

        s6 = make_dll('AAAAAAAAA')
        self.seq6 = s6.nt_sequence
        self.sequence6 = s6

        s7 = make_dll('GTACGATCGATCGATGCTAGC')
        self.seq7 = s7.nt_sequence
        self.sequence7 = s7


# ==========================================
# Tests for Sequence
# ==========================================
class TestCreateKeys(TestSequences):

    def testOmegas(self):
        my_omegas = [0.29327471612351436,
                     0.6550136761581515,
                     1.0699896623909886,
                     1.9817219453273531]
        result = self.sequence7.create_keys(my_omegas)
        expected = (1, 1, 1, 1)
        self.assertEqual(expected, result)

    def testSixOmegaValues(self):
        my_omegas = [0.29327471612351436,
                     0.6550136761581515,
                     1.0699896623909886,
                     1.9817219453273531,
                     1.0699896623909886,
                     1.9817219453273531]
        result = self.sequence7.create_keys(my_omegas)
        expected = (1, 1, 2, 2)
        self.assertEqual(expected, result)

    def testNoOmegas(self):
        my_omegas = []
        result = self.sequence7.create_keys(my_omegas)
        expected = (0, 0, 0, 0)
        self.assertEqual(expected, result)

    def testThreeValues(self):
        my_omegas = [1.9817219453273531, 0.29327471612351436, 0.29327471612351436]
        result = self.sequence7.create_keys(my_omegas)
        expected = (2, 0, 0, 1)
        self.assertEqual(expected, result)


class TestGetSubstitutionRates(unittest.TestCase):

    def testDefaultPi(self):
        random.seed(9001)   # Set seed value to initialize pseudo-random number generator

        seq = Sequence('GTACGATCGATCGATGCTAGC', None, {'+0': [(0, 20)]}, 0.5, None, 0.3)
        nt = seq.nt_sequence.slice_sequence(0, 0)       # First nucleotide is G

        result = seq.get_substitution_rates(nt[0])
        expected = ({'A': 0.04252483383790958,
                     'C': 0.046544550314008,
                     'G': None,
                     'T': 0.08620490462173985},
                    {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 0, 1)})
        self.assertEqual(expected, result)

    def testNonDefaultPi(self):

        random.seed(900)    # Set seed value to initialize pseudo-random number generator

        pi = {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}
        seq = Sequence('TTTTTTCTTTTTTT', None, {'+0': [(0, 11)]}, 0.5, pi, 0.3)
        nt = seq.nt_sequence.slice_sequence(1, 1)

        result = seq.get_substitution_rates(nt[0])
        expected = ({'A': 0.07431457294977574,
                     'C': 0.036659339515439295,
                     'G': 0.07431457294977574,
                     'T': None},
                    {'A': (0, 0, 0, 1), 'C': (1, 0, 0, 0), 'G': (0, 0, 0, 1), 'T': None})
        self.assertEqual(expected, result)


class TestIsTransversion(TestSequences):
    """Since is_transv is a static method, it can be called with or without an instance of a Sequence object"""

    def testSameNucleotide(self):
        result = self.sequence7.is_transv('A', 'A')
        expected = None
        self.assertEqual(expected, result)

    def testTransition(self):
        result = self.sequence7.is_transv('A', 'G')
        expected = False
        self.assertEqual(expected, result)

    def testTransversion(self):
        result = self.sequence7.is_transv('C', 'A')
        expected = True
        self.assertEqual(expected, result)


class TestCodonIterator(TestSequences):
    """Since codon_iterator is a static method, it can be called with or without an instance of a Sequence object"""

    def testForwardStrand(self):
        results = []
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        for codon in self.sequence7.codon_iterator(orf, 0, 20):
            results.append(codon)

        expected =[['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                   ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]
        self.assertEqual(expected, results)

    def testReverseStrand(self):
        results = []
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        for codon in self.sequence7.codon_iterator(orf, 20, 0):
            results.append(codon)

        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'],
                    ['C', 'T', 'A'], ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        self.assertEqual(expected, results)


class TestFindCodons(TestSequences):
    """
    For reverse strand orfs, find_codons reverses the sequence, without complementing the nucleotides.
    is_nonsyn() handles finding the complement of the nucleotides.
    """

    def testFwdOrf(self):
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                    ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'],
                    ['A', 'G', 'C']]
        result = self.sequence7.find_codons('+0', (0, 20))
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.get_state())

    def testReverseOrf(self):
        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'],
                    ['C', 'T', 'A'], ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        result = self.sequence7.find_codons('-0', (20, 0))
        for res in result:
            print(res.nts_in_codon)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-0')
            for pos, nt in enumerate(codon.nts_in_codon):
                self.assertEqual(expected[idx][pos], nt.get_state())


class TestGetFrequencyRates(unittest.TestCase):
    """ get_frequency_rates() is a static method, so it can be called without the instance of a Sequence object"""

    def testSameNucleotide(self):
        result = Sequence.get_frequency_rates('AAAAAAAAA')
        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        self.assertEqual(expected, result)

    def testShortSeq(self):
        result = Sequence.get_frequency_rates('GTACGATCGATCGATGCTAGC')
        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        self.assertEqual(expected, result)

    def testLongerSeq(self):
        s = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCA" \
            "GGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACA" \
            "CCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGG" \
            "AGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTC" \
            "AGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTA" \
            "AGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAA" \
            "ATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG"
        result = Sequence.get_frequency_rates(s)
        expected = {'A': 0.25, 'C': 0.25, 'T': 0.22, 'G': 0.28}
        self.assertEqual(expected, result)



# ==========================================
# Tests for DoubleLinkedList
# ==========================================
class TestInsertNt(unittest.TestCase):

    def testSimpleUse(self):
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


class TestNucleotideAtPos(TestSequences):

    def testFindFirst(self):
        result_nt = self.seq1.nucleotide_at_pos(0)
        expected_state = 'A'
        self.assertEqual(expected_state, result_nt.get_state())

    def testFindLast(self):
        result_nt = self.seq2.nucleotide_at_pos(13)
        expected_state = 'G'
        self.assertEqual(expected_state, result_nt.get_state())

    def testFindMiddle(self):
        result_nt = self.seq3.nucleotide_at_pos(6)
        expected_state = 'C'
        self.assertEqual(expected_state, result_nt.get_state())


class TestSliceSequence(TestSequences):

    def testSliceMiddle(self):
        nts = self.seq7.slice_sequence(1, 4)
        expected = ['T', 'A', 'C', 'G']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSliceFromStart(self):
        nts = self.seq7.slice_sequence(0, 3)
        expected = ['G', 'T', 'A', 'C']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSliceToEnd(self):
        nts = self.seq7.slice_sequence(18, 20)
        expected = ['A', 'G', 'C']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSliceOutOfBounds1(self):
        nts = self.seq7.slice_sequence(21, 23)
        expected = []
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSliceOutOfBounds2(self):
        nts = self.seq7.slice_sequence(-1, -5)
        expected = []
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSliceReverse(self):
        nts = self.seq7.slice_sequence(5, 1)
        expected = ['T', 'A', 'C', 'G', 'A']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)

    def testSameStartSameEnd(self):
        nts = self.seq7.slice_sequence(0, 0)
        expected = ['G']
        result = [nt.get_state() for nt in nts]
        self.assertEqual(expected, result)


# ==========================================
# Tests for Nucleotide
# ==========================================
class TestGetNonsynSubs(TestSequences):

    def testNoORFs(self):
        s = Sequence('ATGCCGTATGC', rcseq=None, sorted_orfs=[], mu=0.5, pi=None, kappa=0.3, )
        s.nt_sequence = DoubleLinkedList()
        for pos, nt in enumerate(s.original_seq):
            s.nt_sequence.insert_nt(nt, pos)

        s.nt_sequence.print_seq()

        nt = s.nt_sequence.nucleotide_at_pos(2)

        expected = {'+0': {'A': None, 'T': None, 'G': None, 'C': None},
                    '+1': {'A': None, 'T': None, 'G': None, 'C': None},
                    '+2': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-0': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-1': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-2': {'A': None, 'T': None, 'G': None, 'C': None}}

        result = nt.get_nonsyn_subs()
        self.assertEqual(expected, result)

    def testOneOrf(self):
        s = Sequence('ATGCCCTGA', rcseq=None, sorted_orfs=None, mu=0.5, pi=None, kappa=0.3)
        seq = DoubleLinkedList()
        for pos, nt in enumerate(s.nt_sequence):
            seq.insert_nt(nt, pos)

        nt = s.nt_sequence.nucleotide_at_pos(0)

        expected = {'+0': {'A': False, 'T': True, 'G': True, 'C': True},
                    '+1': {'A': None,  'T': None, 'G': None, 'C': None},
                    '+2': {'A': None,  'T': None, 'G': None, 'C': None},
                    '-0': {'A': None,  'T': None, 'G': None, 'C': None},
                    '-1': {'A': None,  'T': None, 'G': None, 'C': None},
                    '-2': {'A': None,  'T': None, 'G': None, 'C': None}}

        result = nt.get_nonsyn_subs()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
