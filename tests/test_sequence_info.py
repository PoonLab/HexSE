import unittest

from refact.new_sequence_info import Sequence
from refact.new_sequence_info import DoubleLinkedList
from refact.new_sequence_info import Nucleotide
from refact.new_sequence_info import Codon


def make_dll(original_sequence):
    s = Sequence(original_sequence, rcseq=None, sorted_orfs=None, mu=None, pi=None, kappa=None)
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


class TestIsTransversion(TestSequences):
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

    def testForwardStrand(self):
        results = []
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        for codon in self.sequence7.codon_iterator(orf, 0, 21):
            results.append(codon)

        expected =[['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                   ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]
        self.assertEqual(expected, results)

    def testReverseStrand(self):
        results = []
        orf = ['G', 'T', 'A', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'C', 'G', 'A', 'T', 'G', 'C', 'T', 'A', 'G', 'C']
        for codon in self.sequence7.codon_iterator(orf, 21, 0):
            results.append(codon)

        expected = [['C', 'G', 'A'], ['T', 'C', 'G'], ['T', 'A', 'G'],
                    ['C', 'T', 'A'], ['G', 'C', 'T'], ['A', 'G', 'C'], ['A', 'T', 'G']]
        self.assertEqual(expected, results)


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


# ==========================================
# Tests for Nucleotide
# ==========================================
class TestGetNonsynSubs(TestSequences):

    def testNoORFs(self):
        s = Sequence('ATGCCGTATGC', rcseq=None, sorted_orfs=[], pi=None, kappa=None, mu=None)
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
        s = Sequence('ATGCCCTGA', rcseq=None, sorted_orfs=None, mu=None, pi=None, kappa=None)
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
