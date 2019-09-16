import unittest

from refact.new_sequence_info import Sequence
from refact.new_sequence_info import DoubleLinkedList
from refact.new_sequence_info import Nucleotide
from refact.new_sequence_info import Codon


def make_dll(original_sequence):
    s = Sequence(original_sequence, kappa=None, mu=None)
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
        self.sequence5 = s1


class TestValidSequence(unittest.TestCase):
    """
    Tests valid_sequence from Sequence class
    """

    def testBadInput(self):
        seq = "SGDKLJFAHHH"
        expected = False
        result = Sequence.valid_sequence(seq)
        self.assertEqual(expected, result)

    def testLowerCase(self):
        seq = "atgcatgcatgcatcg"
        expected = True
        result = Sequence.valid_sequence(seq)
        self.assertEqual(expected, result)

    def testShortSequence(self):
        seq = "TGG"
        expected = False
        result = Sequence.valid_sequence(seq)
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        seq = "TGCAGACATAGACAGATCAG"
        expected = True
        result = Sequence.valid_sequence(seq)
        self.assertEqual(expected, result)


class TestValidORFs(unittest.TestCase):
    """
    Tests valid_orfs from Sequence class
    """

    def testBadInput(self):
        s = Sequence("ATGCATGCATGC", kappa=None, mu=None)
        orf = "90, 1010"
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testNotTuple(self):
        s = Sequence("AAAAAAAATGCGTCGA", kappa=None, mu=None)
        orf = (9, 15)
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testCorrectInput(self):
        s = Sequence("ATGCGTAAACGGGCTAGAGCTAGCA", kappa=None, mu=None)
        orfs = [(0, 8)]
        expected = True
        result = s.valid_orfs(orfs)
        self.assertEqual(expected, result)

    def testNotOrfs(self):
        s = Sequence("ATGCGCGCATGACGA", kappa=None, mu=None)
        orf = [(1, 1)]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testInvalidFormat(self):
        s = Sequence("ATGCGCGATGACGA", kappa=None, mu=None)
        orf = [(1, 7, "j")]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testInvalidFormat2(self):
        s = Sequence("ATGTCGATGCATGC", kappa=None, mu=None)
        orf = [(1, 2, 4)]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testNotIntegers(self):
        s = Sequence("ATCGATCGATGC", kappa=None, mu=None)
        orf = [("8", "89")]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)


class TestReverseAndComplement(unittest.TestCase):
    """
    Tests reverse_and_complement from Sequence class
    """

    def testLowerCaseInputSimple(self):
        s = Sequence("atgcatgcatgc", kappa=None, mu=None)
        expected = "GCATGCATGCAT"
        result = s.reverse_and_complement()
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        s = Sequence("ATGAAATTTGGGCCCTAA", kappa=None, mu=None)
        expected = "TTAGGGCCCAAATTTCAT"
        result = s.reverse_and_complement()
        self.assertEqual(expected, result)


class TestGetOpenReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames from Sequence class
    """

    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        s = Sequence("ATGAAATAG", kappa=None, mu=None)
        expected = [(0, 8)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testNoStartNoStop(self):
        s = Sequence("TTTTTTTTTTTTT", kappa=None, mu=None)
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testOnlyStartCodon(self):
        s = Sequence("ATGCCTCTCTCTCTTCTCTC", kappa=None, mu=None)
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testOnlyStopCodon(self):
        s = Sequence("AGAACGTAA", kappa=None, mu=None)
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        # Tests case when the start codon is in the middle of the sequence
        # ORF (forward direction): ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATC TAG
        # ORF (reverse direction): ATG AAG CGA ACA GAT TTT CGT TCA TGA
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC", kappa=None, mu=None)
        expected = [(5, 49), (29, 3)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        # Tests scenario when Start and Stop codons are present, but they are not in the same frame
        # ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATT AG...
        s = Sequence("AAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATTAGGCCTACCC", kappa=None, mu=None)
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testTwoStartSameFrame(self):
        # Tests scenario when two ATGs are present in the same reading frame
        # ORF: ATG CCC ATG CCC TAA TAA
        s = Sequence("ATGCCCATGCCCTAATAA", kappa=None, mu=None)
        expected = [(0, 14)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        # Tests scenario when there are 2 different ORFs
        # ORF1: ATG AAA GTG CAA CAT GGG TAA
        # ORF2: ATG GGT AAA TA
        s = Sequence("ATGAAAGTGCAACATGGGTAAATAG", kappa=None, mu=None)
        expected = [(0, 20), (13, 24)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        s = Sequence("atgaaatag", kappa=None, mu=None)
        expected = [(0, 8)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testStopBeforeStart(self):
        # Tests scenario when a stop codon precedes a start codon
        s = Sequence("TAGATGAAATAG", kappa=None, mu=None)
        expected = [(3, 11)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs1(self):
        # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
        # ORF1: ATG TTT TAA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTAGATGCACTAA", kappa=None, mu=None)
        expected = [(0, 8), (9, 17)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs2(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTGACCCAAAATGCACTAA", kappa=None, mu=None)
        expected = [(0, 8), (15, 23)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs3(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTGACCCAAATGCACTAA", kappa=None, mu=None)
        expected = [(0, 8), (14, 22)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs4(self):
        # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
        # ORF1: ATG GAG TGA
        # ORF2: ATG GAG TGA
        s = Sequence("AATGGAGTGACCCGGGATGGAGTAG", kappa=None, mu=None)
        expected = [(1, 9), (16, 24)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testInternalORF(self):
        # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
        # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
        # ORF2: ATG GCA CAA GTG TAA
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG", kappa=None, mu=None)
        expected = [(0, 23), (5, 19)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)


class TestSortOrfs(unittest.TestCase):
    """
    Tests sort_orfs from Sequence class
    """

    def testSimpleUse(self):
        s = Sequence("ATGAAAGGGTAA", kappa=None, mu=None)
        expected = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(0, 11)])
        self.assertEqual(expected, result)

    def testMultipleOrfs(self):
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG", kappa=None, mu=None)
        expected = {'+0': [(0, 23)], '+1': [], '+2': [(5, 19)], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(0, 23), (5, 19)])
        self.assertEqual(expected, result)

    def testBacktoBackOrfs(self):
        s = Sequence("AATGGAGTGACCCGGGATGGAGTAG", kappa=None, mu=None)
        expected = {'+0': [(1, 9), (16, 24)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(1, 9), (16, 24)])
        self.assertEqual(expected, result)

    def testFwdReverseOrfs(self):
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC", kappa=None, mu=None)
        expected = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        result = s.sort_orfs([(5, 49), (29, 3)])
        self.assertEqual(expected, result)

    def testAllSixOrfs(self):
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC", kappa=None, mu=None)
        expected = {'+0': [(0, 8)], '+1': [(1, 9)], '+2': [(2, 10)],
                    '-0': [(8, 0)], '-1': [(10, 2)], '-2': [(9, 1)]}
        result = s.sort_orfs([(0, 8), (8, 0), (1, 9), (9, 1), (2, 10), (10, 2)])
        self.assertEqual(expected, result)


class TestInsertNt(unittest.TestCase):
    """
    Tests insert_nt from DoubleLinkedList class
    """
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
    """
    Tests nucleotide_at_pos from DoubleLinkedList class
    """
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


class TestGetNonsynSubs(TestSequences):
    """
    Tests get_nonsyn_subs from the Nucleotide class
    """

    def testNoORFs(self):
        s = Sequence('ATGCCGTATGC', kappa=None, mu=None)
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
        s = Sequence('ATGCCCTGA', kappa=None, mu=None)
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
