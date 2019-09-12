import unittest

from refact.new_sequence_info import DoubleLinkedList
from refact.new_sequence_info import Nucleotide
from refact.new_sequence_info import Sequence


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
        s = Sequence("ATGCATGCATGC")
        orf = "90, 1010"
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testNotTuple(self):
        s = Sequence("AAAAAAAATGCGTCGA")
        orf = (9, 15)
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testCorrectInput(self):
        s = Sequence("ATGCGTAAACGGGCTAGAGCTAGCA")
        orfs = [(0, 8)]
        expected = True
        result = s.valid_orfs(orfs)
        self.assertEqual(expected, result)

    def testNotOrfs(self):
        s = Sequence("ATGCGCGCATGACGA")
        orf = [(1, 1)]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testInvalidFormat(self):
        s = Sequence("ATGCGCGATGACGA")
        orf = [(1, 7, "j")]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testInvalidFormat2(self):
        s = Sequence("ATGTCGATGCATGC")
        orf = [(1, 2, 4)]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)

    def testNotIntegers(self):
        s = Sequence("ATCGATCGATGC")
        orf = [("8", "89")]
        expected = False
        result = s.valid_orfs(orf)
        self.assertEqual(expected, result)


class TestReverseAndComplement(unittest.TestCase):
    """
    Tests reverse_and_complement from Sequence class
    """

    def testLowerCaseInputSimple(self):
        s = Sequence("atgcatgcatgc")
        expected = "GCATGCATGCAT"
        result = s.reverse_and_complement()
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        s = Sequence("ATGAAATTTGGGCCCTAA")
        expected = "TTAGGGCCCAAATTTCAT"
        result = s.reverse_and_complement()
        self.assertEqual(expected, result)


class TestGetOpenReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames from Sequence class
    """

    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        s = Sequence("ATGAAATAG")
        expected = [(0, 8)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testNoStartNoStop(self):
        s = Sequence("TTTTTTTTTTTTT")
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testOnlyStartCodon(self):
        s = Sequence("ATGCCTCTCTCTCTTCTCTC")
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testOnlyStopCodon(self):
        s = Sequence("AGAACGTAA")
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        # Tests case when the start codon is in the middle of the sequence
        # ORF (forward direction): ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATC TAG
        # ORF (reverse direction): ATG AAG CGA ACA GAT TTT CGT TCA TGA
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        expected = [(5, 49), (29, 3)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        # Tests scenario when Start and Stop codons are present, but they are not in the same frame
        # ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATT AG...
        s = Sequence("AAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATTAGGCCTACCC")
        expected = []
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testTwoStartSameFrame(self):
        # Tests scenario when two ATGs are present in the same reading frame
        # ORF: ATG CCC ATG CCC TAA TAA
        s = Sequence("ATGCCCATGCCCTAATAA")
        expected = [(0, 14)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        # Tests scenario when there are 2 different ORFs
        # ORF1: ATG AAA GTG CAA CAT GGG TAA
        # ORF2: ATG GGT AAA TA
        s = Sequence("ATGAAAGTGCAACATGGGTAAATAG")
        expected = [(0, 20), (13, 24)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        s = Sequence("atgaaatag")
        expected = [(0, 8)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testStopBeforeStart(self):
        # Tests scenario when a stop codon precedes a start codon
        s = Sequence("TAGATGAAATAG")
        expected = [(3, 11)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs1(self):
        # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
        # ORF1: ATG TTT TAA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTAGATGCACTAA")
        expected = [(0, 8), (9, 17)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs2(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTGACCCAAAATGCACTAA")
        expected = [(0, 8), (15, 23)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs3(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = Sequence("ATGTTTTGACCCAAATGCACTAA")
        expected = [(0, 8), (14, 22)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testBackToBackORFs4(self):
        # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
        # ORF1: ATG GAG TGA
        # ORF2: ATG GAG TGA
        s = Sequence("AATGGAGTGACCCGGGATGGAGTAG")
        expected = [(1, 9), (16, 24)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)

    def testInternalORF(self):
        # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
        # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
        # ORF2: ATG GCA CAA GTG TAA
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG")
        expected = [(0, 23), (5, 19)]
        result = s.get_open_reading_frames()
        self.assertEqual(expected, result)


class TestSortOrfs(unittest.TestCase):
    """
    Tests sort_orfs from Sequence class
    """

    def testSimpleUse(self):
        s = Sequence("ATGAAAGGGTAA")
        expected = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(0, 11)])
        self.assertEqual(expected, result)

    def testMultipleOrfs(self):
        s = Sequence("ATGAGATGGCACAAGTGTAACTAG")
        expected = {'+0': [(0, 23)], '+1': [], '+2': [(5, 19)], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(0, 23), (5, 19)])
        self.assertEqual(expected, result)

    def testBacktoBackOrfs(self):
        s = Sequence("AATGGAGTGACCCGGGATGGAGTAG")
        expected = {'+0': [(1, 9), (16, 24)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = s.sort_orfs([(1, 9), (16, 24)])
        self.assertEqual(expected, result)

    def testFwdReverseOrfs(self):
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        expected = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        result = s.sort_orfs([(5, 49), (29, 3)])
        self.assertEqual(expected, result)

    def testAllSixOrfs(self):
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        expected = {'+0': [(0, 8)], '+1': [(1, 9)], '+2': [(2, 10)],
                    '-0': [(8, 0)], '-1': [(10, 2)], '-2': [(9, 1)]}
        result = s.sort_orfs([(0, 8), (8, 0), (1, 9), (9, 1), (2, 10), (10, 2)])
        self.assertEqual(expected, result)


class TestCreateNtORFDict(unittest.TestCase):
    """
    Tests create_nt_orf_dict from Sequence class
    """

    def testFirstNt(self):
        s = Sequence("ATGAAAGGGTAA")
        expected = {'+0': 0, '+1': None, '+2': None, '-0': None, '-1': None, '-2': None}
        result = s.create_nt_orf_dict(0)
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        s = Sequence('ATGAGATGGCACAAGTGTAACTAG')
        expected = {'+0': 2, '+1': None, '+2': 0, '-0': None, '-1': None, '-2': None}
        result = s.create_nt_orf_dict(5)
        self.assertEqual(expected, result)

    def testLastNt(self):
        s = Sequence('ATGAGATGGCACAAGTGTAACTAG')
        expected = {'+0': 2, '+1': None, '+2': None, '-0': None, '-1': None, '-2': None}
        result = s.create_nt_orf_dict(23)
        self.assertEqual(expected, result)

    def testFwdRevORFs(self):
        s = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        expected = {'+0': 0, '+1': None, '+2': None, '-0': None, '-1': None, '-2': 0}
        result = s.create_nt_orf_dict(5)
        self.assertEqual(expected, result)


class TestInsertNt(unittest.TestCase):
    """
    Tests insert_nt from DoubleLinkedList class
    """

    def testSimpleUse(self):
        seq = DoubleLinkedList()

        # Test adding first nucleotide
        seq.insert_nt('A', 0, {'+0': 0})
        exp_head = Nucleotide('A', 0, {'+0': 0})
        exp_current = exp_head

        self.assertEqual(exp_head.state, seq.head.letter)
        self.assertEqual(exp_head.pos, seq.head.pos)
        self.assertEqual(exp_head.pos_in_codons, seq.head.pos_in_codons)

        self.assertEqual(exp_current.state, seq.current_nt.letter)
        self.assertEqual(exp_current.pos, seq.current_nt.pos)
        self.assertEqual(exp_current.pos_in_codons, seq.current_nt.pos_in_codons)

        # Test adding second nucleotide
        seq.insert_nt('T', 1, {'+0': 1})
        exp_current = Nucleotide('T', 1, {'+0': 1})

        self.assertEqual(exp_head.state, seq.head.letter)
        self.assertEqual(exp_head.pos, seq.head.pos)
        self.assertEqual(exp_head.pos_in_codons, seq.head.pos_in_codons)

        self.assertEqual(exp_current.state, seq.current_nt.letter)
        self.assertEqual(exp_current.pos, seq.current_nt.pos)
        self.assertEqual(exp_current.pos_in_codons, seq.current_nt.pos_in_codons)


class TestGetNonsynSubs(unittest.TestCase):
    """
    Tests get_nonsyn_subs from the Nucleotide class
    """

    def testNoORFs(self):
        s = Sequence('AAAAAAAAAAA')
        seq = DoubleLinkedList()
        for pos, nt in enumerate(s.nt_sequence):
            pos_in_codon = s.create_nt_orf_dict(pos)
            seq.insert_nt(nt, pos, pos_in_codon)

        nt = s.nt_sequence[2]

        expected = {'+0': {'A': None, 'T': None, 'G': None, 'C': None},
                    '+1': {'A': None, 'T': None, 'G': None, 'C': None},
                    '+2': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-0': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-1': {'A': None, 'T': None, 'G': None, 'C': None},
                    '-2': {'A': None, 'T': None, 'G': None, 'C': None}}

        result = nt.get_nonsyn_subs()
        self.assertEqual(expected, result)

    def testOneOrf(self):

        s = Sequence('ATGCCCTGA')
        seq = DoubleLinkedList()
        for pos, nt in enumerate(s.nt_sequence):
            pos_in_codon = s.create_nt_orf_dict(pos)
            seq.insert_nt(nt, pos, pos_in_codon)

        nt = s.nt_sequence[0]

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
