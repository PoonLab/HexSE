import unittest

from src.run_simulation import get_open_reading_frames
from src.run_simulation import reverse_and_complement
from src.run_simulation import sort_orfs
from src.run_simulation import valid_orfs
from src.run_simulation import valid_sequence


class TestValidSequence(unittest.TestCase):
    """
    Tests valid_sequence
    """
    def testBadInput(self):
        seq = "SGDKLJFAHHH"
        expected = False
        result = valid_sequence(seq)
        self.assertEqual(expected, result)

    def testLowerCase(self):
        seq = "atgcatgcatgcatcg"
        expected = True
        result = valid_sequence(seq)
        self.assertEqual(expected, result)

    def testShortSequence(self):
        seq = "TGG"
        expected = False
        result = valid_sequence(seq)
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        seq = "TGCAGACATAGACAGATCAG"
        expected = True
        result = valid_sequence(seq)
        self.assertEqual(expected, result)


class TestValidORFs(unittest.TestCase):
    """
    Tests valid_orfs
    """
    def testCorrectInput(self):
        s = "ATGCGTAAACGGGCTAGAGCTAGCA"
        orfs = [(0, 8)]
        expected = []
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)

    def testNotOrf(self):
        s = "ATGCGCGCATGACGA"
        orfs = [(1, 1)]
        expected = [(1, 1)]
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)

    def testInvalidFormat(self):
        s = "ATGCGCGATGACGA"
        orfs = [(1, 7, "j")]
        expected = [(1, 7, "j")]
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)

    def testInvalidFormat2(self):
        s = "ATGTCGATGCATGC"
        orfs = [(1, 2, 4)]
        expected = [(1, 2, 4)]
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)

    def testOrfTooShort(self):
        s = "ATGTCGATGCATGC"
        orfs = [(0, 3)]
        expected = [(0, 3)]
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)

    def testNotMultipleOfThree(self):
        s = "ATGTCGATGCATGC"
        orfs = [(1, 2)]
        expected = [(1, 2)]
        result = valid_orfs(orfs, s)
        self.assertEqual(expected, result)


class TestReverseAndComplement(unittest.TestCase):
    """
    Tests reverse_and_complement
    """
    def testLowerCaseInputSimple(self):
        s = "atgcatgcatgc"
        expected = "GCATGCATGCAT"
        result = reverse_and_complement(s)
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        s = "ATGAAATTTGGGCCCTAA"
        expected = "TTAGGGCCCAAATTTCAT"
        result = reverse_and_complement(s)
        self.assertEqual(expected, result)


class TestGetOpenReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames
    """
    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        s = "ATGAAATAG"
        expected = [(0, 8)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testNoStartNoStop(self):
        s = "TTTTTTTTTTTTT"
        expected = []
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testOnlyStartCodon(self):
        s = "ATGCCTCTCTCTCTTCTCTC"
        expected = []
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testOnlyStopCodon(self):
        s = "AGAACGTAA"
        expected = []
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        # Tests case when the start codon is in the middle of the sequence
        # ORF (forward direction): ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATC TAG
        # ORF (reverse direction): ATG AAG CGA ACA GAT TTT CGT TCA TGA
        s = "AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC"
        expected = [(5, 49), (29, 3)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        # Tests scenario when Start and Stop codons are present, but they are not in the same frame
        # ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATT AG...
        s = "AAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATTAGGCCTACCC"
        expected = []
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testTwoStartSameFrame(self):
        # Tests scenario when two ATGs are present in the same reading frame
        # ORF: ATG CCC ATG CCC TAA TAA
        s = "ATGCCCATGCCCTAATAA"
        expected = [(0, 14)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        # Tests scenario when there are 2 different ORFs
        # ORF1: ATG AAA GTG CAA CAT GGG TAA
        # ORF2: ATG GGT AAA TA
        s = "ATGAAAGTGCAACATGGGTAAATAG"
        expected = [(0, 20), (13, 24)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        s = "atgaaatag"
        expected = [(0, 8)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testStopBeforeStart(self):
        # Tests scenario when a stop codon precedes a start codon
        s = "TAGATGAAATAG"
        expected = [(3, 11)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs1(self):
        # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
        # ORF1: ATG TTT TAA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTAGATGCACTAA"
        expected = [(0, 8), (9, 17)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs2(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTGACCCAAAATGCACTAA"
        expected = [(0, 8), (15, 23)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs3(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTGACCCAAATGCACTAA"
        expected = [(0, 8), (14, 22)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs4(self):
        # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
        # ORF1: ATG GAG TGA
        # ORF2: ATG GAG TGA
        s = "AATGGAGTGACCCGGGATGGAGTAG"
        expected = [(1, 9), (16, 24)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testInternalORF(self):
        # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
        # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
        # ORF2: ATG GCA CAA GTG TAA
        s = "ATGAGATGGCACAAGTGTAACTAG"
        expected = [(0, 23), (5, 19)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)


class TestSortOrfs(unittest.TestCase):
    """
    Tests sort_orfs
    """
    def testSimpleUse(self):
        expected = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = sort_orfs([(0, 11)])
        self.assertEqual(expected, result)

    def testMultipleOrfs(self):
        expected = {'+0': [(0, 23)], '+1': [], '+2': [(5, 19)], '-0': [], '-1': [], '-2': []}
        result = sort_orfs([(0, 23), (5, 19)])
        self.assertEqual(expected, result)

    def testBacktoBackOrfs(self):
        expected = {'+0': [(1, 9), (16, 24)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        result = sort_orfs([(1, 9), (16, 24)])
        self.assertEqual(expected, result)

    def testFwdReverseOrfs(self):
        expected = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        result = sort_orfs([(5, 49), (29, 3)])
        self.assertEqual(expected, result)

    def testAllSixOrfs(self):
        expected = {'+0': [(0, 8)], '+1': [(1, 9)], '+2': [(2, 10)],
                    '-0': [(8, 0)], '-1': [(10, 2)], '-2': [(9, 1)]}
        result = sort_orfs([(0, 8), (8, 0), (1, 9), (9, 1), (2, 10), (10, 2)])
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
