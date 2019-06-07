import unittest

from sequence_info import Sequence


class TestReverseAndComplement(unittest.TestCase):
    """
    Tests reverse_and_complement
    """
    def testLowerCaseInputSimple(self):
        s2 = Sequence("atgcatgcatgc")
        expected = "GCATGCATGCAT"
        result = Sequence.reverse_and_complement(s2)
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        s3 = Sequence("ATGAAATTTGGGCCCTAA")
        expected = "TTAGGGCCCAAATTTCAT"
        result = Sequence.reverse_and_complement(s3)
        self.assertEqual(expected, result)


class TestValidSequence(unittest.TestCase):
    """
    Tests valid_sequence
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


class TestInvalidORFs(unittest.TestCase):

    def testBadInput(self):
        orfs = "90, 1010"
        expected = False
        result = Sequence.valid_orfs(orfs)
        self.assertEqual(expected, result)

    def testNotTuple(self):
        orfs = (80, 79)
        expected = False
        result = Sequence.valid_orfs(orfs)
        self.assertEqual(expected, result)

    def testCorrectInput(self):
        orfs = [(0, 16)]
        expected = True
        result = Sequence.valid_orfs(orfs)
        self.assertEqual(expected, result)


# class TestGetCodon(unittest.TestCase):
#     """
#     Tests get_codon
#     """
#     def testSimpleUse(self):
#         s1 = Sequence("ATGAAAGGGTAA", [(0, 11)])
#         expected = ("ATG", 0)
#         result = Sequence.get_codon(s1, s1.original_seq, 0, s1.orfs)
#         self.assertEqual(expected, result)
#
#     def testSimpleUse2(self):
#         expected = ("ATG", 1)
#         result = get_codon("ATGAAAGGGTAA", 1, (0, 11))
#         self.assertEqual(expected, result)
#
#     def testSimpleUse3(self):
#         expected = ("GGG", 2)
#         result = get_codon("ATGAAAGGGTAA", 8, (0, 11))
#         self.assertEqual(expected, result)
#
#     def testSamePositionMultipleOrfs1(self):
#         expected = ("CAT", 1)
#         result = get_codon("ATGAAAGTGCAACATGGGTAAATAG", 13, (0, 20))
#         self.assertEqual(expected, result)
#
#     def testSamePositionMultipleOrfs2(self):
#         expected = ("ATG", 0)
#         result = get_codon("ATGAAAGTGCAACATGGGTAAATAG", 13, (13, 24))
#         self.assertEqual(expected, result)
#
#     def testFindStopCodon(self):
#         expected = ("TAA", 2)
#         result = get_codon("ATGAAAGGGTAA", 11, (0, 11))
#         self.assertEqual(expected, result)


class TestGetReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames
    """
    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        s1 = Sequence("ATGAAATAG")
        expected = [(0, 8)]
        result = Sequence.get_reading_frames(s1)
        self.assertEqual(expected, result)

    def testNoStartNoStop(self):
        s2 = Sequence("TTTTTTTTTTTTT")
        expected = []
        result = Sequence.get_reading_frames(s2)
        self.assertEqual(expected, result)

    def testOnlyStartCodon(self):
        s3 = Sequence("ATGCCTCTCTCTCTTCTCTC")
        expected = []
        result = Sequence.get_reading_frames(s3)
        self.assertEqual(expected, result)

    def testOnlyStopCodon(self):
        s4 = Sequence("AGAACGTAA")
        expected = []
        result = Sequence.get_reading_frames(s4)
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        # Tests case when the start codon is in the middle of the sequence
        # ORF (forward direction): ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATC TAG
        # ORF (reverse direction): ATG AAG CGA ACA GAT TTT CGT TCA TGA
        s5 = Sequence("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        expected = [(5, 49), (29, 3)]
        result = Sequence.get_reading_frames(s5)
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        # Tests scenario when Start and Stop codons are present, but they are not in the same frame
        # ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATT AG...
        s6 = Sequence("AAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATTAGGCCTACCC")
        expected = []
        result = Sequence.get_reading_frames(s6)
        self.assertEqual(expected, result)

    # def testTwoStartSameFrame(self):
    #     # Tests scenario when two ATGs are present in the same reading frame
    #     # ORF: ATG CCC ATG CCC TAA TAA
    #     expected = [(0, 14)]
    #     result = get_reading_frames("ATGCCCATGCCCTAATAA")
    #     self.assertEqual(expected, result)
    #
    # def testMultipleORFs(self):
    #     # Tests scenario when there are 2 different ORFs
    #     # ORF1: ATG AAA GTG CAA CAT GGG TAA
    #     # ORF2: ATG GGT AAA TAG
    #     expected = [(0, 20), (13, 24)]
    #     result = get_reading_frames("ATGAAAGTGCAACATGGGTAAATAG")
    #     self.assertEqual(expected, result)
    #
    # def testLowerCaseInput(self):
    #     expected = [(0, 8)]
    #     result = get_reading_frames("atgaaatag")
    #     self.assertEqual(expected, result)
    #
    # def testBadInput(self):
    #     expected = []
    #     result = get_reading_frames("78979GTC")
    #     self.assertEqual(expected, result)
    #
    # def testStopBeforeStart(self):
    #     # Tests scenario when a stop codon precedes a start codon
    #     expected = [(3, 11)]
    #     result = get_reading_frames("TAGATGAAATAG")
    #     self.assertEqual(expected, result)
    #
    # def testBackToBackORFs1(self):
    #     # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
    #     # ORF1: ATG TTT TAA
    #     # ORF2: ATG CAC TAA
    #     expected = [(0, 8), (9, 17)]
    #     result = get_reading_frames("ATGTTTTAGATGCACTAA")
    #     self.assertEqual(expected, result)
    #
    # def testBackToBackORFs2(self):
    #     # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
    #     # ORF1: ATG TTT TGA
    #     # ORF2: ATG CAC TAA
    #     expected = [(0, 8), (15, 23)]
    #     result = get_reading_frames("ATGTTTTGACCCAAAATGCACTAA")
    #     self.assertEqual(expected, result)
    #
    # def testBackToBackORFs3(self):
    #     # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
    #     # ORF1: ATG TTT TGA
    #     # ORF2: ATG CAC TAA
    #     expected = [(0, 8), (14, 22)]
    #     result = get_reading_frames("ATGTTTTGACCCAAATGCACTAA")
    #     self.assertEqual(expected, result)
    #
    # def testBackToBackORFs4(self):
    #     # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
    #     # ORF1: ATG GAG TGA
    #     # ORF2: ATG GAG TGA
    #     expected = [(1, 9), (16, 24)]
    #     result = get_reading_frames("AATGGAGTGACCCGGGATGGAGTAG")
    #     self.assertEqual(expected, result)
    #
    # def testInternalORF(self):
    #     # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
    #     # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
    #     # ORF2: ATG GCA CAA GTG TAA
    #     expected = [(0, 23), (5, 19)]
    #     result = get_reading_frames("ATGAGATGGCACAAGTGTAACTAG")
    #     self.assertEqual(expected, result)
    #
    # def testStartStop(self):
    #     # Tests a scenario when a STOP codon immediately follows a start codon
    #     expected = []
    #     result = get_reading_frames("ATGTAA")
    #     self.assertEqual(expected, result)
    #

