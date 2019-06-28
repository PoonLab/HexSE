import unittest

from src.ovrf_functions import get_frequency_rates
from src.ovrf_functions import draw_omega_values
from src.ovrf_functions import get_reading_frames
from src.ovrf_functions import codon_iterator
from src.ovrf_functions import get_syn_codons
from src.ovrf_functions import get_syn_subs
from src.ovrf_functions import get_codon
from src.ovrf_functions import reverse_and_complement
from src.ovrf_functions import sort_orfs


class TestGetFrequencyRates(unittest.TestCase):
    """
    Tests get_frequency_rates
    """
    def testSimpleUse(self):
        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        result = get_frequency_rates("AAAAAAAAA")
        self.assertEqual(expected, result)

    def testSimpleUse1(self):
        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        result = get_frequency_rates("GTACGATCGATCGATGCTAGC")
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        seq = "tggaagggctaattcactcccaacgaagacaagatatccttgatctgtgg" \
              "atctaccacacacaaggctacttccctgattagcagaactacacaccagg" \
              "gccagggatcagatatccactgacctttggatggtgctacaagctagtac" \
              "cagttgagccagagaagttagaagaagccaacaaaggagagaacaccagc" \
              "ttgttacaccctgtgagcctgcatggaatggatgacccggagagagaagt" \
              "gttagagtggaggtttgacagccgcctagcatttcatcacatggcccgag" \
              "agctgcatccggagtacttcaagaactgctgacatcgagcttgctacaag" \
              "ggactttccgctggggactttccagggaggcgtggcctgggcgggactgg" \
              "ggagtggcgagccctcagatcctgcatataagcagctgctttttgcctgt" \
              "actgggtctctctggttagaccagatctgagcctgggagctctctggcta" \
              "actagggaacccactgcttaagcctcaataaagcttgccttgagtgcttc" \
              "aagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc" \
              "agacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagg" \
              "gacctgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcg"
        expected = {'A': 0.25, 'C': 0.25, 'T': 0.22, 'G': 0.28}
        result = get_frequency_rates(seq.upper())
        self.assertEqual(expected, result)


class TestDrawOmegaValues(unittest.TestCase):
    """
    Tests draw_omega_values
    """
    def testEmptyString(self):
        return None


class TestGetReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames
    """
    def testEmptyString(self):
        expected = []
        result = get_reading_frames("")
        self.assertEqual(expected, result)

    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        expected = [(0, 8)]
        result = get_reading_frames("ATGAAATAG")
        self.assertEqual(expected, result)

    def testNoStartNoStop(self):
        expected = []
        result = get_reading_frames("TTTTTTTTTTTTT")
        self.assertEqual(expected, result)

    def testOnlyStartCodon(self):
        expected = []
        result = get_reading_frames("ATGCCTCTCTCTCTTCTCTC")
        self.assertEqual(expected, result)

    def testOnlyStopCodon(self):
        expected = []
        result = get_reading_frames("AGACTAA")
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        # Tests case when the start codon is in the middle of the sequence
        # ORF: ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATC TAG
        expected = [(5, 49)]
        result = get_reading_frames("AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC")
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        # Tests scenario when Start and Stop codons are present, but they are not in the same frame
        # ATG AAC GAA AAT CTG TTC GCT TCA TTC ATT GCC CCC ACA ATT AG...
        expected = []
        result = get_reading_frames("AAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATTAGGCCTACCC")
        self.assertEqual(expected, result)

    def testTwoStartSameFrame(self):
        # Tests scenario when two ATGs are present in the same reading frame
        # ORF: ATG CCC ATG CCC TAA TAA
        expected = [(0, 14)]
        result = get_reading_frames("ATGCCCATGCCCTAATAA")
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        # Tests scenario when there are 2 different ORFs
        # ORF1: ATG AAA GTG CAA CAT GGG TAA
        # ORF2: ATG GGT AAA TAG
        expected = [(0, 20), (13, 24)]
        result = get_reading_frames("ATGAAAGTGCAACATGGGTAAATAG")
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        expected = [(0, 8)]
        result = get_reading_frames("atgaaatag")
        self.assertEqual(expected, result)

    def testBadInput(self):
        expected = []
        result = get_reading_frames("78979GTC")
        self.assertEqual(expected, result)

    def testStopBeforeStart(self):
        # Tests scenario when a stop codon precedes a start codon
        expected = [(3, 11)]
        result = get_reading_frames("TAGATGAAATAG")
        self.assertEqual(expected, result)

    def testBackToBackORFs1(self):
        # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
        # ORF1: ATG TTT TAA
        # ORF2: ATG CAC TAA
        expected = [(0, 8), (9, 17)]
        result = get_reading_frames("ATGTTTTAGATGCACTAA")
        self.assertEqual(expected, result)

    def testBackToBackORFs2(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        expected = [(0, 8), (15, 23)]
        result = get_reading_frames("ATGTTTTGACCCAAAATGCACTAA")
        self.assertEqual(expected, result)

    def testBackToBackORFs3(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        expected = [(0, 8), (14, 22)]
        result = get_reading_frames("ATGTTTTGACCCAAATGCACTAA")
        self.assertEqual(expected, result)

    def testBackToBackORFs4(self):
        # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
        # ORF1: ATG GAG TGA
        # ORF2: ATG GAG TGA
        expected = [(1, 9), (16, 24)]
        result = get_reading_frames("AATGGAGTGACCCGGGATGGAGTAG")
        self.assertEqual(expected, result)

    def testInternalORF(self):
        # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
        # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
        # ORF2: ATG GCA CAA GTG TAA
        expected = [(0, 23), (5, 19)]
        result = get_reading_frames("ATGAGATGGCACAAGTGTAACTAG")
        self.assertEqual(expected, result)

    def testStartStop(self):
        # Tests a scenario when a STOP codon immediately follows a start codon
        expected = []
        result = get_reading_frames("ATGTAA")
        self.assertEqual(expected, result)


class TestCodonIterator(unittest.TestCase):
    """
    Tests codon_iterator
    """
    def testEmptyString(self):
        expected = []
        result = [i for i in codon_iterator("")]
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        expected = [("ATG", 0), ("SWS", 3), ("-MS", 6), ("TAG", 9)]
        result = [i for i in codon_iterator("ATGSWS-MSTAG")]
        self.assertEqual(expected, result)

    def testSeqLength8(self):
        expected = [("CTA", 0), ("GTC", 3), ("AT", 6)]
        result = [i for i in codon_iterator("CTAGTCAT")]
        self.assertEqual(expected, result)

    def testLowerCase(self):
        expected = [("atg", 0), ("cat", 3), ("gca", 6), ("tgc", 9), ("atg", 12), ("c", 15)]
        result = [i for i in codon_iterator("atgcatgcatgcatgc")]
        self.assertEqual(expected, result)


class TestGetSynCodons(unittest.TestCase):
    """
    Tests get_syn_codons
    """
    def testSimpleUse1(self):
        expected = ["AAA", "AAG"]
        result = get_syn_codons("AAA")
        self.assertEqual(result, expected)

    def testLowerCase(self):
        expected = ["AAA", "AAG"]
        result = get_syn_codons("aaa")
        self.assertEqual(result, expected)

    def testEmptyString(self):
        with self.assertRaises(KeyError):
            get_syn_codons("")

    def testPartialCodon(self):
        with self.assertRaises(KeyError):
            get_syn_codons("AG")

    def testFourCharCodon(self):
        with self.assertRaises(KeyError):
            get_syn_codons(">Seq")

    def testInvalidCodon(self):
        with self.assertRaises(KeyError):
            get_syn_codons("ZXY")


class TestGetSynSubs(unittest.TestCase):
    """
    Tests get_syn_subs
    """
    maxDiff = None

    def testNoSeqNoORFs(self):
        expected = []
        result = get_syn_subs("", [])
        self.assertEqual(expected, result)

    def testNoSeq(self):
        expected = []
        result = get_syn_subs("", [(0, 8), (2, 10)])
        self.assertEqual(expected, result)

    def testNoORFs(self):
        expected = [{'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},

                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},

                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []},
                    {'A': [], 'C': [], 'G': [], 'T': []}]
        result = get_syn_subs("ATGTTTCCC", [])
        self.assertEqual(expected, result)

    def testStartCodon(self):
        expected = [{'A': [0], 'C': [1], 'G': [1], 'T': [1]},   # A
                    {'A': [1], 'C': [1], 'G': [1], 'T': [0]},   # T
                    {'A': [1], 'C': [1], 'G': [0], 'T': [1]}]   # G
        result = get_syn_subs("ATG", [(0, 2)])
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        expected = [{'A': [0], 'C': [1], 'G': [1], 'T': [1]},   # A
                    {'A': [1], 'C': [1], 'G': [1], 'T': [0]},   # T
                    {'A': [1], 'C': [1], 'G': [0], 'T': [1]},   # G

                    {'A': [1], 'C': [1], 'G': [0], 'T': [1]},   # G
                    {'A': [1], 'C': [1], 'G': [0], 'T': [1]},   # G
                    {'A': [0], 'C': [0], 'G': [0], 'T': [0]},   # G

                    {'A': [1], 'C': [1], 'G': [1], 'T': [0]},   # T
                    {'A': [0], 'C': [1], 'G': [1], 'T': [1]},   # A
                    {'A': [0], 'C': [1], 'G': [0], 'T': [1]}]   # G
        result = get_syn_subs("ATGGGGTAG", [(0, 8)])
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        expected = [
                                                                    # ORF (0, 23)           # ORF (5, 19)
            {'A': [0, 0], 'C': [1, 0], 'G': [1, 0], 'T': [1, 0]},   # A
            {'A': [1, 0], 'C': [1, 0], 'G': [1, 0], 'T': [0, 0]},   # T
            {'A': [1, 0], 'C': [1, 0], 'G': [0, 0], 'T': [1, 0]},   # G

            {'A': [0, 0], 'C': [0, 0], 'G': [1, 0], 'T': [1, 0]},   # A
            {'A': [1, 0], 'C': [1, 0], 'G': [0, 0], 'T': [1, 0]},   # G
            {'A': [0, 0], 'C': [1, 1], 'G': [0, 1], 'T': [1, 1]},   # A                     A   5

            {'A': [1, 1], 'C': [1, 1], 'G': [1, 1], 'T': [0, 0]},   # T                     T
            {'A': [1, 1], 'C': [1, 1], 'G': [0, 0], 'T': [1, 1]},   # G                     G
            {'A': [1, 1], 'C': [1, 1], 'G': [0, 0], 'T': [1, 1]},   # G                     G   8

            {'A': [1, 1], 'C': [0, 0], 'G': [1, 1], 'T': [1, 1]},   # C                     C
            {'A': [0, 0], 'C': [1, 0], 'G': [1, 0], 'T': [1, 0]},   # A                     A
            {'A': [1, 1], 'C': [0, 0], 'G': [1, 1], 'T': [0, 1]},   # C                     C   11

            {'A': [0, 0], 'C': [1, 1], 'G': [1, 1], 'T': [1, 1]},   # A                     A
            {'A': [0, 0], 'C': [1, 1], 'G': [1, 0], 'T': [1, 1]},   # A                     A
            {'A': [0, 1], 'C': [1, 1], 'G': [0, 0], 'T': [1, 1]},   # G                     G   14

            {'A': [1, 1], 'C': [1, 1], 'G': [1, 1], 'T': [0, 0]},   # T                     T
            {'A': [1, 0], 'C': [1, 0], 'G': [0, 0], 'T': [1, 0]},   # G                     G
            {'A': [1, 1], 'C': [0, 1], 'G': [1, 1], 'T': [0, 0]},   # T                     T   17

            {'A': [0, 0], 'C': [1, 1], 'G': [1, 0], 'T': [1, 1]},   # A                     A
            {'A': [0, 0], 'C': [1, 1], 'G': [1, 0], 'T': [1, 1]},   # A                     A
            {'A': [1, 0], 'C': [0, 0], 'G': [1, 0], 'T': [0, 0]},   # C

            {'A': [1, 0], 'C': [1, 0], 'G': [1, 0], 'T': [0, 0]},   # T
            {'A': [0, 0], 'C': [1, 0], 'G': [1, 0], 'T': [1, 0]},   # A
            {'A': [0, 0], 'C': [1, 0], 'G': [0, 0], 'T': [1, 0]}]   # G

        result = get_syn_subs("ATGAGATGGCACAAGTGTAACTAG", [(0, 23), (5, 19)])
        self.assertEqual(expected, result)


class TestGetCodon(unittest.TestCase):
    """
    Tests get_codon
    """
    def testSimpleUse(self):
        expected = ("ATG", 0)
        result = get_codon("ATGAAAGGGTAA", 0, (0, 11))
        self.assertEqual(expected, result)

    def testSimpleUse2(self):
        expected = ("ATG", 1)
        result = get_codon("ATGAAAGGGTAA", 1, (0, 11))
        self.assertEqual(expected, result)

    def testSimpleUse3(self):
        expected = ("GGG", 2)
        result = get_codon("ATGAAAGGGTAA", 8, (0, 11))
        self.assertEqual(expected, result)

    def testSamePositionMultipleOrfs1(self):
        expected = ("CAT", 1)
        result = get_codon("ATGAAAGTGCAACATGGGTAAATAG", 13, (0, 20))
        self.assertEqual(expected, result)

    def testSamePositionMultipleOrfs2(self):
        expected = ("ATG", 0)
        result = get_codon("ATGAAAGTGCAACATGGGTAAATAG", 13, (13, 24))
        self.assertEqual(expected, result)

    def testFindStopCodon(self):
        expected = ("TAA", 2)
        result = get_codon("ATGAAAGGGTAA", 11, (0, 11))
        self.assertEqual(expected, result)


class TestReverseAndComplement(unittest.TestCase):
    """
    Tests reverse_and_complement
    """
    def testEmptyString(self):
        expected = ""
        result = reverse_and_complement("")
        self.assertEqual(expected, result)

    def testSimpleUse(self):
        expected = "KSBRA-HKG-DMCNHYG*VWT"
        result = reverse_and_complement("AWB*CRDNGKH-CMD-TYVSM")
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        seq = "awb*crdngkh-cmd-tyvsm"
        expected = "KSBRA-HKG-DMCNHYG*VWT"
        result = reverse_and_complement(seq)
        self.assertEqual(expected, result)

    def testBadInput(self):
        seq = "AtCg-*N..."
        with self.assertRaises(KeyError):
            reverse_and_complement(seq)


class TestSortOrfs(unittest.TestCase):
    """
    Tests sort_orfs
    """
    def testEmptyString(self):
        return None

