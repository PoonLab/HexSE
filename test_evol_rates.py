"""
Unit tests for evol_rates.py
"""
import unittest

from evol_rates import reverse_and_complement
from evol_rates import codon_iterator
from evol_rates import get_reading_frames
from evol_rates import get_syn_codons


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

    def testStartBadInputStop(self):
        # Tests scenario when a start and stop codon are separated by some invalid characters
        expected = []
        result = get_reading_frames("ATG#A>GGGTAA")
        self.assertEqual(expected, result)

    def testBadInput(self):
        expected = []
        result = get_reading_frames("78979GTC")
        self.assertEqual(expected, result)


class TestGetSynCodons(unittest.TestCase):

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

    def testBadInput(self):
        with self.assertRaises(KeyError):
            get_syn_codons(">Seq")
