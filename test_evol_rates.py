import unittest
from src.evol_rates import Rates
from src.sequence_info import Sequence
from src.sequence_info import Nucleotide


class TestDrawOmegaValues(unittest.TestCase):
    pass


class TestGetFrequencyRates(unittest.TestCase):
    """
    Tests get_frequency_rates
    """
    def testSimpleUse(self):
        s = Sequence('AAAAAAAAA', [])
        r = Rates(s, 0.5, None, None, None)
        expected = {'A': 1, 'C': 0, 'T': 0, 'G': 0}
        result = r.get_frequency_rates()
        self.assertEqual(expected, result)

    def testSimpleUse1(self):
        expected = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
        result = Rates.get_frequency_rates("GTACGATCGATCGATGCTAGC")
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
        result = Rates.get_frequency_rates(seq.upper(), seq)
        self.assertEqual(expected, result)


class TestGetSynSubs(unittest.TestCase):
    """
    Tests get_syn_subs
    """
    maxDiff = None

    def testNoSeqNoORFs(self):
        expected = []
        result = Rates.get_syn_subs("", [])
        self.assertEqual(expected, result)

    def testNoSeq(self):
        expected = []
        result = Rates.get_syn_subs("", [(0, 8), (2, 10)])
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
        result = Rates.get_syn_subs("ATGTTTCCC", [])
        self.assertEqual(expected, result)

    def testStartCodon(self):
        expected = [{'A': [0], 'C': [1], 'G': [1], 'T': [1]},   # A
                    {'A': [1], 'C': [1], 'G': [1], 'T': [0]},   # T
                    {'A': [1], 'C': [1], 'G': [0], 'T': [1]}]   # G
        result = Rates.get_syn_subs("ATG", [(0, 2)])
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
        result = Rates.get_syn_subs("ATGGGGTAG", [(0, 8)])
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

        result = Rates.get_syn_subs("ATGAGATGGCACAAGTGTAACTAG", [(0, 23), (5, 19)])
        self.assertEqual(expected, result)
