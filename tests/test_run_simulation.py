import unittest
import os

from src.run_simulation import *

CURR_ABSPATH = os.path.dirname(os.path.abspath(__file__))


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
        orf_locations = {'+': [{'coords': [(0, 9)],
                                'dN_values': [1.42, 0.67, 1.22, 0.74],
                                'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                         '-': []}  # 0-based exclusive indexing
        expected = {'+': [], '-': []}
        result = valid_orfs(orf_locations, len(s))
        self.assertEqual(expected, result)

    def testNotOrf(self):
        s = "ATGCGCGCATGACGA"
        orfs = {'+': [{'coords': [(1, 1)],
                       'dN_values': [1.42, 0.67, 1.22, 0.74],
                       'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': []}
        expected = {'+': [{'coords': [(1, 1)],
                           'dN_values': [1.42, 0.67, 1.22, 0.74],
                           'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                    '-': []}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testMultipleNotOrfs(self):
        s = "ATGCGCGCATGACGA"
        orfs = {'+': [{'coords': [(1, 1)],
                       'dN_values': [[1.42, 0.67, 1.22, 0.74]],
                       'dS_values': [[1.42, 0.67, 1.22, 0.74]]},
                      {'coords': [(2, 4)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(0, 1), (2, 3)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        expected = {'+': [{'coords': [(1, 1)],
                           'dN_values': [[1.42, 0.67, 1.22, 0.74]],
                           'dS_values': [[1.42, 0.67, 1.22, 0.74]]},
                          {'coords': [(2, 4)],
                           'dN_values': [0.56, 0.89, 1.13],
                           'dS_values': [0.56, 0.89, 1.13]}],
                    '-': [{'coords': [(0, 1), (2, 3)],
                           'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                           'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testNotMultipleOfThree(self):
        s = "ATGTCGATGCATGC"
        orfs = {'+': [{'coords': [(1, 11)],
                       'dN_values': [0.12, 0.54, 0.98, 1.26],
                       'dS_values': [0.12, 0.54, 0.98, 1.26]}],
                '-': []}
        expected = {'+': [{'coords': [(1, 11)],
                           'dN_values': [0.12, 0.54, 0.98, 1.26],
                           'dS_values': [0.12, 0.54, 0.98, 1.26]}],
                    '-': []}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testOutsideRangeOfSequence(self):
        s = "ATGTCGATGCATGC"
        orfs = {'+': [{'coords': [(0, 12)],
                       'dN_values': [1.42, 0.67, 1.22, 0.74],
                       'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': [{'coords': [(1, 20)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        expected = {'+': [],
                    '-': [{'coords': [(1, 20)],
                           'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                           'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)


class TestGetOpenReadingFrames(unittest.TestCase):
    """
    Tests get_reading_frames
    """

    def testOneORF(self):
        # Tests one ORF (ATG AAA TAG)
        s = "ATGAAATAG"
        expected = [(0, 9)]
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
        expected = [(5, 50), (30, 3)]
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
        expected = [(0, 15)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testMultipleORFs(self):
        # Tests scenario when there are 2 different ORFs
        # ORF1: ATG AAA GTG CAA CAT GGG TAA
        # ORF2: ATG GGT AAA TA
        s = "ATGAAAGTGCAACATGGGTAAATAG"
        expected = [(0, 21), (13, 25)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testLowerCaseInput(self):
        s = "atgaaatag"
        expected = [(0, 9)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testStopBeforeStart(self):
        # Tests scenario when a stop codon precedes a start codon
        s = "TAGATGAAATAG"
        expected = [(3, 12)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs1(self):
        # Tests a scenario when there are 2 consecutive ORFs in the same reading frame
        # ORF1: ATG TTT TAA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTAGATGCACTAA"
        expected = [(0, 9), (9, 18)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs2(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 6 nucleotides
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTGACCCAAAATGCACTAA"
        expected = [(0, 9), (15, 24)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs3(self):
        # Tests a scenario when there are 2 consecutive ORFs, separated by 5 nucleotide
        # ORF1: ATG TTT TGA
        # ORF2: ATG CAC TAA
        s = "ATGTTTTGACCCAAATGCACTAA"
        expected = [(0, 9), (14, 23)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testBackToBackORFs4(self):
        # Tests a scenario when there are 2 ORFs in the same reading frame (+1), separated by 6 nucleotides
        # ORF1: ATG GAG TGA
        # ORF2: ATG GAG TGA
        s = "AATGGAGTGACCCGGGATGGAGTAG"
        expected = [(1, 10), (16, 25)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)

    def testInternalORF(self):
        # Tests a scenario when there is an ORF encoded in a different reading frame inside of another ORF
        # ORF1: ATG AGA TGG CAC AAG TGT AAC TAG
        # ORF2: ATG GCA CAA GTG TAA
        s = "ATGAGATGGCACAAGTGTAACTAG"
        expected = [(0, 24), (5, 20)]
        result = get_open_reading_frames(s)
        self.assertEqual(expected, result)


class TestSortOrfs(unittest.TestCase):
    """
    Tests sort_orfs
    """

    def testSimpleUse(self):
        expected = {'+0': [{'coords': [(0, 11)],
                            'dN_values': [1.42, 0.67, 1.22, 0.74],
                            'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                    '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        orfs = {'+': [{'coords': [(0, 11)],
                       'dN_values': [1.42, 0.67, 1.22, 0.74],
                       'dS_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testMultipleOrfs(self):
        expected = {'+0': [{'coords': [(0, 23)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [{'coords': [(5, 19)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(0, 23)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]},
                      {'coords': [(5, 19)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testSplicedOrf(self):
        expected = {'+0': [{'coords': [(0, 23), (5, 19)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [{'coords': [(8, 21)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(0, 23), (5, 19)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]},
                      {'coords': [(8, 21)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testBacktoBackOrfs(self):
        expected = {'+0': [{'coords': [(1, 9)],
                            'dN_values': [1.42, 0.67, 1.22, 0.74],
                            'dS_values': [1.42, 0.67, 1.22, 0.74]},
                           {'coords': [(16, 24)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(1, 9)],
                       'dN_values': [1.42, 0.67, 1.22, 0.74],
                       'dS_values': [1.42, 0.67, 1.22, 0.74]},
                      {'coords': [(16, 24)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testFwdReverseOrfs(self):
        expected = {'+0': [{'coords': [(5, 49)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [],
                    '-0': [],
                    '-1': [],
                    '-2': [{'coords': [(3, 29)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        orfs = {'+': [{'coords': [(5, 49)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(3, 29)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testAllSixOrfs(self):
        expected = {'+0': [{'coords': [(0, 8)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '+1': [{'coords': [(1, 9)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '+2': [{'coords': [(2, 10)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '-0': [{'coords': [(0, 8)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-1': [{'coords': [(1, 9)],
                            'dN_values': [0.56, 0.89, 1.13],
                            'dS_values': [0.56, 0.89, 1.13]}],
                    '-2': [{'coords': [(2, 10)],
                            'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                            'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        orfs = {'+': [{'coords': [(0, 8)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]},
                      {'coords': [(1, 9)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]},
                      {'coords': [(2, 10)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(0, 8)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]},
                      {'coords': [(1, 9)],
                       'dN_values': [0.56, 0.89, 1.13],
                       'dS_values': [0.56, 0.89, 1.13]},
                      {'coords': [(2, 10)],
                       'dN_values': [0.35, 0.78, 1.03, 1.43, 1.80],
                       'dS_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)


class TestDiscretization(unittest.TestCase):

    def test_get_rate_values(self):
        alpha = 1.25
        ncat = 4
        expected = [0.18399707613117805, 0.5408613381926227, 1.0316122738771036, 2.2435293117310526]
        result = get_rate_values(alpha, ncat)
        self.assertEqual(expected, result)

    def test_discretize_gamma(self):
        alpha = 1.25
        ncat = 4
        expected = np.array([0.18399708, 0.54086134, 1.03161227, 2.24352931])
        result = discretize_gamma(alpha, ncat)
        self.assertEqual(expected.all(), result.all())


class TestParseInput(unittest.TestCase):

    def testParseGenbank(self):
        in_seq = os.path.join(CURR_ABSPATH, 'fixtures/NC_003977.2_HBV.gb')
        exp_seq = 'AATTCCACAACCTTCCACCAAACTCTGCAAGATCCCAGAGTGAGAGGCCTGTATTTCCCTGCTGGTGGCT'
        exp_orfs = {1: [[(2308, 3182), (0, 1625)], [(2849, 3182), (0, 837)], [(3173, 3182), (0, 837)],
                        [(156, 837)], [(1375, 1840)], [(1815, 2454)], [(1853, 1922)], [(1902, 2454)]],
                    -1: []}

        res_seq, res_orfs = parse_genbank(in_seq)
        self.assertEqual(exp_seq, res_seq[:70])  # First 70 nucleotides of the HBV genome
        self.assertEqual(exp_orfs, res_orfs)

    def testParseFasta(self):
        in_seq = os.path.join(CURR_ABSPATH, 'fixtures/HBV.fasta')
        exp_seq = 'CATTCGGGCTGGGTTTCACCCCACCGCACGGAGGCCTTTTGGGGTGGAGCCCTCAGGCTCAGGGCATACTACAAACTTTGCCAGCAAATCCGCC' \
                  'TCCTGCCTCCACCAATCGCCAGTCAGGAAGGCAGCCTACCCCGCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGG'
        res_seq = parse_fasta(in_seq)
        self.assertEqual(exp_seq, res_seq[3000:])  # Last 182 nucleotides of the HBV genome

    def testCheckOrfs(self):
        #  If user specified ORFs
        in_orfs = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV_orfs.csv')
        exp_orfs = [(2849, 3182), (0, 837), (3173, 3182), (156, 837),
                    (1375, 1840), (1815, 2454), (1853, 1922), (1902, 2454)]
        self.assertEqual(exp_orfs, check_orfs(in_orfs))

        # If user did not specify ORFs
        in_orfs = None
        s = 'ATGAAAGTGCAACATGGGTAAATAG'
        exp_orfs = [(0, 21), (13, 25)]
        self.assertEqual(exp_orfs, check_orfs(in_orfs, s))


if __name__ == '__main__':
    unittest.main()
