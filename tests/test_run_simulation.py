import unittest
import os

from hexse.run_simulation import *

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
        orf_locations = {'+': [{'coords': [(0, 9)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                         '-': []}  # 0-based exclusive indexing
        expected = {'+': [], '-': []}
        result = valid_orfs(orf_locations, len(s))
        self.assertEqual(expected, result)

    def testNotOrf(self):
        s = "ATGCGCGCATGACGA"
        orfs = {'+': [{'coords': [(1, 1)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': []}
        expected = {'+': [{'coords': [(1, 1)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                    '-': []}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testMultipleNotOrfs(self):
        s = "ATGCGCGCATGACGA"
        orfs = {'+': [{'coords': [(1, 1)], 'omega_values': [[1.42, 0.67, 1.22, 0.74]]},
                      {'coords': [(2, 4)], 'omega_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(0, 1), (2, 3)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        expected = {'+': [{'coords': [(1, 1)], 'omega_values': [[1.42, 0.67, 1.22, 0.74]]},
                          {'coords': [(2, 4)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '-': [{'coords': [(0, 1), (2, 3)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testNotMultipleOfThree(self):
        s = "ATGTCGATGCATGC"
        orfs = {'+': [{'coords': [(1, 11)], 'omega_values': [0.12, 0.54, 0.98, 1.26]}],
                '-': []}
        expected = {'+': [{'coords': [(1, 11)], 'omega_values': [0.12, 0.54, 0.98, 1.26]}],
                    '-': []}
        result = valid_orfs(orfs, len(s))
        self.assertEqual(expected, result)

    def testOutsideRangeOfSequence(self):
        s = "ATGTCGATGCATGC"
        orfs = {'+': [{'coords': [(0, 12)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': [{'coords': [(1, 20)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        expected = {'+': [],
                    '-': [{'coords': [(1, 20)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
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
        expected = {'+0': [{'coords': [(0, 11)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                    '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        orfs = {'+': [{'coords': [(0, 11)], 'omega_values': [1.42, 0.67, 1.22, 0.74]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testMultipleOrfs(self):
        expected = {'+0': [{'coords': [(0, 23)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [{'coords': [(5, 19)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(0, 23)], 'omega_values': [0.56, 0.89, 1.13]},
                      {'coords': [(5, 19)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testSplicedOrf(self):
        expected = {'+0': [{'coords': [(0, 23), (5, 19)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [{'coords': [(8, 21)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(0, 23), (5, 19)], 'omega_values': [0.56, 0.89, 1.13]},
                      {'coords': [(8, 21)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testBacktoBackOrfs(self):
        expected = {'+0': [{'coords': [(1, 9)], 'omega_values': [1.42, 0.67, 1.22, 0.74]},
                           {'coords': [(16, 24)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [],
                    '-0': [],
                    '-1': [],
                    '-2': []}
        orfs = {'+': [{'coords': [(1, 9)], 'omega_values': [1.42, 0.67, 1.22, 0.74]},
                      {'coords': [(16, 24)], 'omega_values': [0.56, 0.89, 1.13]}],
                '-': []}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testFwdReverseOrfs(self):
        expected = {'+0': [{'coords': [(5, 49)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '+1': [],
                    '+2': [],
                    '-0': [],
                    '-1': [],
                    '-2': [{'coords': [(3, 29)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        orfs = {'+': [{'coords': [(5, 49)], 'omega_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(3, 29)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)

    def testAllSixOrfs(self):
        expected = {'+0': [{'coords': [(0, 8)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '+1': [{'coords': [(1, 9)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '+2': [{'coords': [(2, 10)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '-0': [{'coords': [(0, 8)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}],
                    '-1': [{'coords': [(1, 9)], 'omega_values': [0.56, 0.89, 1.13]}],
                    '-2': [{'coords': [(2, 10)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        orfs = {'+': [{'coords': [(0, 8)], 'omega_values': [0.56, 0.89, 1.13]},
                      {'coords': [(1, 9)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]},
                      {'coords': [(2, 10)], 'omega_values': [0.56, 0.89, 1.13]}],
                '-': [{'coords': [(0, 8)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]},
                      {'coords': [(1, 9)], 'omega_values': [0.56, 0.89, 1.13]},
                      {'coords': [(2, 10)], 'omega_values': [0.35, 0.78, 1.03, 1.43, 1.80]}]}
        result = sort_orfs(orfs)
        self.assertEqual(expected, result)


class TestCreateValuesDict(unittest.TestCase):

    def test_create_mu_dict(self):
        alpha = 1.15
        ncat = 6
        expected = {'mu1': 0.051710707633483066,
                    'mu2': 0.15181054803756722,
                    'mu3': 0.26809045653750935,
                    'mu4': 0.4186255904232465,
                    'mu5': 0.6442570794470408,
                    'mu6': 1.2255056178040284}

        result = create_values_dict(alpha, ncat, "mu", dist=ss.gamma)
        self.assertEqual(expected, result)

    def test_create_omega_dict(self):
        alpha = 0.75
        ncat = 3
        expected = {'omega1': 0.038415174742770994,
                    'omega2': 0.18912229978125328,
                    'omega3': 0.6724625254759758}
        result = create_values_dict(alpha, ncat, "omega", dist=ss.gamma)
        self.assertEqual(expected, result)


class TestDiscretize(unittest.TestCase):

    def test_discretize_gamma(self):
        alpha = 1.25
        ncat = 4
        expected = np.array([0.09199853806558903, 0.27043066909631136, 0.5158061369385518, 1.1217646558655263])
        result = discretize(alpha, ncat, dist=ss.gamma)
        self.assertEqual(expected.all(), result.all())

    def test_discertize_gamma_str(self):
        alpha = 1.25
        ncat = 4
        expected = np.array([0.09199853806558903, 0.27043066909631136, 0.5158061369385518, 1.1217646558655263])
        result = discretize(alpha, ncat, dist='ss.gamma')
        self.assertEqual(expected.all(), result.all())

    def test_discretize_lognorm(self):
        alpha = 1.25
        ncat = 4
        expected = np.array([0.18399708, 0.54086134, 1.03161227, 2.24352931])
        result = discretize(alpha, ncat, dist=ss.lognorm)
        self.assertEqual(expected.all(), result.all())

    def test_discretize_lognorm_str(self):
        alpha = 1.25
        ncat = 4
        expected = np.array([0.18399708, 0.54086134, 1.03161227, 2.24352931])
        result = discretize(alpha, ncat, dist='ss.lognorm')
        self.assertEqual(expected.all(), result.all())


class TestReadSequence(unittest.TestCase):

    def test_genbank(self):
        in_seq = os.path.join(CURR_ABSPATH, 'fixtures/NC_003977.2_HBV.gb')
        exp_seq = 'AATTCCACAACCTTCCACCAAACTCTGCAAGATCCCAGAGTGAGAGGCCTGTATTTCCCTGCTGGTGGCT'
        res_seq = read_sequence(in_seq)
        self.assertEqual(exp_seq, res_seq[:70])  # First 70 nucleotides of the HBV genome

    def test_fasta(self):
        in_seq = os.path.join(CURR_ABSPATH, 'fixtures/test_seq.fa')
        exp_seq = 'ATGAATGCCTGACTAA'
        res_seq = read_sequence(in_seq)
        self.assertEqual(exp_seq, res_seq)

    def test_invalid_format(self):
        in_seq = os.path.join('fixtures/test_seq_orfs.csv')
        with self.assertRaises(SystemExit) as e:
            res = read_sequence(in_seq)
        self.assertEqual(e.exception.code, 1)


class TestReadOrfs(unittest.TestCase):

    def test_genbank_format(self):
        in_orfs = os.path.join(CURR_ABSPATH, 'fixtures/NC_003977.2_HBV.gb')
        exp_orfs = {'+': [{'coords': [(2308, 3182), (0, 1625)]},
                          {'coords': [(2849, 3182), (0, 837)]},
                          {'coords': [(3173, 3182), (0, 837)]},
                          {'coords': [(156, 837)]},
                          {'coords': [(1375, 1840)]},
                          {'coords': [(1815, 2454)]},
                          {'coords': [(1853, 1922)]},
                          {'coords': [(1902, 2454)]}],
                    '-': []}
        res_orfs = parse_genbank_orfs(in_orfs)
        self.assertEqual(exp_orfs, res_orfs)

    def test_csv_format(self):
        in_orfs = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV_orfs.csv')
        omega_values = {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}
        exp_orfs = {'+': [{'coords': [(2849, 3182)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(0, 837)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(3173, 3182)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(156, 837)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(1375, 1840)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(1815, 2454)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(1853, 1922)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(1902, 2454)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(8312, 10228)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}},
                          {'coords': [(13340, 13348)],
                           'omega_values': {'omega1': 0.56, 'omega2': 0.89, 'omega3': 1.13}}],
                    '-': []}
        res_orfs = parse_orfs_from_csv(in_orfs, omega_values)
        self.assertEqual(exp_orfs, res_orfs)

    def test_yaml_format(self):
        path = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV.yaml')

        with open(path, 'r') as stream:
            settings = yaml.safe_load(stream)

        exp_orfs = {'+': [{'coords': [[2849, 3182]], 'omega_classes': 3, 'omega_shape': 1.5,
                           'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]},
                          {'coords': [[0, 837]], 'omega_classes': 4, 'omega_shape': 1.7,
                           'omega_values': [0.17314944359111226, 0.42077177796287807,
                                            0.7209957283986081, 1.4050830500389533]},
                          {'coords': [[3173, 3182]], 'omega_classes': 5, 'omega_shape': 1.9,
                           'omega_values': [0.1846759438880899, 0.4079605248732454, 0.6343987842679637,
                                            0.9365407795137214, 1.636423967437797]},
                          {'coords': [[156, 837]], 'omega_classes': 6, 'omega_shape': 1.2,
                           'omega_values': [0.05771979106727954, 0.1642586933267217, 0.2853664305928927,
                                            0.44049638575624284, 0.6712493575878171, 1.2609093415549824]},
                          {'coords': [[1375, 1840]], 'omega_classes': 7, 'omega_shape': 1.3,
                           'omega_values': [0.06190206523811568, 0.16401472284710655, 0.2707180969394326,
                                            0.39578942300043607, 0.5562688094260637, 0.7941290708203261,
                                            1.3971778116544598]},
                          {'coords': [[1815, 2454]], 'omega_classes': 3, 'omega_shape': 1.25,
                           'omega_values': [0.1204513150144901, 0.38254280439302735, 0.9970058805432849]},
                          {'coords': [[1853, 1922]], 'omega_classes': 4, 'omega_shape': 1.87,
                           'omega_values': [0.20738045775300373, 0.4790470356727463,
                                            0.7976433287340154, 1.5079291778341823]},
                          {'coords': [[1902, 2454]], 'omega_classes': 5, 'omega_shape': 2.5,
                           'omega_values': [0.3075971532102495, 0.5999467116400771, 0.8730832667988639,
                                            1.2223993729945162, 1.9969734953254998]}],
                    '-': []}

        res_orfs = parse_orfs_from_yaml(settings)
        self.assertEqual(exp_orfs, res_orfs)


class TestHandleStopCodons(unittest.TestCase):

    def test_internal_stop(self):
        seq = "ATGGCCTAATGA"
        exp_count = 2
        res_count = stop_in_seq(seq, start=0, end=12)
        self.assertEqual(exp_count, res_count)

    def test_stop_out_of_frame(self):
        seq = "ATGAAAGCAATGAGGTGA"
        exp_count = 1
        res_count = stop_in_seq(seq, start=0, end=18)
        self.assertEqual(exp_count, res_count)

    def test_no_internal_stop_codons(self):
        expected = 0
        seq = 'ATGGGAGAACGGGCTAGAGCTAGCA'
        result = count_internal_stop_codons(seq, '+', (0, 18))
        self.assertEqual(expected, result)

    def test_internal_stop_codon(self):
        seq = "ATGTGATAA"
        exp = 1
        res = count_internal_stop_codons(seq, '+', (0, 9))
        self.assertEqual(exp, res)


class TestGetParameters(unittest.TestCase):

    def test_get_pi_from_seq(self):
        expected = {'A': 0.44, 'C': 0.0, 'G': 0.22, 'T': 0.33}
        s = 'ATGTGATAA'
        settings = {}
        pi = [None, None, None, None]
        result = get_pi(pi, settings, s)
        self.assertEqual(expected, result)

    def test_get_pi_from_settings(self):
        expected = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        s = 'ATGTGATAA'
        pi = [0.05, 0.05, 0.05, 0.85]   # input file takes priority

        path = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV.yaml')
        with open(path, 'r') as stream:
            settings = yaml.safe_load(stream)

        result = get_pi(pi, settings, s)
        self.assertEqual(expected, result)

    def test_user_spec_pi(self):
        expected = {'A': 0.05, 'C': 0.85, 'G': 0.05, 'T': 0.05}
        s = 'ATGTGATAA'
        pi = [0.05, 0.05, 0.05, 0.85]
        settings = None
        result = get_pi(pi, settings, s)
        self.assertEqual(expected, result)

    def test_invalid_pi(self):
        with self.assertRaises(SystemExit) as e:
            s = 'ATGCTGCATGGCGCA'
            settings = None
            pi = [0, 1, 2, 3]
            get_pi(pi, settings, s)

        self.assertEqual(e.exception.code, 1)

    def test_get_kappa_settings(self):
        expected = 0.3  # Input file takes priority

        path = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV.yaml')
        with open(path, 'r') as stream:
            settings = yaml.safe_load(stream)

        result = get_kappa(0.03, settings)
        self.assertEqual(expected, result)

    def test_user_spec_kappa(self):
        expected = 0.03
        settings = None
        result = get_kappa(0.03, settings)
        self.assertEqual(expected, result)

    def test_get_mu_settings(self):
        expected = 0.0005  # Input file takes priority

        path = os.path.join(CURR_ABSPATH, 'fixtures/test_HBV.yaml')
        with open(path, 'r') as stream:
            settings = yaml.safe_load(stream)

        result = get_global_rate(0.01, settings)
        self.assertEqual(expected, result)

    def test_user_spec_mu(self):
        expected = 0.01
        settings = None
        result = get_global_rate(0.01, settings)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
