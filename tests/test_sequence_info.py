import random
import unittest

from src.sequence_info import Sequence

MAX_DIFF = None
KAPPA = 0.3
GLOBAL_RATE = 1

OMEGA_VALUES_4 = {'omega1':  0.29327471612351436,
                  'omega2': 0.6550136761581515,
                  'omega3': 1.0699896623909886,
                  'omega4': 1.9817219453273531}

OMEGA_VALUES_5 = {'omega1': 0.25496553479835254,
                  'omega2': 0.5492012532921686,
                  'omega3': 0.842647578837806,
                  'omega4': 1.2308533419685328,
                  'omega5': 2.1223322911031457}

MU_VALUES_4 = {'mu1': 0.325970770628452,
               'mu2': 0.7739880200789496,
               'mu3': 1.481877317174674,
               'mu4': 4.351175963585551}

MU_VALUES_3 = {'mu1': 0.3965034306394888,
               'mu2': 1.0832728336432287,
               'mu3': 3.7199827893060657}

EMPTY_EVENT_TREE = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'T': {'class': {}},
                                                'C': {'class': {}},
                                                'G': {'class': {}}}},
                              'T': {'from_nt': {'A': {'class': {}},
                                                'T': None,
                                                'C': {'class': {}},
                                                'G': {'class': {}}}},
                              'C': {'from_nt': {'A': {'class': {}},
                                                'T': {'class': {}},
                                                'C': None,
                                                'G': {'class': {}}}},
                              'G': {'from_nt': {'A': {'class': {}},
                                                'T': {'class': {}},
                                                'C': {'class': {}},
                                                'G': None}}}}


class TestSequence1(unittest.TestCase):
    """
    Sequence: GTACGATCGATCGATGCTAGC
    ORFs:
        (0, 21) (+)
    Notes:
        No START or STOP codons
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9001)     # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, KAPPA, GLOBAL_RATE, pi1, OMEGA_VALUES_4, MU_VALUES_4)

    def testReverseComplement(self):
        s = str(self.sequence1)
        expected = 'GCTAGCATCGATCGATCGTAC'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'CATGCTAGCTAGCTACGATCG'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence1 = self.sequence1.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence1, new_sequence1)
        self.assertEqual(self.sequence1.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence1.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence1.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
        self.assertEqual(self.sequence1.omega_values, new_sequence1.omega_values)
        self.assertEqual(self.sequence1.cat_values, new_sequence1.cat_values)
        self.assertEqual(self.sequence1.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence1.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence1.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq1_nt in enumerate(self.sequence1.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(seq1_nt, new_nt)

            self.assertEqual(seq1_nt.state, new_nt.state)
            self.assertEqual(seq1_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq1_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq1_nt.rates, new_nt.rates)
            self.assertEqual(seq1_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(seq1_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq1_nt.categories_keys, new_nt.categories_keys)
            self.assertEqual(str(seq1_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq1_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq1_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.24, 'C': 0.24, 'G': 0.29, 'T': 0.24}
        result = Sequence.get_frequency_rates(str(self.sequence1))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence1.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)
        exp_sub_rates = {'A': 4.252483383790958e-05, 'C': 4.654455031400801e-05,
                         'G': None,                  'T': 4.654455031400801e-05}
        exp_omegas = {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 1, 0)}
        exp_total_rate = 0.0001356139344659256
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omegas, nt.omega_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence1.nt_sequence[2]  # A
        result = self.sequence1.is_transv(nt.state, 'A')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence1.is_transv(nt.state, 'G')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence1.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

    def testCodonIterator(self):
        orf = self.sequence1.nt_sequence[0: 21]
        result = []
        expected = [['G', 'T', 'A'], ['C', 'G', 'A'], ['T', 'C', 'G'],
                    ['A', 'T', 'C'], ['G', 'A', 'T'], ['G', 'C', 'T'], ['A', 'G', 'C']]

        for codon in self.sequence1.codon_iterator(orf, 0, 21):
            result.append(codon)

        for pos, cdn in enumerate(result):
            for idx, nt in enumerate(cdn):
                self.assertEqual(expected[pos][idx], nt.state)

    def testFindCodons(self):
        # Test forward strand ORF
        expected = ['GTA', 'CGA', 'TCG', 'ATC', 'GAT', 'GCT', 'AGC']
        result = self.sequence1.find_codons('+0', [(0, 21)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStopStartCodon(self):
        nt = self.sequence1.nt_sequence[0]  # G in codon GTA
        expected = False
        result = self.sequence1.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


class TestSequence2(unittest.TestCase):
    """
    Sequence: ATGAATAAACCCGTATGA
    ORFs (indexing relative to forward strand)
        (0, 18) (+) ATG AAT AAA CCC GTA TGA
        (3, 15) (-)
    Notes:
        2 overlapping ORFs
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4001)
        s2 = 'ATGAATAAACCCGTATGA'
        sorted_orfs = {'+0': [[(0, 18)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [[(3, 15)]]}
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, sorted_orfs, KAPPA, GLOBAL_RATE, pi2, OMEGA_VALUES_5, MU_VALUES_3)

    def testReverseComplement(self):
        s = str(self.sequence2)
        expected = 'TCATACGGGTTTATTCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTTATTTGGGCATACT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence2 = self.sequence2.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence2, new_sequence2)
        self.assertEqual(self.sequence2.orfs, new_sequence2.orfs)
        self.assertEqual(self.sequence2.kappa, new_sequence2.kappa)
        self.assertEqual(self.sequence2.global_rate, new_sequence2.global_rate)
        self.assertEqual(self.sequence2.pi, new_sequence2.pi)
        self.assertEqual(self.sequence2.omega_values, new_sequence2.omega_values)
        self.assertEqual(self.sequence2.cat_values, new_sequence2.cat_values)
        self.assertEqual(self.sequence2.is_circular, new_sequence2.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence2.event_tree, new_sequence2.event_tree)
        self.assertCountEqual(self.sequence2.event_tree, new_sequence2.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq2_nt in enumerate(self.sequence2.nt_sequence):
            new_nt = new_sequence2.nt_sequence[pos]
            self.assertIsNot(seq2_nt, new_nt)

            self.assertEqual(seq2_nt.state, new_nt.state)
            self.assertEqual(seq2_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq2_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq2_nt.rates, new_nt.rates)
            self.assertEqual(seq2_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(seq2_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq2_nt.categories_keys, new_nt.categories_keys)
            self.assertEqual(str(seq2_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq2_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq2_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence2.create_event_tree()
        self.assertEqual(expected, result)

    def testGetFrequencyRates(self):
        expected = {'A': 0.44, 'C': 0.17, 'G': 0.17, 'T': 0.22}
        result = Sequence.get_frequency_rates(str(self.sequence2))
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(4001)   # Set seed value to initialize pseudo-random number generator

        # Synonymous mutation in (+) strand, non-synonymous mutation in (-) frame
        nt = self.sequence2.nt_sequence[11]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.2335164156624226, 'C': None, 'G': 0.030341674692605226, 'T': 0.10113891564201742}
        exp_omega_keys = {'A': ['omega4'], 'C': None, 'G': ['omega2'], 'T': ['omega2']}
        exp_cat_keys = {'A': 'mu3', 'G': 'mu2', 'T': 'mu2'}
        exp_total_rate = 0.36499700599704527
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Mutation would destroy a START codon
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omega_keys = {'A': None, 'C': [], 'G': [], 'T': []}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu1'}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Non-synonymous mutation in both (+) and (-) strands
        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.09666062791987462, 'G': 0.32220209306624875, 'T': 0.13672268277658034}
        # FIXME (maybe??): same omega values used, but different keys
        exp_omega_keys = {'A': None, 'C': ['omega4', 'omega2'], 'G': ['omega2', 'omega4'], 'T': ['omega4', 'omega5']}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.5555854037627037
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence2.nt_sequence[9]  # C
        result = self.sequence2.is_transv(nt.state, 'C')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence2.is_transv(nt.state, 'T')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence2.is_transv(nt.state, 'G')
        expected = True
        self.assertEqual(expected, result)

    def testFindCodons(self):

        # Test forward strand ORF
        expected = ['ATG', 'AAT', 'AAA', 'CCC', 'GTA', 'TGA']
        result = self.sequence2.find_codons('+0', [(0, 18)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Test reverse strand ORF
        expected = ['ATG', 'CCC', 'AAA', 'TAA']
        result = self.sequence2.find_codons('-2', [(3, 15)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-2')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence2.nt_sequence[5]  # T in codon AAT (+ strand) and TAA in (- strand)
        expected = True
        result = self.sequence2.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence2.nt_sequence[14]  # A in codon GTA (+ strand) and ATG (- strand)
        expected = True
        result = self.sequence2.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


class TestSequence3(unittest.TestCase):
    """
    Sequence: ATGACGTGGTGA
    ORFs:
        (0, 12) (+)
    Notes:
        Simple ORF
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(555)

        s3 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, OMEGA_VALUES_4, MU_VALUES_3)

    def testReverseComplement(self):
        s = str(self.sequence3)
        expected = 'TCACCACGTCAT'   # Reverse and complement
        result = Sequence.complement(s, rev=True)
        self.assertEqual(expected, result)

        expected = 'TACTGCACCACT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence3 = self.sequence3.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence3, new_sequence3)
        self.assertEqual(self.sequence3.orfs, new_sequence3.orfs)
        self.assertEqual(self.sequence3.kappa, new_sequence3.kappa)
        self.assertEqual(self.sequence3.global_rate, new_sequence3.global_rate)
        self.assertEqual(self.sequence3.pi, new_sequence3.pi)
        self.assertEqual(self.sequence3.omega_values, new_sequence3.omega_values)
        self.assertEqual(self.sequence3.cat_values, new_sequence3.cat_values)
        self.assertEqual(self.sequence3.is_circular, new_sequence3.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence3.event_tree, new_sequence3.event_tree)
        self.assertCountEqual(self.sequence3.event_tree, new_sequence3.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq3_nt in enumerate(self.sequence3.nt_sequence):
            new_nt = new_sequence3.nt_sequence[pos]
            self.assertIsNot(seq3_nt, new_nt)

            self.assertEqual(seq3_nt.state, new_nt.state)
            self.assertEqual(seq3_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq3_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq3_nt.rates, new_nt.rates)
            self.assertEqual(seq3_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(seq3_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq3_nt.categories_keys, new_nt.categories_keys)
            self.assertEqual(str(seq3_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq3_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq3_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.25, 'C': 0.08, 'G': 0.42, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence3))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence3.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):

        # Tests a non-synonymous mutation
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.053216889078518174, 'C': 0.06492879242812255, 'G': 0.008721332329709673, 'T': None}
        exp_omega_keys = {'A': ['omega2'], 'C': ['omega2'], 'G': ['omega1'], 'T': None}
        exp_cat_keys = {'A': 'mu2', 'C': 'mu1', 'G': 'mu1'}
        exp_total_rate = 0.1268670138363504

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation in TGG at position 2
        random.seed(555)
        nt = self.sequence3.nt_sequence[8]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.053456076057738736, 'G': None, 'T': 0.307016589860621}
        exp_omega_keys = {'A': [], 'C': ['omega3'], 'G': None, 'T': ['omega2']}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'T': 'mu3'}
        exp_total_rate = 0.36047266591835975

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.1665314408685853, 'C': 0.13649237703904682, 'G': None, 'T': 0.04995943226057559}
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.3529832501682077

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a mutation that would destroy a STOP codon
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.1665314408685853, 'C': 0.13649237703904682, 'G': None, 'T': 0.04995943226057559}
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu1', 'C': 'mu2', 'T': 'mu1'}
        exp_total_rate = 0.3529832501682077

        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence3.nt_sequence[1]  # T
        result = self.sequence3.is_transv(nt.state, 'T')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence3.is_transv(nt.state, 'C')
        expected = False
        self.assertEqual(expected, result)

        result = self.sequence3.is_transv(nt.state, 'G')
        expected = True
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ACG', 'TGG', 'TGA']
        result = self.sequence3.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence3.nt_sequence[0]      # A in START codon
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[11]      # A in STOP codon (TGA)
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[4]
        expected = False
        result = self.sequence3.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[8]  # Last G in TGG codon
        expected = True     # Introduces a STOP codon
        result = self.sequence3.is_start_stop_codon(nt, 'A')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


class TestSequence4(unittest.TestCase):
    """
    Sequence: ATGATGCCCTAA
    ORFs:
        (0, 12) (+)
    Notes:
        Internal methionine
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(5001)

        s4 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(4000)
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, OMEGA_VALUES_5, MU_VALUES_4)

    def testReverseComplement(self):
        s = str(self.sequence4)
        expected = 'TTAGGGCATCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTACGGGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence4 = self.sequence4.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence4, new_sequence4)
        self.assertEqual(self.sequence4.orfs, new_sequence4.orfs)
        self.assertEqual(self.sequence4.kappa, new_sequence4.kappa)
        self.assertEqual(self.sequence4.global_rate, new_sequence4.global_rate)
        self.assertEqual(self.sequence4.pi, new_sequence4.pi)
        self.assertEqual(self.sequence4.omega_values, new_sequence4.omega_values)
        self.assertEqual(self.sequence4.cat_values, new_sequence4.cat_values)
        self.assertEqual(self.sequence4.is_circular, new_sequence4.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence4.event_tree, new_sequence4.event_tree)
        self.assertCountEqual(self.sequence4.event_tree, new_sequence4.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seq4_nt in enumerate(self.sequence4.nt_sequence):
            new_nt = new_sequence4.nt_sequence[pos]
            self.assertIsNot(seq4_nt, new_nt)

            self.assertEqual(seq4_nt.state, new_nt.state)
            self.assertEqual(seq4_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seq4_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seq4_nt.rates, new_nt.rates)
            self.assertEqual(seq4_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(seq4_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seq4_nt.categories_keys, new_nt.categories_keys)

            self.assertEqual(str(seq4_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seq4_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seq4_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.33, 'C': 0.25, 'G': 0.17, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence4))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.3113585720876178, 'G': 0.0, 'T': 0.06849001095359351}
        exp_omega_keys = {'A': None, 'C': ['omega5'], 'G': [], 'T': ['omega5']}
        exp_cat_keys = {'C': 'mu3', 'G': 'mu4', 'T': 'mu1'}
        exp_total_rate = 0.37984858304121133
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutation in internal methionine
        # ATG (TGT, GTG)
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.019536686677227796, 'G': 0.5420777234203797, 'T': 0.9142294839271567}
        exp_omega_keys = {'A': None, 'C': ['omega1'], 'G': ['omega5'], 'T': ['omega5']}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu2', 'T': 'mu4'}
        exp_total_rate = 1.4758438940247642
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence4.nt_sequence[2]  # G
        result = self.sequence4.is_transv(nt.state, 'G')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence4.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence4.is_transv(nt.state, 'A')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ATG', 'CCC', 'TAA']
        result = self.sequence4.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence4.nt_sequence[3]  # A in internal methione codon
        expected = False
        result = self.sequence4.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence4.nt_sequence[0]  # A in start codon
        expected = True
        result = self.sequence4.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


class TestSequence5(unittest.TestCase):
    """
    Sequence: ATGAATGCCTGACTAA
    ORFs:
        (0, 12) (+)
        (4, 16) (+)
    Notes:
        2 overlapping ORFs in the forward strand
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9991)

        s5 = 'ATGAATGCCTGACTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [[(4, 16)]], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, OMEGA_VALUES_4, MU_VALUES_4)

    def testReverseComplement(self):
        s = str(self.sequence5)
        expected = 'TTAGTCAGGCATTCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTTACGGACTGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence5 = self.sequence5.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence5, new_sequence5)
        self.assertEqual(self.sequence5.orfs, new_sequence5.orfs)
        self.assertEqual(self.sequence5.kappa, new_sequence5.kappa)
        self.assertEqual(self.sequence5.global_rate, new_sequence5.global_rate)
        self.assertEqual(self.sequence5.pi, new_sequence5.pi)
        self.assertEqual(self.sequence5.omega_values, new_sequence5.omega_values)
        self.assertEqual(self.sequence5.cat_values, new_sequence5.cat_values)
        self.assertEqual(self.sequence5.is_circular, new_sequence5.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence5.event_tree, new_sequence5.event_tree)
        self.assertCountEqual(self.sequence5.event_tree, new_sequence5.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, seeq5_nt in enumerate(self.sequence5.nt_sequence):
            new_nt = new_sequence5.nt_sequence[pos]
            self.assertIsNot(seeq5_nt, new_nt)

            self.assertEqual(seeq5_nt.state, new_nt.state)
            self.assertEqual(seeq5_nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(seeq5_nt.complement_state, new_nt.complement_state)
            self.assertEqual(seeq5_nt.rates, new_nt.rates)
            self.assertEqual(seeq5_nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(seeq5_nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(seeq5_nt.categories_keys, new_nt.categories_keys)
            self.assertEqual(str(seeq5_nt.codons), str(new_nt.codons))
            self.assertEqual(len(seeq5_nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(seeq5_nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.38, 'C': 0.19, 'G': 0.19, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence5))
        self.assertEqual(expected, result)

    def testEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_omega_keys = {'A': None, 'C': [], 'G': [], 'T': []}
        exp_cat_keys = {'C': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        self.sequence5.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_cat_keys = {'A': ['omega1'], 'C': ['omega1'], 'G': ['omega4'], 'T': None}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence5.nt_sequence[4]  # A
        result = self.sequence5.is_transv(nt.state, 'A')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence5.is_transv(nt.state, 'C')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence5.is_transv(nt.state, 'G')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        # Check codons in +0 frame
        expected = ['ATG', 'AAT', 'GCC', 'TGA']
        result = self.sequence5.find_codons('+0', [(0, 12)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Check Codons in +1 frame
        expected = ['ATG', 'CCT', 'GAC', 'TAA']
        result = self.sequence5.find_codons('+1', [(4, 16)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+1')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence5.nt_sequence[4]  # A in START codon in ORF (4, 16)
        expected = True
        result = self.sequence5.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

        nt = self.sequence5.nt_sequence[15]     # last A in STOP codon in ORF (4, 16)
        expected = True
        result = self.sequence5.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence5.nt_sequence[3]  # First codon in AAT in ORF (0, 12)
        expected = False
        result = self.sequence5.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


class TestSequence6(unittest.TestCase):
    """
    Sequence: ATGATGGCCCTAA
    ORFs:
        [(0, 5), (6, 13)] (+)
    Notes:
        Spliced ORF (5th nucleotide is excluded)
    """

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4000)

        s6 = 'ATGATGGCCCTAA'
        sorted_orfs = {'+0': [[(0, 5), (6, 13)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, OMEGA_VALUES_5, MU_VALUES_3)

    def testReverseComplement(self):
        s = str(self.sequence6)
        expected = 'TTAGGGCCATCAT'
        result = Sequence.complement(s, rev=True)   # Reverse and complement
        self.assertEqual(expected, result)

        expected = 'TACTACCGGGATT'
        result = Sequence.complement(s, rev=False)  # Complement only
        self.assertEqual(expected, result)

    def testDeepcopy(self):
        # Check that Sequences are different objects with the same attributes
        new_sequence6 = self.sequence6.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence6, new_sequence6)
        self.assertEqual(self.sequence6.orfs, new_sequence6.orfs)
        self.assertEqual(self.sequence6.kappa, new_sequence6.kappa)
        self.assertEqual(self.sequence6.global_rate, new_sequence6.global_rate)
        self.assertEqual(self.sequence6.pi, new_sequence6.pi)
        self.assertEqual(self.sequence6.omega_values, new_sequence6.omega_values)
        self.assertEqual(self.sequence6.cat_values, new_sequence6.cat_values)
        self.assertEqual(self.sequence6.is_circular, new_sequence6.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence6.event_tree, new_sequence6.event_tree)
        self.assertCountEqual(self.sequence6.event_tree, new_sequence6.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence6.nt_sequence):
            new_nt = new_sequence6.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)
            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)
            self.assertEqual(nt.omega_keys, new_nt.omega_keys)
            self.assertEqual(nt.categories_keys, new_nt.categories_keys)
            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.31, 'C': 0.23, 'G': 0.23, 'T': 0.23}
        result = Sequence.get_frequency_rates(str(self.sequence6))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = EMPTY_EVENT_TREE
        result = self.sequence6.create_event_tree()
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': 0.2491527517379426, 'C': 0.07474582552138279, 'G': None, 'T': 0.07474582552138279}
        # No omega keys because nt treated as synonymous mutation
        exp_omega_keys = {'A': [], 'C': [], 'G': None, 'T': []}
        exp_cat_keys = {'A': 'mu2', 'C': 'mu2', 'T': 'mu2'}
        exp_total_rate = 0.3986444027807081
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests internal methionine
        random.seed(4000)
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.12400154884247461, 'G': 1.4194135069676665, 'T': 0.1900007865404335}
        exp_omega_keys = {'A': None, 'C': ['omega4'], 'G': ['omega4'], 'T': ['omega2']}
        exp_cat_keys = {'C': 'mu2', 'G': 'mu3', 'T': 'mu3'}
        exp_total_rate = 1.7334158423505746
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        self.assertEqual(exp_cat_keys, nt.categories_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

    def testIsTransv(self):
        nt = self.sequence6.nt_sequence[6]  # G
        result = self.sequence6.is_transv(nt.state, 'G')
        expected = None
        self.assertEqual(expected, result)

        result = self.sequence6.is_transv(nt.state, 'T')
        expected = True
        self.assertEqual(expected, result)

        result = self.sequence6.is_transv(nt.state, 'A')
        expected = False
        self.assertEqual(expected, result)

    def testFindCodons(self):
        expected = ['ATG', 'ATG', 'CCC', 'TAA']
        result = self.sequence6.find_codons('+0', [(0, 5), (6, 13)])
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStopStartCodon(self):
        nt = self.sequence6.nt_sequence[5]
        expected = False
        result = self.sequence6.is_start_stop_codon(nt, 'G')
        self.assertEqual(expected, result)

        nt = self.sequence6.nt_sequence[2]
        expected = True
        result = self.sequence6.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        pass

    def testMutationRate(self):
        pass


# ==========================================
# Tests for Codon
# ==========================================
class TestCodon(unittest.TestCase):

    def setUp(self):
        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [[(0, 21)]]}, kappa, mu, pi1, dN_values, dS_values)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [[(0, 12)]], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(9001)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values)

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.nt_seq5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values)

        s6 = 'ATGAATGCCTGACTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [[(4, 16)]], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        random.seed(4000)
        self.nt_seq6 = Sequence(s6, sorted_orfs, kappa, mu, pi6, dN_values, dS_values)

    def testNtInPos(self):
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]
        nt = self.nt_seq1.nt_sequence[0]
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[1]
        nt = self.nt_seq1.nt_sequence[3]
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence[4]
        expected = 1
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        nt = self.nt_seq1.nt_sequence[5]
        expected = 2
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        # Test Sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        # CAT GCT AGC TAG CTA CGA TCG
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        exp_codon = ['G', 'T', 'A']
        exp_mutated = ['G', 'T', 'C']
        res_codon, res_mutated = codons[0].mutate_codon(2, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        codons = self.nt_seq1.find_codons('-0', [(0, 12)])
        exp_codon = ['G', 'A', 'T']
        exp_mutated = ['T', 'A', 'T']
        res_codon, res_mutated = codons[0].mutate_codon(0, 'A')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        # Test Sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]                    # GTA = Valine

        # Mutation at wobble position
        expected = False                     # Synonymous
        result = codon.is_nonsyn(2, 'C')     # GTC = Valine
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(0, 'C')    # CTA = Leucine
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATG ACG TGG TGA
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[2]                   # TGG = Tryptophan

        # Mutation at second position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(1, 'A')    # TAG = STOP
        self.assertEqual(expected, result)

        # Mutation at wobble position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(2, 'A')    # TGA = STOP
        self.assertEqual(expected, result)

        # Testing mutation at position 2 in ACG
        codon = codons[1]                   # ACG = Threonine
        expected = False                    # Synonymous
        result = codon.is_nonsyn(2, 'T')    # ACT = Threonine
        self.assertEqual(expected, result)

        # Testing mutation at position 1 in TGA
        expected = False                    # Synonymous
        codon = codons[3]                   # TGA = STOP
        result = codon.is_nonsyn(1, 'A')    # TAA = STOP
        self.assertEqual(expected, result)

    def testIsStop(self):
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.introduces_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.introduces_stop(1, 'A')
        self.assertEqual(expected, result)

        # G to A mutation in middle position (TAA = STOP)
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[3]
        expected = True
        result = codon.introduces_stop(1, 'A')
        self.assertEqual(expected, result)

    def testIsStart(self):
        # Testing sequence 1
        # GTA CGA TCG ATC GAT GCT AGC
        codons = self.nt_seq1.find_codons('+0', [(0, 12)])
        codon = codons[0]               # GTA = Valine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 4
        # ATG ACG TGG TGA
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[0]               # first Methionine
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

        # Testing sequence 5
        # ATG ATG CCC TAA
        codons = self.nt_seq5.find_codons('+0', [(0, 12)])
        codon = codons[1]               # second Methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
