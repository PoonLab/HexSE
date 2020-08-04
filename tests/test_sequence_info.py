import random
import unittest

from src.sequence_info import Sequence

MAX_DIFF = None
DN_VALUES = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
DS_VALUES = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
KAPPA = 0.3
MU = 0.0005


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
        self.sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, KAPPA, MU, pi1, DN_VALUES, DS_VALUES)

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
        self.assertEqual(self.sequence1.mu, new_sequence1.mu)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
        self.assertEqual(self.sequence1.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence1.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence1.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence1.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence1.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence1.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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
        expected = {'A': 0.24, 'C': 0.24, 'G': 0.29, 'T': 0.24}
        result = Sequence.get_frequency_rates(str(self.sequence1))
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.get_substitution_rates(nt)
        exp_sub_rates = {'A': 3.974321933436628e-05, 'C': 0.00015870631784843496,
                         'G': None,                  'T': 0.00015870631784843496}
        exp_dN_values = {'A': (1, 0, 0, 0), 'C': (0, 0, 1, 0), 'G': None, 'T': (0, 0, 1, 0)}
        exp_dS_values = {'A': (0, 0, 1, 0), 'C': (1, 0, 0, 0), 'G': None, 'T': (1, 0, 0, 0)}
        exp_total_rate = 0.0003571558550312362
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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


class TestSequence2(unittest.TestCase):
    """
    Sequence: ATGAATAAACCCGTATGA
    ORFs (indexing relative to forward strand)
        (0, 18) (+)
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
        self.sequence2 = Sequence(s2, sorted_orfs, KAPPA, MU, pi2, DN_VALUES, DS_VALUES)

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
        new_sequence1 = self.sequence2.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence2, new_sequence1)
        self.assertEqual(self.sequence2.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence2.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence2.mu, new_sequence1.mu)
        self.assertEqual(self.sequence2.pi, new_sequence1.pi)
        self.assertEqual(self.sequence2.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence2.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence2.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence2.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence2.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence2.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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
        expected = {'A': 0.44, 'C': 0.17, 'G': 0.17, 'T': 0.22}
        result = Sequence.get_frequency_rates(str(self.sequence2))
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(4001)   # Set seed value to initialize pseudo-random number generator

        nt = self.sequence2.nt_sequence[11]
        self.sequence2.get_substitution_rates(nt)
        exp_sub_rates = {'A': 7.714939618092219e-05, 'C': None,         'G': 7.714939618092219e-05, 'T': 8.5e-05}
        exp_dN_values = {'A': (0, 0, 0, 1),          'C': None,         'G': (0, 0, 0, 1),          'T': (0, 0, 0, 1)}
        exp_dS_values = {'A': (0, 1, 0, 0),          'C': None,         'G': (0, 1, 0, 0),          'T': (0, 0, 0, 1)}
        exp_total_rate = 0.0002392987923618444
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        random.seed(4001)  # Set seed value to initialize pseudo-random number generator
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.get_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_dN_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_dS_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0

        # Substitution rates are 0 because a STOP codon would be introduced
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, MU, pi3, DN_VALUES, DS_VALUES)

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
        new_sequence1 = self.sequence3.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence3, new_sequence1)
        self.assertEqual(self.sequence3.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence3.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence3.mu, new_sequence1.mu)
        self.assertEqual(self.sequence3.pi, new_sequence1.pi)
        self.assertEqual(self.sequence3.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence3.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence3.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence3.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence3.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence3.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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
        expected = {'A': 0.25, 'C': 0.08, 'G': 0.42, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence3))
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.get_substitution_rates(nt)
        exp_sub_rates = {'A': 2.2956308569414034e-05, 'C': 0.000125,     'G': 1.679018660974705e-05, 'T': None}
        exp_dN_values = {'A': (0, 1, 0, 0),           'C': (0, 1, 0, 0), 'G': (1, 0, 0, 0),          'T': None}
        exp_dS_values = {'A': (0, 0, 1, 0),           'C': (0, 1, 0, 0), 'G': (0, 1, 0, 0),          'T': None}
        exp_total_rate = 0.0001647464951791611
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        # Tests a nonsense mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[8]
        self.sequence3.get_substitution_rates(nt)
        exp_sub_rates = {'A': 0.0,  'C': 3.856659839661558e-05, 'G': None, 'T': 6.3e-05}
        exp_dN_values = {'A': None, 'C': (0, 1, 0, 0),          'G': None, 'T': (0, 1, 0, 0)}
        exp_dS_values = {'A': None, 'C': (0, 0, 1, 0),          'G': None, 'T': (0, 1, 0, 0)}
        exp_total_rate = 0.00010156659839661558
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, MU, pi4, DN_VALUES, DS_VALUES)

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
        new_sequence1 = self.sequence4.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence4, new_sequence1)
        self.assertEqual(self.sequence4.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence4.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence4.mu, new_sequence1.mu)
        self.assertEqual(self.sequence4.pi, new_sequence1.pi)
        self.assertEqual(self.sequence4.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence4.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence4.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence4.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence4.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence4.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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
        expected = {'A': 0.33, 'C': 0.25, 'G': 0.17, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence4))
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        random.seed(5001)

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.get_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 2.6726498343139115e-05, 'G': 0.0, 'T': 4.95e-05}
        exp_dN_values = {'A': None, 'C': (0, 0, 1, 0), 'G': None, 'T': (1, 0, 0, 0)}
        exp_dS_values = {'A': None, 'C': (0, 0, 0, 1), 'G': None, 'T': (1, 0, 0, 0)}
        exp_total_rate = 7.622649834313912e-05
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, MU, pi5, DN_VALUES, DS_VALUES)

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
        new_sequence1 = self.sequence5.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence5, new_sequence1)
        self.assertEqual(self.sequence5.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence5.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence5.mu, new_sequence1.mu)
        self.assertEqual(self.sequence5.pi, new_sequence1.pi)
        self.assertEqual(self.sequence5.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence5.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence5.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence5.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence5.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence5.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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
        expected = {'A': 0.38, 'C': 0.19, 'G': 0.19, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence5))
        self.assertEqual(expected, result)

    def testGetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.get_substitution_rates(nt)
        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        exp_dN_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_dS_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, MU, pi6, DN_VALUES, DS_VALUES)

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
        new_sequence1 = self.sequence6.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence6, new_sequence1)
        self.assertEqual(self.sequence6.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence6.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence6.mu, new_sequence1.mu)
        self.assertEqual(self.sequence6.pi, new_sequence1.pi)
        self.assertEqual(self.sequence6.dN_values, new_sequence1.dN_values)
        self.assertEqual(self.sequence6.dS_values, new_sequence1.dS_values)
        self.assertEqual(self.sequence6.is_circular, new_sequence1.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence6.event_tree, new_sequence1.event_tree)
        self.assertCountEqual(self.sequence6.event_tree, new_sequence1.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence6.nt_sequence):
            new_nt = new_sequence1.nt_sequence[pos]
            self.assertIsNot(nt, new_nt)

            self.assertEqual(nt.state, new_nt.state)
            self.assertEqual(nt.pos_in_seq, new_nt.pos_in_seq)
            self.assertEqual(nt.complement_state, new_nt.complement_state)
            self.assertEqual(nt.rates, new_nt.rates)
            self.assertEqual(nt.mutation_rate, new_nt.mutation_rate)

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

    def testGetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.get_substitution_rates(nt)
        exp_sub_rates = {'A': 0.000115, 'C': 3.45e-05, 'G': None, 'T': 3.45e-05}
        # No dN and dS values because nt treated as synonymous mutation
        exp_dN_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_dS_values = {'A': None, 'C': None, 'G': None, 'T': None}
        exp_total_rate = 0.000184
        self.assertEqual(exp_sub_rates, nt.rates)
        self.assertEqual(exp_dN_values, nt.dN_keys)
        self.assertEqual(exp_dS_values, nt.dS_keys)
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


class TestSequenceInfo(unittest.TestCase):

    def testCreateEventTree(self):
        result = self.sequence4.create_event_tree()
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.25},
                              'C': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': None,
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1}},
                                    'stationary_frequency': 0.08},
                              'G': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'G': None,
                                                'T': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3}},
                                    'stationary_frequency': 0.42},
                              'T': {'from_nt': {'A': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'C': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': False,
                                                      'kappa': 1},
                                                'G': {'is_nonsyn': {
                                                        'dN': {},
                                                        'dS': {}},
                                                      'is_syn': [],
                                                      'is_trv': True,
                                                      'kappa': 0.3},
                                                'T': None},
                                    'stationary_frequency': 0.25}}}
        self.assertEqual(expected, result)

    def testNtsOnTips(self):
        # Get the Nucleotide objects
        a3 = self.sequence4.nt_sequence[3]
        c4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        t6 = self.sequence4.nt_sequence[6]
        g7 = self.sequence4.nt_sequence[7]
        g8 = self.sequence4.nt_sequence[8]
        t9 = self.sequence4.nt_sequence[9]
        g10 = self.sequence4.nt_sequence[10]
        a11 = self.sequence4.nt_sequence[11]

        self.sequence4.nt_in_event_tree(g5)
        result = self.sequence4.get_nts_on_tips()

        expected = \
            {'to_nt':
                 {'A': {'stationary_frequency': 0.25,
                        'from_nt': {'A': None,
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6, t9]},
                                                        'dS': {(0, 0, 0, 1): [t6],
                                                               (1, 0, 0, 0): [t9]}},
                                          'is_syn': [],
                                          'nts_in_subs': {t6: None, t9: None},
                                          'number_of_events': 0},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 1, 0): [c4]},
                                                        'dS': {(0, 0, 1, 0): [c4]}},
                                          'is_syn': [],
                                          'nts_in_subs': {c4: None},
                                          'number_of_events': 0},
                                    'G': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'dN': {}, 'dS': {}},
                                          'is_syn': [g5],
                                          'nts_in_subs': {g5: None},
                                          'number_of_events': 1}},
                        'events_for_nt': 1},
                  'T': {'stationary_frequency': 0.25,
                        'from_nt': {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3, a11]},
                                                        'dS': {(1, 0, 0, 0): [a3],
                                                               (0, 0, 1, 0): [a11]}},
                                          'is_syn': [],
                                          'nts_in_subs': {a3: None, a11: None},
                                          'number_of_events': 0},
                                    'T': None,
                                    'C': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'dN': {(0, 0, 0, 1): [c4]},
                                                        'dS': {(0, 1, 0, 0): [c4]}},
                                          'is_syn': [],
                                          'nts_in_subs': {c4: None},
                                          'number_of_events': 0},
                                    'G': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 1, 0, 0): [g7, g8],
                                                               (1, 0, 0, 0): [g10]},
                                                        'dS': {(0, 0, 0, 1): [g7, g10],
                                                               (0, 0, 1, 0): [g8]}},
                                          'is_syn': [g5],
                                          'nts_in_subs': {g7: None, g8: None, g10: None, g5: None},
                                          'number_of_events': 1}},
                        'events_for_nt': 1},
                  'C': {'stationary_frequency': 0.08,
                        'from_nt': {'A': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(1, 0, 0, 0): [a3],
                                                               (0, 0, 1, 0): [a11]},
                                                        'dS': {(0, 0, 1, 0): [a3],
                                                               (1, 0, 0, 0): [a11]}},
                                          'is_syn': [],
                                          'nts_in_subs': {a3: None, a11: None},
                                          'number_of_events': 0},
                                    'T': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'dN': {(0, 1, 0, 0): [t6],
                                                               (0, 0, 0, 1): [t9]},
                                                        'dS': {(0, 1, 0, 0): [t6],
                                                               (1, 0, 0, 0): [t9]}},
                                          'is_syn': [],
                                          'nts_in_subs': {t6: None, t9: None},
                                          'number_of_events': 0},
                                    'C': None,
                                    'G': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 1, 0): [g7],
                                                               (0, 0, 0, 1): [g8, g10]},
                                                        'dS': {(0, 1, 0, 0): [g7],
                                                               (0, 0, 1, 0): [g8, g10]}},
                                          'is_syn': [g5],
                                          'nts_in_subs': {g7: None, g8: None, g10: None, g5: None},
                                          'number_of_events': 1}},
                        'events_for_nt': 1},
                  'G': {'stationary_frequency': 0.42,
                        'from_nt': {'A': {'is_trv': False,
                                          'kappa': 1,
                                          'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3],
                                                               (0, 1, 0, 0): [a11]},
                                                        'dS': {(1, 0, 0, 0): [a3],
                                                               (0, 0, 1, 0): [a11]}},
                                          'is_syn': [],
                                          'nts_in_subs': {a3: None, a11: None},
                                          'number_of_events': 0},
                                    'T': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6],
                                                               (0, 1, 0, 0): [t9]},
                                                        'dS': {(0, 1, 0, 0): [t6],
                                                               (1, 0, 0, 0): [t9]}},
                                          'is_syn': [],
                                          'nts_in_subs': {t6: None, t9: None},
                                          'number_of_events': 0},
                                    'C': {'is_trv': True,
                                          'kappa': 0.3,
                                          'is_nonsyn': {'dN': {(0, 0, 1, 0): [c4]},
                                                        'dS': {(1, 0, 0, 0): [c4]}},
                                          'is_syn': [],
                                          'nts_in_subs': {c4: None},
                                          'number_of_events': 0},
                                    'G': None},
                        'events_for_nt': 0}},
             'total_events': 3}
        self.assertEqual(expected, result)


# ==========================================
# Tests for Nucleotide
# ==========================================
class TestNucleotide(unittest.TestCase):

    def setUp(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        self.nt_seq1 = Sequence(s1, {'+0': [[(0, 21)]]}, kappa, mu, pi1, dN_values, dS_values).nt_sequence

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [[(5, 50)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [[(30, 3)]]}
        pi3 = Sequence.get_frequency_rates(s3)
        self.nt_seq3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, dN_values, dS_values).nt_sequence

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.nt_seq4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values).nt_sequence

        s5 = 'ATGATGCCCTAA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        random.seed(4000)
        self.nt_seq5 = Sequence(s5, sorted_orfs, kappa, mu, pi5, dN_values, dS_values).nt_sequence

    def testGetMutationRate(self):

        seq1_expected_rates = [0.0003571558550312362,               # g0
                               0.0005827273266040434,               # t1
                               0.000192,                            # a2
                               7.2e-05,                             # c3
                               0.00034766669196033067,              # g4
                               0.000192,                            # a5
                               0.00025618838497572276,              # t6
                               0.0010541264093863016,               # c7
                               0.000232,                            # g8
                               0.0003079830809958123,               # a9
                               0.0002408033470208297,               # t10
                               0.0002873431595987048,               # c11
                               0.0001337127084216265,               # g12
                               0.0008629850479040508,               # a13
                               0.00016147550956710228,              # t14
                               0.00012183126917068458,              # g15
                               0.0004350559820278692,               # c16
                               0.000192,                            # t17
                               0.00016586728204163575,              # a18
                               0.0006069065088774146,               # g19
                               0.00021652268848481484]              # c20
        for pos, exp_seq1_rate in enumerate(seq1_expected_rates):
            nt = self.nt_seq1[pos]
            self.assertEqual(exp_seq1_rate, nt.mutation_rate)

        seq2_expected_rates = [0.0005891339975452602,               # t0
                               0.000744,                            # t1
                               0.0005656952120995977,               # t2
                               0.0006931397654833255,               # t3
                               0.004506792415479317,                # t4
                               0.0005785558496056398,               # t5
                               0.00010811765387774733,              # t6
                               0.0016808290584385383,               # t7
                               0.000744,                            # t8
                               0.00028363540788326586,              # t9
                               0.0006008687472494076,               # t10
                               0.0008628671796935054,               # t11
                               0.000744,                            # t12
                               0.000744]                            # t13
        for pos, exp_seq2_rate in enumerate(seq2_expected_rates):
            nt = self.nt_seq2[pos]
            self.assertEqual(exp_seq2_rate, nt.mutation_rate)

        seq3_expected_rates = [0.00022400000000000002,              # a0
                               0.00022400000000000002,              # a1
                               0.00022400000000000002,              # t2
                               0.00022400000000000002,              # t3
                               0.000264,                            # c4
                               0,                                   # a5
                               0,                                   # t6
                               0,                                   # g7
                               0.00010471857777732271,              # a8
                               0.0006302898914741413,               # a9
                               0.00022182549714269624,              # c10
                               6.266381714279864e-05,               # g11
                               0.0002583615486840646,               # a12
                               0.0003090695937097542,               # a13
                               0.0004090398495395135,               # a14
                               0.00034485756178951926,              # a15
                               0.00022400000000000002,              # t16
                               0.00026400417380813366,              # c17
                               0.0002545348481202625,               # t18
                               9.6e-05,                             # g19
                               0.00010471857777732271,              # t20
                               0.00019966332957355537,              # t21
                               0.0003956207184376176,               # c22
                               0.00023764224407273933,              # g23
                               0.0006377200406261789,               # c24
                               0.00022400000000000002,              # t25
                               0.00018821557333319682,              # t26
                               0.0004992019752883201,               # c27
                               0.00022400000000000002,              # a28
                               0.00016970355199247907,              # t29
                               0.0003773596763284332,               # t30
                               0.000264,                            # c31
                               0.00020467702889720896,              # a32
                               0.0010376810049062904,               # t33
                               0.00018821557333319685,              # t34
                               0.00014488998560907796,              # g35
                               0.00024122649834313912,              # c36
                               0.000264,                            # c37
                               0.00018589300994994543,              # c38
                               0.0008604108919944202,               # c39
                               0.000264,                            # c40
                               0.00016970355199247907,              # a41
                               0.00029517211496089,                 # c42
                               0.00022400000000000002,              # a43
                               0.0005288819356437898,               # a44
                               0.0004300188373987396,               # t45
                               0.0005489824183629611,               # c46
                               0.00030241338586062196,              # t47
                               0.0002597879680797651,               # a48
                               2.7718726670232408e-05,              # g49
                               9.6e-05,                             # g50
                               0.000264,                            # c51
                               0.000264,                            # c52
                               0.00022400000000000002,              # t53
                               0.00022400000000000002,              # a54
                               0.000264,                            # c55
                               0.000264,                            # c56
                               0.000264]                            # c57
        for pos, exp_seq3_rate in enumerate(seq3_expected_rates):
            nt = self.nt_seq3[pos]
            self.assertEqual(exp_seq3_rate, nt.mutation_rate)

        seq4_expected_rates = [0,                                   # a0
                               0,                                   # t1
                               0,                                   # g2
                               0.0002814065924469639,               # a3
                               7.808147629798601e-05,               # c4
                               0.00033600000000000004,              # g5
                               0.00019227093666957,                 # t6
                               9.120751350437504e-05,               # g7
                               0.0001339496956925102,               # g8
                               0.00018274734722965082,              # t9
                               0.0001414794858434976,               # g10
                               0.00025206864734794677]              # a11
        for pos, exp_seq4_rate in enumerate(seq4_expected_rates):
            nt = self.nt_seq4[pos]
            self.assertEqual(exp_seq4_rate, nt.mutation_rate)

        seq5_expected_rates = [0,                                   # a0
                               0,                                   # t1
                               0,                                   # g2
                               0.0006403806519537574,               # a3
                               0.0005152115627357426,               # t4
                               0.00035554475079678744,              # g5
                               0.0003177555392993143,               # c6
                               0.00014555683758674118,              # c7
                               0.00019999999999999998,              # c8
                               0.00014249115743216948,              # t9
                               5.345299668627823e-05,               # a10
                               7.980232731162653e-05]               # a11
        for pos, exp_seq6_rate in enumerate(seq5_expected_rates):
            nt = self.nt_seq5[pos]
            self.assertEqual(exp_seq6_rate, nt.mutation_rate)


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
        result = codon.is_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.is_stop(1, 'A')
        self.assertEqual(expected, result)

        # G to A mutation in middle position (TAA = STOP)
        codons = self.nt_seq4.find_codons('+0', [(0, 12)])
        codon = codons[3]
        expected = True
        result = codon.is_stop(1, 'A')
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
