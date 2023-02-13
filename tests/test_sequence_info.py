import random
import unittest
import sys
import numpy as np

from hexse.sequence_info import Sequence, Nucleotide

MAX_DIFF = None
CAT_VALUES = {'mu1': 0.051710707633483066, 'mu2': 0.15181054803756722, 'mu3': 0.26809045653750935,
              'mu4': 0.4186255904232465, 'mu5': 0.6442570794470408, 'mu6': 1.2255056178040284}
KAPPA = 0.3
GLOBAL_RATE = 0.0005


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

        orfs = {'+0': [{'coords': [[0, 21]], 
                        'omega_classes': 3, 'omega_shape': 1.5, 
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404], 
                        'orf_map': np.array([1])}], 
                '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

        s1 = 'GTACGATCGATCGATGCTAGC'
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, orfs, KAPPA, GLOBAL_RATE, pi1, CAT_VALUES)
        orf = {'coords': [[0, 21]],
               'omega_classes': 3, 'omega_shape': 1.5,
               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}
        self.seq1_codons = self.sequence1.find_codons('+0', orf)

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

        # self.assertEqual(self.sequence1.orfs, new_sequence1.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence1.orfs[strand], new_sequence1.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence1.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence1.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence1.pi, new_sequence1.pi)
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
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.24, 'C': 0.24, 'G': 0.29, 'T': 0.24}
        result = Sequence.get_frequency_rates(str(self.sequence1))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None, 
                                                'C': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'G': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'T': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}}}, 
                              'C': {'from_nt': {'A': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'C': None, 
                                                'G': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'T': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}}}, 
                              'G': {'from_nt': {'A': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'C': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'G': None, 
                                                'T': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}}}, 
                              'T': {'from_nt': {'A': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'C': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'G': {
                                                    'mu1': {(1,): {}}, 
                                                    'mu2': {(1,): {}}, 
                                                    'mu3': {(1,): {}}, 
                                                    'mu4': {(1,): {}}, 
                                                    'mu5': {(1,): {}}, 
                                                    'mu6': {(1,): {}}}, 
                                                'T': None}}}}
        result = self.sequence1.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)

        exp_sub_rates = {'A': 1.280932279322075e-06,
                         'C': 1.9922704712789297e-06,
                         'G': None, 
                         'T': 1.9922704712789297e-06}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omegas = {'A': (0.1708353283825978,), 'C': (0.1708353283825978,), 'G': None, 'T': (0.1708353283825978,)}
        self.assertEqual(exp_omegas, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1', 'C': 'mu3', 'T': 'mu3'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 5.265473221879935e-06
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {(0.1708353283825978,): {'value': 0.1708353283825978, 'nt_events': 1},
                            (-1,): {'value': 1, 'nt_events': 1},
                            (0.4810288100937172,): {'value': 0.4810288100937172, 'nt_events': 1}, 
                            (1.1481358615121404,): {'value': 1.1481358615121404, 'nt_events': 1}}
        self.assertEqual(exp_total_omegas, self.sequence1.total_omegas)

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
        orf = {'coords': [[0, 21]],
               'omega_classes': 3, 'omega_shape': 1.5,
               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}
        result = self.sequence1.find_codons('+0', orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStopStartCodon(self):
        nt = self.sequence1.nt_sequence[0]  # G in codon GTA
        expected = False
        result = self.sequence1.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None,
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}},
                                                      'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}, ((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 0, 1, 0),): {'prob': 0.010258648679344635, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}, ((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': None,
                                                                       'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}, ((0, 0, 0, 1),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}, ((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}, ((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 0, 0, 1),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 0, 0, 1),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'C': None,
                                                                       'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 1, 0, 0),): {'prob': 0.0018007182628609634, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 0, 0, 1),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}, ((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}, ((1, 0, 0, 0),): {'prob': 0.0015264218065529106, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}, ((0, 0, 1, 0),): {'prob': 0.010258648679344635, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'G': None}}}}

        result = self.sequence1.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(9001)

        g0 = self.sequence1.nt_sequence[0]
        t1 = self.sequence1.nt_sequence[1]
        a2 = self.sequence1.nt_sequence[2]
        c3 = self.sequence1.nt_sequence[3]
        g4 = self.sequence1.nt_sequence[4]
        a5 = self.sequence1.nt_sequence[5]
        t6 = self.sequence1.nt_sequence[6]
        c7 = self.sequence1.nt_sequence[7]
        g8 = self.sequence1.nt_sequence[8]
        a9 = self.sequence1.nt_sequence[9]
        t10 = self.sequence1.nt_sequence[10]
        c11 = self.sequence1.nt_sequence[11]
        g12 = self.sequence1.nt_sequence[12]
        a13 = self.sequence1.nt_sequence[13]
        t14 = self.sequence1.nt_sequence[14]
        g15 = self.sequence1.nt_sequence[15]
        c16 = self.sequence1.nt_sequence[16]
        t17 = self.sequence1.nt_sequence[17]
        a18 = self.sequence1.nt_sequence[18]
        g19 = self.sequence1.nt_sequence[19]
        c20 = self.sequence1.nt_sequence[20]

        exp_event_tree = \
                        {"to_nt": {"A": {"from_nt": {"A": None,
                                                                "C": {"mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu5": {(1,): {(1.1481358615121404,): [c20],
                                                                                    "nt_events": 1, "region_weight": 1.1481358615121404,},
                                                                                    "nt_events": 1,},
                                                                        "mu6": {(1,): {(-1,): [c3, c11], (0.4810288100937172,): [c16],
                                                                                    "nt_events": 3, "region_weight": 2.4810288100937172,},
                                                                                    "nt_events": 3,},
                                                                                    "nt_events": 4,},
                                                                    "G": {
                                                                        "mu1": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [g15],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [g4],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu4": {
                                                                            (1,): {(-1,): [g8], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [g0, g0],
                                                                                (1.1481358615121404,): [g12, g19],
                                                                                "nt_events": 3,
                                                                                "region_weight": 2.4671070514068787,
                                                                            },
                                                                            "nt_events": 3,
                                                                        },
                                                                        "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "nt_events": 6,
                                                                    },
                                                                    "T": {
                                                                        "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [t10],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [t14],
                                                                                "nt_events": 1,
                                                                                "region_weight": 1.1481358615121404,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {(-1,): [t17], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [t1],
                                                                                (0.4810288100937172,): [t6],
                                                                                "nt_events": 2,
                                                                                "region_weight": 0.651864138476315,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "nt_events": 5,
                                                                    },
                                                                },
                                                                "nt_events": 15,
                                                            },
                                                            "C": {
                                                                "from_nt": {
                                                                    "A": {
                                                                        "mu1": {
                                                                            (1,): {(-1,): [a5], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [a13, a18],
                                                                                "nt_events": 2,
                                                                                "region_weight": 2.296271723024281,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [a9],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {(-1,): [a2], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "nt_events": 5,
                                                                    },
                                                                    "C": None,
                                                                    "G": {
                                                                        "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {
                                                                            (1,): {(-1,): [g8], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [g12, g19],
                                                                                "nt_events": 2,
                                                                                "region_weight": 2.296271723024281,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu5": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [g0, g0],
                                                                                (0.4810288100937172,): [g4],
                                                                                "nt_events": 2,
                                                                                "region_weight": 0.651864138476315,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu6": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [g15],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "nt_events": 6,
                                                                    },
                                                                    "T": {
                                                                        "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [t6],
                                                                                (0.1708353283825978,): [t10],
                                                                                "nt_events": 2,
                                                                                "region_weight": 0.651864138476315,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu6": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [t1],
                                                                                (-1,): [t14, t17],
                                                                                "nt_events": 3,
                                                                                "region_weight": 2.170835328382598,
                                                                            },
                                                                            "nt_events": 3,
                                                                        },
                                                                        "nt_events": 5,
                                                                    },
                                                                },
                                                                "nt_events": 16,
                                                            },
                                                            "G": {
                                                                "from_nt": {
                                                                    "A": {
                                                                        "mu1": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [a18],
                                                                                "nt_events": 1,
                                                                                "region_weight": 1.1481358615121404,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu2": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [a9],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [a13],
                                                                                "nt_events": 1,
                                                                                "region_weight": 1.1481358615121404,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {(-1,): [a5], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {
                                                                            (1,): {(-1,): [a2], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "nt_events": 5,
                                                                    },
                                                                    "C": {
                                                                        "mu1": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [c16],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [c3],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [c7],
                                                                                (0.1708353283825978,): [c11],
                                                                                "nt_events": 2,
                                                                                "region_weight": 0.651864138476315,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu6": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [c20],
                                                                                "nt_events": 1,
                                                                                "region_weight": 1.1481358615121404,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "nt_events": 5,
                                                                    },
                                                                    "G": None,
                                                                    "T": {
                                                                        "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [t1],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [t10],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [t6],
                                                                                (1.1481358615121404,): [t14],
                                                                                (-1,): [t17],
                                                                                "nt_events": 3,
                                                                                "region_weight": 2.6291646716058574,
                                                                            },
                                                                            "nt_events": 3,
                                                                        },
                                                                        "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "nt_events": 5,
                                                                    },
                                                                },
                                                                "nt_events": 15,
                                                            },
                                                            "T": {
                                                                "from_nt": {
                                                                    "A": {
                                                                        "mu1": {
                                                                            (1,): {
                                                                                (1.1481358615121404,): [a13],
                                                                                "nt_events": 1,
                                                                                "region_weight": 1.1481358615121404,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [a9],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.1708353283825978,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (-1,): [a5],
                                                                                (1.1481358615121404,): [a18],
                                                                                "nt_events": 2,
                                                                                "region_weight": 2.1481358615121406,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {(-1,): [a2], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "nt_events": 5,
                                                                    },
                                                                    "C": {
                                                                        "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [c16],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu5": {
                                                                            (1,): {(-1,): [c11], "nt_events": 1, "region_weight": 1},
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [c7],
                                                                                (-1,): [c20],
                                                                                "nt_events": 2,
                                                                                "region_weight": 1.4810288100937172,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "nt_events": 4,
                                                                    },
                                                                    "G": {
                                                                        "mu1": {
                                                                            (1,): {
                                                                                (0.1708353283825978,): [g0, g0],
                                                                                (1.1481358615121404,): [g19],
                                                                                "nt_events": 2,
                                                                                "region_weight": 1.3189711898947383,
                                                                            },
                                                                            "nt_events": 2,
                                                                        },
                                                                        "mu2": {
                                                                            (1,): {
                                                                                (-1,): [g8],
                                                                                (1.1481358615121404,): [g12],
                                                                                (0.4810288100937172,): [g15],
                                                                                "nt_events": 3,
                                                                                "region_weight": 2.629164671605858,
                                                                            },
                                                                            "nt_events": 3,
                                                                        },
                                                                        "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "mu5": {
                                                                            (1,): {
                                                                                (0.4810288100937172,): [g4],
                                                                                "nt_events": 1,
                                                                                "region_weight": 0.4810288100937172,
                                                                            },
                                                                            "nt_events": 1,
                                                                        },
                                                                        "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                                        "nt_events": 6,
                                                                    },
                                                                    "T": None,
                                                                },
                                                                "nt_events": 15,
                                                        },
                                                    }
                                                }

        res_omega_key = self.sequence1.nt_in_event_tree(g0)
        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 61
        res_count = self.sequence1.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 10)
        nt.rates = {'A': None, 'C': 0.3132846693781597, 'G': 1.0442822312605322, 'T': 0.05572713744568437}
        nt.get_mutation_rate()
        exp_mutation_rate = 1.4132940380843764
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': 1.595897111087155e-05, 'C': 4.787691333261464e-06, 'G': None, 'T': 3.8427968379662253e-07},
            {'A': 7.536947567481705e-06, 'C': 2.512315855827235e-05, 'G': 1.648775562437735e-06, 'T': None},
            {'A': None, 'C': 2.3193254860093468e-05, 'G': 0.00014706067413648342, 'T': 2.3193254860093468e-05},
            {'A': 4.4118202240945026e-05, 'C': None, 'G': 4.642532399005903e-06, 'T': None},
            {'A': 1.869908882932933e-05, 'C': 1.348092040995756e-05, 'G': None, 'T': 1.348092040995756e-05},
            {'A': None, 'C': 1.8615854748053904e-06, 'G': 7.73108495336449e-05, 'T': 1.5070521255236873e-05},
            {'A': 2.1222126327435754e-05, 'C': 8.763029673862326e-06, 'G': 1.1156623787551084e-05, 'T': None},
            {'A': None, 'C': None, 'G': 7.2493549068986665e-06, 'T': 7.074042109145251e-05},
            {'A': 6.070071061137074e-05, 'C': 6.603758839634174e-06, 'G': None, 'T': 6.603758839634174e-06},
            {'A': None, 'C': 2.574577447535311e-06, 'G': 3.112152579112792e-06, 'T': 1.648775562437735e-06},
            {'A': 9.336457737338377e-07, 'C': 3.112152579112792e-06, 'G': 2.574577447535311e-06, 'T': None},
            {'A': 4.4118202240945026e-05, 'C': None, 'G': 2.574577447535311e-06, 'T': 7.73108495336449e-05},
            {'A': 0.00010725572525720244, 'C': 1.338948562667462e-05, 'G': None, 'T': 7.582012344561795e-06},
            {'A': None, 'C': 1.1080953622075548e-05, 'G': 5.7676686349394705e-05, 'T': 2.137353042894174e-06},
            {'A': 1.7303005904818413e-05, 'C': 0.00014706067413648342, 'G': 2.6629007650064053e-05, 'T': None},
            {'A': 3.606779323495576e-06, 'C': 2.5643402645651533e-05, 'G': None, 'T': 3.1765982567750934e-06},
            {'A': 2.1222126327435754e-05, 'C': None, 'G': 8.954762458333844e-07, 'T': 2.4164516356328888e-05},
            {'A': 2.3193254860093468e-05, 'C': 0.00014706067413648342, 'G': 2.3193254860093468e-05, 'T': None},
            {'A': None, 'C': 1.1080953622075548e-05, 'G': 7.1245101429805794e-06, 'T': 1.7303005904818413e-05},
            {'A': 0.00010725572525720244, 'C': 1.338948562667462e-05, 'G': None, 'T': 2.5826349268304602e-06},
            {'A': 2.6629007650064053e-05, 'C': None, 'G': 5.0653690138274265e-05, 'T': 0.00014706067413648342}
            ]

        exp_total_rates = [2.1130942127929635e-05,
                           3.430888168819179e-05,
                           0.00019344718385667036,
                           4.876073463995093e-05,
                           4.5660929649244447e-05,
                           9.424295626368716e-05,
                           4.1141779788849165e-05,
                           7.798977599835118e-05,
                           7.390822829063909e-05,
                           7.3355055890858385e-06,
                           6.62037580038194e-06,
                           0.00012400362922212525,
                           0.00012822722322843884,
                           7.089499301436444e-05,
                           0.0001909926876913659,
                           3.2426780225922204e-05,
                           4.6282118929598024e-05,
                           0.00019344718385667036,
                           3.550846966987454e-05,
                           0.00012322784581070753,
                           0.00022434337192482173]

        for pos, nt in enumerate(self.sequence1.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq1_codons[0]     # GTA
        nt = self.sequence1.nt_sequence[1]  # T
        expected = 1
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['G', 'T', 'A']
        exp_mutated = ['G', 'T', 'C']
        res_codon, res_mutated = self.seq1_codons[0].mutate_codon(2, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        exp_codon = ['G', 'A', 'T']
        exp_mutated = ['T', 'A', 'T']
        res_codon, res_mutated = self.seq1_codons[4].mutate_codon(0, 'T')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq1_codons[0]          # GTA = Valine

        # Mutation at wobble position
        expected = False                     # Synonymous
        result = codon.is_nonsyn(2, 'C')     # GTC = Valine
        self.assertEqual(expected, result)

        # Mutation at first position
        expected = True                     # Non-synonymous
        result = codon.is_nonsyn(0, 'C')    # CTA = Leucine
        self.assertEqual(expected, result)

    def testIsStop(self):
        codon = self.seq1_codons[2]   # TCG = Serine

        # T to G mutation at first position (GCG = Alanine)
        expected = False
        result = codon.creates_stop(0, 'G')
        self.assertEqual(expected, result)

        # C to A mutation in middle position (TAG = STOP)
        expected = True
        result = codon.creates_stop(1, 'A')
        self.assertEqual(expected, result)

    def testIsStart(self):
        codon = self.seq1_codons[0]               # GTA = Valine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq1_codons[2]     # TCG
        expected = True
        result = codon.creates_stop(1, 'A')     # TCG --> TAG
        self.assertEqual(expected, result)

        result = codon.creates_stop(1, 'G')     # TCG -> TGG
        expected = False
        self.assertEqual(expected, result)

    def testGetCodons(self):
        expected = self.seq1_codons
        result = self.sequence1.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("GTACGATCGATCGATGCTAGC")]
        for pos, nt in enumerate(nts):
            result = self.sequence1.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("GTACGATCGATCGATGCTAGC")]
        for pos, nt in enumerate(nts):
            result = self.sequence1.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])


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

        orfs = {'+0': [{'coords': [[0, 21]],
                        'omega_classes': 3, 'omega_shape': 1.5,
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                        'orf_map': np.array([1, 0])}],
                '+1': [], '+2': [], '-0': [], '-1': [],
                '-2': [{'coords': [[3, 15]],
                        'omega_classes': 4, 'omega_shape': 1.25,
                        'omega_values': [0.09199853806558903, 0.27043066909631136,
                                         0.5158061369385518, 1.1217646558655263], 
                        'orf_map': np.array([0, 1])}]}

        s2 = 'ATGAATAAACCCGTATGA'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, orfs, KAPPA, GLOBAL_RATE, pi2, CAT_VALUES)
        plus_orf = {'coords': [[0, 21]],
                    'omega_classes': 3, 'omega_shape': 1.5,
                    'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}
        minus_orf = {'coords': [[3, 15]],
                     'omega_classes': 4, 'omega_shape': 1.25,
                     'omega_values': [0.09199853806558903, 0.27043066909631136,
                                      0.5158061369385518, 1.1217646558655263]}
        self.plus_0_codons = self.sequence2.find_codons('+0', plus_orf)
        self.minus_2_codons = self.sequence2.find_codons('-2', minus_orf)

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

        # self.assertEqual(self.sequence2.orfs, new_sequence2.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence2.orfs[strand], new_sequence2.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence2.kappa, new_sequence2.kappa)
        self.assertEqual(self.sequence2.global_rate, new_sequence2.global_rate)
        self.assertEqual(self.sequence2.pi, new_sequence2.pi)
        self.assertEqual(self.sequence2.is_circular, new_sequence2.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence2.event_tree, new_sequence2.event_tree)
        self.assertCountEqual(self.sequence2.event_tree, new_sequence2.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence2.nt_sequence):
            new_nt = new_sequence2.nt_sequence[pos]
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
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testCreateEventTree(self):
        expected = \
                    {
                        "to_nt": {
                            "A": {
                                "from_nt": {
                                    "A": None,
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                }
                            },
                            "C": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "C": None,
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                }
                            },
                            "G": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "G": None,
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                }
                            },
                            "T": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}},
                                    },
                                    "T": None,
                                }
                            },
                        }
                    }

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

        exp_sub_rates = {'A': 1.1974784480129713e-05, 'C': None, 'G': 3.5055586634238734e-05, 'T': 1.4475135109970076e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (-1, 1.1217646558655263), 'C': None, 'G': (-1, 1.1217646558655263), 'T': (-1, 1.1217646558655263)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'G': 'mu6', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 6.150550622433852e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_total_omegas = {
                                (0.1708353283825978, 0.27043066909631136): {
                                    "value": 0.04619911215979399,
                                    "nt_events": 1,
                                },
                                (0.1708353283825978, -1): {"value": 0.1708353283825978, "nt_events": 1},
                                (-1, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                                (0.4810288100937172, -1): {"value": 0.4810288100937172, "nt_events": 1},
                                (0.4810288100937172, 1.1217646558655263): {
                                    "value": 0.5396011176161822,
                                    "nt_events": 1,
                                },
                                (-1, 1.1217646558655263): {"value": 1.1217646558655263, "nt_events": 1},
                            }

        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation would destroy START CODON
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        # Tuples are of length 1 because the nucleotide is part of only 1 ORF
        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation in STOP codon on reverse strand ORF
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        # Tuples are of length 2 because the nucleotide is part of 2 ORFs
        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Non-synonymous mutation in both (+) and (-) strands
        random.seed(4001)
        nt = self.sequence2.nt_sequence[10]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': 5.760216329605508e-06, 'C': None, 'G': 1.6862747125805073e-05, 'T': 6.962957017894693e-06}
        self.assertEqual(exp_sub_rates, nt.rates) # AssertionError

        exp_omega_keys = {'A': (0.4810288100937172, 1.1217646558655263), 'C': None, 'G': (0.4810288100937172, 1.1217646558655263), 'T': (0.4810288100937172, 1.1217646558655263)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'G': 'mu6', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys) # AssertionError

        exp_total_rate = 2.9585920473305272e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate) # AssertionError

        exp_total_omegas = {
                                (0.1708353283825978, 0.27043066909631136): {
                                    "value": 0.04619911215979399,
                                    "nt_events": 1,
                                },
                                (0.1708353283825978, -1): {"value": 0.1708353283825978, "nt_events": 1},
                                (-1, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                                (0.4810288100937172, -1): {"value": 0.4810288100937172, "nt_events": 1},
                                (0.4810288100937172, 1.1217646558655263): {
                                    "value": 0.5396011176161822,
                                    "nt_events": 1,
                                },
                                (-1, 1.1217646558655263): {"value": 1.1217646558655263, "nt_events": 1},
                            }

        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

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
        result = self.sequence2.find_codons('+0', {'coords': [[0, 21]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Test reverse strand ORF
        expected = ['ATG', 'CCC', 'AAA', 'TAA']
        result = self.sequence2.find_codons('-2', {'coords': [[3, 15]],
                                                   'omega_classes': 4, 'omega_shape': 1.25,
                                                   'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                                    0.5158061369385518, 1.1217646558655263]})
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

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None,
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0), (0, 0, 0, 0, 1)): {'prob': 0.07156927976490862, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): {'prob': 6.513801256693088e-05, 'number_of_events': 0}, ((0, 0, 0, 1), (0, 1, 0, 0, 0)): {'prob': 4.2429606811696455e-09, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): {'prob': 2.5645398616067244e-06, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': None,
                                                'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0), (0, 0, 0, 1, 0)): {'prob': 0.00464913259374619, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 1, 0), (0, 0, 0, 0, 1)): {'prob': 0.17082398175834337, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 0, 1), (1, 0, 0, 0, 0)): {'prob': 0.01368788931269199, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): {'prob': 6.513801256693088e-05, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): {'prob': 2.5645398616067244e-06, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): {'prob': 2.5645398616067244e-06, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': None,
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0), (0, 0, 0, 0, 1)): {'prob': 0.17082398175834337, 'number_of_events': 0}, ((0, 0, 0, 1), (0, 0, 1, 0, 0)): {'prob': 0.0010060810470185998, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 0, 1, 0), (0, 0, 1, 0, 0)): {'prob': 0.000451987256652962, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((1, 0, 0, 0), (0, 0, 0, 1, 0)): {'prob': 0.00464913259374619, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 1, 0, 0), (0, 0, 0, 0, 1)): {'prob': 0.07156927976490862, 'number_of_events': 0}, ((0, 0, 0, 1), (0, 0, 1, 0, 0)): {'prob': 0.00042151280621773425, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': None}}}}
        result = self.sequence2.create_probability_tree()
        print(result)
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(4001)

        # All other nucleotides are part of START or STOP codons
        a6 = self.sequence2.nt_sequence[6]
        a7 = self.sequence2.nt_sequence[7]
        a8 = self.sequence2.nt_sequence[8]
        c9 = self.sequence2.nt_sequence[9]
        c10 = self.sequence2.nt_sequence[10]
        c11 = self.sequence2.nt_sequence[11]

        exp_event_tree = \
                          {
                            "to_nt": {
                                "A": {
                                    "from_nt": {
                                        "A": None,
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.4810288100937172, 1.1217646558655263): [c10],
                                                    "nt_events": 1,
                                                    "region_weight": 0.5396011176161822,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.4810288100937172, -1): [c9],
                                                    (-1, 1.1217646558655263): [c11, c11],
                                                    "nt_events": 2,
                                                    "region_weight": 1.6027934659592435,
                                                },
                                                "nt_events": 2,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 3,
                                },
                                "C": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.27043066909631136): [a8],
                                                    "nt_events": 1,
                                                    "region_weight": 0.04619911215979399,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.27043066909631136): [a7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.04619911215979399,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.27043066909631136): [a6],
                                                    "nt_events": 1,
                                                    "region_weight": 0.04619911215979399,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 3,
                                        },
                                        "C": None,
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 3,
                                },
                                "G": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 0.27043066909631136): [a8],
                                                    "nt_events": 1,
                                                    "region_weight": 0.27043066909631136,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, -1): [a6],
                                                    "nt_events": 1,
                                                    "region_weight": 0.1708353283825978,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.27043066909631136): [a7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.04619911215979399,
                                                },
                                                "nt_events": 1,
                                            },
                                            "nt_events": 3,
                                        },
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.4810288100937172, -1): [c9],
                                                    (0.4810288100937172, 1.1217646558655263): [c10],
                                                    "nt_events": 2,
                                                    "region_weight": 1.0206299277098994,
                                                },
                                                "nt_events": 2,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 1.1217646558655263): [c11, c11],
                                                    "nt_events": 1,
                                                    "region_weight": 1.1217646558655263,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": None,
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 6,
                                },
                                "T": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.27043066909631136): [a7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.04619911215979399,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 1,
                                        },
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.4810288100937172, -1): [c9],
                                                    "nt_events": 1,
                                                    "region_weight": 0.4810288100937172,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.4810288100937172, 1.1217646558655263): [c10],
                                                    "nt_events": 1,
                                                    "region_weight": 0.5396011176161822,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 1.1217646558655263): [c11, c11],
                                                    "nt_events": 1,
                                                    "region_weight": 1.1217646558655263,
                                                },
                                                "nt_events": 1,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": None,
                                    },
                                    "nt_events": 4,
                                },
                            }
                        }

        res_omega_key = self.sequence2.nt_in_event_tree(c11)
        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 16
        res_count = self.sequence2.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('T', 2)
        nt.rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        nt.get_mutation_rate()
        exp_mutation_rate = 0
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': 1.9644309348256e-06, 'G': 1.5733528846049123e-05, 'T': None},
            {'A': None, 'C': 1.276448619924954e-06, 'G': 1.2455799727664908e-05, 'T': 8.17445710656087e-07},
            {'A': None, 'C': 4.6289182736460586e-07, 'G': 3.0765154782486437e-06, 'T': None},
            {'A': 1.5032339481933658e-05, 'C': None, 'G': 6.342956741319806e-07, 'T': 1.7116532419066297e-05},
            {'A': 7.115304686096532e-07, 'C': None, 'G': 7.115304686096532e-07, 'T': 2.9549556408649665e-05},
            {'A': 3.5055586634238734e-05, 'C': None, 'G': 1.8428972935878388e-05, 'T': 0.00011685195544746245},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [
            0,
            0,
            0,
            0,
            0,
            0,
            1.7697959780874724e-05,
            1.4549694058245949e-05,
            3.5394073056132494e-06,
            3.2783167575131933e-05,
            3.097261734586897e-05,
            0.00017033651501757957,
            0,
            0,
            0,
            0,
            0,
            0,
        ]

        for pos, nt in enumerate(self.sequence2.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        plus_0_codon = self.plus_0_codons[4]    # GTA
        nt = self.sequence2.nt_sequence[12]     # G
        expected = 0
        result = plus_0_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        minus_2_codon = self.minus_2_codons[0]  # ATG
        expected = 2
        result = minus_2_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['C', 'C', 'C']
        exp_mutated = ['C', 'C', 'G']
        res_codon, res_mutated = self.plus_0_codons[3].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        res_codon, res_mutated = self.minus_2_codons[1].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        plus_0_codon = self.plus_0_codons[2]     # AAA
        expected = False                         # Synonymous
        result = plus_0_codon.is_nonsyn(2, 'G')
        self.assertEqual(expected, result)

        expected = True
        result = plus_0_codon.is_nonsyn(2, 'C')  # Lys --> Asn
        self.assertEqual(expected, result)

    def testIsStop(self):
        codon = self.plus_0_codons[5]
        exp = True
        result = codon.is_stop()
        self.assertEqual(exp, result)

        codon = self.minus_2_codons[2]
        exp = False
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        plus_0_codon = self.plus_0_codons[0]    # ATG
        expected = True
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

        plus_0_codon = self.plus_0_codons[4]    # GTA
        expected = False
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

        minus_2_codon = self.minus_2_codons[0]  # ATG
        expected = True
        result = minus_2_codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.minus_2_codons[3]
        exp = True
        result = codon.creates_stop(2, 'G')
        self.assertEqual(exp, result)

        codon = self.plus_0_codons[1]
        exp = False
        result = codon.creates_stop(2, 'A')
        self.assertEqual(exp, result)

    def testGetCodons(self):
        expected = self.plus_0_codons
        result = self.sequence2.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGAATAAACCCGTATGA")]
        for pos, nt in enumerate(nts):
            result = self.sequence2.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGAATAAACCCGTATGA")]
        for pos, nt in enumerate(nts):
            result = self.sequence2.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])


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
        sorted_orfs = {'+0': [{'coords': [[0, 12]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                               'orf_map': np.array([1])}],
                        '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, CAT_VALUES)
        self.seq3_codons = self.sequence3.find_codons('+0', {'coords': [[0, 12]],
                                                             'omega_classes': 3, 'omega_shape': 1.5,
                                                             'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})

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

        # self.assertEqual(self.sequence3.orfs, new_sequence3.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence3.orfs[strand], new_sequence3.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence3.kappa, new_sequence3.kappa)
        self.assertEqual(self.sequence3.global_rate, new_sequence3.global_rate)
        self.assertEqual(self.sequence3.pi, new_sequence3.pi)
        self.assertEqual(self.sequence3.is_circular, new_sequence3.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence3.event_tree, new_sequence3.event_tree)
        self.assertCountEqual(self.sequence3.event_tree, new_sequence3.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence3.nt_sequence):
            new_nt = new_sequence3.nt_sequence[pos]
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
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.25, 'C': 0.08, 'G': 0.42, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence3))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None, 'C': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'G': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'T': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}}}, 'C': {'from_nt': {'A': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'C': None, 'G': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'T': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}}}, 'G': {'from_nt': {'A': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'C': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'G': None, 'T': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}}}, 'T': {'from_nt': {'A': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'C': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'G': {'mu1': {(1,): {}}, 'mu2': {(1,): {}}, 'mu3': {(1,): {}}, 'mu4': {(1,): {}}, 'mu5': {(1,): {}}, 'mu6': {(1,): {}}}, 'T': None}}}}
        result = self.sequence3.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': None, 'C': None, 'G': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {(0.4810288100937172,): {'value': 0.4810288100937172, 'nt_events': 1}, (-1,): {'value': 1, 'nt_events': 1}, (0.1708353283825978,): {'value': 0.1708353283825978, 'nt_events': 1}}
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)

        exp_sub_rates = {'A': 3.1880215087889116e-05, 'C': 1.688969876186309e-05, 'G': None, 'T': 9.564064526366735e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (-1,), 'C': (-1,), 'G': None, 'T': (-1,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu2', 'C': 'mu3', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 5.833397837611894e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a mutation that would destroy a STOP codon
        random.seed(555)
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': None, 'C': None, 'G': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

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
        result = self.sequence3.find_codons('+0', {'coords': [[0, 12]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

    def testIsStartStopCodon(self):
        nt = self.sequence3.nt_sequence[0]  # A in START codon
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[11]  # A in STOP codon (TGA)
        expected = True
        result = self.sequence3.is_start_stop_codon(nt, 'T')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[4]
        expected = False
        result = self.sequence3.is_start_stop_codon(nt, 'C')
        self.assertEqual(expected, result)

        nt = self.sequence3.nt_sequence[8]  # Last G in TGG codon
        expected = True  # Introduces a STOP codon
        result = self.sequence3.is_start_stop_codon(nt, 'A')
        self.assertEqual(expected, result)

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None,
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': None,
                                                'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': None,
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': None}}}}

        result = self.sequence3.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(555)

        a3 = self.sequence3.nt_sequence[3]
        c4 = self.sequence3.nt_sequence[4]
        g5 = self.sequence3.nt_sequence[5]
        t6 = self.sequence3.nt_sequence[6]
        g7 = self.sequence3.nt_sequence[7]
        g8 = self.sequence3.nt_sequence[8]

        exp_event_tree = \
                            {
                                "to_nt": {
                                    "A": {
                                        "from_nt": {
                                            "A": None,
                                            "C": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "G": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {
                                                    (1,): {(-1,): [g5], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.1708353283825978,): [t6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 3,
                                    },
                                    "C": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "C": None,
                                            "G": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.1708353283825978,): [g7],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {
                                                    (1,): {
                                                        (-1,): [g5],
                                                        (0.1708353283825978,): [g8],
                                                        "nt_events": 2,
                                                        "region_weight": 1.1708353283825979,
                                                    },
                                                    "nt_events": 2,
                                                },
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 3,
                                            },
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {
                                                    (1,): {
                                                        (0.1708353283825978,): [t6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 5,
                                    },
                                    "G": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "G": None,
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {
                                                    (1,): {
                                                        (0.1708353283825978,): [t6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 3,
                                    },
                                    "T": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "G": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.1708353283825978,): [g8],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu2": {
                                                    (1,): {
                                                        (0.1708353283825978,): [g7],
                                                        "nt_events": 1,
                                                        "region_weight": 0.1708353283825978,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {
                                                    (1,): {(-1,): [g5], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 3,
                                            },
                                            "T": None,
                                        },
                                        "nt_events": 5,
                                    },
                                }
                            }


        res_omega_key = self.sequence3.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 16
        res_count = self.sequence3.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('C', 4)
        nt.rates = {'A': 0.0584793504496421, 'C': None, 'G': 0.0076247167865047885, 'T': 0.06286076399166855}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.12896483122781544
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': 1.162148311203238e-05, 'G': 7.368793863692968e-05, 'T': 9.327877560764421e-07},
            {'A': 8.763029673862327e-07, 'C': None, 'G': 1.5475107996686342e-06, 'T': 8.054838785442963e-06},
            {'A': 3.1880215087889116e-05, 'C': 4.058819600516357e-05, 'G': None, 'T': 7.72068539216538e-05},
            {'A': 1.717474544205974e-06, 'C': 1.3757733716268578e-05, 'G': 2.6818515078492823e-06, 'T': None},
            {'A': None, 'C': 5.565429903261429e-07, 'G': None, 'T': 1.6338801040342159e-06},
            {'A': None, 'C': 6.933897792999362e-06, 'G': None, 'T': 5.565429903261429e-07},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None}
        ]

        exp_total_rates = [0,
            0,
            0,
            8.624220950503851e-05,
            1.047865255249783e-05,
            0.00014967526501470646,
            1.815705976832383e-05,
            2.1904230943603587e-06,
            7.490440783325505e-06,
            0,
            0,
            0]

        for pos, nt in enumerate(self.sequence3.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq3_codons[2]  # TGG
        nt = self.sequence3.nt_sequence[6]  # T
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['T', 'G', 'A']
        exp_mutated = ['T', 'C', 'A']
        res_codon, res_mutated = self.seq3_codons[3].mutate_codon(1, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq3_codons[0]  # ATG
        expected = True
        result = codon.is_nonsyn(2, 'T')  # ATT
        self.assertEqual(expected, result)

        codon = self.seq3_codons[3]  # TGA
        expected = False
        result = codon.is_nonsyn(1, 'A')
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = False
        codon = self.seq3_codons[0]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = True
        codon = self.seq3_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq3_codons[2]  # TGG
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq3_codons[2]  # TGG
        exp = True
        result = codon.creates_stop(1, 'A')  # Creates TAG
        self.assertEqual(exp, result)

        exp = False
        result = codon.creates_stop(1, 'C')  # Creates TCG
        self.assertEqual(exp, result)

    def testGetCodons(self):
        expected = self.seq3_codons
        result = self.sequence3.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGACGTGGTGA")]
        for pos, nt in enumerate(nts):
            result = self.sequence3.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGACGTGGTGA")]
        for pos, nt in enumerate(nts):
            result = self.sequence3.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])


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
        sorted_orfs = {'+0': [{'coords': [[0, 12]],
                               'omega_classes': 3, 'omega_shape': 1.5, 
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                               'orf_map': np.array([1])}],
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(4000)
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, CAT_VALUES)
        self.seq4_codons = self.sequence4.find_codons('+0', {'coords': [[0, 12]],
                                                             'omega_classes': 3, 'omega_shape': 1.5,
                                                             'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})

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

        # self.assertEqual(self.sequence4.orfs, new_sequence4.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence4.orfs[strand], new_sequence4.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence4.kappa, new_sequence4.kappa)
        self.assertEqual(self.sequence4.global_rate, new_sequence4.global_rate)
        self.assertEqual(self.sequence4.pi, new_sequence4.pi)
        self.assertEqual(self.sequence4.is_circular, new_sequence4.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence4.event_tree, new_sequence4.event_tree)
        self.assertCountEqual(self.sequence4.event_tree, new_sequence4.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence4.nt_sequence):
            new_nt = new_sequence4.nt_sequence[pos]
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
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.33, 'C': 0.25, 'G': 0.17, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence4))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = \
                    {
                        "to_nt": {
                            "A": {
                                "from_nt": {
                                    "A": None,
                                    "C": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "G": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "T": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                }
                            },
                            "C": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "C": None,
                                    "G": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "T": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                }
                            },
                            "G": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "C": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "G": None,
                                    "T": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                }
                            },
                            "T": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "C": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "G": {
                                        "mu1": {(1,): {}},
                                        "mu2": {(1,): {}},
                                        "mu3": {(1,): {}},
                                        "mu4": {(1,): {}},
                                        "mu5": {(1,): {}},
                                        "mu6": {(1,): {}},
                                    },
                                    "T": None,
                                }
                            },
                        }
                    }


        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)

        # Only 3 possible omegas because there is only 1 ORF
        exp_total_omegas = {(0.4810288100937172,): {'value': 0.4810288100937172, 'nt_events': 1}, (-1,): {'value': 1, 'nt_events': 1}}

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence4.total_omegas)

        # Tests mutation in internal methionine
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 1.534035770788274e-05, 'G': 2.1278273495443722e-05, 'T': 2.9180423700224156e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': (0.4810288100937172,), 'G': (0.4810288100937172,), 'T': (0.4810288100937172,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu5', 'G': 'mu3', 'T': 'mu6'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 6.579905490355062e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence4.total_omegas)

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
        result = self.sequence4.find_codons('+0', {'coords': [[0, 12]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
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

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None,
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': None,
                                                                       'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': None,
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': None}}}}

        result = self.sequence4.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(5001)

        a3 = self.sequence4.nt_sequence[3]
        t4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        c6 = self.sequence4.nt_sequence[6]
        c7 = self.sequence4.nt_sequence[7]
        c8 = self.sequence4.nt_sequence[8]

        exp_event_tree = \
                            {
                                "to_nt": {
                                    "A": {
                                        "from_nt": {
                                            "A": None,
                                            "C": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c7],
                                                        (-1,): [c8],
                                                        "nt_events": 2,
                                                        "region_weight": 1.4810288100937172,
                                                    },
                                                    "nt_events": 2,
                                                },
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 3,
                                            },
                                            "G": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g5],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 5,
                                    },
                                    "C": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "C": None,
                                            "G": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g5],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 3,
                                    },
                                    "G": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c6, c7],
                                                        "nt_events": 2,
                                                        "region_weight": 0.9620576201874343,
                                                    },
                                                    "nt_events": 2,
                                                },
                                                "mu4": {
                                                    (1,): {(-1,): [c8], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 3,
                                            },
                                            "G": None,
                                            "T": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 5,
                                    },
                                    "T": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3, a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c7],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {
                                                    (1,): {(-1,): [c8], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 3,
                                            },
                                            "G": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g5],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    "nt_events": 1,
                                                },
                                                "mu2": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu3": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu4": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu5": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "mu6": {(1,): {"nt_events": 0, "region_weight": 0}, "nt_events": 0},
                                                "nt_events": 1,
                                            },
                                            "T": None,
                                        },
                                        "nt_events": 5,
                                    },
                                }
                            }

        res_omega_key = self.sequence4.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 18
        res_count = self.sequence4.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 10)
        nt.rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        nt.get_mutation_rate()
        exp_mutation_rate = 0
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': 3.6147497404682092e-06, 'G': 5.1134525692942473e-05, 'T': 2.9180423700224156e-05},
            {'A': 2.2106381591078907e-05, 'C': 2.5171371204509256e-05, 'G': 4.8359712489644815e-06, 'T': None},
            {'A': 1.7116532419066297e-05, 'C': 1.5032339481933658e-05, 'G': None, 'T': 6.342956741319806e-07},
            {'A': 4.8359712489644815e-06, 'C': None, 'G': 4.8359712489644815e-06, 'T': 1.611990416321494e-05},
            {'A': 9.327877560764421e-07, 'C': None, 'G': 4.8359712489644815e-06, 'T': 3.109292520254807e-06},
            {'A': 1.939151536255615e-06, 'C': None, 'G': 1.569845964087174e-05, 'T': 8.05321349308801e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           8.392969913363484e-05,
                           5.211372404455265e-05,
                           3.278316757513194e-05,
                           2.57918466611439e-05,
                           8.87805152529573e-06,
                           9.816974610800746e-05,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence4.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq4_codons[3]         # TAA
        nt = self.sequence4.nt_sequence[9]  # T
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['A', 'T', 'G']
        exp_mutated = ['T', 'T', 'G']
        res_codon, res_mutated = self.seq4_codons[1].mutate_codon(0, 'T')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq4_codons[1]         # ATG (Met)
        expected = True
        result = codon.is_nonsyn(2, 'A')    # ATA (Ile)
        self.assertEqual(expected, result)

        codon = self.seq4_codons[2]         # CCC (Pro)
        expected = False
        result = codon.is_nonsyn(2, 'A')    # CCA (Pro)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.seq4_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq4_codons[2]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq4_codons[1]     # internal methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

        codon = self.seq4_codons[0]     # First methionine in ORF
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        exp = True
        codon = self.seq4_codons[3]
        result = codon.creates_stop(2, 'G')  # TAA to TAG (STOP to STOP)
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq4_codons[1]
        result = codon.creates_stop(2, 'A')  # Creates ATA
        self.assertEqual(exp, result)

    def testGetCodons(self):
        expected = self.seq4_codons
        result = self.sequence4.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGATGCCCTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence4.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGATGCCCTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence4.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])


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
        sorted_orfs = {'+0': [{'coords': [[0, 12]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                               'orf_map': np.array([1, 0])}],
                        '+1': [{'coords': [[4, 16]],
                                'omega_classes': 4, 'omega_shape': 1.25,
                                'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                 0.5158061369385518, 1.1217646558655263],
                                'orf_map': np.array([0, 1])}],
                        '+2': [], '-0': [], '-1': [], '-2': []}

        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, CAT_VALUES)

        self.plus_0_codons = self.sequence5.find_codons('+0', {'coords': [[0, 12]],
                                                               'omega_classes': 3, 'omega_shape': 1.5,
                                                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.plus_1_codons = self.sequence5.find_codons('+1', {'coords': [[4, 16]],
                                                               'omega_classes': 4, 'omega_shape': 1.25,
                                                               'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                                                0.5158061369385518, 1.1217646558655263]})

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

        # self.assertEqual(self.sequence5.orfs, new_sequence5.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence5.orfs[strand], new_sequence5.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence5.kappa, new_sequence5.kappa)
        self.assertEqual(self.sequence5.global_rate, new_sequence5.global_rate)
        self.assertEqual(self.sequence5.pi, new_sequence5.pi)
        self.assertEqual(self.sequence5.is_circular, new_sequence5.is_circular)

        # Event trees reference different Nucleotides, but the Nucleotides have the same states
        self.assertNotEqual(self.sequence5.event_tree, new_sequence5.event_tree)
        self.assertCountEqual(self.sequence5.event_tree, new_sequence5.event_tree)

        # Check that Nucleotides are different objects with the same attributes
        for pos, nt in enumerate(self.sequence5.nt_sequence):
            new_nt = new_sequence5.nt_sequence[pos]
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
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.38, 'C': 0.19, 'G': 0.19, 'T': 0.25}
        result = Sequence.get_frequency_rates(str(self.sequence5))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = \
                    {
                        "to_nt": {
                            "A": {
                                "from_nt": {
                                    "A": None,
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                }
                            },
                            "C": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "C": None,
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                }
                            },
                            "G": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "G": None,
                                    "T": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                }
                            },
                            "T": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "C": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "G": {
                                        "mu1": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu2": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu3": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu4": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu5": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                        "mu6": {(1, 0): {}, (1, 1): {}, (0, 1): {}},
                                    },
                                    "T": None,
                                }
                            },
                        }
                    }

        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)

        exp_all_omegas = {
                            (0.1708353283825978, None): {"value": 0.1708353283825978, "nt_events": 1},
                            (0.1708353283825978, 0.09199853806558903): {
                                "value": 0.015716600461153824,
                                "nt_events": 1,
                            },
                            (-1, 0.09199853806558903): {"value": 0.09199853806558903, "nt_events": 1},
                            (None, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                            (None, -1): {"value": 1, "nt_events": 1},
                        }


        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        exp_sub_rates = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': None, 'G': None, 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': None, 'C': None, 'G': None}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_all_omegas = {
                            (0.1708353283825978, None): {"value": 0.1708353283825978, "nt_events": 1},
                            (0.1708353283825978, 0.09199853806558903): {
                                "value": 0.015716600461153824,
                                "nt_events": 1,
                            },
                            (-1, 0.09199853806558903): {"value": 0.09199853806558903, "nt_events": 1},
                            (None, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                            (None, -1): {"value": 1, "nt_events": 1},
                        }

        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)
        # Tests a mutation that is synonymous in one frame and non-synonymous in the other
        random.seed(9991)
        nt = self.sequence5.nt_sequence[8]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 1.355833208815998e-07, 'C': None, 'G': 1.355833208815998e-07, 'T': 4.51944402938666e-07}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (-1, 0.09199853806558903), 'C': None, 'G': (-1, 0.09199853806558903), 'T': (-1, 0.09199853806558903)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 7.231110447018656e-07
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_all_omegas = {
                            (0.1708353283825978, None): {"value": 0.1708353283825978, "nt_events": 1},
                            (0.1708353283825978, 0.09199853806558903): {
                                "value": 0.015716600461153824,
                                "nt_events": 1,
                            },
                            (-1, 0.09199853806558903): {"value": 0.09199853806558903, "nt_events": 1},
                            (None, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                            (None, -1): {"value": 1, "nt_events": 1},
                        }

        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests a nucleotide involved in multiple non-synonymous codons
        random.seed(9991)
        nt = self.sequence5.nt_sequence[7]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 2.3162421146011233e-08, 'C': None, 'G': 2.3162421146011233e-08, 'T': 7.72080704867041e-08}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (0.1708353283825978, 0.09199853806558903), 'C': None, 'G': (0.1708353283825978, 0.09199853806558903), 'T': (0.1708353283825978, 0.09199853806558903)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 1.2353291277872657e-07
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_all_omegas = {
                            (0.1708353283825978, None): {"value": 0.1708353283825978, "nt_events": 1},
                            (0.1708353283825978, 0.09199853806558903): {
                                "value": 0.015716600461153824,
                                "nt_events": 1,
                            },
                            (-1, 0.09199853806558903): {"value": 0.09199853806558903, "nt_events": 1},
                            (None, 0.27043066909631136): {"value": 0.27043066909631136, "nt_events": 1},
                            (None, -1): {"value": 1, "nt_events": 1},
                        }

        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

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
        result = self.sequence5.find_codons('+0', {'coords': [[0, 12]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Check Codons in +1 frame
        expected = ['ATG', 'CCT', 'GAC', 'TAA']
        result = self.sequence5.find_codons('+1', {'coords': [[4, 16]],
                                                   'omega_classes': 4, 'omega_shape': 1.25,
                                                   'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                                    0.5158061369385518, 1.1217646558655263]})
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

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None, 'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                      'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                      'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                      'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                      'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                      'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0), (0, 0, 0, 1, 0)): {'prob': 0.003693863037341954, 'number_of_events': 0}, ((0, 0, 0, 1), (0, 1, 0, 0, 0)): {'prob': 1.9117194029041184e-07, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0, 0),): {'prob': 0.07107960659988767, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {((1, 0, 0, 0),): {'prob': 0.02354161198404715, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': None,
                                                'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): {'prob': 0.07107960659988767, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 1, 0, 0), (0, 0, 0, 1, 0)): {'prob': 0.003693863037341954, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((1, 0, 0, 0),): {'prob': 0.02354161198404715, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},

                                                'C': None,
                                                'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0,
                                    'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 0, 1, 0),): {'prob': 0.15821650716270594, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): {'prob': 0.07107960659988767, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 0, 1, 0), (0, 0, 0, 1, 0)): {'prob': 0.003693863037341954, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0, 0),): {'prob': 0.037266143605810084, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': None}}}}

        result = self.sequence5.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(9991)

        # Only 4 possible substitution sites because the
        a3 = self.sequence5.nt_sequence[3]
        c7 = self.sequence5.nt_sequence[7]
        c8 = self.sequence5.nt_sequence[8]
        c12 = self.sequence5.nt_sequence[12]

        exp_event_tree = \
                        {
                            "to_nt": {
                                "A": {
                                    "from_nt": {
                                        "A": None,
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 0.09199853806558903): [c8],
                                                    "nt_events": 1,
                                                    "region_weight": 0.09199853806558903,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {
                                                    (None, 0.27043066909631136): [c12],
                                                    "nt_events": 1,
                                                    "region_weight": 0.27043066909631136,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.09199853806558903): [c7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.015716600461153824,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 3,
                                },
                                "C": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {
                                                    (0.1708353283825978, None): [a3, a3],
                                                    "nt_events": 1,
                                                    "region_weight": 0.1708353283825978,
                                                },
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "nt_events": 1,
                                        },
                                        "C": None,
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 1,
                                },
                                "G": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {
                                                    (0.1708353283825978, None): [a3, a3],
                                                    "nt_events": 1,
                                                    "region_weight": 0.1708353283825978,
                                                },
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 1,
                                        },
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {
                                                    (None, 0.27043066909631136): [c12],
                                                    "nt_events": 1,
                                                    "region_weight": 0.27043066909631136,
                                                },
                                                "nt_events": 1,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.09199853806558903): [c7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.015716600461153824,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 0.09199853806558903): [c8],
                                                    "nt_events": 1,
                                                    "region_weight": 0.09199853806558903,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": None,
                                        "T": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                    },
                                    "nt_events": 4,
                                },
                                "T": {
                                    "from_nt": {
                                        "A": {
                                            "mu1": {
                                                (1, 0): {
                                                    (0.1708353283825978, None): [a3, a3],
                                                    "nt_events": 1,
                                                    "region_weight": 0.1708353283825978,
                                                },
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 1,
                                        },
                                        "C": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (0.1708353283825978, 0.09199853806558903): [c7],
                                                    "nt_events": 1,
                                                    "region_weight": 0.015716600461153824,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {(None, -1): [c12], "nt_events": 1, "region_weight": 1},
                                                "nt_events": 1,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {
                                                    (-1, 0.09199853806558903): [c8],
                                                    "nt_events": 1,
                                                    "region_weight": 0.09199853806558903,
                                                },
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 1,
                                            },
                                            "nt_events": 3,
                                        },
                                        "G": {
                                            "mu1": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu2": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu3": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu4": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu5": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "mu6": {
                                                (1, 0): {"nt_events": 0, "region_weight": 0},
                                                (1, 1): {"nt_events": 0, "region_weight": 0},
                                                (0, 1): {"nt_events": 0, "region_weight": 0},
                                                "nt_events": 0,
                                            },
                                            "nt_events": 0,
                                        },
                                        "T": None,
                                    },
                                    "nt_events": 4,
                                },
                            }
                        }


        res_omega_key = self.sequence5.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 12
        res_count = self.sequence5.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('A', 3)
        nt.rates = {'A': None, 'C': 0.04954407504576764, 'G': 0.3688469654724257, 'T': 0.1106540896417277}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.529045130159921
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': 1.1933500315179366e-05, 'G': 4.927574916928588e-06, 'T': 5.035388960093674e-07},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': 5.489322915008737e-07, 'C': None, 'G': 1.2008391189617825e-07, 'T': 2.2666484428304817e-07},
            {'A': 7.029220070174351e-07, 'C': None, 'G': 1.6892152192547133e-06, 'T': 1.0710748896768025e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': 4.965465882349613e-06, 'C': None, 'G': 1.1700455003277953e-06, 'T': 2.546859337106339e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [
            0,
            0,
            0,
            1.7364614128117323e-05,
            0,
            0,
            0,
            8.956810476801001e-07,
            1.3102886123040173e-05,
            0,
            0,
            0,
            3.1604104753740795e-05,
            0,
            0,
            0,
        ]

        for pos, nt in enumerate(self.sequence5.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        plus_0_codon = self.plus_0_codons[2]    # GCC
        nt = self.sequence5.nt_sequence[7]      # First C
        expected = 1
        result = plus_0_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        plus_1_codon = self.plus_1_codons[1]    # CCT
        expected = 0
        result = plus_1_codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['G', 'C', 'C']
        exp_mutated = ['G', 'C', 'G']
        res_codon, res_mutated = self.plus_0_codons[2].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

        exp_codon = ['G', 'A', 'C']
        exp_mutated = ['G', 'A', 'G']
        res_codon, res_mutated = self.plus_1_codons[2].mutate_codon(2, 'G')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        plus_0_codon = self.plus_0_codons[1]     # AAT (Asn)
        expected = True
        result = plus_0_codon.is_nonsyn(0, 'C')  # CAT (His)
        self.assertEqual(expected, result)

        plus_1_codon = self.plus_1_codons[1]     # CCT (Pro)
        expected = False
        result = plus_1_codon.is_nonsyn(2, 'C')  # CCT (Pro)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.plus_1_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.plus_0_codons[0]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        plus_1_codon = self.plus_1_codons[0]    # ATG
        expected = True
        result = plus_1_codon.is_start()
        self.assertEqual(expected, result)

        plus_0_codon = self.plus_0_codons[3]    # TGA
        expected = False
        result = plus_0_codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.plus_0_codons[3]
        exp = True
        result = codon.creates_stop(1, 'A')
        self.assertEqual(exp, result)

        codon = self.plus_1_codons[2]
        exp = False
        result = codon.creates_stop(0, 'T')
        self.assertEqual(exp, result)

    def testGetCodons(self):
        expected = self.plus_0_codons
        result = self.sequence5.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGAATGCCTGACTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence5.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGAATGCCTGACTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence5.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])
            

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
        sorted_orfs = {'+0': [{'coords': [(0, 5), (6, 13)],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                               'orf_map': np.array([1])}], 
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, CAT_VALUES)
        self.seq6_codons = self.sequence6.find_codons('+0', {'coords': [(0, 5), (6, 13)],
                                                             'omega_classes': 3, 'omega_shape': 1.5,
                                                             'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})

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

        # self.assertEqual(self.sequence6.orfs, new_sequence6.orfs) # This assertion fails if 'orf_map' np.array has more than one element
        strands = ['+0', '+1', '+2', '-0', '-1', '-2']
        for strand in strands:
            for orf, new_orf in zip(self.sequence6.orfs[strand], new_sequence6.orfs[strand]):
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())

        self.assertEqual(self.sequence6.kappa, new_sequence6.kappa)
        self.assertEqual(self.sequence6.global_rate, new_sequence6.global_rate)
        self.assertEqual(self.sequence6.pi, new_sequence6.pi)
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

            self.assertEqual(str(nt.codons), str(new_nt.codons))
            self.assertEqual(len(nt.codons), len(new_nt.codons))

            # Check that Codons are different objects with the same attributes
            for i, codon in enumerate(nt.codons):
                new_codon = new_nt.codons[i]
                self.assertIsNot(codon, new_codon)
                # self.assertEqual(codon.orf, new_codon.orf) # This assertion fails if 'orf_map' np.array has more than one element

                orf, new_orf = codon.orf, new_codon.orf
                self.assertEqual(orf['coords'], new_orf['coords'])
                self.assertEqual(orf['omega_classes'], new_orf['omega_classes'])
                self.assertEqual(orf['omega_shape'], new_orf['omega_shape'])
                self.assertEqual(orf['omega_values'], new_orf['omega_values'])
                self.assertTrue((orf['orf_map'] == new_orf['orf_map']).all())
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.31, 'C': 0.23, 'G': 0.23, 'T': 0.23}
        result = Sequence.get_frequency_rates(str(self.sequence6))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = {
                        "to_nt": {
                            "A": {
                                "from_nt": {
                                    "A": None,
                                    "C": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "G": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "T": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                }
                            },
                            "C": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "C": None,
                                    "G": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "T": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                }
                            },
                            "G": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "C": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "G": None,
                                    "T": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                }
                            },
                            "T": {
                                "from_nt": {
                                    "A": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "C": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "G": {
                                        "mu1": {(1,): {}, (0,): {}},
                                        "mu2": {(1,): {}, (0,): {}},
                                        "mu3": {(1,): {}, (0,): {}},
                                        "mu4": {(1,): {}, (0,): {}},
                                        "mu5": {(1,): {}, (0,): {}},
                                        "mu6": {(1,): {}, (0,): {}},
                                    },
                                    "T": None,
                                }
                            },
                        }
                    }

        result = self.sequence6.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': 4.814194289867335e-05, 'C': 1.4442582869602003e-05, 'G': None, 'T': 1.4442582869602003e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        # No omega keys because nt treated as synonymous mutation
        exp_omega_keys = {'A': (None,), 'C': (None,), 'G': None, 'T': (None,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'C': 'mu4', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 7.702710863787735e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {(0.4810288100937172,): {'value': 0.4810288100937172, 'nt_events': 1}, (None,): {'value': 1, 'nt_events': 1}, (-1,): {'value': 1, 'nt_events': 1}}
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

        # Tests internal methionine
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 2.7411913172937843e-05, 'G': 1.1318913328738838e-05, 'T': 1.4410639058920151e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': (0.4810288100937172,), 'G': (0.4810288100937172,), 'T': (0.4810288100937172,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu6', 'G': 'mu2', 'T': 'mu5'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 5.314146556059683e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

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
        result = self.sequence6.find_codons('+0', {'coords': [(0, 5), (6, 13)],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
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

    @unittest.skip('create_probability_tree is no longer defined')
    def testCreateProbabilityTree(self):
        expected = {'to_nt': {'A': {'number_of_events': 0,
                                    'from_nt': {'A': None,
                                                'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}, ((1, 0, 0, 0),): {'prob': 0.0006395174037506402, 'number_of_events': 0}, ((0, 0, 0, 1),): {'prob': 0.0006395174037506402, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                           'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                           'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                'G': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                             'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                             'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                             'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                             'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                             'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'T': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': None,
                                                                       'C': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'C': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'C': None,
                                                                       'G': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0}}},
                              'G': {'number_of_events': 0, 'from_nt': {'A': {'prob': 0.625, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                    'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                    'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'T': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {((1, 0, 0, 0),): {'prob': 0.061012617279750776, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'C': {'prob': 0.18749999999999997, 'cat': {'mu1': {'prob': 0.01873576363611503, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu2': {'prob': 0.05500382175507591, 'omega': {((0, 0, 0, 1),): {'prob': 1, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu3': {'prob': 0.09713422338727759, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu4': {'prob': 0.15167593856558387, 'omega': {}, 'number_of_events': 0},
                                                                                                                  'mu5': {'prob': 0.23342647807042782, 'omega': {((0, 1, 0, 0),): {'prob': 0.17179600360560732, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                                                                  'mu6': {'prob': 0.4440237745855198, 'omega': {((0, 0, 1, 0),): {'prob': 0.4100485219703123, 'number_of_events': 0}}, 'number_of_events': 0}}, 'number_of_events': 0},
                                                                       'G': None}}}}

        result = self.sequence6.create_probability_tree()
        self.assertEqual(expected, result)

    def testNtInEventTree(self):
        random.seed(4000)

        a3 = self.sequence6.nt_sequence[3]
        t4 = self.sequence6.nt_sequence[4]
        g5 = self.sequence6.nt_sequence[5]
        g6 = self.sequence6.nt_sequence[6]
        c7 = self.sequence6.nt_sequence[7]
        c8 = self.sequence6.nt_sequence[8]
        c9 = self.sequence6.nt_sequence[9]

        exp_event_tree = \
                            {
                                "to_nt": {
                                    "A": {
                                        "from_nt": {
                                            "A": None,
                                            "C": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c7, c8],
                                                        "nt_events": 2,
                                                        "region_weight": 0.9620576201874343,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 2,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {(-1,): [c9], "nt_events": 1, "region_weight": 1},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 3,
                                            },
                                            "G": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {(None,): [g5, g5], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 2,
                                            },
                                            "T": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 6,
                                    },
                                    "C": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 1,
                                            },
                                            "C": None,
                                            "G": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {(None,): [g5, g5], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 2,
                                            },
                                            "T": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 4,
                                    },
                                    "G": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c7],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c8],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu5": {
                                                    (1,): {(-1,): [c9], "nt_events": 1, "region_weight": 1},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 3,
                                            },
                                            "G": None,
                                            "T": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [t4],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 1,
                                            },
                                        },
                                        "nt_events": 5,
                                    },
                                    "T": {
                                        "from_nt": {
                                            "A": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {
                                                        (0.4810288100937172,): [a3],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "nt_events": 1,
                                            },
                                            "C": {
                                                "mu1": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c7],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {(-1,): [c9], "nt_events": 1, "region_weight": 1},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {
                                                        (0.4810288100937172,): [c8],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 3,
                                            },
                                            "G": {
                                                "mu1": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {(None,): [g5, g5], "nt_events": 1, "region_weight": 1},
                                                    "nt_events": 1,
                                                },
                                                "mu2": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu3": {
                                                    (1,): {
                                                        (0.4810288100937172,): [g6],
                                                        "nt_events": 1,
                                                        "region_weight": 0.4810288100937172,
                                                    },
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 1,
                                                },
                                                "mu4": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu5": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "mu6": {
                                                    (1,): {"nt_events": 0, "region_weight": 0},
                                                    (0,): {"nt_events": 0, "region_weight": 0},
                                                    "nt_events": 0,
                                                },
                                                "nt_events": 2,
                                            },
                                            "T": None,
                                        },
                                        "nt_events": 6,
                                    },
                                }
                            }

        res_omega_key = self.sequence6.nt_in_event_tree(g5)
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

        # No new omega keys are created initially
        self.assertEqual(None, res_omega_key)

    @unittest.skip('count_nts_on_event_tree is no longer defined')
    def testCountNtsOnEventTree(self):
        exp_count = 18
        res_count = self.sequence6.count_nts_on_event_tree()
        self.assertEqual(exp_count, res_count)

    def testGetMutationRate(self):
        nt = Nucleotide('C', 10)
        nt.rates = {'A': 0.25667881246211854, 'C': None, 'G': 0.02735873671412473, 'T': 0.2491527517379426}
        nt.get_mutation_rate()
        exp_mutation_rate = 0.5331903009141858
        self.assertEqual(exp_mutation_rate, nt.mutation_rate)

    def testCheckMutationRates(self):
        exp_sub_rates = [
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': 3.395673998621651e-06, 'G': 4.803546352973384e-05, 'T': 2.7411913172937843e-05},
            {'A': 2.0337871063792594e-05, 'C': 2.315766150814852e-05, 'G': 4.449093549047323e-06, 'T': None},
            {'A': 4.814194289867335e-05, 'C': 4.227994381423898e-05, 'G': None, 'T': 1.7840194133551658e-06},
            {'A': 1.4830311830157746e-05, 'C': 4.449093549047323e-06, 'G': None, 'T': 4.449093549047323e-06},
            {'A': 8.581647355903268e-07, 'C': None, 'G': 4.449093549047323e-06, 'T': 2.8605491186344227e-06},
            {'A': 8.581647355903268e-07, 'C': None, 'G': 6.9472984524445545e-06, 'T': 3.5639214876899296e-05},
            {'A': 4.227994381423898e-05, 'C': None, 'G': 2.2226869240922907e-05, 'T': 3.083040250181358e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           7.884305070129333e-05,
                           4.794462612098844e-05,
                           9.22059061262675e-05,
                           2.3728498928252395e-05,
                           8.167807403272074e-06,
                           4.344467806493418e-05,
                           9.533721555697547e-05,
                           0,
                           0,
                           0]

        for pos, nt in enumerate(self.sequence6.nt_sequence):
            self.assertEqual(exp_sub_rates[pos], nt.rates)
            self.assertEqual(exp_total_rates[pos], nt.mutation_rate)

    def testNtInPos(self):
        codon = self.seq6_codons[0]         # ATG
        nt = self.sequence6.nt_sequence[0]  # A
        expected = 0
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

        codon = self.seq6_codons[3]          # TAA
        nt = self.sequence6.nt_sequence[12]  # A
        expected = 2
        result = codon.nt_in_pos(nt)
        self.assertEqual(expected, result)

    def testMutateCodon(self):
        exp_codon = ['A', 'T', 'G']
        exp_mutated = ['C', 'T', 'G']
        res_codon, res_mutated = self.seq6_codons[1].mutate_codon(0, 'C')
        self.assertEqual(exp_codon, res_codon)
        self.assertEqual(exp_mutated, res_mutated)

    def testIsNonSyn(self):
        codon = self.seq6_codons[3]             # TAA (STOP)
        expected = True
        result = codon.is_nonsyn(1, 'C')        # TAC (Tyr)
        self.assertEqual(expected, result)

        expected = False
        result = codon.is_nonsyn(1, 'G')        # TGA (STOP)
        self.assertEqual(expected, result)

    def testIsStop(self):
        exp = True
        codon = self.seq6_codons[3]
        result = codon.is_stop()
        self.assertEqual(exp, result)

        exp = False
        codon = self.seq6_codons[2]
        result = codon.is_stop()
        self.assertEqual(exp, result)

    def testIsStart(self):
        codon = self.seq6_codons[0]     # First methionine
        expected = True
        result = codon.is_start()
        self.assertEqual(expected, result)

        codon = self.seq6_codons[1]     # Internal methionine
        expected = False
        result = codon.is_start()
        self.assertEqual(expected, result)

    def testCreatesStop(self):
        codon = self.seq6_codons[3]
        exp = True
        result = codon.creates_stop(2, 'G')
        self.assertEqual(exp, result)

        exp = False
        result = codon.creates_stop(0, 'G')
        self.assertEqual(exp, result)

    def testGetCodons(self):
        expected = self.seq6_codons
        result = self.sequence6.get_codons()
        self.assertEqual(expected, result)

    def testGetRightNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGATGGCCCTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence6.get_right_nt(pos)
            self.assertEqual(result, nts[pos + 1])

    def testGetLeftNT(self):
        nts = [Nucleotide(nt, pos) for pos, nt in enumerate("ATGATGGCCCTAA")]
        for pos, nt in enumerate(nts):
            result = self.sequence6.get_left_nt(pos)
            self.assertEqual(result, nts[pos - 1])


def main(out = sys.stderr, verbosity = 2):
    loader = unittest.TestLoader()
  
    suite = loader.loadTestsFromModule(sys.modules[__name__])
    unittest.TextTestRunner(out, verbosity = verbosity).run(suite)

if __name__ == '__main__':
    with open('testing.out', 'w') as f:
        main(f)
