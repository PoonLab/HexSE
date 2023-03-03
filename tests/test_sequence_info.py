import random
import unittest
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
                        'dn_values': [0.42584203488769556, 1.0711311227395655, 1.7848172815920647,
                             	      2.780153100609863, 5.1880564601470684],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([1])}], 
                '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

        s1 = 'GTACGATCGATCGATGCTAGC'
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, orfs, KAPPA, GLOBAL_RATE, pi1, CAT_VALUES)
        self.circular_sequence1 = Sequence(s1, orfs, KAPPA, GLOBAL_RATE, pi1, CAT_VALUES, circular=True)
        self.orf = {'coords': [[0, 21]],
                    'omega_classes': 3, 'omega_shape': 1.5,
                    'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                    'dn_values': [0.42584203488769556, 1.0711311227395655, 1.7848172815920647,
                                  2.780153100609863, 5.1880564601470684],
                    'ds_values': [0.6137056388801096, 3.386294361119891],
                    'orf_map': np.array([1])}
        self.seq1_codons = self.sequence1.find_codons('+0', self.orf)

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

        exp_sub_rates = {'A': 9.42914477978867e-07,
                         'C': 1.4665417538020733e-06,
                         'G': None,
                         'T': 1.4665417538020733e-06}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0.12575458287886826,), 'C': (0.12575458287886826,), 'G': None,
                          'T': (0.12575458287886826,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'A': 'mu1',
                        'C': 'mu3',
                        'T': 'mu3'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 3.875997985583014e-06
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {(0.12575458287886826,): {'value': 0.12575458287886826, 'nt_events': 1},
                            (3.386294361119891,): {'value': 3.386294361119891, 'nt_events': 1},
                            (0.6137056388801096,): {'value': 0.6137056388801096, 'nt_events': 1},
                            (2.9082628030744515,): {'value': 2.9082628030744515, 'nt_events': 1},
                            (1.532074860270361,): {'value': 1.532074860270361, 'nt_events': 1},
                            (8.45365616912733,): {'value': 8.45365616912733, 'nt_events': 1},
                            (0.5270709191984724,): {'value': 0.5270709191984724, 'nt_events': 1}}
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
        result = self.sequence1.find_codons('+0', self.orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence1.get_codons()
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.532074860270361,
                                                  (1.532074860270361,): [c20]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [c3]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [c11]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 2.9082628030744515,
                                                  (2.9082628030744515,): [c16]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 4},
                             'G': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [g8]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 2.9082628030744515,
                                                  (2.9082628030744515,): [g4]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 2.0591457794688335,
                                                  (0.5270709191984724,): [g12],
                                                  (1.532074860270361,): [g19]}},
                                   'mu6': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 3.03401738595332,
                                                  (0.12575458287886826,): [g0,
                                                                           g0],
                                                  (2.9082628030744515,): [g15]}},
                                   'nt_events': 6},
                             'T': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [t14]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [t10]}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [t17]}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.12575458287886826,
                                                  (0.12575458287886826,): [t1]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.532074860270361,
                                                  (1.532074860270361,): [t6]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 5}},
                 'nt_events': 15},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.532074860270361,
                                                  (1.532074860270361,): [a18]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [a2]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [a13]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [a9]}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [a5]}},
                                   'nt_events': 5},
                             'C': None,
                             'G': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 3,
                                           (1,): {'nt_events': 3,
                                                  'region_weight': 5.942280189027771,
                                                  (0.12575458287886826,): [g0,
                                                                           g0],
                                                  (2.9082628030744515,): [g4,
                                                                          g15]}},
                                   'mu3': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 2.0591457794688335,
                                                  (0.5270709191984724,): [g12],
                                                  (1.532074860270361,): [g19]}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [g8]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 6},
                             'T': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 2.1457804991504705,
                                                  (0.6137056388801096,): [t17],
                                                  (1.532074860270361,): [t6]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [t14]}},
                                   'mu6': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 8.579410752006197,
                                                  (0.12575458287886826,): [t1],
                                                  (8.45365616912733,): [t10]}},
                                   'nt_events': 5}},
                 'nt_events': 16},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [a9]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [a2]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [a5]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [a13]}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.532074860270361,
                                                  (1.532074860270361,): [a18]}},
                                   'nt_events': 5},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 3,
                                           (1,): {'nt_events': 3,
                                                  'region_weight': 5.972412523615174,
                                                  (1.532074860270361,): [c7,
                                                                         c20],
                                                  (2.9082628030744515,): [c3]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [c11]}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 2.9082628030744515,
                                                  (2.9082628030744515,): [c16]}},
                                   'nt_events': 5},
                             'G': None,
                             'T': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [t10]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.532074860270361,
                                                  (1.532074860270361,): [t6]}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [t17]}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.12575458287886826,
                                                  (0.12575458287886826,): [t1]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [t14]}},
                                   'nt_events': 5}},
                 'nt_events': 15},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [a13]}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [a2]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [a5]}},
                                   'mu6': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 9.98573102939769,
                                                  (1.532074860270361,): [a18],
                                                  (8.45365616912733,): [a9]}},
                                   'nt_events': 5},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 4.918369221390252,
                                                  (1.532074860270361,): [c7],
                                                  (3.386294361119891,): [c20]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 3.521968441954561,
                                                  (0.6137056388801096,): [c11],
                                                  (2.9082628030744515,): [c16]}},
                                   'nt_events': 4},
                             'G': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.5270709191984724,
                                                  (0.5270709191984724,): [g12]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 2.9082628030744515,
                                                  (2.9082628030744515,): [g15]}},
                                   'mu4': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 4.440337663344812,
                                                  (1.532074860270361,): [g19],
                                                  (2.9082628030744515,): [g4]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.12575458287886826,
                                                  (0.12575458287886826,): [g0,
                                                                           g0]}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [g8]}},
                                   'nt_events': 6},
                             'T': None},
                 'nt_events': 15}}}

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
                {'A': 2.2346377428485024e-05, 'C': 8.304529383108346e-07, 'G': None, 'T': 3.524295192728992e-06},
                {'A': 1.895187114219431e-06, 'C': 1.8493553733918643e-05, 'G': 1.895187114219431e-06, 'T': None},
                {'A': None, 'C': 1.8506707300012832e-05, 'G': 6.168902433337609e-05, 'T': 5.103322114574608e-05},
                {'A': 3.354011617396852e-06, 'C': None, 'G': 2.806839009386231e-05, 'T': None},
                {'A': 6.401822064594056e-05, 'C': 1.920546619378217e-05, 'G': None, 'T': 5.296008563737085e-05},
                {'A': None, 'C': 2.707558949252105e-05, 'G': 3.082954625067139e-05, 'T': 1.423383129162287e-05},
                {'A': 3.5533802698992575e-05, 'C': 4.9288157848742635e-05, 'G': 8.37306447020002e-06, 'T': None},
                {'A': None, 'C': None, 'G': 1.4786447354622791e-05, 'T': 7.696388915439513e-05},
                {'A': 2.5390613261973534e-05, 'C': 6.166514221777651e-05, 'G': None, 'T': 0.00018052164019443805},
                {'A': None, 'C': 0.00019606780203017157, 'G': 5.2457345111488055e-05, 'T': 0.0003729601125449721},
                {'A': 4.620075033442972e-05, 'C': 0.0012432003751499069, 'G': 1.5737203533446416e-05, 'T': None},
                {'A': 9.248863875201417e-06, 'C': None, 'G': 0.00019606780203017157, 'T': 9.025196497507017e-05},
                {'A': 4.9237529804319885e-05, 'C': 6.146666725966998e-06, 'G': None, 'T': 1.1856016439080772e-06},
                {'A': None, 'C': 7.943233490797814e-06, 'G': 4.074830052771301e-05, 'T': 5.086896600800274e-06},
                {'A': 9.811875673722017e-07, 'C': 0.00026179729382917006, 'G': 2.3253421408519e-05, 'T': None},
                {'A': 0.0005167933984666436, 'C': 1.920546619378217e-05, 'G': None, 'T': 3.391597136341696e-05},
                {'A': 6.745208039183558e-05, 'C': None, 'G': 0.00012830732651585632, 'T': 0.00042769108838618775},
                {'A': 5.923030496652448e-06, 'C': 1.9743434988841493e-05, 'G': 5.923030496652448e-06, 'T': None},
                {'A': None, 'C': 2.8520883061938023e-06, 'G': 0.00022530796177891793, 'T': 6.759238853367538e-05},
                {'A': 0.0001431222608709423, 'C': 1.7866957220169205e-05, 'G': None, 'T': 2.7899409818468235e-05},
                {'A': 2.8520883061938023e-06, 'C': None, 'G': 1.4786447354622791e-05, 'T': 0.00017011073715248696}
            ]

        exp_total_rates = [2.6701125559524848e-05,
                           2.2283927962357505e-05,
                           0.000131228952779135,
                           3.142240171125916e-05,
                           0.00013618377247709358,
                           7.213896703481531e-05,
                           9.319502501793524e-05,
                           9.175033650901792e-05,
                           0.0002675773956741881,
                           0.0006214852596866318,
                           0.001305138329017783,
                           0.00029556863088044317,
                           5.656979817419496e-05,
                           5.3778430619311096e-05,
                           0.0002860319028050613,
                           0.0005699148360238426,
                           0.0006234504952938797,
                           3.158949598214639e-05,
                           0.0002957524386187871,
                           0.00018888862790957976,
                           0.00018774927281330354,]

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

    def testGetComplementState(self):
        nts = self.sequence1.nt_sequence
        self.assertEqual(nts, self.sequence1.get_sequence()) # testing get_sequence() also
        complements = "CATGCTAGCTAGCTACGATCG"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)
        
    def testGetRightLeftNT(self):
        seq1 = ''.join(nt.state for nt in self.sequence1.nt_sequence)

        for pos, nt in enumerate(seq1[:-1]):
            result = self.sequence1.get_right_nt(pos)
            self.assertEqual(result.state, seq1[pos + 1])

        for pos, nt in enumerate(seq1[1:], 1):
            result = self.sequence1.get_left_nt(pos)
            self.assertEqual(result.state, seq1[pos - 1])

        circular_seq1 = ''.join(nt.state for nt in self.circular_sequence1.nt_sequence)

        for pos, nt in enumerate(circular_seq1):
            result = self.circular_sequence1.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence1[pos + 1])

        for pos, nt in enumerate(circular_seq1):
            result = self.circular_sequence1.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq1[pos - 1])

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
                        'dn_values': [0.13695378264465718, 0.4767518562354524,
                                      0.9999999999999997, 2.3862943611198904],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([1, 0])}],
                '+1': [], '+2': [], '-0': [], '-1': [],
                '-2': [{'coords': [[3, 15]],
                        'omega_classes': 4, 'omega_shape': 1.25,
                        'omega_values': [0.13695378264465718, 0.4767518562354524,
                                         0.9999999999999997, 2.3862943611198904], 
                        'dn_values': [1.8965767845633317, 6.103423215436677],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([0, 1])}]}

        s2 = 'ATGAATAAACCCGTATGA'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, orfs, KAPPA, GLOBAL_RATE, pi2, CAT_VALUES)
        self.circular_sequence2 = Sequence(s2, orfs, KAPPA, GLOBAL_RATE, pi2, CAT_VALUES, circular=True)
        plus_orf = {'coords': [[0, 21]],
                        'omega_classes': 3, 'omega_shape': 1.5,
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                        'dn_values': [0.13695378264465718, 0.4767518562354524,
                                      0.9999999999999997, 2.3862943611198904],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([1, 0])}
        minus_orf = {'coords': [[3, 15]],
                        'omega_classes': 4, 'omega_shape': 1.25,
                        'omega_values': [0.13695378264465718, 0.4767518562354524,
                                      0.9999999999999997, 2.3862943611198904], 
                        'dn_values': [1.8965767845633317, 6.103423215436677],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([0, 1])}
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
                   {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'G': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'T': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}}}},
                              'C': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'C': None,
                                                'G': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'T': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}}}},
                              'G': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                                (1, 1): {}}},
                                                'C': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'G': None,
                                                'T': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}}}},
                              'T': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'C': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'G': {'mu1': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu2': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu3': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu4': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu5': {(1, 0): {},
                                                              (1, 1): {}},
                                                      'mu6': {(1, 0): {},
                                                              (1, 1): {}}},
                                                'T': None}}}}
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

        exp_sub_rates = {'A': 3.6692034228945255e-06, 'C': None, 'G': 1.0741410727129097e-05, 'T': 4.435337886915156e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (0.6137056388801096, 0.5600744006011661), 'C': None, 'G': (0.6137056388801096, 0.5600744006011661), 'T': (0.6137056388801096, 0.5600744006011661)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'G': 'mu6', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 1.884595203693878e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_total_omegas = {(0.7046919454251792, 3.090368842013911): {'value': 2.177758031360141,
                                                                      'nt_events': 1},
                            (0.7046919454251792, 0.6137056388801096): {'value': 0.432473420580827,
                                                                       'nt_events': 1},
                            (3.386294361119891, 3.090368842013911): {'value': 10.464898583492312,
                                                                     'nt_events': 1},
                            (3.888337029906392, 3.386294361119891): {'value': 13.167053758505679,
                                                                     'nt_events': 1},
                            (3.888337029906392, 0.5600744006011661): {'value': 2.177758031360141,
                                                                      'nt_events': 1},
                            (0.6137056388801096, 0.5600744006011661): {'value': 0.3437208178413331,
                                                                       'nt_events': 1}}

        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation would destroy START CODON
        # random.seed(4001)
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

        exp_sub_rates = {'A': 2.32474636627662e-05, 'C': None, 'G': 6.80557948920672e-05, 'T': 2.8101564419889702e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (3.888337029906392, 0.5600744006011661), 'C': None, 'G': (3.888337029906392, 0.5600744006011661), 'T': (3.888337029906392, 0.5600744006011661)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'G': 'mu6', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0.0001194048229747231
        self.assertEqual(exp_total_rate, nt.mutation_rate)
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
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                        'dn_values': [0.13695378264465718, 0.4767518562354524,
                                      0.9999999999999997, 2.3862943611198904],
                        'ds_values': [0.6137056388801096, 3.386294361119891],
                        'orf_map': np.array([1, 0])})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence2.get_codons()[:6]
        self.assertEqual(len(expected), len(result))
        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Test reverse strand ORF
        expected = ['ATG', 'CCC', 'AAA', 'TAA']
        result = self.sequence2.find_codons('-2', {'coords': [[3, 15]],
                                            'omega_classes': 4, 'omega_shape': 1.25,
                                            'omega_values': [0.13695378264465718, 0.4767518562354524,
                                                        0.9999999999999997, 2.3862943611198904], 
                                            'dn_values': [1.8965767845633317, 6.103423215436677],
                                            'ds_values': [0.6137056388801096, 3.386294361119891],
                                            'orf_map': np.array([0, 1])})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '-2')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence2.get_codons()[6:]
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
                   {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.3437208178413331,
                                                    (0.6137056388801096, 0.5600744006011661): [c11,
                                                                                               c11]}},
                                   'mu3': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 13.167053758505679,
                                                    (3.888337029906392, 3.386294361119891): [c9]}},
                                   'mu4': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 2.177758031360141,
                                                    (3.888337029906392, 0.5600744006011661): [c10]}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 3},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 2,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 2,
                                                    'region_weight': 4.355516062720282,
                                                    (0.7046919454251792, 3.090368842013911): [a7,
                                                                                              a8]}},
                                   'mu2': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 2.177758031360141,
                                                    (0.7046919454251792, 3.090368842013911): [a6]}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'C': None,
                             'G': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 3},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 2.177758031360141,
                                                    (0.7046919454251792, 3.090368842013911): [a7]}},
                                   'mu3': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 10.464898583492312,
                                                    (3.386294361119891, 3.090368842013911): [a8]}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.432473420580827,
                                                    (0.7046919454251792, 0.6137056388801096): [a6]}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'C': {'mu1': {'nt_events': 2,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 2,
                                                    'region_weight': 15.344811789865819,
                                                    (3.888337029906392, 0.5600744006011661): [c10],
                                                    (3.888337029906392, 3.386294361119891): [c9]}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.3437208178413331,
                                                    (0.6137056388801096, 0.5600744006011661): [c11,
                                                                                               c11]}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'G': None,
                             'T': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 6},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 2.177758031360141,
                                                    (0.7046919454251792, 3.090368842013911): [a7]}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 13.167053758505679,
                                                    (3.888337029906392, 3.386294361119891): [c9]}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 2.177758031360141,
                                                    (3.888337029906392, 0.5600744006011661): [c10]}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.3437208178413331,
                                                    (0.6137056388801096, 0.5600744006011661): [c11,
                                                                                               c11]}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': None},
                 'nt_events': 4}}}
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
            {'A': None, 'C': 2.182003825544377e-05, 'G': 6.129729383401259e-05, 'T': None},
            {'A': None, 'C': 7.432484984504838e-06, 'G': 7.273346085147923e-05, 'T': 3.8533185560467476e-05},
            {'A': None, 'C': 7.432484984504838e-06, 'G': 0.0006172186765507815, 'T': None},
            {'A': 9.001401706097875e-05, 'C': None, 'G': 1.7362380516161258e-05, 'T': 5.787460172053753e-05},
            {'A': 2.32474636627662e-05, 'C': None, 'G': 2.871641925831414e-06, 'T': 4.962607231272327e-05},
            {'A': 1.3306013660745468e-06, 'C': None, 'G': 5.646836541313942e-06, 'T': 1.8822788471046477e-05},
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
            8.311733208945635e-05,
            0.00011869913139645155,
            0.0006246511615352864,
            0.00016525099929767752,
            7.574517790132089e-05,
            2.5800226378434966e-05,
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

    def testGetComplementState(self):
        nts = self.sequence2.nt_sequence
        self.assertEqual(nts, self.sequence2.get_sequence()) # testing get_sequence() also
        complements = "TACTTATTTGGGCATACT"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)

    def testGetRightLeftNT(self):
        seq2 = ''.join(nt.state for nt in self.sequence2.nt_sequence)
        for pos, nt in enumerate(seq2[:-1]):
            result = self.sequence2.get_right_nt(pos)
            self.assertEqual(result.state, seq2[pos + 1])
        
        for pos, nt in enumerate(seq2[1:], 1):
            result = self.sequence2.get_left_nt(pos)
            self.assertEqual(result.state, seq2[pos - 1])

        circular_seq2 = ''.join(nt.state for nt in self.circular_sequence2.nt_sequence)

        for pos, nt in enumerate(circular_seq2):
            result = self.circular_sequence2.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence2[pos + 1])

        for pos, nt in enumerate(circular_seq2):
            result = self.circular_sequence2.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq2[pos - 1])

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
                               'dn_values': [0.42584203488769556, 1.0711311227395655, 1.7848172815920647,
                                             2.780153100609863, 5.1880564601470684],
                               'ds_values': [0.6137056388801096, 3.386294361119891],
                               'orf_map': np.array([1])}],
                        '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, CAT_VALUES)
        self.circular_sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, CAT_VALUES, circular=True)
        self.seq3_codons = self.sequence3.find_codons('+0', {'coords': [[0, 12]],
                                                             'omega_classes': 3, 'omega_shape': 1.5,
                                                             'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                                                             'dn_values': [0.42584203488769556,
                                                                           1.0711311227395655,
                                                                           1.7848172815920647,
                                                                           2.780153100609863,
                                                                           5.1880564601470684],
                                                             'ds_values': [0.6137056388801096, 3.386294361119891],
                                                             'orf_map': np.array([1])})

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
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'C': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': None, 'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'G': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': None, 'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'T': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': None}}}}
        
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

        exp_total_omegas = {(1.7453499770576755,): {'value': 1.7453499770576755, 'nt_events': 1}, (0.6137056388801096,): {'value': 0.6137056388801096, 'nt_events': 1}, (8.45365616912733,): {'value': 8.45365616912733, 'nt_events': 1}}
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)
        exp_sub_rates = {'A': 1.95650677681483e-05, 'C': 1.0365303369141784e-05, 'G': None, 'T': 5.86952033044449e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': (0.6137056388801096,), 'C': (0.6137056388801096,), 'G': None, 'T': (0.6137056388801096,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu2', 'C': 'mu3', 'T': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 3.579989146773457e-05
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
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404],
                                                   'dn_values': [0.42584203488769556,
                                                                 1.0711311227395655,
                                                                 1.7848172815920647,
                                                                 2.780153100609863,
                                                                 5.1880564601470684],
                                                   'ds_values': [0.6137056388801096, 3.386294361119891],
                                                   'orf_map': np.array([1])})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence3.get_codons()
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
                   {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [c4]}},
                                   'nt_events': 1},
                             'G': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [g5]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'T': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [t6]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 3},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [a3,
                                                                          a3]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': None,
                             'G': {'mu1': {'nt_events': 2,
                                           (1,): {'nt_events': 2,
                                                  'region_weight': 9.067361808007439,
                                                  (0.6137056388801096,): [g5],
                                                  (8.45365616912733,): [g8]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [g7]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'T': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [t6]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 5},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [a3,
                                                                          a3]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [c4]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'G': None,
                             'T': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [t6]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 3},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [a3,
                                                                          a3]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.7453499770576755,
                                                  (1.7453499770576755,): [c4]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'G': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [g7]}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.6137056388801096,
                                                  (0.6137056388801096,): [g5]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 8.45365616912733,
                                                  (8.45365616912733,): [g8]}},
                                   'nt_events': 3},
                             'T': None},
                 'nt_events': 5}}}
        
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
            {'A': None, 'C': 2.739930617403634e-05, 'G': 3.312031706681014e-05, 'T': 4.216702795620516e-05},
            {'A': 2.566723442301976e-05, 'C': None, 'G': 5.614940066005484e-06, 'T': 4.4978163153285514e-05},
            {'A': 5.395170593867494e-05, 'C': 1.999314630504405e-06, 'G': None, 'T': 5.86952033044449e-06},
            {'A': 0.00020423729378142873, 'C': 5.4643067824466724e-05, 'G': 0.0001327093801888518, 'T': None},
            {'A': None, 'C': 0.00022295175871727107, 'G': None, 'T': 2.754010618353123e-05},
            {'A': None, 'C': 2.754010618353123e-05, 'G': None, 'T': 0.0006526801969537012},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [
            0,
            0,
            0,
            0.00010268665119705163,
            7.626033764231076e-05,
            6.182054089962384e-05,
            0.0003915897417947473,
            0.0002504918649008023,
            0.0006802203031372324,
            0,
            0,
            0,
        ]

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

    def testGetComplementState(self):
        nts = self.sequence3.nt_sequence
        self.assertEqual(nts, self.sequence3.get_sequence()) # testing get_sequence() also
        complements = "TACTGCACCACT"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)

    def testGetRightLeftNT(self):
        seq3 = ''.join(nt.state for nt in self.sequence3.nt_sequence)
        for pos, nt in enumerate(seq3[:-1]):
            result = self.sequence3.get_right_nt(pos)
            self.assertEqual(result.state, seq3[pos + 1])

        for pos, nt in enumerate(seq3[1:], 1):
            result = self.sequence3.get_left_nt(pos)
            self.assertEqual(result.state, seq3[pos - 1])


        circular_seq3 = ''.join(nt.state for nt in self.circular_sequence3.nt_sequence)

        for pos, nt in enumerate(circular_seq3):
            result = self.circular_sequence3.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence3[pos + 1])

        for pos, nt in enumerate(circular_seq3):
            result = self.circular_sequence3.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq3[pos - 1])

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
                               'dn_values': [0.3329677267246186, 1.0887942245237032, 2.8982380487141928],
                               'ds_values': [0.6137056388801096, 3.386294361119891],
                               'orf_map': np.array([1])}],
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(4000)
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, CAT_VALUES)
        self.circular_sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, CAT_VALUES, circular=True)
        self.plus_0_orf = sorted_orfs['+0'][0]
        self.seq4_codons = self.sequence4.find_codons('+0', self.plus_0_orf)

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
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'C': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': None,
                                                'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'G': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': None,
                                                'T': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}}}},
                                'T': {'from_nt': {'A': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'C': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'G': {'mu1': {(1,): {}},
                                                        'mu2': {(1,): {}},
                                                        'mu3': {(1,): {}},
                                                        'mu4': {(1,): {}},
                                                        'mu5': {(1,): {}},
                                                        'mu6': {(1,): {}}},
                                                'T': None}}}}

        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)

        # Only 3 possible omegas because there is only 1 ORF
        exp_total_omegas = {(1.774131041895811,): {'value': 1.774131041895811, 'nt_events': 1}, (0.8558730398605124,): {'value': 0.8558730398605124, 'nt_events': 1}, (3.386294361119891,): {'value': 3.386294361119891, 'nt_events': 1}}

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

        exp_sub_rates = {'A': None, 'C': 5.657832593860247e-05, 'G': 7.847855416157049e-05, 'T': 0.0001076232741489096}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (1.774131041895811,), 'G': (1.774131041895811,), 'T': (1.774131041895811,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu5', 'G': 'mu3', 'T': 'mu6'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.00024268015424908256
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
        result = self.sequence4.find_codons('+0', self.plus_0_orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence4.get_codons()
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
            {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c8]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c6]}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g5]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'T': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 5},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3,
                                                                         a3]}},
                                   'nt_events': 1},
                             'C': None,
                             'G': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g5]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'T': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 3},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3,
                                                                         a3]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c8]}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c6]}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'G': None,
                             'T': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 5},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3,
                                                                         a3]}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c8]}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c6]}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 1,
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g5]}},
                                   'mu2': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'T': None},
                 'nt_events': 5}}}
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
            {'A': None, 'C': 0.0001076232741489096, 'G': 1.5137359315808562e-05, 'T': 2.3543566248471145e-05},
            {'A': 1.7836035036720565e-05, 'C': 5.945345012240189e-05, 'G': 3.4403089354110365e-06, 'T': None},
            {'A': 4.0428346083233286e-05, 'C': 2.3394100760795048e-06, 'G': None, 'T': 2.3394100760795048e-06},
            {'A': 1.3435888373960466e-05, 'C': None, 'G': 2.0677584938924793e-05, 'T': 0.00013110965230950865},
            {'A': 2.0677584938924793e-05, 'C': None, 'G': 8.604427274788151e-06, 'T': 1.6241319403975378e-05},
            {'A': 6.566537912579363e-06, 'C': None, 'G': 1.927782010418003e-05, 'T': 0.00017719868453384057},
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None}, 
        ]

        exp_total_rates = [       
            0,
            0,
            0,
            0.0001463041997131893,
            8.072979409453349e-05,
            4.5107166235392293e-05,
            0.0001652231256223939,
            4.5523331617688324e-05,
            0.00020304304255059995,
            0,
            0,
            0,]

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

    def testGetComplementState(self):
        nts = self.sequence4.nt_sequence
        self.assertEqual(nts, self.sequence4.get_sequence()) # testing get_sequence() also
        complements = "TACTACGGGATT"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)

    def testGetRightLeftNT(self):
        seq4 = ''.join(nt.state for nt in self.sequence4.nt_sequence)
        for pos, nt in enumerate(seq4[:-1]):
            result = self.sequence4.get_right_nt(pos)
            self.assertEqual(result.state, seq4[pos + 1])

        for pos, nt in enumerate(seq4[1:], 1):
            result = self.sequence4.get_left_nt(pos)
            self.assertEqual(result.state, seq4[pos - 1])

        circular_seq4 = ''.join(nt.state for nt in self.circular_sequence4.nt_sequence)

        for pos, nt in enumerate(circular_seq4):
            result = self.circular_sequence4.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence4[pos + 1])

        for pos, nt in enumerate(circular_seq4):
            result = self.circular_sequence4.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq4[pos - 1])

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
                               'dn_values': [0.13695378264465718, 0.4767518562354524,
                                             0.9999999999999997, 2.3862943611198904],
                               'ds_values': [0.6137056388801096, 3.386294361119891],
                               'orf_map': np.array([1, 0])}],
                        '+1': [{'coords': [[4, 16]],
                                'omega_classes': 4, 'omega_shape': 1.25,
                                'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                 0.5158061369385518, 1.1217646558655263],
                                'dn_values': [0.13695378264465718, 0.4767518562354524,
                                              0.9999999999999997, 2.3862943611198904],
                                'ds_values': [0.6137056388801096, 3.386294361119891],
                                'orf_map': np.array([0, 1])}],
                        '+2': [], '-0': [], '-1': [], '-2': []}

        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, CAT_VALUES)
        self.circular_sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, CAT_VALUES, circular=True)
        self.plus_0_orf = sorted_orfs['+0'][0]
        self.plus_1_orf = sorted_orfs['+1'][0]
        self.plus_0_codons = self.sequence5.find_codons('+0', self.plus_0_orf)
        self.plus_1_codons = self.sequence5.find_codons('+1', self.plus_1_orf)

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
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'G': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'T': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}}}},
                                'C': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'C': None,
                                                'G': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'T': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}}}},
                                'G': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'C': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'G': None,
                                                'T': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}}}},
                                'T': {'from_nt': {'A': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'C': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'G': {'mu1': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu2': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu3': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu4': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu5': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}},
                                                        'mu6': {(1, 0): {},
                                                                (1, 1): {},
                                                                (0, 1): {}}},
                                                'T': None}}}}
        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)

        exp_total_omegas = {(0.04044355511945654, None): {'value': 0.04044355511945654, 'nt_events': 1}, (3.888337029906392, 1.6294456766354635): {'value': 6.335833962682549, 'nt_events': 1}, (0.6137056388801096, 1.6294456766354635): {'value': 0.9999999999999997, 'nt_events': 1}, (None, 1.6294456766354635): {'value': 1.6294456766354635, 'nt_events': 1}, (None, 0.6137056388801096): {'value': 0.6137056388801096, 'nt_events': 1}}

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
        self.assertEqual(exp_total_omegas, self.sequence5.total_omegas)

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
        self.assertEqual(exp_total_omegas, self.sequence5.total_omegas)

        # Tests a mutation that is synonymous in one frame and non-synonymous in the other
        random.seed(9991)
        nt = self.sequence5.nt_sequence[8]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 1.473755167554267e-06, 'C': None, 'G': 1.473755167554267e-06, 'T': 4.91251722518089e-06}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (0.6137056388801096, 1.6294456766354635), 'C': None, 'G': (0.6137056388801096, 1.6294456766354635), 'T': (0.6137056388801096, 1.6294456766354635)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 7.860027560289425e-06
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence5.total_omegas)

        # Tests a nucleotide involved in multiple non-synonymous codons
        random.seed(9991)
        nt = self.sequence5.nt_sequence[7]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 9.337468043269239e-06, 'C': None, 'G': 9.337468043269239e-06, 'T': 3.1124893477564126e-05}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': (3.888337029906392, 1.6294456766354635), 'C': None, 'G': (3.888337029906392, 1.6294456766354635), 'T': (3.888337029906392, 1.6294456766354635)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 4.979982956410261e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence5.total_omegas)

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
        result = self.sequence5.find_codons('+0', self.plus_0_orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence5.get_codons()[:4]
        self.assertEqual(len(expected), len(result))
        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Check Codons in +1 frame
        expected = ['ATG', 'CCT', 'GAC', 'TAA']
        result = self.sequence5.find_codons('+1', self.plus_1_orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+1')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence5.get_codons()[4:]
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0, 1): {'nt_events': 1,
                                                    'region_weight': 1.6294456766354635,
                                                    (None, 1.6294456766354635): [c12]},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 6.335833962682549,
                                                    (3.888337029906392, 1.6294456766354635): [c7]}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.9999999999999997,
                                                    (0.6137056388801096, 1.6294456766354635): [c8]}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 3},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 1,
                                                    'region_weight': 0.04044355511945654,
                                                    (0.04044355511945654, None): [a3,
                                                                                  a3]},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 1},
                             'C': None,
                             'G': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 1},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 1,
                                                    'region_weight': 0.04044355511945654,
                                                    (0.04044355511945654, None): [a3,
                                                                                  a3]},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 6.335833962682549,
                                                    (3.888337029906392, 1.6294456766354635): [c7]}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 2,
                                           (0, 1): {'nt_events': 1,
                                                    'region_weight': 1.6294456766354635,
                                                    (None, 1.6294456766354635): [c12]},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.9999999999999997,
                                                    (0.6137056388801096, 1.6294456766354635): [c8]}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 3},
                             'G': None,
                             'T': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0}},
                 'nt_events': 4},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 1,
                                                    'region_weight': 0.04044355511945654,
                                                    (0.04044355511945654, None): [a3,
                                                                                  a3]},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 0.9999999999999997,
                                                    (0.6137056388801096, 1.6294456766354635): [c8]}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 2,
                                           (0, 1): {'nt_events': 1,
                                                    'region_weight': 0.6137056388801096,
                                                    (None, 0.6137056388801096): [c12]},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 1,
                                                    'region_weight': 6.335833962682549,
                                                    (3.888337029906392, 1.6294456766354635): [c7]}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0, 1): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 0): {'nt_events': 0,
                                                    'region_weight': 0},
                                           (1, 1): {'nt_events': 0,
                                                    'region_weight': 0}},
                                   'nt_events': 0},
                             'T': None},
                 'nt_events': 4}}}

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
            {'A': None, 'C': 9.650503070197635e-07, 'G': 3.2168343567325454e-06, 'T': 3.4996622123451855e-07},
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': 7.559175365337173e-05, 'C': None, 'G': 2.741262314527127e-05, 'T': 0.0007376370109003972},
            {'A': 3.4926910107414796e-05, 'C': None, 'G': 1.1930829327062522e-05, 'T': 1.4422002063568881e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': 1.244990680753976e-05, 'C': None, 'G': 1.944063826565763e-05, 'T': 7.144947227193056e-05},
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [
            0,
            0,
            0,
            4.531850884986827e-06,
            0,
            0,
            0,
            0.0008406413876990402,
            6.12797414980462e-05,
            0,
            0,
            0,
            0.00010334001734512795,
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

    def testGetComplementState(self):
        nts = self.sequence5.nt_sequence
        self.assertEqual(nts, self.sequence5.get_sequence()) # testing get_sequence() also
        complements = "TACTTACGGACTGATT"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)

    def testGetRightNT(self):
        seq5 = ''.join(nt.state for nt in self.sequence5.nt_sequence)
        for pos, nt in enumerate(seq5[:-1]):
            result = self.sequence5.get_right_nt(pos)
            self.assertEqual(result.state, seq5[pos + 1])

        for pos, nt in enumerate(seq5[1:], 1):
            result = self.sequence5.get_left_nt(pos)
            self.assertEqual(result.state, seq5[pos - 1])

        circular_seq5 = ''.join(nt.state for nt in self.circular_sequence5.nt_sequence)

        for pos, nt in enumerate(circular_seq5):
            result = self.circular_sequence5.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence5[pos + 1])

        for pos, nt in enumerate(circular_seq5):
            result = self.circular_sequence5.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq5[pos - 1])           

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
                               'dn_values': [0.3329677267246186,
                                             1.0887942245237032,
                                             2.8982380487141928],
                               'ds_values': [0.6137056388801096, 3.386294361119891],
                               'orf_map': np.array([1])}], 
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, CAT_VALUES)
        self.circular_sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, CAT_VALUES, circular=True)
        self.plus_0_orf = sorted_orfs['+0'][0]
        self.seq6_codons = self.sequence6.find_codons('+0', self.plus_0_orf)

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
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'G': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'T': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}}}},
                                'C': {'from_nt': {'A': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'C': None,
                                                'G': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'T': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}}}},
                                'G': {'from_nt': {'A': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'C': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'G': None,
                                                'T': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}}}},
                                'T': {'from_nt': {'A': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'C': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'G': {'mu1': {(1,): {},
                                                                (0,): {}},
                                                        'mu2': {(1,): {},
                                                                (0,): {}},
                                                        'mu3': {(1,): {},
                                                                (0,): {}},
                                                        'mu4': {(1,): {},
                                                                (0,): {}},
                                                        'mu5': {(1,): {},
                                                                (0,): {}},
                                                        'mu6': {(1,): {},
                                                                (0,): {}}},
                                                'T': None}}}}
        
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
        exp_total_omegas = {(1.774131041895811,): {'value': 1.774131041895811, 'nt_events': 1}, (None,): {'value': 1, 'nt_events': 1}, (0.8558730398605124,): {'value': 0.8558730398605124, 'nt_events': 1}, (3.386294361119891,): {'value': 3.386294361119891, 'nt_events': 1}}
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

        # Tests internal methionine
        random.seed(4000)
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 3.4535394452935715e-05, 'G': 0.00011511798150978571, 'T': 3.4535394452935715e-05}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omega_keys = {'A': None, 'C': (1.774131041895811,), 'G': (1.774131041895811,), 'T': (1.774131041895811,)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)
        exp_cat_keys = {'C': 'mu4', 'G': 'mu4', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 0.00018418877041565714
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
        result = self.sequence6.find_codons('+0', self.plus_0_orf)
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        result = self.sequence6.get_codons()
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'mu1': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c8]}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c9]}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 1,
                                                  'region_weight': 1,
                                                  (None,): [g5, g5]},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g6]}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 2},
                             'T': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 6},
           'C': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3]}},
                                   'nt_events': 1},
                             'C': None,
                             'G': {'mu1': {'nt_events': 1,
                                           (0,): {'nt_events': 1,
                                                  'region_weight': 1,
                                                  (None,): [g5, g5]},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g6]}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 2},
                             'T': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 4},
           'G': {'from_nt': {'A': {'mu1': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3]}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c8]}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu4': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c9]}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'G': None,
                             'T': {'mu1': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [t4]}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1}},
                 'nt_events': 5},
           'T': {'from_nt': {'A': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [a3]}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 1},
                             'C': {'mu1': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c7]}},
                                   'mu3': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 3.386294361119891,
                                                  (3.386294361119891,): [c9]}},
                                   'mu4': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 0.8558730398605124,
                                                  (0.8558730398605124,): [c8]}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'nt_events': 3},
                             'G': {'mu1': {'nt_events': 1,
                                           (0,): {'nt_events': 1,
                                                  'region_weight': 1,
                                                  (None,): [g5, g5]},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu2': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu3': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu4': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu5': {'nt_events': 0,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 0,
                                                  'region_weight': 0}},
                                   'mu6': {'nt_events': 1,
                                           (0,): {'nt_events': 0,
                                                  'region_weight': 0},
                                           (1,): {'nt_events': 1,
                                                  'region_weight': 1.774131041895811,
                                                  (1.774131041895811,): [g6]}},
                                   'nt_events': 2},
                             'T': None},
                 'nt_events': 6}}}

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
            {'A': None, 'C': 0.00010110065147321811, 'G': 1.4219943599698952e-05, 'T': 2.21166834455335e-05},
            {'A': 1.640915223378292e-05, 'C': 5.469717411260974e-05, 'G': 3.1650842205781537e-06, 'T': None},
            {'A': 3.083040250181358e-05, 'C': 1.7840194133551658e-06, 'G': None, 'T': 1.7840194133551658e-06},
            {'A': 8.541011531371199e-05, 'C': 3.943337868448051e-05, 'G': None, 'T': 7.501016077045216e-05},
            {'A': 1.902337814381081e-05, 'C': None, 'G': 7.9160730928051e-06, 'T': 1.4942013851657346e-05},
            {'A': 1.5268941184784538e-06, 'C': None, 'G': 4.482604155497203e-06, 'T': 4.120339101347877e-05},
            {'A': 3.1320245442884365e-05, 'C': None, 'G': 4.8906836931339996e-05, 'T': 0.0001044008181429479},
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},        
            {'A': None, 'C': None, 'G': None, 'T': None},
        ]

        exp_total_rates = [       
                0,
                0,
                0,
                0.00013743727851845056,
                7.427141056697082e-05,
                3.4398441328523915e-05,
                0.00019985365476864466,
                4.188146508827326e-05,
                4.721288928745443e-05,
                0.00018462790051717226,
                0,
                0,
                0,
            ]

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

    def testGetComplementState(self):
        nts = self.sequence6.nt_sequence
        self.assertEqual(nts, self.sequence6.get_sequence()) # testing get_sequence() also
        complements = "TACTACCGGGATT"
        for pos, nt in enumerate(nts):
            result = nt.get_complement_state()
            expected = complements[pos]
            self.assertEqual(result, expected)

    def testGetRightNT(self):
        seq6 = ''.join(nt.state for nt in self.sequence6.nt_sequence)
        for pos, nt in enumerate(seq6[:-1]):
            result = self.sequence6.get_right_nt(pos)
            self.assertEqual(result.state, seq6[pos + 1])

        for pos, nt in enumerate(seq6[1:], 1):
            result = self.sequence6.get_left_nt(pos)
            self.assertEqual(result.state, seq6[pos - 1])

        circular_seq6 = ''.join(nt.state for nt in self.circular_sequence6.nt_sequence)

        for pos, nt in enumerate(circular_seq6):
            result = self.circular_sequence6.get_right_nt(pos)
            self.assertEqual(result.state, self.circular_sequence6[pos + 1])

        for pos, nt in enumerate(circular_seq6):
            result = self.circular_sequence6.get_left_nt(pos)
            self.assertEqual(result.state, circular_seq6[pos - 1])


if __name__ == '__main__':
    unittest.main()
