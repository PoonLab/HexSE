import random
import unittest

from src.sequence_info import Sequence, Nucleotide

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

        orfs = {'+0': [{'coords': [[(0, 21)]],
                        'omega_classes': 3, 'omega_shape': 1.5,
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

        s1 = 'GTACGATCGATCGATGCTAGC'
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, orfs, KAPPA, GLOBAL_RATE, pi1, CAT_VALUES)
        orf = {'coords': [[(0, 21)]],
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
        self.assertEqual(self.sequence1.orfs, new_sequence1.orfs)
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
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testGetFrequencyRates(self):
        expected = {'A': 0.24, 'C': 0.24, 'G': 0.29, 'T': 0.24}
        result = Sequence.get_frequency_rates(str(self.sequence1))
        self.assertEqual(expected, result)

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': None}}}}
        result = self.sequence1.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(9001)  # Set seed value to initialize pseudo-random number generator

        nt = self.sequence1.nt_sequence[0]  # First nucleotide is G
        self.sequence1.set_substitution_rates(nt)

        exp_sub_rates = {'A': 6.640901570929766e-06,
                         'C': 1.082033797048673e-06,
                         'G': None,
                         'T': 1.3389485626674619e-05}
        self.assertEqual(exp_sub_rates, nt.rates)
        exp_omegas = {'A': ((1, 0, 0, 0),), 'C': ((0, 1, 0, 0),), 'G': None, 'T': ((0, 0, 1, 0),)}
        self.assertEqual(exp_omegas, nt.omega_keys)
        exp_cat_keys = {'A': 'mu3', 'C': 'mu1', 'T': 'mu3'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)
        exp_total_rate = 2.1112420994653056e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {((0, 0, 1, 0),): 1.1481358615121404,
                            ((0, 1, 0, 0),): 0.4810288100937172,
                            ((1, 0, 0, 0),): 0.1708353283825978}
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
        orf = {'coords': [[(0, 21)]],
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
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1),): [c11],
                                                                                 ((0, 1, 0, 0),): [c16]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [c3]},
                                                                         'mu6': {((0, 1, 0, 0),): [c20]}}},
                                                      'G': {'category': {'mu1': {((0, 1, 0, 0),): [g15]},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [g12,
                                                                                                   g19],
                                                                                 ((1, 0, 0, 0),): [g0]},
                                                                         'mu4': {((0, 0, 0, 1),): [g8]},
                                                                         'mu5': {((0, 0, 1, 0),): [g4]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 0, 0, 1),): [t17]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 0, 1, 0),): [t6]},
                                                                         'mu5': {((0, 0, 1, 0),): [t1,
                                                                                                   t10]},
                                                                         'mu6': {((0, 0, 1, 0),): [t14]}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 1, 0),): [a18]},
                                                                         'mu2': {((0, 0, 0, 1),): [a5],
                                                                                 ((0, 1, 0, 0),): [a9]},
                                                                         'mu3': {((0, 0, 0, 1),): [a2]},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [a13]},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {((0, 1, 0, 0),): [g0]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [g19],
                                                                                 ((1, 0, 0, 0),): [g12]},
                                                                         'mu5': {((0, 0, 0, 1),): [g8],
                                                                                 ((1, 0, 0, 0),): [g4]},
                                                                         'mu6': {((0, 0, 1, 0),): [g15]}}},
                                                      'T': {'category': {'mu1': {((0, 0, 1, 0),): [t1]},
                                                                         'mu2': {((0, 0, 0, 1),): [t17],
                                                                                 ((1, 0, 0, 0),): [t6]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1),): [t14],
                                                                                 ((1, 0, 0, 0),): [t10]}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [a9]},
                                                                         'mu5': {((0, 0, 0, 1),): [a2,
                                                                                                   a5],
                                                                                 ((1, 0, 0, 0),): [a18]},
                                                                         'mu6': {((0, 0, 1, 0),): [a13]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 1, 0),): [c20],
                                                                                 ((1, 0, 0, 0),): [c7]},
                                                                         'mu3': {((0, 0, 1, 0),): [c3]},
                                                                         'mu4': {((0, 1, 0, 0),): [c16]},
                                                                         'mu5': {((0, 1, 0, 0),): [c11]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [t6]},
                                                                         'mu3': {((0, 1, 0, 0),): [t14]},
                                                                         'mu4': {((0, 0, 0, 1),): [t17]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [t1],
                                                                                 ((1, 0, 0, 0),): [t10]}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {((1, 0, 0, 0),): [a18]},
                                                                         'mu3': {((1, 0, 0, 0),): [a13]},
                                                                         'mu4': {((0, 0, 0, 1),): [a5],
                                                                                 ((1, 0, 0, 0),): [a9]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1),): [a2]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1),): [c20]},
                                                                         'mu4': {((0, 0, 1, 0),): [c16]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1),): [c11],
                                                                                 ((0, 1, 0, 0),): [c7]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [g0,
                                                                                                   g12,
                                                                                                   g15]},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [g8],
                                                                                 ((1, 0, 0, 0),): [g19]},
                                                                         'mu6': {((0, 1, 0, 0),): [g4]}}},
                                                      'T': None}}}}

        res_omega_key = self.sequence1.nt_in_event_tree(g0)
        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': 6.640901570929766e-06, 'C': 1.082033797048673e-06, 'G': None, 'T': 1.3389485626674619e-05},
            {'A': 2.6629007650064053e-05, 'C': 7.1245101429805794e-06, 'G': 5.065369013827426e-05, 'T': None},
            {'A': None, 'C': 9.651256435350337e-06, 'G': 7.73108495336449e-05, 'T': 4.4118202240945026e-05},
            {'A': 2.3193254860093468e-05, 'C': None, 'G': 1.1080953622075548e-05, 'T': 0.0},
            {'A': 0.00010725572525720243, 'C': 4.787691333261464e-06, 'G': None, 'T': 2.5643402645651533e-05},
            {'A': None, 'C': 5.46517972935242e-06, 'G': 7.73108495336449e-05, 'T': 1.5070521255236873e-05},
            {'A': 1.730300590481841e-05, 'C': 3.112152579112792e-06, 'G': 2.628908902158698e-06, 'T': None},
            {'A': 0.0, 'C': None, 'G': 9.336457737338375e-07, 'T': 7.074042109145251e-05},
            {'A': 6.070071061137074e-05, 'C': 2.8025182955946275e-05, 'G': None, 'T': 2.8025182955946275e-05},
            {'A': None, 'C': 2.628908902158698e-06, 'G': 2.4164516356328888e-05, 'T': 2.574577447535311e-06},
            {'A': 2.6629007650064053e-05, 'C': 2.512315855827235e-05, 'G': 7.536947567481705e-06, 'T': None},
            {'A': 1.8615854748053904e-06, 'C': None, 'G': 1.1156623787551084e-05, 'T': 0.00014706067413648342},
            {'A': 4.463161875558207e-05, 'C': 3.1109477491051676e-06, 'G': None, 'T': 1.3389485626674619e-05},
            {'A': None, 'C': 3.96222731028535e-06, 'G': 0.00016884563379424754, 'T': 1.648775562437735e-06},
            {'A': 5.065369013827426e-05, 'C': 0.00014706067413648342, 'G': 4.642532399005903e-06, 'T': None},
            {'A': 3.606779323495576e-06, 'C': 6.120654225041473e-05, 'G': None, 'T': 1.3389485626674619e-05},
            {'A': 8.954762458333844e-07, 'C': None, 'G': 7.2493549068986665e-06, 'T': 5.7676686349394705e-05},
            {'A': 1.8615854748053904e-06, 'C': 1.8217265764508066e-05, 'G': 1.5070521255236873e-05, 'T': None},
            {'A': None, 'C': 2.137353042894174e-06, 'G': 1.3207424367617833e-05, 'T': 9.336457737338375e-07},
            {'A': 4.463161875558207e-05, 'C': 8.759637179169221e-06, 'G': None, 'T': 4.787691333261464e-06},
            {'A': 2.1222126327435754e-05, 'C': None, 'G': 6.274768836878727e-06, 'T': 3.217085478450112e-05},
        ]

        exp_total_rates = [2.1112420994653056e-05,
                           8.440720793131888e-05,
                           0.00013108030820994026,
                           3.427420848216902e-05,
                           0.00013768681923611542,
                           9.78465505182342e-05,
                           2.30440673860899e-05,
                           7.167406686518635e-05,
                           0.00011675107652326329,
                           2.9368002706022896e-05,
                           5.928911377581811e-05,
                           0.00016007888339883989,
                           6.113205213136186e-05,
                           0.00017445663666697063,
                           0.00020235689667376356,
                           7.820280720058493e-05,
                           6.582151750212676e-05,
                           3.5149372494550326e-05,
                           1.6278423184245844e-05,
                           5.817894726801275e-05,
                           5.96677499488156e-05]

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

        orfs = {'+0': [{'coords': [[(0, 21)]],
                        'omega_classes': 3, 'omega_shape': 1.5,
                        'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                '+1': [], '+2': [], '-0': [], '-1': [],
                '-2': [{'coords': [[(3, 15)]],
                        'omega_classes': 4, 'omega_shape': 1.25,
                        'omega_values': [0.09199853806558903, 0.27043066909631136,
                                         0.5158061369385518, 1.1217646558655263]}]}
        s2 = 'ATGAATAAACCCGTATGA'
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, orfs, KAPPA, GLOBAL_RATE, pi2, CAT_VALUES)
        plus_orf = {'coords': [[(0, 21)]],
                    'omega_classes': 3, 'omega_shape': 1.5,
                    'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}
        minus_orf = {'coords': [[(3, 15)]],
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
        self.assertEqual(self.sequence2.orfs, new_sequence2.orfs)
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
                self.assertEqual(codon.orf, new_codon.orf)
                self.assertEqual(codon.frame, new_codon.frame)
                self.assertEqual(str(codon.nts_in_codon), str(new_codon.nts_in_codon))

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
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

        exp_sub_rates = {'A': 3.5055586634238734e-05, 'C': None, 'G': 2.886834562234422e-06, 'T': 9.62278187411474e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 1), (0, 0, 0, 1, 0)),
                          'C': None,
                          'G': ((0, 0, 0, 1), (0, 1, 0, 0, 0)),
                          'T': ((0, 0, 0, 1), (0, 1, 0, 0, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu6', 'G': 'mu4', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 4.7565203070587895e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_total_omegas = {((0, 0, 0, 1), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                            ((0, 0, 0, 1), (1, 0, 0, 0, 0)): 0.09199853806558903,
                            ((0, 0, 1, 0), (0, 0, 0, 0, 1)): 1.1481358615121404,
                            ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 0, 1, 0), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((0, 1, 0, 0), (0, 0, 0, 0, 1)): 0.4810288100937172,
                            ((0, 1, 0, 0), (0, 1, 0, 0, 0)): 0.27043066909631136,
                            ((0, 1, 0, 0), (1, 0, 0, 0, 0)): 0.09199853806558903,
                            ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((1, 0, 0, 0), (1, 0, 0, 0, 0)): 0.09199853806558903}
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation would destroy START CODON
        nt = self.sequence2.nt_sequence[0]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)

        # Tuples are of length 1 because the nucleotide is part of only 1 ORF
        exp_omega_keys = {'A': None,
                          'C': ((0, 0, 0, 0),),
                          'G': ((0, 0, 0, 0),),
                          'T': ((0, 0, 0, 0),)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu4', 'G': 'mu6', 'T': 'mu6'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Mutation in STOP codon on reverse strand ORF
        nt = self.sequence2.nt_sequence[4]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)

        # Tuples are of length 2 because the nucleotide is part of 2 ORFs
        exp_omega_keys = {'A': None,
                          'C': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'G': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'T': ((0, 0, 0, 0), (0, 0, 0, 0, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu4', 'G': 'mu6', 'T': 'mu5'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0.0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence2.total_omegas)

        # Non-synonymous mutation in both (+) and (-) strands
        nt = self.sequence2.nt_sequence[10]
        self.sequence2.set_substitution_rates(nt)

        exp_sub_rates = {'A': 1.4476499901684232e-06, 'C': None, 'G': 1.6862747125805073e-05, 'T': 5.717766785346913e-07}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((1, 0, 0, 0), (0, 0, 1, 0, 0)),
                          'C': None,
                          'G': ((0, 1, 0, 0), (0, 0, 0, 1, 0)),
                          'T': ((0, 1, 0, 0), (0, 1, 0, 0, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu5', 'G': 'mu6', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 1.8882173794508187e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_total_omegas = {((0, 0, 0, 1), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                            ((0, 0, 0, 1), (1, 0, 0, 0, 0)): 0.09199853806558903,
                            ((0, 0, 1, 0), (0, 0, 0, 0, 1)): 1.1481358615121404,
                            ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 0, 1, 0), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((0, 1, 0, 0), (0, 0, 0, 0, 1)): 0.4810288100937172,
                            ((0, 1, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 1, 0, 0), (0, 1, 0, 0, 0)): 0.27043066909631136,
                            ((0, 1, 0, 0), (1, 0, 0, 0, 0)): 0.09199853806558903,
                            ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((1, 0, 0, 0), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((1, 0, 0, 0), (1, 0, 0, 0, 0)): 0.09199853806558903}
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
        result = self.sequence2.find_codons('+0', {'coords': [[(0, 21)]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Test reverse strand ORF
        expected = ['ATG', 'CCC', 'AAA', 'TAA']
        result = self.sequence2.find_codons('-2', {'coords': [[(3, 15)]],
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {((0, 1, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): [c11],
                                                                                 ((0, 1, 0, 0), (0, 1, 0, 0, 0)): [c10]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): [a6]},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): [a8]},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [a8],
                                                                                 ((0, 0, 1, 0), (0, 0, 0, 0, 1)): [a6]},
                                                                         'mu2': {((0, 0, 1, 0), (0, 0, 1, 0, 0)): [a7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((1, 0, 0, 0), (0, 0, 0, 1, 0)): [c10]},
                                                                         'mu2': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [c11],
                                                                                 ((0, 1, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 1, 0), (0, 0, 0, 1, 0)): [c10]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 1, 0), (0, 0, 0, 0, 1)): [c9]},
                                                                         'mu6': {((0, 0, 0, 1), (1, 0, 0, 0, 0)): [c11]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}

        res_omega_key = self.sequence2.nt_in_event_tree(c11)
        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': None, 'C': 4.439677739210545e-07, 'G': 1.3061601928797729e-05, 'T': 0.0},
            {'A': None, 'C': 6.682850533879727e-07, 'G': 1.977900389625616e-05, 'T': 1.574724181334861e-07},
            {'A': None, 'C': 7.830280216562374e-07, 'G': 5.867994075412873e-06, 'T': 0.0},
            {'A': 6.342956741319806e-07, 'C': None, 'G': 1.8621438056957443e-06, 'T': 6.287404584042902e-05},
            {'A': 4.065205624182217e-06, 'C': None, 'G': 2.5269700838806675e-07, 'T': 5.661016765365336e-06},
            {'A': 8.451064757202812e-06, 'C': None, 'G': 1.996772714409441e-06, 'T': 9.583301644476654e-06},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        ]

        exp_total_rates = [
            0,
            0,
            0,
            0,
            0,
            0,
            1.3505569702718784e-05,
            2.060476136777762e-05,
            6.65102209706911e-06,
            6.537048532025675e-05,
            9.97891939793562e-06,
            2.0031139116088907e-05,
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
        sorted_orfs = {'+0': [{'coords': [[(0, 12)]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, sorted_orfs, KAPPA, GLOBAL_RATE, pi3, CAT_VALUES)
        self.seq3_codons = self.sequence3.find_codons('+0', {'coords': [[(0, 12)]],
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
        new_sequence1 = self.sequence3.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence3, new_sequence1)
        self.assertEqual(self.sequence3.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence3.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence3.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence3.pi, new_sequence1.pi)
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

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': None}}}}
        result = self.sequence3.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(555)    # Set seed value to initialize pseudo-random number generator
        nt = self.sequence3.nt_sequence[9]
        self.sequence3.set_substitution_rates(nt)

        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 0),), 'C': ((0, 0, 0, 0),), 'G': ((0, 0, 0, 0),), 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu2', 'C': 'mu3', 'G': 'mu2'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {((0, 0, 1, 0),): 1.1481358615121404,
                            ((0, 1, 0, 0),): 0.4810288100937172,
                            ((1, 0, 0, 0),): 0.1708353283825978}
        self.assertEqual(exp_total_omegas, self.sequence3.total_omegas)

        # Tests a synonymous mutation
        random.seed(555)
        nt = self.sequence3.nt_sequence[5]
        self.sequence3.set_substitution_rates(nt)

        exp_sub_rates = {'A': 3.1880215087889116e-05,
                         'C': 1.688969876186309e-05,
                         'G': None,
                         'T': 9.564064526366735e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 1),), 'C': ((0, 0, 0, 1),), 'G': None, 'T': ((0, 0, 0, 1),)}  # Synonymous
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

        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 0),), 'C': ((0, 0, 0, 0),), 'G': ((0, 0, 0, 0),), 'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu2', 'C': 'mu3', 'G': 'mu2'}
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
        result = self.sequence3.find_codons('+0', {'coords': [[(0, 12)]],
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {((0, 1, 0, 0),): [c4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 0, 0, 1),): [g5]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 0, 1, 0),): [t6]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [a3]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {((0, 0, 0, 1),): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [g7]},
                                                                         'mu6': {((0, 0, 1, 0),): [g8]}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [t6]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [c4]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 1, 0),): [t6]},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [a3]},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [c4]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {((0, 0, 0, 1),): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [g7]},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g8]},
                                                                         'mu6': {}}},
                                                      'T': None}}}}

        res_omega_key = self.sequence3.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 2.7384467730819763e-06, 'G': 0.0001758808685356745, 'T': 1.1621483112032378e-05},
            {'A': 2.984920819444615e-07, 'C': None, 'G': 1.3207424367617832e-06, 'T': 2.8606416083725682e-06},
            {'A': 8.791137398888177e-05, 'C': 3.257774580909433e-06, 'G': None, 'T': 3.257774580909433e-06},
            {'A': 2.2264094196814308e-06, 'C': 5.724915147353247e-06, 'G': 2.7738549635483385e-05, 'T': None},
            {'A': 0.0, 'C': 1.9524091628214396e-05, 'G': None, 'T': 8.124431698260329e-06},
            {'A': 0.0, 'C': 8.864395774197996e-05, 'G': None, 'T': 6.933897792999363e-06},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           0.00019024079842078887,
                           4.4798761270788134e-06,
                           9.442692315070064e-05,
                           3.5689874202518064e-05,
                           2.7648523326474725e-05,
                           9.557785553497932e-05,
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
        sorted_orfs = {'+0': [{'coords': [[(0, 12)]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        random.seed(4000)
        self.sequence4 = Sequence(s4, sorted_orfs, KAPPA, GLOBAL_RATE, pi4, CAT_VALUES)
        self.seq4_codons = self.sequence4.find_codons('+0', {'coords': [[(0, 12)]],
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
        new_sequence1 = self.sequence4.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence4, new_sequence1)
        self.assertEqual(self.sequence4.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence4.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence4.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence4.pi, new_sequence1.pi)
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

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': None}}}}
        result = self.sequence4.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        random.seed(5001)

        # Only 3 possible omegas because there is only 1 ORF
        exp_total_omegas = {((0, 0, 1, 0),): 1.1481358615121404,
                            ((0, 1, 0, 0),): 0.4810288100937172,
                            ((1, 0, 0, 0),): 0.1708353283825978}

        # Tests nucleotide involved in a stop codon
        nt = self.sequence4.nt_sequence[11]
        self.sequence4.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': ((0, 0, 0, 0),), 'G': ((0, 0, 0, 0),), 'T': ((0, 0, 0, 0),)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu5', 'G': 'mu3', 'T': 'mu6'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_total_omegas, self.sequence4.total_omegas)

        # Tests mutation in internal methionine
        nt = self.sequence4.nt_sequence[3]
        self.sequence4.set_substitution_rates(nt)

        exp_sub_rates = {'A': None,
                         'C': 2.918042370022416e-05,
                         'G': 0.00012204961839612692,
                         'T': 4.37283778113398e-07}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': ((0, 1, 0, 0),), 'G': ((0, 0, 1, 0),), 'T': ((1, 0, 0, 0),)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu6', 'G': 'mu5', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0.00015166732587446445
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
        result = self.sequence4.find_codons('+0', {'coords': [[(0, 12)]],
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

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {((0, 0, 1, 0),): [c7]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [c6]},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [g5]}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [a3]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [g5]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 1, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 1, 0, 0),): [c7]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [c6]},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [t4]},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [c7]},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [c6]},
                                                                         'mu6': {((0, 0, 0, 1),): [c8]}}},
                                                      'G': {'category': {'mu1': {((1, 0, 0, 0),): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}

        res_omega_key = self.sequence4.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 9.967862996985666e-06, 'G': 9.796201446598298e-06, 'T': 6.3834820486331165e-06},
            {'A': 9.32787756076442e-07, 'C': 3.1092925202548074e-06, 'G': 2.6818515078492823e-06, 'T': None},
            {'A': 0.00011959899060425867, 'C': 7.849008815636847e-06, 'G': None, 'T': 2.2526740084629595e-07},
            {'A': 2.6818515078492823e-06, 'C': None, 'G': 7.551411361352776e-06, 'T': 3.873827704010793e-05},
            {'A': 2.2264094196814308e-06, 'C': None, 'G': 9.32787756076442e-07, 'T': 1.611990416321494e-05},
            {'A': 2.415964047926403e-05, 'C': None, 'G': 2.415964047926403e-05, 'T': 0.00015318820222550355},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           2.614754649221708e-05,
                           6.723931784180533e-06,
                           0.00012767326682074182,
                           4.8971539909309985e-05,
                           1.9279101338972813e-05,
                           0.0002015074831840316,
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
        sorted_orfs = {'+0': [{'coords': [[(0, 12)]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                       '+1': [{'coords': [[(4, 16)]],
                               'omega_classes': 4, 'omega_shape': 1.25,
                               'omega_values': [0.09199853806558903, 0.27043066909631136,
                                                0.5158061369385518, 1.1217646558655263]}],
                       '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, sorted_orfs, KAPPA, GLOBAL_RATE, pi5, CAT_VALUES)

        self.plus_0_codons = self.sequence5.find_codons('+0', {'coords': [[(0, 12)]],
                                                               'omega_classes': 3, 'omega_shape': 1.5,
                                                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.plus_1_codons = self.sequence5.find_codons('+1', {'coords': [[(4, 16)]],
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
        new_sequence1 = self.sequence5.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence5, new_sequence1)
        self.assertEqual(self.sequence5.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence5.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence5.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence5.pi, new_sequence1.pi)
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

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': None}}}}
        result = self.sequence5.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        random.seed(9991)

        exp_all_omegas = {((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                          ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                          ((0, 0, 1, 0),): 1.1481358615121404,
                          ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 0, 1, 0, 0),): 0.5158061369385518,
                          ((0, 1, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 1, 0, 0, 0),): 0.27043066909631136,
                          ((1, 0, 0, 0),): 0.1708353283825978,
                          ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263}

        # Tests a nucleotide that is involved in multiple codons, one of which is a start codon
        nt = self.sequence5.nt_sequence[4]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None,
                          'C': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'G': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'T': ((0, 0, 0, 0), (0, 0, 0, 0, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu1', 'G': 'mu1', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests mutations that would destroy a stop codon in +1 frame
        random.seed(9991)
        nt = self.sequence5.nt_sequence[9]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'C': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'G': ((0, 0, 0, 0), (0, 0, 0, 0, 0)),
                          'T': None}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'C': 'mu1', 'G': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_total_omegas = {((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                            ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                            ((0, 0, 1, 0),): 1.1481358615121404,
                            ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 0, 1, 0, 0),): 0.5158061369385518,
                            ((0, 1, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                            ((0, 1, 0, 0, 0),): 0.27043066909631136,
                            ((1, 0, 0, 0),): 0.1708353283825978,
                            ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263}
        self.assertEqual(exp_total_omegas, self.sequence5.total_omegas)

        # Tests a mutation that is synonymous in one frame and non-synonymous in the other
        random.seed(9991)
        nt = self.sequence5.nt_sequence[8]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 1.355833208815998e-07,
                         'C': None,
                         'G': 1.0976188559998076e-06,
                         'T': 5.510688194538514e-06}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((0, 0, 0, 1), (1, 0, 0, 0, 0)),
                          'C': None,
                          'G': ((0, 0, 0, 1), (1, 0, 0, 0, 0)),
                          'T': ((0, 0, 0, 1), (0, 0, 0, 1, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu4', 'T': 'mu1'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 6.7438903714199215e-06
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_all_omegas = {((0, 0, 0, 1), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                          ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                          ((0, 0, 0, 1), (1, 0, 0, 0, 0)): 0.09199853806558903,
                          ((0, 0, 1, 0),): 1.1481358615121404,
                          ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 0, 1, 0, 0),): 0.5158061369385518,
                          ((0, 1, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 1, 0, 0, 0),): 0.27043066909631136,
                          ((1, 0, 0, 0),): 0.1708353283825978,
                          ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263}
        self.assertEqual(exp_all_omegas, self.sequence5.total_omegas)

        # Tests a nucleotide involved in multiple non-synonymous codons
        random.seed(9991)
        nt = self.sequence5.nt_sequence[7]
        self.sequence5.set_substitution_rates(nt)

        exp_sub_rates = {'A': 2.316242114601123e-08,
                         'C': None,
                         'G': 7.952399355049068e-07,
                         'T': 1.0258954496212767e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': ((1, 0, 0, 0), (1, 0, 0, 0, 0)),
                          'C': None,
                          'G': ((0, 1, 0, 0), (0, 0, 0, 1, 0)),
                          'T': ((1, 0, 0, 0), (0, 0, 1, 0, 0))}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu1', 'G': 'mu1', 'T': 'mu6'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 1.1077356852863685e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)

        exp_all_omegas = {((0, 0, 0, 1), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 0, 0, 1), (0, 0, 1, 0, 0)): 0.5158061369385518,
                          ((0, 0, 0, 1), (0, 1, 0, 0, 0)): 0.27043066909631136,
                          ((0, 0, 0, 1), (1, 0, 0, 0, 0)): 0.09199853806558903,
                          ((0, 0, 1, 0),): 1.1481358615121404,
                          ((0, 0, 1, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 0, 1, 0, 0),): 0.5158061369385518,
                          ((0, 1, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((0, 1, 0, 0, 0),): 0.27043066909631136,
                          ((1, 0, 0, 0),): 0.1708353283825978,
                          ((1, 0, 0, 0), (0, 0, 0, 1, 0)): 1.1217646558655263,
                          ((1, 0, 0, 0), (0, 0, 1, 0, 0)): 0.5158061369385518,
                          ((1, 0, 0, 0), (1, 0, 0, 0, 0)): 0.09199853806558903}
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
        result = self.sequence5.find_codons('+0', {'coords': [[(0, 12)]],
                                                   'omega_classes': 3, 'omega_shape': 1.5,
                                                   'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]})
        self.assertEqual(len(expected), len(result))

        for idx, codon in enumerate(result):
            self.assertEqual(codon.frame, '+0')
            self.assertEqual(expected[idx], str(codon))

        # Check Codons in +1 frame
        expected = ['ATG', 'CCT', 'GAC', 'TAA']
        result = self.sequence5.find_codons('+1', {'coords': [[(4, 16)]],
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
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): [c8],
                                                                                 ((1, 0, 0, 0), (0, 0, 0, 1, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0, 0),): [c12]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {((1, 0, 0, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [c8]},
                                                                         'mu2': {((0, 0, 1, 0), (0, 0, 0, 1, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0, 0),): [c12]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 0, 1),): [c12]},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [c8]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0), (0, 0, 0, 1, 0)): [c7]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        res_omega_key = self.sequence5.nt_in_event_tree(a3)
        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 5.035388960093674e-07, 'G': 5.8482810783176505e-05, 'T': 2.6105613071930803e-06},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': 8.291369071409994e-07, 'C': None, 'G': 5.572394341408454e-06, 'T': 6.282199909613655e-05},
            {'A': 1.170045500327795e-06, 'C': None, 'G': 7.601719597693948e-07, 'T': 1.3136856759987016e-05},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 1.8015514577705695e-05, 'C': None, 'G': 4.965465882349613e-06, 'T': 4.912517225180892e-06},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        ]

        exp_total_rates = [
            0,
            0,
            0,
            6.159691098637895e-05,
            0,
            0,
            0,
            6.922353034468601e-05,
            1.5067074220084205e-05,
            0,
            0,
            0,
            2.78934976852362e-05,
            0,
            0,
            0
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
        sorted_orfs = {'+0': [{'coords': [[(0, 5), (6, 13)]],
                               'omega_classes': 3, 'omega_shape': 1.5,
                               'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                       '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, sorted_orfs, KAPPA, GLOBAL_RATE, pi6, CAT_VALUES)
        self.seq6_codons = self.sequence6.find_codons('+0', {'coords': [[(0, 5), (6, 13)]],
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
        new_sequence1 = self.sequence6.__deepcopy__(memodict={})
        self.assertIsNot(self.sequence6, new_sequence1)
        self.assertEqual(self.sequence6.orfs, new_sequence1.orfs)
        self.assertEqual(self.sequence6.kappa, new_sequence1.kappa)
        self.assertEqual(self.sequence6.global_rate, new_sequence1.global_rate)
        self.assertEqual(self.sequence6.pi, new_sequence1.pi)
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

    def testCreateEventTree(self):
        expected = {'to_nt': {'A': {'from_nt': {'A': None,
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': None,
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': None,
                                                'T': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}}}},
                              'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'C': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'G': {'category': {'mu1': {},
                                                                   'mu2': {},
                                                                   'mu3': {},
                                                                   'mu4': {},
                                                                   'mu5': {},
                                                                   'mu6': {}}},
                                                'T': None}}}}
        result = self.sequence6.create_event_tree()
        self.assertEqual(expected, result)

    def testSetSubstitutionRates(self):
        # Tests nucleotide not involved in an ORF
        random.seed(4000)
        nt = self.sequence6.nt_sequence[5]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': 4.814194289867335e-05,
                         'C': 1.4442582869602003e-05,
                         'G': None,
                         'T': 1.4442582869602003e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        # No omega keys because nt treated as synonymous mutation
        exp_omega_keys = {'A': (), 'C': (), 'G': None, 'T': ()}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'A': 'mu4', 'C': 'mu4', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 7.702710863787735e-05
        self.assertEqual(exp_total_rate, nt.mutation_rate)
        exp_total_omegas = {((0, 0, 1, 0),): 1.1481358615121404,
                            ((0, 1, 0, 0),): 0.4810288100937172,
                            ((1, 0, 0, 0),): 0.1708353283825978}
        self.assertEqual(exp_total_omegas, self.sequence6.total_omegas)

        # Tests internal methionine
        nt = self.sequence6.nt_sequence[3]
        self.sequence6.set_substitution_rates(nt)

        exp_sub_rates = {'A': None,
                         'C': 8.104909747635022e-06,
                         'G': 0.0002180922769842364,
                         'T': 2.2349715960390447e-05}
        self.assertEqual(exp_sub_rates, nt.rates)

        exp_omega_keys = {'A': None, 'C': ((0, 0, 1, 0),), 'G': ((0, 0, 1, 0),), 'T': ((0, 0, 1, 0),)}
        self.assertEqual(exp_omega_keys, nt.omega_keys)

        exp_cat_keys = {'C': 'mu2', 'G': 'mu6', 'T': 'mu4'}
        self.assertEqual(exp_cat_keys, nt.cat_keys)

        exp_total_rate = 0.00024854690269226186
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
        result = self.sequence6.find_codons('+0', {'coords': [[(0, 5), (6, 13)]],
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
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1),): [c9],
                                                                                 ((0, 1, 0, 0),): [c7],
                                                                                 ((1, 0, 0, 0),): [c8]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [g6]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [a3]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {((1, 0, 0, 0),): [g6]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 1, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 0, 1),): [c9]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [c8]},
                                                                         'mu6': {((0, 0, 1, 0),): [c7]}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [t4]},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [c7]},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [c9]},
                                                                         'mu6': {((0, 0, 1, 0),): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [g6]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}


        res_omega_key = self.sequence6.nt_in_event_tree(g5)
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

        # No new omega keys are created initially
        self.assertEqual({}, res_omega_key)

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
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': 0.0, 'C': 0.0, 'G': None, 'T': 0.0},
            {'A': None, 'C': 9.363750088077444e-06, 'G': 9.202492268016581e-06, 'T': 5.9966043487159576e-06},
            {'A': 8.581647355903267e-07, 'C': 2.8605491186344227e-06, 'G': 2.4673033872213396e-06, 'T': None},
            {'A': 7.40895641364097e-05, 'C': 4.227994381423898e-05, 'G': None, 'T': 2.2226869240922907e-05},
            {'A': 8.397903437451395e-06, 'C': 8.947438664949276e-07, 'G': None, 'T': 4.4490935490473226e-06},
            {'A': 4.4490935490473226e-06, 'C': None, 'G': 4.854311971584616e-05, 'T': 5.266921935564986e-06},
            {'A': 1.5800765806694958e-06, 'C': None, 'G': 1.0691764463069788e-05, 'T': 0.00016181039905282054},
            {'A': 9.249120750544072e-06, 'C': None, 'G': 5.237463907296068e-06, 'T': 7.40895641364097e-05},
            {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': None},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0},
            {'A': None, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        ]

        exp_total_rates = [0,
                           0,
                           0,
                           2.4562846704809983e-05,
                           6.186017241446089e-06,
                           0.0001385963771915716,
                           1.3741740852993645e-05,
                           5.8259135200458473e-05,
                           0.00017408224009655984,
                           8.857614879424983e-05,
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


if __name__ == '__main__':
    unittest.main()
