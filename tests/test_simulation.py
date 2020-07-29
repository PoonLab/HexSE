import os
import random
import unittest

import numpy as np
from Bio import Phylo

from src.sequence_info import Sequence
from src.simulation import SimulateOnBranch
from src.simulation import SimulateOnTree

TEST_TREE = os.path.join(os.path.dirname(__file__), 'fixtures/test_tree.txt')


# ==========================================
# Tests for SimulateOnBranch
# ==========================================
class TestSimulateOnBranch(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        random.seed(9001)  # Set seed for pseudo-random number generator
        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.5
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, kappa, mu, pi1, dN_values, dS_values)
        branch_length = 0.1
        self.sim_on_branch1 = SimulateOnBranch(sequence1, branch_length)

        random.seed(9001)
        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        sequence2 = Sequence(s2, {'+0': [[(0, 12)]]}, kappa, mu, pi2, dN_values, dS_values)
        self.sim_on_branch2 = SimulateOnBranch(sequence2, branch_length)

        random.seed(9001)
        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [[(5, 50)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [[(30, 3)]]}
        pi3 = Sequence.get_frequency_rates(s3)
        sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, dN_values, dS_values)
        self.sim_on_branch3 = SimulateOnBranch(sequence3, branch_length)

        random.seed(9001)
        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [[(0, 12)]], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, dN_values, dS_values)
        self.sim_on_branch4 = SimulateOnBranch(sequence4, branch_length)

    def testGetSubstitution(self):

        random.seed(9001)       # Set seed for pseudo-random number generator
        # GTA CGA TCG ATC GAT GCT AGC
        # Expected GCT --> GAT substitution
        expected = (self.sim_on_branch1.sequence.nt_sequence[16], 'A')
        result = self.sim_on_branch1.get_substitution()
        self.assertEqual(expected, result)

        # Expected TTT --> TCT substitution
        expected = (self.sim_on_branch2.sequence.nt_sequence[1], 'C')
        result = self.sim_on_branch2.get_substitution()
        self.assertEqual(expected, result)

        # AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC
        # Expected TGT --> TGC substitution
        expected = (self.sim_on_branch3.sequence.nt_sequence[29], 'C')
        result = self.sim_on_branch3.get_substitution()
        self.assertEqual(expected, result)

        # ATG ACG TGG TGA
        # expected TGA -> TTA substitution
        expected = (self.sim_on_branch4.sequence.nt_sequence[10], 'T')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

        # ATG ACG TGG TGA
        # expected ACG -> ACA
        expected = (self.sim_on_branch4.sequence.nt_sequence[5], 'A')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 7.018351252068179
        result = self.sim_on_branch1.sum_rates()
        self.assertEqual(expected, result)

        expected = 17.430491845377023
        result = self.sim_on_branch2.sum_rates()
        self.assertEqual(expected, result)

        expected = 16.793266966449128
        result = self.sim_on_branch3.sum_rates()
        self.assertEqual(expected, result)

        expected = 3.2295309842774507
        result = self.sim_on_branch4.sum_rates()
        self.assertEqual(expected, result)

    def testSelectWeightedValues(self):

        # Selecting the to_nucleotide for seq1
        random.seed(9001)
        seq1_dict = self.sim_on_branch1.event_tree['to_nt']
        seq1_num_events = self.sim_on_branch1.event_tree['total_events']

        res_seq1_to_nt = self.sim_on_branch1.select_weighted_values(seq1_dict, seq1_num_events,
                                                                    'events_for_nt', 'stationary_frequency')
        exp_seq1_to_nt = 'A'
        self.assertEqual(exp_seq1_to_nt, res_seq1_to_nt)

        # Selecting the possible from_nucleotide for seq1
        random.seed(9001)
        seq1_from_dict = self.sim_on_branch1.event_tree['to_nt'][exp_seq1_to_nt]['from_nt']
        events_for_nt = self.sim_on_branch1.event_tree['to_nt'][exp_seq1_to_nt]['events_for_nt']

        res_seq1_from_nt = self.sim_on_branch1.select_weighted_values(seq1_from_dict, events_for_nt,
                                                                      'number_of_events', 'kappa')
        exp_seq1_from_nt = 'T'
        self.assertEqual(exp_seq1_from_nt, res_seq1_from_nt)

        # Selecting the to_nucleotide for seq4
        random.seed(9001)
        seq4_dict = self.sim_on_branch4.event_tree['to_nt']
        seq4_num_events = self.sim_on_branch4.event_tree['total_events']

        res_seq4_to_nt = self.sim_on_branch4.select_weighted_values(seq4_dict, seq4_num_events,
                                                                    'events_for_nt', 'stationary_frequency')
        exp_seq4_to_nt = 'A'
        self.assertEqual(exp_seq4_to_nt, res_seq4_to_nt)

        # Selecting the possible from_nucleotide for seq1
        random.seed(9001)
        seq4_from_dict = self.sim_on_branch4.event_tree['to_nt'][exp_seq4_to_nt]['from_nt']
        events_for_nt = self.sim_on_branch4.event_tree['to_nt'][exp_seq4_to_nt]['events_for_nt']

        res_seq4_from_nt = self.sim_on_branch1.select_weighted_values(seq4_from_dict, events_for_nt,
                                                                      'number_of_events', 'kappa')
        exp_seq4_from_nt = 'G'
        self.assertEqual(exp_seq4_from_nt, res_seq4_from_nt)

    def testWeightedRandomChoice(self):

        # Select to_nucleotide for seq4
        # ATG ACG TGG TGA
        random.seed(9001)
        np.random.seed(9001)

        seq4_to_nt_sum = 0.19333333333333333
        seq4_to_nt = {'A': 0.08333333333333333,
                      'T': 0.08333333333333333,
                      'C': 0.026666666666666665,
                      'G': 0.0}
        exp_to_nt = 'A'
        seq4_result = self.sim_on_branch4.weighted_random_choice(seq4_to_nt, seq4_to_nt_sum)
        self.assertEqual(exp_to_nt, seq4_result)

        # Selecting the possible from_nucleotide for seq4
        seq4_from_nt_sum = 0.2488235294117647
        seq4_from_nt = {'A': None,
                        'T': 0.075,
                        'C': 0.15,
                        'G': 0.25}
        exp_from_nt = 'T'
        seq4_result = self.sim_on_branch4.weighted_random_choice(seq4_from_nt, seq4_from_nt_sum)
        self.assertEqual(exp_from_nt, seq4_result)

        # Selecting the final from_nucleotide for seq1
        seq4_candidate_nts = self.sim_on_branch4.event_tree['to_nt'][exp_to_nt]['from_nt'][exp_from_nt]['nts_in_subs']
        seq4_rates_list = [nt.mutation_rate for nt in seq4_candidate_nts]
        seq4_candidate_dict = dict(zip(seq4_candidate_nts, seq4_rates_list))

        exp_final_nt = 't9'
        result_final_nt = self.sim_on_branch1.weighted_random_choice(seq4_candidate_dict, sum(seq4_rates_list))
        self.assertEqual(repr(result_final_nt), exp_final_nt)

    def testRemoveNt(self):

        random.seed(9001)
        # Testing sequence 2
        # Nucleotides in the start codon (a0, t1, g2) do not appear in the event tree
        a3 = self.sim_on_branch4.sequence.nt_sequence[3]
        c4 = self.sim_on_branch4.sequence.nt_sequence[4]
        g5 = self.sim_on_branch4.sequence.nt_sequence[5]
        t6 = self.sim_on_branch4.sequence.nt_sequence[6]
        g7 = self.sim_on_branch4.sequence.nt_sequence[7]
        g8 = self.sim_on_branch4.sequence.nt_sequence[8]
        t9 = self.sim_on_branch4.sequence.nt_sequence[9]
        g10 = self.sim_on_branch4.sequence.nt_sequence[10]
        a11 = self.sim_on_branch4.sequence.nt_sequence[11]

        # Remove nucleotide all from the sequence 4 event tree
        self.sim_on_branch4.remove_nt(a11)

        expected_event_tree = \
            {'to_nt':
                {'A': {'stationary_frequency': 0.25,
                       'from_nt':
                           {'A': None,
                            'T': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6, t9]},
                                                'dS': {(0, 0, 0, 1): [t6],
                                                       (1, 0, 0, 0): [t9]}},
                                  'is_syn': [],
                                  'nts_in_subs': {t9: None, t6: None},
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
                                  'is_nonsyn': {'dN': {},
                                                'dS': {}},
                                  'is_syn': [g5],
                                  'nts_in_subs': {g5: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 1},
                 'T': {'stationary_frequency': 0.25,
                       'from_nt':
                           {'A': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3]},
                                                'dS': {(1, 0, 0, 0): [a3]}},
                                  'is_syn': [],
                                  'nts_in_subs': {a3: None},
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
                                  'nts_in_subs': {g7: None, g5: None, g10: None, g8: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 1},
                 'C': {'stationary_frequency': 0.08,
                       'from_nt':
                           {'A': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(1, 0, 0, 0): [a3]},
                                                'dS': {(0, 0, 1, 0): [a3]}},
                                  'is_syn': [],
                                  'nts_in_subs': {a3: None},
                                  'number_of_events': 0},
                            'T': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(0, 1, 0, 0): [t6],
                                                       (0, 0, 0, 1): [t9]},
                                                'dS': {(0, 1, 0, 0): [t6],
                                                       (1, 0, 0, 0): [t9]}},
                                  'is_syn': [],
                                  'nts_in_subs': {t9: None, t6: None},
                                  'number_of_events': 0},
                            'C': None,
                            'G': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [g7],
                                                       (0, 0, 0, 1): [g8, g10]},
                                                'dS': {(0, 1, 0, 0): [g7],
                                                       (0, 0, 1, 0): [g8, g10]}},
                                  'is_syn': [g5],
                                  'nts_in_subs': {g7: None, g5: None, g10: None, g8: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 1},
                 'G': {'stationary_frequency': 0.42,
                       'from_nt':
                           {'A': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [a3]},
                                                'dS': {(1, 0, 0, 0): [a3]}},
                                  'is_syn': [],
                                  'nts_in_subs': {a3: None},
                                  'number_of_events': 0},
                            'T': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6],
                                                       (0, 1, 0, 0): [t9]},
                                                'dS': {(0, 1, 0, 0): [t6],
                                                       (1, 0, 0, 0): [t9]}},
                                  'is_syn': [],
                                  'nts_in_subs': {t9: None, t6: None},
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

        result_event_tree = self.sim_on_branch4.event_tree
        self.assertEqual(expected_event_tree, result_event_tree)

        # Testing sequence 1
        g0 = self.sim_on_branch1.sequence.nt_sequence[0]
        t1 = self.sim_on_branch1.sequence.nt_sequence[1]
        a2 = self.sim_on_branch1.sequence.nt_sequence[2]
        c3 = self.sim_on_branch1.sequence.nt_sequence[3]
        g4 = self.sim_on_branch1.sequence.nt_sequence[4]
        a5 = self.sim_on_branch1.sequence.nt_sequence[5]
        t6 = self.sim_on_branch1.sequence.nt_sequence[6]
        c7 = self.sim_on_branch1.sequence.nt_sequence[7]
        g8 = self.sim_on_branch1.sequence.nt_sequence[8]
        a9 = self.sim_on_branch1.sequence.nt_sequence[9]
        t10 = self.sim_on_branch1.sequence.nt_sequence[10]
        c11 = self.sim_on_branch1.sequence.nt_sequence[11]
        g12 = self.sim_on_branch1.sequence.nt_sequence[12]
        a13 = self.sim_on_branch1.sequence.nt_sequence[13]
        t14 = self.sim_on_branch1.sequence.nt_sequence[14]
        g15 = self.sim_on_branch1.sequence.nt_sequence[15]
        c16 = self.sim_on_branch1.sequence.nt_sequence[16]
        t17 = self.sim_on_branch1.sequence.nt_sequence[17]
        a18 = self.sim_on_branch1.sequence.nt_sequence[18]
        g19 = self.sim_on_branch1.sequence.nt_sequence[19]
        c20 = self.sim_on_branch1.sequence.nt_sequence[20]

        # GTACGATCGATCGATGCTAGC
        # Remove nucleotide c7 from the sequence 1 event tree
        random.seed(9001)
        self.sim_on_branch1.remove_nt(c7)

        expected_event_tree = \
            {'to_nt':
                {'A': {'stationary_frequency': 0.24,
                       'from_nt':
                           {'A': None,
                            'T': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [t1, t10],
                                                       (0, 1, 0, 0): [t6, t14]},
                                                'dS': {(0, 0, 1, 0): [t1, t14],
                                                       (0, 0, 0, 1): [t6],
                                                       (1, 0, 0, 0): [t10]}},
                                  'is_syn': [t17],
                                  'nts_in_subs': {t17: None, t10: None, t6: None, t14: None, t1: None},
                                  'number_of_events': 1},
                            'C': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [c16],
                                                       (0, 1, 0, 0): [c20]},
                                                'dS': {(0, 0, 0, 1): [c16],
                                                       (1, 0, 0, 0): [c20]}},
                                  'is_syn': [c3, c11],
                                  'nts_in_subs': {c3: None, c20: None, c11: None, c16: None},
                                  'number_of_events': 2},
                            'G': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(1, 0, 0, 0): [g0],
                                                       (0, 1, 0, 0): [g4, g15],
                                                       (0, 0, 1, 0): [g12],
                                                       (0, 0, 0, 1): [g19]},
                                                'dS': {(0, 0, 1, 0): [g0, g15],
                                                       (0, 1, 0, 0): [g4, g19],
                                                       (0, 0, 0, 1): [g12]}},
                                  'is_syn': [g8],
                                  'nts_in_subs': {g19: None, g4: None, g12: None, g0: None, g15: None, g8: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 4},
                 'T': {'stationary_frequency': 0.24,
                       'from_nt':
                           {'A': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(1, 0, 0, 0): [a9, a13, a18]},
                                                'dS': {(0, 0, 0, 1): [a9],
                                                       (1, 0, 0, 0): [a13],
                                                       (0, 0, 1, 0): [a18]}},
                                  'is_syn': [a2, a5],
                                  'nts_in_subs': {a2: None, a18: None, a9: None, a5: None, a13: None},
                                  'number_of_events': 2},
                            'T': None,
                            'C': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [c16]},
                                                'dS': {(0, 1, 0, 0): [c16]}},
                                  'is_syn': [c11, c20],
                                  'nts_in_subs': {c20: None, c11: None, c16: None},
                                  'number_of_events': 2},
                            'G': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [g0, g4],
                                                       (0, 0, 0, 1): [g12],
                                                       (1, 0, 0, 0): [g15],
                                                       (0, 1, 0, 0): [g19]},
                                                'dS': {(1, 0, 0, 0): [g0, g19],
                                                       (0, 1, 0, 0): [g4],
                                                       (0, 0, 0, 1): [g12, g15]}},
                                  'is_syn': [g8],
                                  'nts_in_subs': {g19: None, g4: None, g12: None, g0: None, g15: None, g8: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 5},
                 'C': {'stationary_frequency': 0.24,
                       'from_nt':
                           {'A': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 1, 0, 0): [a9, a18],
                                                       (1, 0, 0, 0): [a13]},
                                                'dS': {(1, 0, 0, 0): [a9],
                                                       (0, 1, 0, 0): [a13, a18]}},
                                  'is_syn': [a2, a5],
                                  'nts_in_subs': {a2: None, a18: None, a9: None, a5: None, a13: None},
                                  'number_of_events': 2},
                            'T': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [t1],
                                                       (0, 0, 0, 1): [t6],
                                                       (0, 1, 0, 0): [t10]},
                                                'dS': {(1, 0, 0, 0): [t1],
                                                       (0, 0, 1, 0): [t6, t10]}},
                                  'is_syn': [t14, t17],
                                  'nts_in_subs': {t17: None, t10: None, t6: None, t14: None, t1: None},
                                  'number_of_events': 2},
                            'C': None,
                            'G': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 1, 0): [g0, g19],
                                                       (0, 0, 0, 1): [g4],
                                                       (1, 0, 0, 0): [g12],
                                                       (0, 1, 0, 0): [g15]},
                                                'dS': {(1, 0, 0, 0): [g0],
                                                       (0, 1, 0, 0): [g4, g19],
                                                       (0, 0, 1, 0): [g12, g15]}},
                                  'is_syn': [g8],
                                  'nts_in_subs': {g19: None, g4: None, g12: None, g0: None, g15: None, g8: None},
                                  'number_of_events': 1}},
                       'events_for_nt': 5},
                 'G': {'stationary_frequency': 0.29,
                       'from_nt':
                           {'A': {'is_trv': False,
                                  'kappa': 1,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [a9, a13],
                                                       (0, 1, 0, 0): [a18]},
                                                'dS': {(0, 0, 1, 0): [a9],
                                                       (1, 0, 0, 0): [a13],
                                                       (0, 1, 0, 0): [a18]}},
                                  'is_syn': [a2, a5],
                                  'nts_in_subs': {a2: None, a18: None, a9: None, a5: None, a13: None},
                                  'number_of_events': 2},
                            'T': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [t1],
                                                       (0, 1, 0, 0): [t6],
                                                       (0, 0, 1, 0): [t10, t14]},
                                                'dS': {(0, 1, 0, 0): [t1],
                                                       (0, 0, 1, 0): [t6, t10],
                                                       (0, 0, 0, 1): [t14]}},
                                  'is_syn': [t17],
                                  'nts_in_subs': {t17: None, t10: None, t6: None, t14: None, t1: None},
                                  'number_of_events': 1},
                            'C': {'is_trv': True,
                                  'kappa': 0.3,
                                  'is_nonsyn': {'dN': {(0, 0, 0, 1): [c3],
                                                       (0, 0, 1, 0): [c11],
                                                       (0, 1, 0, 0): [c16],
                                                       (1, 0, 0, 0): [c20]},
                                                'dS': {(0, 0, 0, 1): [c3],
                                                       (1, 0, 0, 0): [c11],
                                                       (0, 1, 0, 0): [c16, c20]}},
                                  'is_syn': [],
                                  'nts_in_subs': {c3: None, c20: None, c11: None, c16: None},
                                  'number_of_events': 0},
                            'G': None},
                       'events_for_nt': 3}},
             'total_events': 17}

        result_event_tree = self.sim_on_branch1.event_tree
        self.assertEqual(expected_event_tree, result_event_tree)

    def testUpdateNt(self):

        random.seed(9001)
        # Testing sequence 4: ATG ACG TGG TGA
        # Mutate the G in the 5th position (ACG --> ACA)
        nt_to_mutate = self.sim_on_branch4.sequence.nt_sequence[5]
        selected_mutation = self.sim_on_branch4.get_substitution()
        new_state = str(selected_mutation[1])       # A
        self.sim_on_branch4.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': None, 'C': 0.0375, 'G': 0.125, 'T': 0.0375}

        self.assertEqual('A', nt_to_mutate.state)
        self.assertEqual('T', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testMutateOnBranch(self):
        np.random.seed(9001)    # Used to draw waiting times
        random.seed(9001)

        # Original sequence = ATGACGTGGTGA
        # ATG --> ATA substitution (g5 replaced with a5)
        self.sim_on_branch4.mutate_on_branch()
        exp_seq = 'ATGACATGGTGA'
        exp_res = ''.join(str(nt) for nt in self.sim_on_branch4.sequence.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        # a0, t1, g2 are not in event tree
        a3 = self.sim_on_branch4.sequence.nt_sequence[3]
        c4 = self.sim_on_branch4.sequence.nt_sequence[4]
        a5 = self.sim_on_branch4.sequence.nt_sequence[5]
        t6 = self.sim_on_branch4.sequence.nt_sequence[6]
        g7 = self.sim_on_branch4.sequence.nt_sequence[7]
        g8 = self.sim_on_branch4.sequence.nt_sequence[8]
        t9 = self.sim_on_branch4.sequence.nt_sequence[9]
        g10 = self.sim_on_branch4.sequence.nt_sequence[10]
        a11 = self.sim_on_branch4.sequence.nt_sequence[11]

        exp_event_tree = \
            {'to_nt': {'A': {'stationary_frequency': 0.25,
                             'from_nt':
                                 {'A': None,
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
                                        'is_nonsyn': {'dN': {(0, 1, 0, 0): [c4]},
                                                      'dS': {(0, 0, 0, 1): [c4]}},
                                        'is_syn': [],
                                        'nts_in_subs': {c4: None},
                                        'number_of_events': 0},
                                  'G': {'is_trv': False,
                                        'kappa': 1,
                                        'is_nonsyn': {'dN': {},
                                                      'dS': {}},
                                        'is_syn': [],
                                        'nts_in_subs': {},
                                        'number_of_events': 0}},
                             'events_for_nt': 0},
                       'T': {'stationary_frequency': 0.25,
                             'from_nt':
                                 {'A': {'is_trv': True,
                                        'kappa': 0.3,
                                        'is_nonsyn': {'dN': {(0, 0, 1, 0): [a11],
                                                             (1, 0, 0, 0): [a3]},
                                                      'dS': {(0, 0, 1, 0): [a11],
                                                             (0, 0, 0, 1): [a3]}},
                                        'is_syn': [a5],
                                        'nts_in_subs': {a11: None, a5: None, a3: None},
                                        'number_of_events': 1},
                                  'T': None,
                                  'C': {'is_trv': False,
                                        'kappa': 1,
                                        'is_nonsyn': {'dN': {(0, 1, 0, 0): [c4]},
                                                      'dS': {(0, 0, 0, 1): [c4]}},
                                        'is_syn': [],
                                        'nts_in_subs': {c4: None},
                                        'number_of_events': 0},
                                  'G': {'is_trv': True,
                                        'kappa': 0.3,
                                        'is_nonsyn': {'dN': {(0, 1, 0, 0): [g7, g8],
                                                             (1, 0, 0, 0): [g10]},
                                                      'dS': {(0, 0, 0, 1): [g7, g10],
                                                             (0, 0, 1, 0): [g8]}},
                                        'is_syn': [],
                                        'nts_in_subs': {g7: None, g10: None, g8: None},
                                        'number_of_events': 0}},
                             'events_for_nt': 1},
                       'C': {'stationary_frequency': 0.08,
                             'from_nt':
                                 {'A': {'is_trv': True,
                                        'kappa': 0.3,
                                        'is_nonsyn': {'dN': {(0, 0, 1, 0): [a11],
                                                             (1, 0, 0, 0): [a3]},
                                                      'dS': {(1, 0, 0, 0): [a11],
                                                             (0, 0, 1, 0): [a3]}},
                                        'is_syn': [a5],
                                        'nts_in_subs': {a11: None, a5: None, a3: None},
                                        'number_of_events': 1},
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
                                        'is_syn': [],
                                        'nts_in_subs': {g7: None, g10: None, g8: None},
                                        'number_of_events': 0}},
                             'events_for_nt': 1},
                       'G': {'stationary_frequency': 0.42,
                             'from_nt':
                                 {'A': {'is_trv': False,
                                        'kappa': 1,
                                        'is_nonsyn': {'dN': {(0, 1, 0, 0): [a11],
                                                             (0, 0, 1, 0): [a3]},
                                                      'dS': {(0, 0, 1, 0): [a11, a3]}},
                                        'is_syn': [a5],
                                        'nts_in_subs': {a11: None, a5: None, a3: None},
                                        'number_of_events': 1},
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
                                        'is_nonsyn': {'dN': {(0, 0, 0, 1): [c4]},
                                                      'dS': {(0, 1, 0, 0): [c4]}},
                                        'is_syn': [],
                                        'nts_in_subs': {c4: None},
                                        'number_of_events': 0},
                                  'G': None},
                             'events_for_nt': 1}},
             'total_events': 3}

        res_event_tree = self.sim_on_branch4.event_tree
        self.assertEqual(exp_event_tree, res_event_tree)

    def testUpdateNucleotideOnTree(self):

        # ATG ACG TGG TGA
        np.random.seed(9001)
        random.seed(9001)

        a3 = self.sim_on_branch4.sequence.nt_sequence[3]
        c4 = self.sim_on_branch4.sequence.nt_sequence[4]
        g5 = self.sim_on_branch4.sequence.nt_sequence[5]
        t6 = self.sim_on_branch4.sequence.nt_sequence[6]
        g7 = self.sim_on_branch4.sequence.nt_sequence[7]
        g8 = self.sim_on_branch4.sequence.nt_sequence[8]
        t9 = self.sim_on_branch4.sequence.nt_sequence[9]
        g10 = self.sim_on_branch4.sequence.nt_sequence[10]
        a11 = self.sim_on_branch4.sequence.nt_sequence[11]

        # Change nucleotide 'A' at position 3 to T
        # a3 should no longer appear in the tree, while t3 should appear in the 'from_nt'-T branches of the tree
        self.sim_on_branch4.remove_nt(a3)
        self.sim_on_branch4.update_nucleotide(a3, 'T')
        self.sim_on_branch4.update_nt_on_tree(a3, 'T')
        t3 = self.sim_on_branch4.sequence.nt_sequence[3]  # Get mutated nucleotide
        self.sim_on_branch4.sequence.get_nts_on_tips()

        expected_event_tree = \
            {'to_nt': {'A': {'stationary_frequency': 0.25,
                             'from_nt': {'A': None,
                                         'T': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6, t9],
                                                                    (1, 0, 0, 0): [t3]},
                                                             'dS': {(0, 0, 0, 1): [t6],
                                                                    (1, 0, 0, 0): [t9],
                                                                    (0, 0, 1, 0): [t3]}},
                                               'is_syn': [],
                                               'nts_in_subs': {t9: None, t3: None, t6: None},
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
                                               'is_nonsyn': {'dN': {},
                                                             'dS': {}},
                                               'is_syn': [g5],
                                               'nts_in_subs': {g5: None},
                                               'number_of_events': 1}},
                             'events_for_nt': 1},
                       'T': {'stationary_frequency': 0.25,
                             'from_nt': {'A': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [a11]},
                                                             'dS': {(0, 0, 1, 0): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': {a11: None},
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
                                               'nts_in_subs': {g7: None, g5: None, g10: None, g8: None},
                                               'number_of_events': 1}},
                             'events_for_nt': 1},
                       'C': {'stationary_frequency': 0.08,
                             'from_nt': {'A': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [a11]},
                                                             'dS': {(1, 0, 0, 0): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': {a11: None},
                                               'number_of_events': 0},
                                         'T': {'is_trv': False,
                                               'kappa': 1,
                                               'is_nonsyn': {'dN': {(0, 1, 0, 0): [t6],
                                                                    (0, 0, 0, 1): [t9],
                                                                    (0, 0, 1, 0): [t3]},
                                                             'dS': {(0, 1, 0, 0): [t6],
                                                                    (1, 0, 0, 0): [t9, t3]}},
                                               'is_syn': [],
                                               'nts_in_subs': {t9: None, t3: None, t6: None},
                                               'number_of_events': 0},
                                         'C': None,
                                         'G': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 1, 0): [g7],
                                                                    (0, 0, 0, 1): [g8, g10]},
                                                             'dS': {(0, 1, 0, 0): [g7],
                                                                    (0, 0, 1, 0): [g8, g10]}},
                                               'is_syn': [g5],
                                               'nts_in_subs': {g7: None, g5: None, g10: None, g8: None},
                                               'number_of_events': 1}},
                             'events_for_nt': 1},
                       'G': {'stationary_frequency': 0.42,
                             'from_nt': {'A': {'is_trv': False,
                                               'kappa': 1,
                                               'is_nonsyn': {'dN': {(0, 1, 0, 0): [a11]},
                                                             'dS': {(0, 0, 1, 0): [a11]}},
                                               'is_syn': [],
                                               'nts_in_subs': {a11: None},
                                               'number_of_events': 0},
                                         'T': {'is_trv': True,
                                               'kappa': 0.3,
                                               'is_nonsyn': {'dN': {(0, 0, 0, 1): [t6],
                                                                    (0, 1, 0, 0): [t9],
                                                                    (0, 0, 1, 0): [t3]},
                                                             'dS': {(0, 1, 0, 0): [t6],
                                                                    (1, 0, 0, 0): [t9, t3]}},
                                               'is_syn': [],
                                               'nts_in_subs': {t9: None, t3: None, t6: None},
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

        result = self.sim_on_branch4.event_tree
        self.assertEqual(expected_event_tree, result)


# ==========================================
# Tests for SimulateOnTree
# ==========================================
class TestSimulateOnTree(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        random.seed(9001)  # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        dN_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        dS_values = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [[(0, 21)]]}, kappa, mu, pi1, dN_values, dS_values)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree1 = SimulateOnTree(sequence1, phylo_tree)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree1.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree1.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree1.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)


if __name__ == '__main__':
    unittest.main()
