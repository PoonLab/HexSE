import unittest
import random
import os
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
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas)
        branch_length = 0.1
        self.sim_on_branch1 = SimulateOnBranch(sequence1, branch_length)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        sequence2 = Sequence(s2, {'+0': [(0, 12)]}, kappa, mu, pi2, omegas)
        self.sim_on_branch2 = SimulateOnBranch(sequence2, branch_length)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 50)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(30, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas)
        self.sim_on_branch3 = SimulateOnBranch(sequence3, branch_length)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 12)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)
        self.sim_on_branch4 = SimulateOnBranch(sequence4, branch_length)

    def testGetSubstitution(self):
        random.seed(9001)       # Set seed for pseudo-random number generator

        # GTA CGA TCG ATC GAT GCT AGC
        # Expected GCT --> GAT substitution
        expected = (self.sim_on_branch1.sequence.nt_sequence[16], 'A')
        result = self.sim_on_branch1.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch2.sequence.nt_sequence[1], 'C')
        result = self.sim_on_branch2.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch3.sequence.nt_sequence[29], 'C')
        result = self.sim_on_branch3.get_substitution()
        self.assertEqual(expected, result)

        # ATG ACG TGG TGA
        # expected TGG -> TTG substitution
        expected = (self.sim_on_branch4.sequence.nt_sequence[7], 'T')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 4.25294227720946
        result = self.sim_on_branch1.sum_rates()
        self.assertEqual(expected, result)

        expected = 10.771830262362425
        result = self.sim_on_branch2.sum_rates()
        self.assertEqual(expected, result)

        expected = 12.48467014546256
        result = self.sim_on_branch3.sum_rates()
        self.assertEqual(expected, result)

        expected = 2.716842288378799
        result = self.sim_on_branch4.sum_rates()
        self.assertEqual(expected, result)

    def testSelectWeightedValues(self):

        # Selecting the to_nucleotide for seq1
        random.seed(9001)
        seq1_dict = self.sim_on_branch1.event_tree['to_nt']
        seq1_num_events = self.sim_on_branch1.event_tree['total_events']

        seq1_to_result = self.sim_on_branch1.select_weighted_values(seq1_dict, seq1_num_events,
                                                                    'events_for_nt', 'stationary_frequency')
        expected_seq1_to = 'A'
        self.assertEqual(expected_seq1_to, seq1_to_result)

        # Selecting the possible from_nucleotide for seq1
        random.seed(9001)
        seq1_from_dict = self.sim_on_branch1.event_tree['to_nt'][expected_seq1_to]['from_nt']
        events_for_nt = self.sim_on_branch1.event_tree['to_nt'][expected_seq1_to]['events_for_nt']

        seq1_from_result = self.sim_on_branch1.select_weighted_values(seq1_from_dict, events_for_nt,
                                                                      'number_of_events', 'kappa')
        expected_seq1_from = 'T'
        self.assertEqual(expected_seq1_from, seq1_from_result)

    def testWeightedRandomChoice(self):

        # Select to_nucleotide for seq1
        random.seed(9001)
        seq1_to_nt_sum = 0.2488235294117647
        seq1_to_nt = {'A': 0.05647058823529411,
                      'T': 0.07058823529411765,
                      'C': 0.07058823529411765,
                      'G': 0.051176470588235295}
        exp_to_nt = 'A'
        seq1_result = self.sim_on_branch1.weighted_random_choice(seq1_to_nt, seq1_to_nt_sum)
        self.assertEqual(exp_to_nt, seq1_result)

        # Selecting the possible from_nucleotide for seq1
        random.seed(9001)
        seq1_from_nt_sum = 0.2488235294117647
        seq1_from_nt = {'A': None,
                        'T': 0.075,
                        'C': 0.15,
                        'G': 0.25}
        exp_from_nt = 'T'
        seq1_result = self.sim_on_branch1.weighted_random_choice(seq1_from_nt, seq1_from_nt_sum)
        self.assertEqual(exp_from_nt, seq1_result)

        # Selecting the final from_nucleotide for seq1
        random.seed(9001)
        seq1_candidate_nts = self.sim_on_branch1.event_tree['to_nt'][exp_to_nt]['from_nt'][exp_from_nt]['nts_in_subs']
        seq1_rates_list = [nt.mutation_rate for nt in seq1_candidate_nts]
        seq1_candidate_dict = dict(zip(seq1_candidate_nts, seq1_rates_list))

        expected_final_nt = 't17'
        result_final_nt = self.sim_on_branch1.weighted_random_choice(seq1_candidate_dict, sum(seq1_rates_list))
        self.assertEqual(expected_final_nt, repr(result_final_nt))

    def testRemoveNt(self):

        # Testing sequence 2
        a0 = self.sim_on_branch4.sequence.nt_sequence[0]
        t1 = self.sim_on_branch4.sequence.nt_sequence[1]
        g2 = self.sim_on_branch4.sequence.nt_sequence[2]
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
        random.seed(9001)
        self.sim_on_branch4.remove_nt(a11)

        expected_event_tree = \
            {'to_nt': {'A': {'events_for_nt': 1,
                             'from_nt': {
                                 'A': None,
                                 'C': {'is_nonsyn': {(1, 0, 0, 0): [c4]},
                                       'is_syn': [],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [c4],
                                       'number_of_events': 0},
                                 'G': {'is_nonsyn': {(0, 0, 1, 0): [g2]},
                                       'is_syn': [g5],
                                       'is_trv': False,
                                       'kappa': 1,
                                       'nts_in_subs': [g5, g2],
                                       'number_of_events': 1},
                                 'T': {'is_nonsyn': {(0, 0, 0, 1): [t9],
                                                     (0, 0, 1, 0): [t1],
                                                     (0, 1, 0, 0): [t6]},
                                       'is_syn': [],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [t1, t6, t9],
                                       'number_of_events': 0}},
                             'stationary_frequency': 0.25},

                       'C': {'events_for_nt': 1,
                             'from_nt': {'A': {'is_nonsyn': {  # (0, 0, 0, 1): [],      # a11 was here
                                                             (0, 0, 1, 0): [a3],
                                                             (1, 0, 0, 0): [a0]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [a0, a3],
                                               'number_of_events': 0},
                                         'C': None,
                                         'G': {'is_nonsyn': {(0, 0, 0, 1): [g2, g8],
                                                             (0, 0, 1, 0): [g7],
                                                             (0, 1, 0, 0): [g10]},
                                               'is_syn': [g5],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [g5, g2, g8, g7, g10],
                                               'number_of_events': 1},
                                         'T': {'is_nonsyn': {(0, 0, 1, 0): [t1, t9],
                                                             (0, 1, 0, 0): [t6]},
                                               'is_syn': [],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [t1, t9, t6],
                                               'number_of_events': 0}},
                             'stationary_frequency': 0.08},

                       'G': {'events_for_nt': 0,
                             'from_nt': {'A': {'is_nonsyn': {(0, 0, 0, 1): [a0],
                                                             (0, 0, 1, 0): [a3]},
                                                             #  (1, 0, 0, 0): []},     # a11 was here
                                               'is_syn': [],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [a0, a3],
                                               'number_of_events': 0},
                                         'C': {'is_nonsyn': {(0, 0, 1, 0): [c4]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': None,
                                         'T': {'is_nonsyn': {(0, 0, 1, 0): [t6, t9],
                                                             (0, 1, 0, 0): [t1]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [t1, t6, t9],
                                               'number_of_events': 0}},
                             'stationary_frequency': 0.42},

                       'T': {'events_for_nt': 1,
                             'from_nt': {'A': {'is_nonsyn': {  # (0, 0, 0, 1): [],   a11 was here
                                                             (0, 1, 0, 0): [a0, a3]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [a0, a3],
                                               'number_of_events': 0},
                                         'C': {'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                               'is_syn': [],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': {'is_nonsyn': {(0, 0, 0, 1): [g2, g8, g10],
                                                             (0, 0, 1, 0): [g7]},
                                               'is_syn': [g5],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [g5, g2, g8, g10, g7],
                                               'number_of_events': 1},
                                         'T': None},
                             'stationary_frequency': 0.25}},
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
            {'to_nt': {'A': {'events_for_nt': 4,
                             'from_nt': {
                                 'A': None,
                                 'C': {'is_nonsyn': {(0, 0, 1, 0): [c16, c20]},
                                       'is_syn': [c3, c11],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [c3, c11, c16, c20],
                                       'number_of_events': 2},
                                 'G': {'is_nonsyn': {(0, 0, 0, 1): [g12],
                                                     (0, 0, 1, 0): [g4, g19],
                                                     (0, 1, 0, 0): [g15],
                                                     (1, 0, 0, 0): [g0]},
                                       'is_syn': [g8],
                                       'is_trv': False,
                                       'kappa': 1,
                                       'nts_in_subs': [g8, g0, g4, g19, g12, g15],
                                       'number_of_events': 1},
                                 'T': {'is_nonsyn': {(0, 0, 0, 1): [t6, t14],
                                                     (0, 0, 1, 0): [t10],
                                                     (1, 0, 0, 0): [t1]},
                                       'is_syn': [t17],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [t17, t1, t6, t14, t10],
                                       'number_of_events': 1}},
                                 'stationary_frequency': 0.24},
                       'C': {'events_for_nt': 5,
                             'from_nt': {
                                 'A': {'is_nonsyn': {(0, 0, 1, 0): [a13, a18],
                                                     (0, 1, 0, 0): [a9]},
                                       'is_syn': [a2, a5],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [a2, a5, a9, a13, a18],
                                       'number_of_events': 2},
                                 'C': None,
                                 'G': {'is_nonsyn': {(0, 0, 1, 0): [g0, g4, g12, g19],
                                                     (1, 0, 0, 0): [g15]},
                                       'is_syn': [g8],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [g8, g0, g4, g12, g19, g15],
                                       'number_of_events': 1},
                                 'T': {'is_nonsyn': {(0, 0, 1, 0): [t1],
                                                     (0, 1, 0, 0): [t6, t10]},
                                       'is_syn': [t14, t17],
                                       'is_trv': False,
                                       'kappa': 1,
                                       'nts_in_subs': [t14, t17, t1, t6, t10],
                                       'number_of_events': 2}},
                                 'stationary_frequency': 0.24},
                       'G': {'events_for_nt': 3,
                             'from_nt': {
                                 'A': {'is_nonsyn': {(0, 0, 0, 1): [a9, a13],
                                                     (1, 0, 0, 0): [a18]},
                                       'is_syn': [a2, a5],
                                       'is_trv': False,
                                       'kappa': 1,
                                       'nts_in_subs': [a2, a5, a9, a13, a18],
                                       'number_of_events': 2},
                                 'C': {'is_nonsyn': {(0, 0, 0, 1): [c11],
                                                     (0, 0, 1, 0): [c3],
                                                     (1, 0, 0, 0): [c16, c20]},
                                       'is_syn': [],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [c3, c11, c16, c20],
                                       'number_of_events': 0},
                                 'G': None,
                                 'T': {'is_nonsyn': {(0, 0, 0, 1): [t6],
                                                     (0, 1, 0, 0): [t10],
                                                     (1, 0, 0, 0): [t1, t14]},
                                       'is_syn': [t17],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [t17, t1, t14, t6, t10],
                                       'number_of_events': 1}},
                                 'stationary_frequency': 0.29},
                       'T': {'events_for_nt': 5,
                             'from_nt': {
                                 'A': {'is_nonsyn': {(0, 1, 0, 0): [a9, a18],
                                                     (1, 0, 0, 0): [a13]},
                                       'is_syn': [a2, a5],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [a2, a5, a9, a18, a13],
                                       'number_of_events': 2},
                                 'C': {'is_nonsyn': {(0, 0, 0, 1): [c16]},
                                                     # (0, 1, 0, 0): []},
                                       'is_syn': [c11, c20],
                                       'is_trv': False,
                                       'kappa': 1,
                                       'nts_in_subs': [c11, c20, c16],
                                       'number_of_events': 2},
                                 'G': {'is_nonsyn': {(0, 0, 0, 1): [g15],
                                                     (0, 0, 1, 0): [g0, g19],
                                                     (0, 1, 0, 0): [g12],
                                                     (1, 0, 0, 0): [g4]},
                                       'is_syn': [g8],
                                       'is_trv': True,
                                       'kappa': 0.3,
                                       'nts_in_subs': [g8, g0, g19, g4, g12, g15],
                                       'number_of_events': 1},
                                 'T': None},
                                 'stationary_frequency': 0.24}},
             'total_events': 17}

        result_event_tree = self.sim_on_branch1.event_tree
        print(result_event_tree)
        self.assertEqual(expected_event_tree, result_event_tree)

    def testUpdateNt(self):

        # Testing sequence 4: ATGACGTGGTGA
        random.seed(9001)

        # Mutate the G in the 5th position
        nt_to_mutate = self.sim_on_branch4.sequence.nt_sequence[5]
        selected_mutation = self.sim_on_branch4.get_substitution()
        new_state = str(selected_mutation[1])
        self.sim_on_branch4.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': None,
                          'C': 0.0375,
                          'G': 0.125,
                          'T': 0.0375}
        expected_omegas = {'A': None, 'C': (0, 0, 0, 0), 'G': (0, 0, 0, 0), 'T': (0, 0, 0, 0)}

        self.assertEqual('A', nt_to_mutate.state)
        self.assertEqual('T', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)
        self.assertEqual(expected_omegas, nt_to_mutate.my_omegas)

    def testMutateOnBranch(self):

        np.random.seed(9001)    # Used to draw waiting times
        random.seed(9001)

        # Original sequence = ATGACGTGGTGA
        expected = 'ATAACGTGGTGA'
        self.sim_on_branch4.mutate_on_branch()
        result = ''.join(str(nt) for nt in self.sim_on_branch4.sequence.nt_sequence)

        self.assertEqual(expected, result)

    def testUpdateNucleotideOnTree(self):

        np.random.seed(9001)
        random.seed(9001)

        a0 = self.sim_on_branch4.sequence.nt_sequence[0]
        t1 = self.sim_on_branch4.sequence.nt_sequence[1]
        g2 = self.sim_on_branch4.sequence.nt_sequence[2]
        a3 = self.sim_on_branch4.sequence.nt_sequence[3]
        c4 = self.sim_on_branch4.sequence.nt_sequence[4]
        g5 = self.sim_on_branch4.sequence.nt_sequence[5]
        t6 = self.sim_on_branch4.sequence.nt_sequence[6]
        g7 = self.sim_on_branch4.sequence.nt_sequence[7]
        g8 = self.sim_on_branch4.sequence.nt_sequence[8]
        t9 = self.sim_on_branch4.sequence.nt_sequence[9]
        g10 = self.sim_on_branch4.sequence.nt_sequence[10]
        a11 = self.sim_on_branch4.sequence.nt_sequence[11]

        # Change nucleotide 'A' at position 0 to T
        # a0 should no longer appear in the tree, while t0 should appear in the 'from_nt'-T branches of the tree
        # TODO: verify synonymous mutations count 
        self.sim_on_branch4.update_nucleotide(a0, 'T')
        self.sim_on_branch4.update_nt_on_tree(a0, 'T')

        t0 = self.sim_on_branch4.sequence.nt_sequence[0]  # Get mutated nucleotide

        expected_event_tree = \
            {'to_nt': {'A': {'events_for_nt': 1,
                             'from_nt': {'A': None,
                                         'C': {'is_nonsyn': {(1, 0, 0, 0): [c4]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': {'is_nonsyn': {(0, 0, 1, 0): [g2]},
                                               'is_syn': [g5],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [g5, g2],
                                               'number_of_events': 1},
                                         'T': {'is_nonsyn': {(0, 0, 0, 1): [t9],
                                                             (0, 0, 1, 0): [t1],
                                                             (0, 1, 0, 0): [t6],
                                                             (1, 0, 0, 0): [t0]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [t1, t6, t9, t0],
                                               'number_of_events': 0}},
                             'stationary_frequency': 0.25},
                       'C': {'events_for_nt': 1,
                             'from_nt': {'A': {'is_nonsyn': {(0, 0, 0, 1): [a11],
                                                             (0, 0, 1, 0): [a3]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [a3, a11],
                                               'number_of_events': 0},
                                         'C': None,
                                         'G': {'is_nonsyn': {(0, 0, 0, 1): [g2, g8],
                                                             (0, 0, 1, 0): [g7],
                                                             (0, 1, 0, 0): [g10]},
                                               'is_syn': [g5],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [g5, g2, g8, g7, g10],
                                               'number_of_events': 1},
                                         'T': {'is_nonsyn': {(0, 0, 1, 0): [t1, t9],
                                                             (0, 1, 0, 0): [t6]},
                                               'is_syn': [t0],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [t1, t9, t6, t0],
                                               'number_of_events': 1}},
                             'stationary_frequency': 0.08},
                       'G': {'events_for_nt': 0,
                             'from_nt': {'A': {'is_nonsyn': {(0, 0, 1, 0): [a3],
                                                             (1, 0, 0, 0): [a11]},
                                               'is_syn': [],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [a3, a11],
                                               'number_of_events': 0},
                                         'C': {'is_nonsyn': {(0, 0, 1, 0): [c4]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': None,
                                         'T': {'is_nonsyn': {(0, 0, 1, 0): [t6, t9, t0],
                                                             (0, 1, 0, 0): [t1]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [t1, t6, t9, t0],
                                               'number_of_events': 0}},
                             'stationary_frequency': 0.42},
                       'T': {'events_for_nt': 1,
                             'from_nt': {'A': {'is_nonsyn': {(0, 0, 0, 1): [a11],
                                                             (0, 1, 0, 0): [a3]},
                                               'is_syn': [],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [a3, a11],
                                               'number_of_events': 0},
                                         'C': {'is_nonsyn': {(0, 0, 0, 1): [c4]},
                                               'is_syn': [],
                                               'is_trv': False,
                                               'kappa': 1,
                                               'nts_in_subs': [c4],
                                               'number_of_events': 0},
                                         'G': {'is_nonsyn': {(0, 0, 0, 1): [g2, g8, g10],
                                                             (0, 0, 1, 0): [g7]},
                                               'is_syn': [g5],
                                               'is_trv': True,
                                               'kappa': 0.3,
                                               'nts_in_subs': [g5, g2, g8, g10, g7],
                                               'number_of_events': 1},
                                         'T': None},
                             'stationary_frequency': 0.25}},
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
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [(0, 21)]}, kappa, mu, pi1, omegas)

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
