import unittest
import random
import os

from Bio import Phylo
from refact.new_sequence_info import Sequence
from refact.simulation_refactore import SimulateOnBranch
from refact.simulation_refactore import SimulateOnTree

TEST_TREE = os.path.abspath('../refact/test_tree.txt')


# ==========================================
# Tests for SimulateOnBranch
# ==========================================
class TestSimulateOnBranch(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        random.seed(9001)  # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)
        branch_length = 10
        self.sim_on_branch1 = SimulateOnBranch(sequence1, branch_length)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        sequence2 = Sequence(s2, {'+0': [(0, 11)]}, kappa, mu, pi2, omegas)
        self.sim_on_branch2 = SimulateOnBranch(sequence2, branch_length)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas)
        self.sim_on_branch3 = SimulateOnBranch(sequence3, branch_length)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)
        self.sim_on_branch4 = SimulateOnBranch(sequence4, branch_length)

    def testGetSubstitution(self):
        random.seed(9001)       # Set seed for pseudo-random number generator

        expected = (self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(7), 'A')
        result = self.sim_on_branch1.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch2.sequence.nt_sequence.nucleotide_at_pos(0), 'C')
        result = self.sim_on_branch2.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch3.sequence.nt_sequence.nucleotide_at_pos(6), 'C')
        result = self.sim_on_branch3.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch4.sequence.nt_sequence.nucleotide_at_pos(8), 'A')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.004637967436947257
        result = self.sim_on_branch1.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.00954128577225131
        result = self.sim_on_branch2.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.013399271035090315
        result = self.sim_on_branch3.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.0030337488758443035
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

        g0 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(0)
        t1 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(1)
        a2 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(2)
        c3 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(3)
        g4 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(4)
        a5 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(5)
        t6 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(6)
        c7 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(7)
        g8 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(8)
        a9 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(9)
        t10 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(10)
        c11 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(11)
        g12 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(12)
        a13 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(13)
        t14 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(14)
        g15 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(15)
        c16 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(16)
        t17 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(17)
        a18 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(18)
        g19 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(19)
        c20 = self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(20)

        # Remove nucleotide t17 from the sequence 1 event tree
        random.seed(9001)
        self.sim_on_branch1.remove_nt(t17)

        expected_event_tree = \
            {'to_nt': {
                'A': {'stationary_frequency': 0.24,
                      'from_nt': {
                          'A': None,
                          'T': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(1, 0, 0, 0): [t1],
                                            (0, 1, 0, 0): [t6, t10, t14]},
                              'is_syn': [],
                              'nts_in_subs': [t1, t6, t10, t14],
                              'number_of_events': 0},
                          'C': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 1, 0, 0): [c7],
                                            (0, 0, 0, 1): [c16],
                                            (0, 0, 1, 0): [c20]},
                              'is_syn': [c3, c11],
                              'nts_in_subs': [c3, c11, c7, c16, c20],
                              'number_of_events': 2},
                          'G': {
                              'is_trv': False,
                              'kappa': 1,
                              'is_nonsyn': {(1, 0, 0, 0): [g0],
                                            (0, 0, 1, 0): [g4, g19],
                                            (0, 1, 0, 0): [g12],
                                            (0, 0, 0, 1): [g15]},
                              'is_syn': [g8],
                              'nts_in_subs': [g8, g0, g4, g19, g12, g15],
                              'number_of_events': 1}},
                      'events_for_nt': 3},

                'T': {'stationary_frequency': 0.24,
                      'from_nt': {
                          'A': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 1, 0, 0): [a9],
                                            (1, 0, 0, 0): [a13],
                                            (0, 0, 1, 0): [a18]},
                              'is_syn': [a2, a5],
                              'nts_in_subs': [a2, a5, a9, a13, a18],
                              'number_of_events': 2},
                          'T': None,
                          'C': {
                              'is_trv': False,
                              'kappa': 1,
                              'is_nonsyn': {(0, 0, 1, 0): [c3],
                                            (0, 0, 0, 1): [c7],
                                            (1, 0, 0, 0): [c16]},
                              'is_syn': [c11, c20],
                              'nts_in_subs': [c11, c20, c3, c7, c16],
                              'number_of_events': 2},
                          'G': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 0, 1, 0): [g0],
                                            (0, 0, 0, 1): [g4, g12],
                                            (1, 0, 0, 0): [g15, g19]},
                              'is_syn': [g8],
                              'nts_in_subs': [g8, g0, g4, g12, g15, g19],
                              'number_of_events': 1}},
                      'events_for_nt': 5},

                'C': {'stationary_frequency': 0.24,
                      'from_nt': {
                          'A': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 1, 0, 0): [a9, a18],
                                            (1, 0, 0, 0): [a13]},
                              'is_syn': [a2, a5],
                              'nts_in_subs': [a2, a5, a9, a18, a13],
                              'number_of_events': 2},
                          'T': {
                              'is_trv': False,
                              'kappa': 1,
                              'is_nonsyn': {(0, 0, 1, 0): [t1],
                                            (0, 0, 0, 1): [t6, t10]},
                              'is_syn': [t14],
                              'nts_in_subs': [t14, t1, t6, t10],
                              'number_of_events': 1},
                          'C': None,
                          'G': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 0, 1, 0): [g0, g12, g15, g19],
                                            (1, 0, 0, 0): [g4]},
                              'is_syn': [g8],
                              'nts_in_subs': [g8, g0, g12, g15, g19, g4],
                              'number_of_events': 1}},
                      'events_for_nt': 4},

                'G': {'stationary_frequency': 0.29,
                      'from_nt': {
                          'A': {
                              'is_trv': False,
                              'kappa': 1,
                              'is_nonsyn': {(0, 0, 1, 0): [a9, a18],
                                            (0, 0, 0, 1): [a13]},
                              'is_syn': [a2, a5],
                              'nts_in_subs': [a2, a5, a9, a18, a13],
                              'number_of_events': 2},
                          'T': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(1, 0, 0, 0): [t1, t14],
                                            (0, 0, 0, 1): [t6, t10]},
                              'is_syn': [],
                              'nts_in_subs': [t1, t14, t6, t10],
                              'number_of_events': 0},
                          'C': {
                              'is_trv': True,
                              'kappa': 0.3,
                              'is_nonsyn': {(0, 0, 1, 0): [c3, c11, c16],
                                            (0, 1, 0, 0): [c7],
                                            (0, 0, 0, 1): [c20]},
                              'is_syn': [],
                              'nts_in_subs': [c3, c11, c16, c7, c20],
                              'number_of_events': 0},
                          'G': None},
                      'events_for_nt': 2}},
                'total_events': 14}

        result_event_tree = self.sim_on_branch1.event_tree

        self.assertEqual(expected_event_tree, result_event_tree)


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
        sequence1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)

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
