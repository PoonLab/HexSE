import os
import random
import unittest

import numpy as np
from Bio import Phylo

from hexse.sequence_info import Sequence
from hexse.simulation import SimulateOnBranch
from hexse.simulation import SimulateOnTree

TEST_TREE = os.path.join(os.path.dirname(__file__), 'fixtures/test_tree.txt')

KAPPA = 0.3
GLOBAL_RATE = 0.0005
CAT_VALUES = {'mu1': 0.051710707633483066, 'mu2': 0.15181054803756722, 'mu3': 0.26809045653750935,
              'mu4': 0.4186255904232465, 'mu5': 0.6442570794470408, 'mu6': 1.2255056178040284}
MAX_DIFF = None


class Seq1(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9001)  # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        orfs1 = {'+0': [{'coords': [[0, 21]],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi1 = Sequence.get_frequency_rates(s1)
        self.sequence1 = Sequence(s1, orfs1, KAPPA, GLOBAL_RATE, pi1, CAT_VALUES)

        branch_length = 0.1
        self.sim_on_branch1 = SimulateOnBranch(self.sequence1, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree1 = SimulateOnTree(self.sequence1, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence1.probability_tree['to_nt']['G']['from_nt']
        expected = 'A'
        result = self.sim_on_branch1.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence1.probability_tree['to_nt']['G']['from_nt']['A']['cat']
        expected = 'mu5'
        result = self.sim_on_branch1.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        random.seed(9001)  # Set seed for pseudo-random number generator
        # GTA CGA TCG ATC GAT GCT AGC
        # Expected GCT --> GAT substitution
        expected = (self.sim_on_branch1.sequence.nt_sequence[20], 'A')
        result = self.sim_on_branch1.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        # Select to_nucleotide
        random.seed(9001)
        np.random.seed(9001)

        seq1_to_nt_sum = 15.389999999999999
        seq1_to_nt = {'A': 3.5999999999999996, 'C': 3.84, 'T': 3.5999999999999996, 'G': 4.35}

        expected = 'A'
        result = self.sim_on_branch1.weighted_random_choice(seq1_to_nt, seq1_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.0017178571290980837
        result = self.sim_on_branch1.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch1.mutate_on_branch()
        exp_seq = 'GTACGATCGATCGATGCTAGC'
        exp_res = ''.join(str(nt) for nt in self.sequence1.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

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

        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

    def testRemoveNt(self):
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

        self.sim_on_branch1.remove_nt(g0)

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
                                                                                 ((1, 0, 0, 0),): []},
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
                                                      'G': {'category': {'mu1': {((0, 1, 0, 0),): []},
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
                                                                         'mu3': {((0, 0, 1, 0),): [g12,
                                                                                                   g15]},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [g8],
                                                                                 ((1, 0, 0, 0),): [g19]},
                                                                         'mu6': {((0, 1, 0, 0),): [g4]}}},
                                                      'T': None}}}}

        self.assertEqual(exp_event_tree, self.sequence1.event_tree)

    def testUpdateNt(self):
        random.seed(9001)
        nt_to_mutate = self.sequence1.nt_sequence[5]
        selected_mutation = self.sim_on_branch1.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch1.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': None,
                          'C': 4.4118202240945026e-05,
                          'G': 0.00014706067413648342,
                          'T': 9.651256435350337e-06}

        self.assertEqual('A', nt_to_mutate.state)
        self.assertEqual('T', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

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

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree1.traverse_tree()

        exp_sequences = ['GTACGATCGATCGATGCTAGC', 'GTACGATCGATCGATGCTAGC', 'GTACGATCGATCGATGCTAGC',
                         'GTACGATCGATCGATGCTAGC', 'GTACGATCGATCGATGCTAGC', 'GTACGATCGATCGATGCTAGC']

        res_sequences = []
        for clade in self.sim_on_tree1.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nGTACGATCGATCGATGCTAGC\n>B \nGTACGATCGATCGATGCTAGC\n>C \nGTACGATCGATCGATGCTAGC\n>D \nGTACGATCGATCGATGCTAGC\n'
        outfile = './seq1_output.txt'
        self.sim_on_tree1.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


class Seq2(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4002)

        s2 = 'ATGAATAAACCCGTATGA'
        orfs2 = {'+0': [{'coords': [[0, 18]],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [], '+2': [], '-0': [], '-1': [],
                 '-2': [{'coords': [[3, 15]],
                         'omega_classes': 4, 'omega_shape': 1.25,
                         'omega_values': [0.09199853806558903, 0.27043066909631136,
                                          0.5158061369385518, 1.1217646558655263]}]}
        pi2 = Sequence.get_frequency_rates(s2)
        self.sequence2 = Sequence(s2, orfs2, KAPPA, GLOBAL_RATE, pi2, CAT_VALUES)

        branch_length = 0.5
        self.sim_on_branch2 = SimulateOnBranch(self.sequence2, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree2 = SimulateOnTree(self.sequence2, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence2.probability_tree['to_nt']['T']['from_nt']
        expected = 'C'
        result = self.sim_on_branch2.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence2.probability_tree['to_nt']['G']['from_nt']['A']['cat']
        expected = 'mu4'
        result = self.sim_on_branch2.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        random.seed(4002)
        expected = (self.sim_on_branch2.sequence.nt_sequence[6], 'G')
        result = self.sim_on_branch2.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        random.seed(4002)
        np.random.seed(4002)

        seq2_to_nt_sum = 3.73
        seq2_to_nt = {'A': 1.32, 'C': 0.51, 'T': 0.88, 'G': 1.02}

        expected = 'G'
        result = self.sim_on_branch2.weighted_random_choice(seq2_to_nt, seq2_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.0001019775055701344
        result = self.sim_on_branch2.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch2.mutate_on_branch()
        exp_seq = 'ATGAATAAACCCGTATGA'
        exp_res = ''.join(str(nt) for nt in self.sequence2.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        a6 = self.sequence2.nt_sequence[6]
        a7 = self.sequence2.nt_sequence[7]
        a8 = self.sequence2.nt_sequence[8]
        c9 = self.sequence2.nt_sequence[9]
        c10 = self.sequence2.nt_sequence[10]
        c11 = self.sequence2.nt_sequence[11]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                             'C': {'category': {'mu1': {((1, 0, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                'mu2': {},
                                                'mu3': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): [c10]},
                                                'mu4': {},
                                                'mu5': {},
                                                'mu6': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): [c11]}}},
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
                                                'mu2': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a8]},
                                                'mu3': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                'mu4': {},
                                                'mu5': {},
                                                'mu6': {((1, 0, 0, 0), (0, 0, 1, 0, 0)): [a6]}}},
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
           'G': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [a8]},
                                                'mu2': {},
                                                'mu3': {},
                                                'mu4': {((1, 0, 0, 0), (0, 0, 0, 0, 1)): [a6],
                                                        ((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                'mu5': {},
                                                'mu6': {}}},
                             'C': {'category': {'mu1': {},
                                                'mu2': {},
                                                'mu3': {((0, 0, 0, 1), (1, 0, 0, 0, 0)): [c11]},
                                                'mu4': {((0, 0, 1, 0), (0, 0, 0, 0, 1)): [c9],
                                                        ((0, 1, 0, 0), (0, 0, 0, 1, 0)): [c10]},
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
                                                'mu5': {((1, 0, 0, 0), (0, 0, 0, 1, 0)): [a7]},
                                                'mu6': {}}},
                             'C': {'category': {'mu1': {},
                                                'mu2': {},
                                                'mu3': {((0, 1, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                'mu4': {},
                                                'mu5': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): [c11],
                                                        ((0, 1, 0, 0), (1, 0, 0, 0, 0)): [c10]},
                                                'mu6': {}}},
                             'G': {'category': {'mu1': {},
                                                'mu2': {},
                                                'mu3': {},
                                                'mu4': {},
                                                'mu5': {},
                                                'mu6': {}}},
                             'T': None}}}}

        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

    def testRemoveNt(self):
        random.seed(9001)

        a6 = self.sequence2.nt_sequence[6]
        a7 = self.sequence2.nt_sequence[7]
        a8 = self.sequence2.nt_sequence[8]
        c9 = self.sequence2.nt_sequence[9]
        c10 = self.sequence2.nt_sequence[10]
        c11 = self.sequence2.nt_sequence[11]

        self.sim_on_branch2.remove_nt(c11)

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {((1, 0, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): [c10]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): []}}},
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
                                                                         'mu2': {((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a8]},
                                                                         'mu3': {((0, 1, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((1, 0, 0, 0), (0, 0, 1, 0, 0)): [a6]}}},
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
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [a8]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0), (0, 0, 0, 0, 1)): [a6],
                                                                                 ((1, 0, 0, 0), (1, 0, 0, 0, 0)): [a7]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1), (1, 0, 0, 0, 0)): []},
                                                                         'mu4': {((0, 0, 1, 0), (0, 0, 0, 0, 1)): [c9],
                                                                                 ((0, 1, 0, 0), (0, 0, 0, 1, 0)): [c10]},
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
                                                                         'mu5': {((1, 0, 0, 0), (0, 0, 0, 1, 0)): [a7]},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0), (0, 0, 0, 0, 1)): [c9]},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1), (0, 1, 0, 0, 0)): [],
                                                                                 ((0, 1, 0, 0), (1, 0, 0, 0, 0)): [c10]},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence2.event_tree)

    def testUpdateNt(self):
        random.seed(4002)
        nt_to_mutate = self.sequence2.nt_sequence[9]
        selected_mutation = self.sim_on_branch2.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch2.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': 9.355258927062633e-06,
                          'C': 1.5139584053833733e-06,
                          'G': None,
                          'T': 5.3386711936328745e-06}

        self.assertEqual('G', nt_to_mutate.state)
        self.assertEqual('C', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree2.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree2.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree2.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree2.traverse_tree()

        exp_sequences = ['ATGAATAAACCCGTATGA', 'ATGAATAAACCCGTATGA', 'ATGAATAAACCCGTATGA', 'ATGAATAAACCCGTATGA',
                         'ATGAATAAACCCGTATGA', 'ATGAATAAACCCGTATGA']

        res_sequences = []
        for clade in self.sim_on_tree2.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nATGAATAAACCCGTATGA\n>B \nATGAATAAACCCGTATGA\n>C \nATGAATAAACCCGTATGA\n>D \nATGAATAAACCCGTATGA\n'
        outfile = './seq2_output.txt'
        self.sim_on_tree2.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


class Seq3(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(987)

        s3 = 'ATGACGTGGTGA'
        orfs3 = {'+0': [{'coords': [[0, 12]],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi3 = Sequence.get_frequency_rates(s3)
        self.sequence3 = Sequence(s3, orfs3, KAPPA, GLOBAL_RATE, pi3, CAT_VALUES)

        branch_length = 0.1
        self.sim_on_branch3 = SimulateOnBranch(self.sequence3, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree3 = SimulateOnTree(self.sequence3, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence3.probability_tree['to_nt']['G']['from_nt']
        expected = 'T'
        result = self.sim_on_branch3.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence3.probability_tree['to_nt']['G']['from_nt']['A']['cat']
        expected = 'mu6'
        result = self.sim_on_branch3.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        expected = (self.sim_on_branch3.sequence.nt_sequence[6], 'G')
        result = self.sim_on_branch3.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        seq3_to_nt_sum = 3.66
        seq3_to_nt = {'A': 0.75, 'C': 0.4, 'T': 1.25, 'G': 1.26}

        expected = 'G'
        result = self.sim_on_branch3.weighted_random_choice(seq3_to_nt, seq3_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.00043717228977664685
        result = self.sim_on_branch3.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch3.mutate_on_branch()
        exp_seq = 'ATGACGTGGTGA'
        exp_res = ''.join(str(nt) for nt in self.sequence3.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        a3 = self.sequence3.nt_sequence[3]
        c4 = self.sequence3.nt_sequence[4]
        g5 = self.sequence3.nt_sequence[5]
        t6 = self.sequence3.nt_sequence[6]
        g7 = self.sequence3.nt_sequence[7]
        g8 = self.sequence3.nt_sequence[8]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 1, 0),): [c4]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {((0, 0, 0, 1),): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((1, 0, 0, 0),): [t6]},
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
                                                                         'mu3': {((0, 1, 0, 0),): [g8]},
                                                                         'mu4': {((0, 0, 1, 0),): [g7]},
                                                                         'mu5': {((0, 0, 0, 1),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 0, 1, 0),): [t6]},
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
                                                                         'mu6': {((0, 0, 1, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c4]}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [t6]},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [c4]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {((1, 0, 0, 0),): [g7]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1),): [g5],
                                                                                 ((1, 0, 0, 0),): [g8]}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

    def testRemoveNt(self):
        random.seed(9001)

        a3 = self.sequence3.nt_sequence[3]
        c4 = self.sequence3.nt_sequence[4]
        g5 = self.sequence3.nt_sequence[5]
        t6 = self.sequence3.nt_sequence[6]
        g7 = self.sequence3.nt_sequence[7]
        g8 = self.sequence3.nt_sequence[8]

        self.sim_on_branch3.remove_nt(a3)

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 1, 0),): [c4]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {((0, 0, 0, 1),): [g5]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((1, 0, 0, 0),): [t6]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {((1, 0, 0, 0),): []},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [g8]},
                                                                         'mu4': {((0, 0, 1, 0),): [g7]},
                                                                         'mu5': {((0, 0, 0, 1),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 0, 1, 0),): [t6]},
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
                                                                         'mu6': {((0, 0, 1, 0),): []}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c4]}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [t6]},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): []}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [c4]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': {'category': {'mu1': {((1, 0, 0, 0),): [g7]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 1),): [g5],
                                                                                 ((1, 0, 0, 0),): [g8]}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence3.event_tree)

    def testUpdateNt(self):
        random.seed(987)
        nt_to_mutate = self.sequence3.nt_sequence[6]
        selected_mutation = self.sim_on_branch3.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch3.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': 1.717474544205974e-06,
                          'C': 3.2418256032424914e-06,
                          'G': 7.85098704946011e-06,
                          'T': None}

        self.assertEqual('T', nt_to_mutate.state)
        self.assertEqual('A', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree3.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree3.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree3.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree3.traverse_tree()

        exp_sequences = ['ATGACGTGGTGA', 'ATGACGTGGTGA', 'ATGACGTGGTGA', 'ATGACGTGGTGA', 'ATGACGTGGTGA', 'ATGACGTGGTGA']

        res_sequences = []
        for clade in self.sim_on_tree3.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nATGACGTGGTGA\n>B \nATGACGTGGTGA\n>C \nATGACGTGGTGA\n>D \nATGACGTGGTGA\n'
        outfile = './seq3_output.txt'
        self.sim_on_tree3.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


class Seq4(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(9001)

        s4 = 'ATGATGCCCTAA'
        orfs4 = {'+0': [{'coords': [[0, 12]],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        self.sequence4 = Sequence(s4, orfs4, KAPPA, GLOBAL_RATE, pi4, CAT_VALUES)

        branch_length = 0.0075
        self.sim_on_branch4 = SimulateOnBranch(self.sequence4, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree4 = SimulateOnTree(self.sequence4, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence4.probability_tree['to_nt']['T']['from_nt']
        expected = 'C'
        result = self.sim_on_branch4.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence4.probability_tree['to_nt']['T']['from_nt']['C']['cat']
        expected = 'mu5'
        result = self.sim_on_branch4.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        random.seed(9001)
        expected = (self.sim_on_branch4.sequence.nt_sequence[6], 'A')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        seq4_to_nt_sum = 3.66
        seq4_to_nt = {'A': 0.75, 'C': 0.4, 'T': 1.25, 'G': 1.26}

        expected = 'G'
        result = self.sim_on_branch4.weighted_random_choice(seq4_to_nt, seq4_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.0004851692506965592
        result = self.sim_on_branch4.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch4.mutate_on_branch()
        exp_seq = 'ATGATGCCCTAA'
        exp_res = ''.join(str(nt) for nt in self.sequence4.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        a3 = self.sequence4.nt_sequence[3]
        t4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        c6 = self.sequence4.nt_sequence[6]
        c7 = self.sequence4.nt_sequence[7]
        c8 = self.sequence4.nt_sequence[8]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [c7]},
                                                                         'mu4': {((0, 0, 0, 1),): [c8]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c6]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [t4]}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((1, 0, 0, 0),): [a3]}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [g5]}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [t4]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [c6],
                                                                                 ((1, 0, 0, 0),): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 1, 0),): [t4]},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [c6]},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {((0, 1, 0, 0),): [c7]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

    def testRemoveNt(self):
        random.seed(9001)

        a3 = self.sequence4.nt_sequence[3]
        t4 = self.sequence4.nt_sequence[4]
        g5 = self.sequence4.nt_sequence[5]
        c6 = self.sequence4.nt_sequence[6]
        c7 = self.sequence4.nt_sequence[7]
        c8 = self.sequence4.nt_sequence[8]

        self.sim_on_branch4.remove_nt(c7)

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): []},
                                                                         'mu4': {((0, 0, 0, 1),): [c8]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c6]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [t4]}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((1, 0, 0, 0),): [a3]}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [g5]}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 1, 0, 0),): [t4]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 1, 0),): [a3]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [c6],
                                                                                 ((1, 0, 0, 0),): []},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 1, 0),): [t4]},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 1, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((1, 0, 0, 0),): [c6]},
                                                                         'mu5': {((0, 0, 0, 1),): [c8]},
                                                                         'mu6': {((0, 1, 0, 0),): []}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g5]},
                                                                         'mu6': {}}},
                                                      'T': None}}}}

        self.assertEqual(exp_event_tree, self.sequence4.event_tree)

    def testUpdateNt(self):
        random.seed(9001)

        nt_to_mutate = self.sequence4.nt_sequence[6]
        selected_mutation = self.sim_on_branch4.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch4.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': None,
                          'C': 6.96488239401271e-05,
                          'G': 5.1134525692942473e-05,
                          'T': 3.661488551883807e-05}

        self.assertEqual('A', nt_to_mutate.state)
        self.assertEqual('T', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree4.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree4.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree4.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree4.traverse_tree()

        exp_sequences = ['ATGATGCCCTAA', 'ATGATGCCCTAA', 'ATGATGCCCTAA', 'ATGATGCCCTAA', 'ATGATGCCCTAA', 'ATGATGCCCTAA']

        res_sequences = []
        for clade in self.sim_on_tree4.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nATGATGCCCTAA\n>B \nATGATGCCCTAA\n>C \nATGATGCCCTAA\n>D \nATGATGCCCTAA\n'
        outfile = './seq4_output.txt'
        self.sim_on_tree4.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


class Seq5(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(121)

        s5 = 'ATGAATGCCTGACTAA'
        orfs5 = {'+0': [{'coords': [[0, 12]],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [{'coords': [[4, 16]],
                         'omega_classes': 4, 'omega_shape': 1.25,
                         'omega_values': [0.09199853806558903, 0.27043066909631136,
                                          0.5158061369385518, 1.1217646558655263]}],
                 '+2': [], '-0': [], '-1': [], '-2': []}
        pi5 = Sequence.get_frequency_rates(s5)
        self.sequence5 = Sequence(s5, orfs5, KAPPA, GLOBAL_RATE, pi5, CAT_VALUES)

        branch_length = 0.00057
        self.sim_on_branch5 = SimulateOnBranch(self.sequence5, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree5 = SimulateOnTree(self.sequence5, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence5.probability_tree['to_nt']['A']['from_nt']
        expected = 'C'
        result = self.sim_on_branch5.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence5.probability_tree['to_nt']['A']['from_nt']['C']['cat']
        expected = 'mu5'
        result = self.sim_on_branch5.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        random.seed(121)
        expected = (self.sim_on_branch5.sequence.nt_sequence[8], 'A')
        result = self.sim_on_branch5.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        seq5_to_nt_sum = 3.09
        seq5_to_nt = {'A': 1.1400000000000001, 'C': 0.19, 'T': 1.0, 'G': 0.76}

        expected = 'A'
        result = self.sim_on_branch5.weighted_random_choice(seq5_to_nt, seq5_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.00019842775070282476
        result = self.sim_on_branch5.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch5.mutate_on_branch()
        exp_seq = 'ATGAATGCCTGACTAA'
        exp_res = ''.join(str(nt) for nt in self.sequence5.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        a3 = self.sequence5.nt_sequence[3]
        c7 = self.sequence5.nt_sequence[7]
        c8 = self.sequence5.nt_sequence[8]
        c12 = self.sequence5.nt_sequence[12]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [c8],
                                                                                 ((1, 0, 0, 0, 0),): [c12]},
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
                                                                         'mu3': {((0, 1, 0, 0),): [a3]},
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
                                                                         'mu6': {((1, 0, 0, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {((0, 1, 0, 0, 0),): [c12]},
                                                                         'mu2': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [c8]},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {((0, 1, 0, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [c8]},
                                                                         'mu2': {((1, 0, 0, 0), (0, 0, 1, 0, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 0, 1),): [c12]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}

        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

    def testRemoveNt(self):
        random.seed(9001)

        a3 = self.sequence5.nt_sequence[3]
        c7 = self.sequence5.nt_sequence[7]
        c8 = self.sequence5.nt_sequence[8]
        c12 = self.sequence5.nt_sequence[12]

        self.sim_on_branch5.remove_nt(c12)

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [c8],
                                                                                 ((1, 0, 0, 0, 0),): []},
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
                                                                         'mu3': {((0, 1, 0, 0),): [a3]},
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
                                                                         'mu6': {((1, 0, 0, 0),): [a3]}}},
                                                      'C': {'category': {'mu1': {((0, 1, 0, 0, 0),): []},
                                                                         'mu2': {((0, 0, 0, 1), (0, 0, 0, 1, 0)): [c8]},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0), (0, 1, 0, 0, 0)): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {((0, 1, 0, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1), (0, 0, 1, 0, 0)): [c8]},
                                                                         'mu2': {((1, 0, 0, 0), (0, 0, 1, 0, 0)): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 0, 0, 0, 1),): []}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence5.event_tree)

    def testUpdateNt(self):
        random.seed(121)
        nt_to_mutate = self.sequence5.nt_sequence[8]
        selected_mutation = self.sim_on_branch5.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch5.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': None,
                          'C': 7.88211405599221e-06,
                          'G': 6.296871779876605e-05,
                          'T': 3.306412916723108e-06}

        self.assertEqual('A', nt_to_mutate.state)
        self.assertEqual('T', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree5.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree5.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree5.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree5.traverse_tree()

        exp_sequences = ['ATGAATGCCTGACTAA', 'ATGAATGCCTGACTAA', 'ATGAATGCCTGACTAA', 'ATGAATGCCTGACTAA',
                         'ATGAATGCCTGACTAA', 'ATGAATGCCTGACTAA']

        res_sequences = []
        for clade in self.sim_on_tree5.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nATGAATGCCTGACTAA\n>B \nATGAATGCCTGACTAA\n>C \nATGAATGCCTGACTAA\n>D \nATGAATGCCTGACTAA\n'
        outfile = './seq5_output.txt'
        self.sim_on_tree5.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


class Seq6(unittest.TestCase):

    def setUp(self):
        self.maxDiff = MAX_DIFF
        random.seed(4675)

        s6 = 'ATGATGGCCCTAA'
        orfs6 = {'+0': [{'coords': [(0, 5), (6, 13)],
                         'omega_classes': 3, 'omega_shape': 1.5,
                         'omega_values': [0.1708353283825978, 0.4810288100937172, 1.1481358615121404]}],
                 '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi6 = Sequence.get_frequency_rates(s6)
        self.sequence6 = Sequence(s6, orfs6, KAPPA, GLOBAL_RATE, pi6, CAT_VALUES)

        branch_length = 0.00087
        self.sim_on_branch6 = SimulateOnBranch(self.sequence6, branch_length)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree6 = SimulateOnTree(self.sequence6, phylo_tree)

    def testSelectKey(self):
        from_tree = self.sequence6.probability_tree['to_nt']['C']['from_nt']
        expected = 'T'
        result = self.sim_on_branch6.select_key(from_tree)
        self.assertEqual(expected, result)

        cat_tree = self.sequence6.probability_tree['to_nt']['C']['from_nt']['T']['cat']
        expected = 'mu3'
        result = self.sim_on_branch6.select_key(cat_tree)
        self.assertEqual(expected, result)

    def testGetSubstitution(self):
        random.seed(4675)
        expected = (self.sim_on_branch6.sequence.nt_sequence[8], 'T')
        result = self.sim_on_branch6.get_substitution()
        self.assertEqual(expected, result)

    def testWeightedRandomChoice(self):
        seq6_to_nt_sum = 4.540000000000001
        seq6_to_nt = {'A': 1.55, 'C': 0.6900000000000001, 'T': 1.1500000000000001, 'G': 1.1500000000000001}

        expected = 'C'
        result = self.sim_on_branch6.weighted_random_choice(seq6_to_nt, seq6_to_nt_sum)
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.0003336861302652134
        result = self.sim_on_branch6.sum_rates()
        self.assertEqual(expected, result)

    def testMutateOnBranch(self):
        random.seed(9991)
        np.random.seed(9991)
        self.sim_on_branch6.mutate_on_branch()
        exp_seq = 'ATGATGGCCCTAA'
        exp_res = ''.join(str(nt) for nt in self.sequence6.nt_sequence)
        self.assertEqual(exp_seq, exp_res)

        a3 = self.sequence6.nt_sequence[3]
        t4 = self.sequence6.nt_sequence[4]
        g6 = self.sequence6.nt_sequence[6]
        c7 = self.sequence6.nt_sequence[7]
        c8 = self.sequence6.nt_sequence[8]
        c9 = self.sequence6.nt_sequence[9]

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1),): [c9]},
                                                                         'mu4': {((1, 0, 0, 0),): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [g6]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): [t4]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 1, 0),): [a3]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [g6]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [t4]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((1, 0, 0, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1),): [c9]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {((1, 0, 0, 0),): [c8]}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): [t4]},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 1, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((1, 0, 0, 0),): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 0, 0, 1),): [c9]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g6]},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

    def testRemoveNt(self):
        random.seed(9001)

        a3 = self.sequence6.nt_sequence[3]
        t4 = self.sequence6.nt_sequence[4]
        g6 = self.sequence6.nt_sequence[6]
        c7 = self.sequence6.nt_sequence[7]
        c8 = self.sequence6.nt_sequence[8]
        c9 = self.sequence6.nt_sequence[9]

        self.sim_on_branch6.remove_nt(t4)

        exp_event_tree = {'to_nt': {'A': {'from_nt': {'A': None,
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((0, 0, 0, 1),): [c9]},
                                                                         'mu4': {((1, 0, 0, 0),): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((0, 1, 0, 0),): [g6]},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {((0, 1, 0, 0),): []},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'C': {'from_nt': {'A': {'category': {'mu1': {},
                                                                         'mu2': {((0, 0, 1, 0),): [a3]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': None,
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {((0, 1, 0, 0),): [g6]},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): []},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'G': {'from_nt': {'A': {'category': {'mu1': {((1, 0, 0, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {((0, 0, 0, 1),): [c9]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 1, 0, 0),): [c7]},
                                                                         'mu5': {},
                                                                         'mu6': {((1, 0, 0, 0),): [c8]}}},
                                                      'G': None,
                                                      'T': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {((1, 0, 0, 0),): []},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}}}},
                                    'T': {'from_nt': {'A': {'category': {'mu1': {((0, 0, 1, 0),): [a3]},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {},
                                                                         'mu6': {}}},
                                                      'C': {'category': {'mu1': {},
                                                                         'mu2': {((1, 0, 0, 0),): [c7]},
                                                                         'mu3': {},
                                                                         'mu4': {((0, 0, 0, 1),): [c9]},
                                                                         'mu5': {},
                                                                         'mu6': {((0, 1, 0, 0),): [c8]}}},
                                                      'G': {'category': {'mu1': {},
                                                                         'mu2': {},
                                                                         'mu3': {},
                                                                         'mu4': {},
                                                                         'mu5': {((1, 0, 0, 0),): [g6]},
                                                                         'mu6': {}}},
                                                      'T': None}}}}
        self.assertEqual(exp_event_tree, self.sequence6.event_tree)

    def testUpdateNt(self):
        random.seed(4675)
        nt_to_mutate = self.sequence6.nt_sequence[8]
        selected_mutation = self.sim_on_branch6.get_substitution()
        new_state = str(selected_mutation[1])  # A
        self.sim_on_branch6.update_nucleotide(nt_to_mutate, new_state)

        expected_rates = {'A': 2.5519465664644717e-05,
                          'C': 1.0159118077381973e-06,
                          'G': 7.2229080855033005e-06,
                          'T': None}

        self.assertEqual('T', nt_to_mutate.state)
        self.assertEqual('A', nt_to_mutate.get_complement_state())
        self.assertEqual(expected_rates, nt_to_mutate.rates)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree6.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree6.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree6.get_parent_clade(child_clade)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)

    def testTraverseTree(self):
        random.seed(9001)
        np.random.seed(9001)
        self.sim_on_tree6.traverse_tree()

        exp_sequences = ['ATGATGGCCCTAA', 'ATGATGGCCCTAA', 'ATGATGGCCCTAA',
                         'ATGATGGCCCTAA', 'ATGATGGCCCTAA', 'ATGATGGCCCTAA']

        res_sequences = []
        for clade in self.sim_on_tree6.phylo_tree.find_clades(order='level'):
            res_sequences.append(str(clade.sequence))
        self.assertEqual(exp_sequences, res_sequences)

    def testGetAlignment(self):
        expected = '>A \nATGATGGCCCTAA\n>B \nATGATGGCCCTAA\n>C \nATGATGGCCCTAA\n>D \nATGATGGCCCTAA\n'
        outfile = './seq6_output.txt'
        self.sim_on_tree6.get_alignment(outfile)
        with open(outfile) as result_handle:
            result = result_handle.read()
        os.remove(outfile)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()
