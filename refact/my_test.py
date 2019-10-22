# Testing
import random

from new_sequence_info import Sequence
from run_simulate import reverse_and_complement
from run_simulate import sort_orfs
from run_simulate import valid_orfs
from simulation_refactore import SimulateOnBranch, SimulateOnTree
from Bio import Phylo

tree = Phylo.read('test_tree.txt','newick', rooted=True)
orfs = [(0,5), (9,1), (7,18)]
sorted_orfs = sort_orfs(orfs)
ori_seq = 'GTACGGGGATCGATAAACGATG'
rcseq = reverse_and_complement(ori_seq)
mu = 0.0005
pi = None
kappa = 0.3

sequence = Sequence(ori_seq, rcseq, sorted_orfs, mu, pi, kappa)
sim_tree = SimulateOnTree(sequence, tree)
print(sim_tree.get_alignment())