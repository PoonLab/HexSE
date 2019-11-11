# Testing
import random

from new_sequence_info import Sequence
from run_simulate import reverse_and_complement
from run_simulate import sort_orfs
from run_simulate import valid_orfs
from simulation_refactore import SimulateOnBranch, SimulateOnTree
from Bio import Phylo, SeqIO
import copy

tree = Phylo.read('test_tree.txt','newick', rooted=True)
orfs = [(0,5), (9,1), (7,18)]
sorted_orfs = sort_orfs(orfs)
ori_seq = 'GTACGGGGATCGATAAACGATG'
print(type(ori_seq))
rcseq = reverse_and_complement(ori_seq)
mu = 0.0005
pi = None
kappa = 0.3

sequence = Sequence(ori_seq, rcseq, sorted_orfs, mu, pi, kappa)
sim_tree = SimulateOnTree(sequence, tree)
print(sim_tree.get_alignment())



######
# TEST with most realistic data
# tree = Phylo.read('../HBV/abayomi_tree.nwk', 'newick', rooted = True)
# seq = str([seq_record.seq for seq_record in SeqIO.parse('../HBV/HBV.fasta' , "fasta")][0])
# print(seq)
# print(type(seq))
# sequence = Sequence(seq, rcseq, sorted_orfs, mu, pi, kappa)
# #print(getattr(sequence, 'event_tree'))
# event_tree = getattr(sequence, 'event_tree')
# print(event_tree)
#copy_tree = copy.deepcopy(event_tree)
#copy_seq = copy.deepcopy(sequence)
# sim_tree = SimulateOnTree(sequence, tree)
# print(sim_tree.get_alignment())
