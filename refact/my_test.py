# Testing
import random

from new_sequence_info import Sequence
from run_simulate import reverse_and_complement
from run_simulate import sort_orfs
from run_simulate import valid_orfs
from simulation_refactore import Simulate

orfs = [(0,5), (9,1), (7,18)]
sorted_orfs = sort_orfs(orfs)
ori_seq = 'GTACGGGGATCGATAAACGATGCTAGC'
rcseq = reverse_and_complement(ori_seq)
mu = 0.5
pi = None
kappa = 0.3

sequence = Sequence(ori_seq, rcseq, sorted_orfs, mu, pi, kappa)
event_tree = sequence.event_tree
#print(sequence.nt_sequence)
#print(event_tree)

simulation = Simulate(sequence)
substitution = simulation.get_substitution()
print(substitution)
my_nt = substitution[0]
print(simulation.sum_rates())
simulation.mutate_on_branch(15)

#print(my_nt.rates.values())
#print(simulation)
# print("My nucleotide: {}, rates for nucleotide: {}, to state: {} \n".format(substitution[0].get_state(), substitution[0].rates, substitution[1]))
#
