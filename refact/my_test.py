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
#print(event_tree, "\n")

simulation = Simulate(sequence)
#substitution = simulation.get_substitution()
#print(substitution)
#my_nt = substitution[0]
# print("My Nucleotide information: {}, {}, {}, {}, RATES AND OMEGAS: {}, {}".format(my_nt.codons, my_nt.pos_in_seq, my_nt.state, my_nt.complement_state, my_nt.my_omegas, my_nt.rates))
# print("New subs rates", sequence.get_substitution_rates(my_nt))
# print("Local Event Tree: {}".format(event_tree['to_nt'][substitution[1]]['from_nt'][substitution[0].state]))
# print("Keys for nucleotide: {}".format(my_nt.my_omegas))
# simulation.find_and_remove_nt(my_nt)

simulation.mutate_on_branch(15)