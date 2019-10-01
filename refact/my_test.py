# Testing
import random

from new_sequence_info import Sequence
from run_simulate import reverse_and_complement
from run_simulate import sort_orfs
from run_simulate import valid_orfs
NUCLEOTIDES = ['A', 'C', 'G', 'T']


# orfs = [(0,20), (18,1), (16,2), (7,18)]
# sorted_orfs = sort_orfs(orfs)
# ori_seq = 'GTACGATCGATCGATGCTAGC'
# rcseq = reverse_and_complement(ori_seq)
# mu = 0.5
# pi = None
# kappa = 0.3



orfs = [(0,5), (9,1), (7,18)]
sorted_orfs = sort_orfs(orfs)
ori_seq = 'GTACGGGGATCGATAAACGATGCTAGC'
rcseq = reverse_and_complement(ori_seq)
mu = 0.5
pi = None
kappa = 0.3

sequence = Sequence(ori_seq, rcseq, sorted_orfs, mu, pi, kappa)
event_tree = sequence.event_tree

#print(my_tree['to_nt']['A']['from_nt']['C']['number_of_events'])
#print(my_tree['total_events'])
#print(my_tree['to_nt']['A']['events_for_nt'])

# Set parameters
percentage = random.uniform(0,1)
total = event_tree['total_events']
limit = total*percentage
sum = 0

# Draw event
my_tree = event_tree['to_nt']
iter_object = iter(my_tree.items())

result = None
while sum < limit:
    try:
        (key, to_nt) = next(iter_object)
        result = to_nt
        number = to_nt['events_for_nt']
        sum += number
        print("limit: {}, sum: {} ".format(limit, sum))
    except StopIteration:
        break

print(result)