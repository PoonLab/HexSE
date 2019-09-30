# Testing
import sys
print(sys.path)

from new_sequence_info import Sequence
from run_simulate import reverse_and_complement
from run_simulate import sort_orfs
from run_simulate import valid_orfs


orfs = [(0,20), (18,1), (16,2), (7,18)]
sorted_orfs = sort_orfs(orfs)
ori_seq = 'GTACGATCGATCGATGCTAGC'
rcseq = reverse_and_complement(ori_seq)
mu = 0.5
pi = None
kappa = 0.3

sequence = Sequence(ori_seq, rcseq, sorted_orfs, mu, pi, kappa)