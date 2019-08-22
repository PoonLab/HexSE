from  new_sequence_info import Sequence

sequence = 'ACGAGCTAGCTAGCTATAGCT'
test = [(0,8),(12,4)]

seq = Sequence(sequence, test)
for nt in seq:
    print(seq.orfs)
