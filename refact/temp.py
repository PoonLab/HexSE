from new_sequence_info import Nucleotide
from new_sequence_info import Sequence
from new_sequence_info import DoubleLinkedList
l1 = 'A'
l2 = 'T'
l3 = 'G'
pos1 = 1 
pos2 = 2
pos3 = 3

#nt1 = Nucleotide(l1, pos1)
#nt2 = Nucleotide(l2, pos2)
#nt3 = Nucleotide(l3, pos3)

#print(nt2.left_nt)
#nt1.set_right_nt(nt2)
#nt2.set_left_nt(nt1)
#nt2.set_right_nt(nt3)
#nt3.set_left_nt(nt2)
#print(nt1.right_nt.get_letter())
#print(nt1.right_nt.right_nt.get_letter())

# Testing DoubleLinkedList
#seq = DoubleLinkedList()
#seq.insert_nt('A', 0)
#seq.insert_nt('G', 1)
#seq.insert_nt('T', 2)

s = 'AAACCCGGGTTTATATCTATCTAGAGCTAGATAGCTA'
orfs = [(0,8), (8,0), (12,4), (9,17), (11,3)]
mu = 0.3
kappa = 0.5
pi = 0.3
omega = 1
seq = Sequence(s,mu, kappa, orfs) 
print(seq.orfs)
print(seq.create_nt_orf_dict(4))

print(seq.get_sequence())   # double linked list
print(seq.get_sequence().get_head())   # first nucleotide
print(seq.get_sequence().get_head().get_pos_in_codons()) # get orf dictionary for this specific nt
print(seq.get_sequence().print_seq())
