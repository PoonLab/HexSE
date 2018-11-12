import random

nucleotides = ["A","G","C","T"]


# Create a string with 100 random selected nucleotides
nucleotide_sequence=' '
for i in range(0,6):
    nucleotide_sequence+= random.choice(nucleotides)
print(nucleotide_sequence)

# Mutate a nucleotide sequence in a random selected position with a randomly chosen nucleotide
def mutate(sequence):
     sequence_list=list(sequence)
     mutation_site=random.randint (0,len(sequence_list)-1)
     sequence_list[mutation_site]=random.choice(nucleotides)
     return ''.join(sequence_list)

m_sequence=mutate(nucleotide_sequence)

print ("This is the original sequence", nucleotide_sequence)
print ("This is the mutated sequence ", m_sequence)