**Simulator description**

1. Sequence info script:
*Sequence class*
- Check that sequence is a string of nucleotides longer than 9
- Check that orfs is a list of tuples that are divisible by three and are inside the sequence lenght
- Generates the reverse complement of seq
- Gets the START and STOP codons for each open reading frame (ORF) both in Forward and Reverse strands (This is only used in case that the used does not specify the orfs for the sequence)
- Sort ORFs according to plus zero ORF (first of the list), and classified them as (+0, +1, +2, -0, -1, -2). The ouput is a dictionary with six keys that access the reading frames for each option. 

*Nucleotide class*
- Store information related to each nucleotide (state, position in sequence, pointer to nucleotide on the left, pointer to nucleotide on the right)

*Double linked list class*
- Lists together objects of class Nucleotide
- Default initialization with empty head node
- Method: Slice sequence (to partition the sequence according to reading frames)
- Method: nucleotide at position (retrieve the nucleotide at a specific position in sequence)

*Codon class*
- Finds positition of the Nucleotide in the condon
- Checks if a mutation on that nucleotide is syn or nonsyn

**EvenTree ckass**
- Dictionary with all posible mutations that may occur in the sequence
-- First level is 'to_nt' and second level is 'from_nt', this last one contains a list of all nucleotides that can be part of that category (ie., A can change to T, C or G, so every A on sequence is going to be in the tips with the path event_tree['to_nt'][C]['from_nt']['A'], event_tree['to_nt'][G]['from_nt']['A'] ,event_tree['to_nt'][T]['from_nt']['A'])

