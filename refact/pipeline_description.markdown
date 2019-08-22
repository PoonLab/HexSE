**Simulator description**

1. Sequence info script:
*Sequence class*
- Check that sequence is a string of nucleotides longer than 9
- Check that orfs is a list of tuples that are divisible by three and are inside the sequence lenght
- Generates the reverse complement of seq
- Gets the START and STOP codons for each open reading frame (ORF) both in Forward and Reverse strands (This is only used in case that the used does not specify the orfs for the sequence)
- Sort ORFs according to plus zero ORF (first of the list), and classified them as (+0, +1, +2, -0, -1, -2). The ouput is a dictionary with six keys that access the reading frames for each option. 

*Nucleotide class*
- Check in which orf the nucleotide is involve
- Retrieve the codon and the position in the codon for each orf

