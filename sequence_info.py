# Store information of sequence
from ovrf_functions import sort_orfs
from ovrf_functions import get_codon

class Sequence(list):
    """
    List of nucleotides in seq and its orfs.
    :param original_seq: Nucleotide sequence as string
    :param sorted_orfs: list of ORFs classified according to their reading frame shift relative to the first orf (+0, +1, +2, -0, -1, -2)
                        (output of sort_orfs)
    :param self.codon: for each nucleotide, list of tuples (for each ORF) indicating:
                        - the codon the nucleotide is part of
                        - the nucleotide's specific position within the codon


    """

    def __init__(self, original_seq, sorted_orfs):
        super(Sequence, self).__init__()
        self.original_seq = original_seq
        self.orfs = sorted_orfs
        self.codon = []

        for position in range(len(original_seq)):
            self.append(Nucleotide(original_seq[position], position, sorted_orfs))

            sequence = list(original_seq)
            local_codon = []
            #Get codon for every nucleotide given reading frames
            for orf in self.orfs:
                if type(orf) == tuple:
                    out = get_codon(sequence, position, orf)
                    # Note: if nt is not in orf, then codon will be an empty string
                    local_codon.append(out)
                else:
                    local_codon.append(0)
            self.codon.append(local_codon)


class Nucleotide(str):
    """
    Nucleotide information (letter, position in seq and orfs to which it belongs to)
    """

    # method '__new__': https://stackoverflow.com/a/30045261/638425
    def __new__(self, letter, *args, **kwargs):
        return super(Nucleotide, self).__new__(self, letter.upper())

    def __init__(self, letter, position, sorted_orfs):
        """
        To which orfs does the nucleotide belong to
        :param letter: Nucleotide A, C, G or T
        :param position: position in <seq>
        :param sorted_orfs: list of orf as tuples defined by user
        """
        self.position = position
        in_orf = [False,] * 6
        for i in range(len(sorted_orfs)):
            orf = sorted_orfs[i]
            if type(orf) == tuple:
                if orf[0] < orf[1]: # positive strand
                    if self.position in range(orf[0],orf[1]+1):
                        in_orf[i] = True
                elif orf[0] > orf[1]: # negative strand
                    if self.position in range(orf[1],orf[0]+1):
                        in_orf[i] = True
        self.in_orf = in_orf
