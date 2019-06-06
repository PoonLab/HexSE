# Store information of sequence
from ovrf_functions import CODON_DICT
from ovrf_functions import COMPLEMENT_DICT

class Sequence(list):
    """
    List of nucleotides in seq and its orfs.
    :param original_seq: Nucleotide sequence as string
    :param sorted_orfs: list of ORFs classified according to their reading frame shift relative
                        to the first orf (+0, +1, +2, -0, -1, -2) (output of sort_orfs)
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
            # Get codon for every nucleotide given reading frames
            for orf in self.orfs:
                if type(orf) == tuple:
                    if position in range(min(orf[0], orf[1]), max(orf[0], orf[1])+1):
                    # if there is and orf and nucleotide is in that orf
                        out = self.get_codon(sequence, position, orf)
                        local_codon.append(out)
                    else:
                        local_codon.append(0)

                else:
                    local_codon.append(0)
            self.codon.append(local_codon)

    def get_codon(self, seq, position, orf):
        """
        Get codon sequence, and position of my_nt in the codon
        :param seq: parental sequence as list of nucleotides
        :param position: position of the nucleotide in <seq>
        :param orf: tuple indicating first and last nucleotide of an open reading frame
        :return codon: tuple with nucleotide triplet and position of the nucleotide in the codon
        """

        if orf[1] > orf[0]:  # positive strand
            my_orf = ''.join(seq[orf[0]:orf[1] + 1])
            position_in_orf = position - orf[0]
        else:  # negative strand
            sub_seq = seq[orf[1]:orf[0] + 1]
            my_orf = self.reverse_and_complement(sub_seq)
            position_in_orf = orf[0] - position

        try:
            position_in_orf < 0
        except ValueError as e:
            raise ValueError("Invalid position: {}".format(position_in_orf))

        if position_in_orf % 3 == 0:
            position_in_codon = 0
            codon = my_orf[position_in_orf:position_in_orf + 3]
        elif position_in_orf % 3 == 1:
            position_in_codon = 1
            codon = my_orf[position_in_orf - 1:position_in_orf + 2]
        else:
            position_in_codon = 2
            codon = my_orf[position_in_orf - 2:position_in_orf + 1]

        return codon, position_in_codon

    def reverse_and_complement(self, seq):
        """
        Generates the reverse complement of a DNA sequence
        :param seq: the DNA sequence
        :return: the reverse complement of the sequence
        """
        string_seq = ''.join(seq)
        rseq = reversed(string_seq.upper())
        rcseq = ''
        for i in rseq:  # reverse order
            try:
                rcseq += COMPLEMENT_DICT[i]
            except KeyError as e:
                raise KeyError("Invalid character '{}' in sequence".format(i))

        return rcseq


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
        super(Nucleotide, self).__init__()
        self.position = position
        in_orf = [False,] * 6
        for i in range(len(sorted_orfs)):
            orf = sorted_orfs[i]
            if type(orf) == tuple:
                if orf[0] < orf[1]:  # positive strand
                    if self.position in range(orf[0], orf[1]+1):
                        in_orf[i] = True
                elif orf[0] > orf[1]:  # negative strand
                    if self.position in range(orf[1], orf[0]+1):
                        in_orf[i] = True
        self.in_orf = in_orf
