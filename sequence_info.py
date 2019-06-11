import re

NUCLEOTIDES = ['A', 'C', 'G', 'T']

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'W', 'R': 'Y', 'K': 'M', 'Y': 'R',
                   'S': 'S', 'M': 'K', 'B': 'V', 'D': 'H',
                   'H': 'D', 'V': 'B', '*': '*', 'N': 'N',
                   '-': '-'}

CODON_DICT = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
              'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
              'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
              '---': '-', 'XXX': '?'}


class Sequence:
    """
    Stores information about the sequence
    :ivar original_seq: the nucleotide sequence as a string
    :ivar orfs: list of ORFs classified according to their reading frame shift relative
                        to the first orf (+0, +1, +2, -0, -1, -2)
    :ivar codons: for each nucleotide, a list of tuples (for each ORF) indicating:
                        - the codon the nucleotide is part of
                        - the nucleotide's position within the codon
    """

    def __init__(self, original_seq, orfs=None):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons .
        :param original_seq: Nucleotide sequence as string
        :param orfs: <option> A list of tuples containing the start and end positions of open reading frames.
                        If no reading frames are specified, get_reading_frames will find orfs
        :raises ValueError: if the sequence contains invalid characters or if it contains no ORFs
        """

        if not Sequence.valid_sequence(original_seq):
            raise ValueError("Invalid sequence: {}".format(original_seq))

        self.original_seq = original_seq
        self.orfs = {}
        self.codons = []

        if orfs is not None:
            if not Sequence.valid_orfs(self, orfs):
                raise ValueError("Invalid orf: {} "
                                 "\nOpen reading frames must of the following format: [(0, 8), ...] where 0  "
                                 "is the start position of the orf, and 8 is the end position.".format(orfs))
            else:
                self.orfs = self.sort_orfs(orfs)

        else:
            unsorted_orfs = self.get_open_reading_frames()
            if unsorted_orfs:
                self.orfs = self.sort_orfs(unsorted_orfs)

                # Set codons based on the ORFs in the original sequence
                for pos, nt in enumerate(self.original_seq):
                    local_codon = []
                    # Get the codon each nucleotide in the ORF belongs to
                    for orf in self.orfs:
                        out = self.get_codon(pos, orf)
                        local_codon.append(out)
                    self.codons.append(local_codon)

    @staticmethod
    def valid_sequence(seq):
        """
        Verifies that the length of the input sequence is valid and the sequence is composed of only nucleotides.
        A valid sequence is assumed to be composed of a START codon, at least one amino acid codon, and a STOP codon.
        :return is_valid: true if the sequence is valid, false otherwise
        """
        is_valid = len(seq) >= 9 and all(pos in NUCLEOTIDES for pos in seq.upper())
        return is_valid

    @staticmethod
    def valid_orfs(self, orfs):
        """
        Verifies that the input orfs are a list of tuples containing the start and end positions of orfs.
        Example of valid input: [(1, 9), (27, 13)]
        Example of invalid input: (1, 9), (27, 13) is not valid
        :param orfs: the list of open reading frames
        :return: true if the orfs are valid, false otherwise
        """
        if not type(orfs) == list:
            return False
        for orf in orfs:
            # Check that each orfs is a tuple of length 2
            if type(orf) is not tuple or len(orf) is not 2:
                return False
            # Check that the orf range is valid
            if orf[0] == orf[1]:
                return False
            # Check that the start and end positions are integers
            if type(orf[0]) is not int or type(orf[1]) is not int:
                return False
            # Check that the start and stop positions are in te range of the sequence
            if 0 > orf[0] or len(self.original_seq) < orf[0] or 0 > orf[1] or len(self.original_seq) < orf[1]:
                return False

        return True

    def reverse_and_complement(self):
        """
        Generates the reverse complement of a DNA sequence
        :return rcseq: the reverse complement of the sequence
        """
        rseq = reversed(self.original_seq.upper())
        rcseq = ''
        for i in rseq:  # reverse order
            rcseq += COMPLEMENT_DICT[i]
        return rcseq

    def get_open_reading_frames(self):
        """
        Creates a list with tuples containing the first and last position
        of forwards are reverse reading frames in sequence according to START and STOP codons.
        Open reading frames are indexed relative to the forward strand.
        :return reading_frames: a list of tuples containing the index of the first nucleotide of
                    the START codon and the index of the last nucleotide of the STOP codon
        """
        start_codon = re.compile('ATG', flags=re.IGNORECASE)
        stop = re.compile('(TAG)|(TAA)|(TGA)', flags=re.IGNORECASE)
        fwd_reading_frames = []

        # Record positions of all potential START codons in the forward (positive) reading frame
        fwd_start_positions = [match.start() for match in start_codon.finditer(self.original_seq)]

        # Find open forward open reading frames
        for position in fwd_start_positions:
            frame = position % 3

            internal_met = False
            # If the ATG codon is an internal methionine and not an initiation codon
            for orf in reversed(fwd_reading_frames):

                # If the START codon and the potential START codon are in the same reading frame
                # and the existing ORF ends before the potential ORF, stop searching
                if orf[0] % 3 == frame and orf[1] < position:
                    break

                # If the potential START codon is between the range of the START and STOP codons,
                # and it is in the same frame, the codon is an internal methionine
                if orf[0] < position < orf[1] and orf[0] % 3 == frame:
                    internal_met = True
                    break

            # If the ATG is a START codon and not simply methionine
            if not internal_met:
                for match in stop.finditer(self.original_seq, position):
                    orf_length = match.end() - position
                    # Find a stop codon and ensure ORF length is sufficient in the forward strand
                    if match.start() % 3 == frame and orf_length >= 8:
                        # Get the positions in the sequence for the first and last nt of the RF
                        orf = (position, match.end() - 1)
                        fwd_reading_frames.append(orf)
                        break

        # Forward (positive) reading frames of the reverse complement of the original
        # sequence is equivalent to reverse (negative) reading frames of the original sequence
        rcseq = self.reverse_and_complement()
        rev_reading_frames = []

        # Record positions of all potential START codons in the reverse (negative) reading frame
        rev_start_positions = [match.start() for match in start_codon.finditer(rcseq)]

        # Find reverse open reading frames
        for position in rev_start_positions:
            frame = position % 3

            internal_met = False
            # If the ATG codon is an internal methionine and not an initiation codon
            for orf in reversed(rev_reading_frames):

                # If the START codon and the potential START codon are in the same reading frame
                # and the existing ORF ends before the potential ORF, stop searching
                if orf[0] % 3 == frame and orf[1] < position:
                    break

                # If the potential START codon is between the range of the START and STOP codons,
                # and it is in the same frame, the codon is an internal methionine
                if orf[0] < position < orf[1] and orf[0] % 3 == frame:
                    internal_met = True
                    break

            # If the ATG is a START codon and not simply methionine
            if not internal_met:
                for match in stop.finditer(rcseq, position):
                    orf_length = match.end() - position
                    # Find a stop codon and ensure ORF length is sufficient in the forward strand
                    if match.start() % 3 == frame and orf_length >= 8:
                        # Get the positions in the sequence for the first and last nt of the RF
                        orf = (len(rcseq) - 1 - position, len(rcseq) - match.end())
                        rev_reading_frames.append(orf)
                        break

        return fwd_reading_frames + rev_reading_frames

    def sort_orfs(self, unsorted_orfs):
        """
        Store orfs in position according to plus zero orf (first of the list).
        They will be classified as (+0, +1, +2, -0, -1, -2)
        :return orf_position: list of ORFs classified according to their shift relative to
                    the plus zero reading frame  (+0, +1, +2, -0, -1, -2)
        """
        sorted_orfs = {}
        if unsorted_orfs:
            plus_zero_orf = unsorted_orfs[0]

            for orf in unsorted_orfs:
                difference = abs(orf[0] - plus_zero_orf[0]) % 3
                if orf[0] < orf[1]:     # positive strand
                    if difference == 0:
                        sorted_orfs[orf] = '+0'
                    elif difference == 1:     # plus one
                        sorted_orfs[orf] = '+1'
                    elif difference == 2:   # plus two
                        sorted_orfs[orf] = '+2'

                elif orf[0] > orf[1]:   # negative strand
                    if difference == 0 or difference == plus_zero_orf[0] % 3:
                        sorted_orfs[orf] = '-0'
                    elif difference == 1:
                        sorted_orfs[orf] = '-1'
                    elif difference == 2:
                        sorted_orfs[orf] = '-2'
        return sorted_orfs

    def get_codon(self, position, orf):
        """
        Get codon sequence, and position of my_nt in the codon
        :param position: position of the nucleotide in the sequence
        :param orf: tuple indicating first and last nucleotide of an open reading frame
        :return codon, position_in_codon: tuple with nucleotide triplet and position of the nucleotide in the codon
        """
        if orf[1] > orf[0]:  # positive strand
            my_orf = self.original_seq[orf[0]:orf[1] + 1]
            position_in_orf = position - orf[0]
        else:  # negative strand
            rseq = self.reverse_and_complement()
            my_orf = rseq[orf[1]:orf[0] + 1]
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

    def nt_in_orfs(self, position):
        """
        Checks which open reading frames the nucleotide is part of
        :return: True if the nucleotide is part of the orf, False otherwise
        """
        if 0 > position or len(self.original_seq) < position:
            raise ValueError("Invalid position {} \n "
                             "Position must be between 0 and {}".format(position, len(self.original_seq)))
        in_orf = []
        for orf in self.orfs:
            if orf[0] < orf[1]:  # positive strand
                if orf[0] <= position <= orf[1]:
                    in_orf.append(orf)
            elif orf[0] > orf[1]:  # negative strand
                if orf[0] >= position >= orf[0]:
                    in_orf.append(orf)
        return in_orf
