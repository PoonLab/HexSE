# Store information of sequence

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
    :ivar sequence: list of Nucleotide objects
    :ivar codons: for each nucleotide, a list of tuples (for each ORF) indicating:
                        - the codon the nucleotide is part of
                        - the nucleotide's position within the codon
    """
    def __init__(self, original_seq, orfs=None):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param original_seq: A list of Nucleotides
        :param orfs: <option> A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        :raises ValueError: If the sequence contains invalid characters or if it contains no ORFs
        """

        # Check if sequence is valid
        if not Sequence.valid_sequence(original_seq):
            raise ValueError("Invalid sequence: {}".format(original_seq))

        self.original_seq = original_seq
        self.orfs = {}
        self.sequence = []
        self.codons = []

        # Check if ORFs are valid and sort ORFs by reading frame
        if orfs is not None and Sequence.valid_orfs(self, orfs):
            self.orfs = self.sort_orfs(orfs)
        else:
            unsorted_orfs = self.get_reading_frames()
            if unsorted_orfs:
                self.orfs = self.sort_orfs(unsorted_orfs)

        # Create Nucleotides
        for pos, nt in enumerate(self.original_seq):
            for reading_frame in self.orfs:
                orf_list = self.orfs[reading_frame]
                self.sequence.append(Nucleotide(pos, nt, orf_list, self.original_seq))
                sequence = list(original_seq)
                local_codon = []

            # Get codon for every nucleotide given reading frames
            for reading_frame in self.orfs:
                orfs_list = self.orfs[reading_frame]
                for orf in orfs_list:
                    if pos in range(min(orf[0], orf[1]), max(orf[0], orf[1]) + 1):
                        # if there is an ORF and nucleotide is in that ORF
                        out = self.get_codon(pos, orf)
                        local_codon.append(out)
                    else:
                        local_codon.append(0)

        # Set codons based on the ORFs in the original sequence
        for pos, nt in enumerate(self.original_seq):
            local_codon = []
            # Get the codon each nucleotide in the ORF belongs to
            for reading_frame in self.orfs:
                orf_list = self.orfs[reading_frame]
                for orf in orf_list:
                    out = self.get_codon(pos, orf)
                    local_codon.append(out)
            self.codons.append(local_codon)


    @staticmethod
    def valid_sequence(seq):
        """
        Verifies that the length of the input sequence is valid and the sequence is composed of only nucleotides.
        A valid sequence is assumed to be composed of a START codon, at least one amino acid codon, and a STOP codon.
        :return is_valid: <True> if the sequence is valid, <False> otherwise
        """
        is_valid = len(seq) >= 9 and all(pos in NUCLEOTIDES for pos in seq.upper())
        return is_valid

    def valid_orfs(self, orfs):
        """
        Verifies that the input ORFs are a list of tuples containing the start and end positions of ORFs.
        Example of valid input: [(1, 9), (27, 13)]
        Example of invalid input: (1, 9), (27, 13)
        :param orfs: The list of open reading frames
        :return: <True> if the ORFs are valid, <False> otherwise
        """
        if not type(orfs) == list:
            print("Invalid format: {} \nOpen reading frames must of the following format: [(0, 8), ...] where 0  "
                  "is the start position of the orf, and 8 is the end position.".format(orfs))
            return False

        for orf in orfs:
            # Check that each orfs is a tuple of length 2
            if type(orf) is not tuple or len(orf) is not 2:
                print("Invalid orf: {} \nOpen reading frames must of the following format: [(0, 8), ...] where 0  "
                      "is the start position of the orf, and 8 is the end position.".format(orf))
                return False

            # Check that the ORF range is valid
            if orf[0] == orf[1]:
                print("Invalid orf: {}".format(orf))
                return False

            # Check that the start and end positions are integers
            if type(orf[0]) is not int or type(orf[1]) is not int:
                print("Invalid orf: {} \nStart and end positions must be integers.".format(orf))
                return False

            # Check that the start and stop positions are in te range of the sequence
            if 0 > orf[0] or len(self.original_seq) < orf[0] or 0 > orf[1] or len(self.original_seq) < orf[1]:
                print("Invalid orf: {} \nPositions must be between 0 and {}".format(orf, len(self.original_seq)))
                return False

        return True

    def reverse_and_complement(self):
        """
        Generates the reverse complement of a DNA sequence
        :return rcseq: The reverse complement of the sequence
        """
        rseq = reversed(self.original_seq.upper())
        rcseq = ''
        for i in rseq:  # reverse order
            rcseq += COMPLEMENT_DICT[i]
        return rcseq

    def get_reading_frames(self):
        """
        Gets positions of the START and STOP codons for each open reading frame in the forward and reverse directions.
        Positions of the START and STOP codons indexed relative to the forward strand.
        :return reading_frames: a list of tuples containing the index of the first nucleotide of
                    the START codon and the index of the last nucleotide of the STOP codon
        """
        start_codon = re.compile('ATG', flags=re.IGNORECASE)
        stop = re.compile('(TAG)|(TAA)|(TGA)', flags=re.IGNORECASE)
        reading_frames = []

        # Record positions of all potential START codons in the forward (positive) reading frame
        fwd_start_positions = [match.start() for match in start_codon.finditer(self.original_seq)]

        # Find open forward open reading frames
        for position in fwd_start_positions:
            frame = position % 3

            internal_met = False
            # If the ATG codon is an internal methionine and not an initiation codon
            for orf in reversed(reading_frames):

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
                        reading_frames.append(orf)
                        break

        # Forward (positive) reading frames of the reverse complement of the original
        # sequence is equivalent to reverse (negative) reading frames of the original sequence
        rcseq = self.reverse_and_complement()

        # Record positions of all potential START codons in the reverse (negative) reading frame
        rev_start_positions = [match.start() for match in start_codon.finditer(rcseq)]

        # Find reverse open reading frames
        for position in rev_start_positions:
            frame = position % 3

            internal_met = False
            # If the ATG codon is an internal methionine and not an initiation codon
            for orf in reversed(reading_frames):

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
                        reading_frames.append(orf)
                        break

        return reading_frames

    @staticmethod
    def sort_orfs(unsorted_orfs):
        """
        Store ORFs in position according to plus zero ORF (first of the list).
        They will be classified as (+0, +1, +2, -0, -1, -2)
        :return orf_position: List of ORFs classified according to their shift relative to
                    the plus zero reading frame  (+0, +1, +2, -0, -1, -2)
        """
        plus_zero, plus_one, plus_two, minus_zero, minus_one, minus_two = [], [], [], [], [], []

        if unsorted_orfs:
            plus_zero_orf = unsorted_orfs[0]

            for orf in unsorted_orfs:
                difference = abs(orf[0] - plus_zero_orf[0]) % 3
                if orf[0] < orf[1]:     # positive strand
                    if difference == 0:
                        plus_zero.append(orf)
                    elif difference == 1:     # plus one
                        plus_one.append(orf)
                    elif difference == 2:   # plus two
                        plus_two.append(orf)

                elif orf[0] > orf[1]:   # negative strand
                    if difference == 0 or difference == plus_zero_orf[0] % 3:
                        minus_zero.append(orf)
                    elif difference == 1:
                        minus_one.append(orf)
                    elif difference == 2:
                        minus_two.append(orf)

        sorted_orfs = {'+0': plus_zero, '+1': plus_one, '+2': plus_two,
                       '-0': minus_zero, '-1': minus_one, '-2': minus_two}

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
            codon = my_orf[position_in_orf: position_in_orf + 3]
        elif position_in_orf % 3 == 1:
            position_in_codon = 1
            codon = my_orf[position_in_orf - 1: position_in_orf + 2]
        else:
            position_in_codon = 2
            codon = my_orf[position_in_orf - 2: position_in_orf + 1]

        return codon, position_in_codon

    def nt_in_orfs(self, position):
        """
        Checks which open reading frames the nucleotide is part of
        :return: True if the nucleotide is part of the ORF, False otherwise
        """
        if 0 > position or len(self.original_seq) < position:
            raise ValueError("Invalid position {} \n "
                             "Position must be between 0 and {}".format(position, len(self.original_seq)))
        in_orf = []
        for reading_frame in self.orfs:
            orfs_list = self.orfs[reading_frame]
            for orf in orfs_list:
                if orf[0] < orf[1]:  # positive strand
                    if orf[0] <= position <= orf[1]:
                        in_orf.append(orf)
                elif orf[0] > orf[1]:  # negative strand
                    if orf[0] >= position >= orf[0]:
                        in_orf.append(orf)
        return in_orf


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    """

    def __init__(self, letter, position, sorted_orfs):
        """
        To which orfs does the nucleotide belong to
        :param letter: Nucleotide A, C, G or T
        :param position: Position of the nucleotide in the sequence
        :param sorted_orfs: <option> A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        """

        self.letter = letter
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



