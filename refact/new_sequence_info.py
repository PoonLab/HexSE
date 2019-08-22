# store sequence information
import re
import sys

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
    Store inputs and create sequence objects
    :ivar original_seq: the nucleotide sequence as a string
    :ivar orfs: list of ORFs (as tuples) classified according to their reading frame shift relative
                        to the first orf (+0, +1, +2, -0, -1, -2)
    :ivar nt_sequence: list of Nucleotide objects
    """

    def __init__(self, original_seq, unsorted_orfs = None):
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
            print("Invalid sequence: {}".format(original_seq))
            sys.exit(0)

        self.original_seq = original_seq
        self.orfs = {}
        self.nt_sequence = []
        self.rcseq = self.reverse_and_complement(original_seq)

        # Check if the ORFs the user specified are valid
        if unsorted_orfs is not None:
            if not Sequence.valid_orfs(self, unsorted_orfs):
                sys.exit(0)

            # Since ORFs are valid, sort the ORFs by reading frame
            self.orfs = self.sort_orfs(unsorted_orfs)

        # If the user did not specify ORFs
        else:
            unsorted_orfs = self.get_open_reading_frames()
            self.orfs = self.sort_orfs(unsorted_orfs)

        # Create Nucleotides
        for pos, nt in enumerate(self.original_seq):
            self.nt_sequence.append(Nucleotide(nt, pos, self.original_seq, self.orfs))


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

            # Check that the start and stop positions are in the range of the sequence
            if 0 > orf[0] or len(self.original_seq) < orf[0] or 0 > orf[1] or len(self.original_seq) < orf[1]:
                print("Invalid orf: {} \nPositions must be between 0 and {}".format(orf, len(self.original_seq)))
                return False

            # Check that the ORF is composed of codons
            if orf[1] > orf[0]:     # Forward strand
                if (orf[1] - orf[0]) % 3 != 2:
                    print("Invalid orf: {}\n ORFs must be composed of codons".format(orf))
                    return False
            if orf[0] > orf[1]:  # Reverse strand
                if (orf[0] - orf[1]) % 3 != 2:
                    print("Invalid orf: {}\n ORFs must be composed of codons".format(orf))
                    return False

        return True

    def reverse_and_complement(self, my_region=None):
        """
        Generates the reverse complement of a DNA sequence
        :param: <option> A sub-sequence of the original sequence
        :return rcseq: The reverse complement of the sequence
        """
        if my_region is None:
            seq = self.original_seq
        else:
            seq = my_region

        rseq = reversed(seq.upper())
        rcseq = ''
        for i in rseq:  # reverse order
            rcseq += COMPLEMENT_DICT[i]
        return rcseq

    def get_open_reading_frames(self):
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


class Nucleotide:

    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    :ivar seq: nucleotide sequence of origin
    :ivar letter: the nucleotide base (A, T, G, C)
    :ivar position: the position of the nucleotide relative to the start of the sequence
    :ivar in_orfs: A list of 6 Boolean values, indicating the reading frames the Nucleotide is part of
                 If the Nucleotide is in True if the nucleotide is in the frame, False otherwise
    :ivar codons: A list of tuples that stores:
                    - the codon each Nucleotide is part of
                    - the Nucleotides's position in the codon
    """

    def __init__(self, letter, position, seq, orfs):
        """
        :param letter: Nucleotide A, C, G or T
        :param position: Position of the nucleotide in the sequence
        :param sorted_orfs: <option> A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        """
        self.seq = seq
        self.letter = letter
        self.position = position
        self.seq_orfs = orfs
        self.orfs_for_nt = self.get_orfs_for_nt(self.seq_orfs)


    def get_orfs_for_nt(self, seq_orfs):
        """
        Store the information of the nucleotide for which it makes part of
        """
        orfs_for_nt =  {'+0':[], '+1': [], '+2': [],
                        '-0': [], '-1': [], '-2': []}

        # Find the orfs in which nt is involved
        for key, values in seq_orfs:
            if values:
                for i in values:
                    from_to = values[i]
                    max = max(from_to[0], from_to[1])
                    min = min(from_to[0], from_to[1])
                    # If position in orf
                    if pos in range(min, max):
                        # Store codon information
                        orfs_for_nt[key].append(get_codon(from_to[0], from_to[1]))
                        print(orfs_for_nt)

    def get_codon(self, initial, final):
        """
        Get codon sequence, and position of my_nt in the codon
        :param position: position of the nucleotide in the sequence
        :param orf: tuple indicating first and last nucleotide of an open reading frame
        :return codon, position_in_codon: tuple with nucleotide triplet and position of the nucleotide in the codon
        """

        if initial > final:  # positive strand
            my_orf = self.original_seq[initial:final + 1]
        else:  # negative strand
            my_orf = self.rcseq[initial: final + 1]

        position_in_orf = abs(initial - self.position)
        position_in_codon = position_in_orf % 3

        if position_in_codon == 0:
            codon = my_orf[position_in_orf: position_in_orf + 3]
        elif position_in_codon == 1:
            codon = my_orf[position_in_orf - 1: position_in_orf + 2]
        else:
            codon = my_orf[position_in_orf - 2: position_in_orf + 1]

        return codon, position_in_codon