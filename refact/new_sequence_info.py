# Store sequence information

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
    :ivar unsorted_orfs: list of ORFs (as tuples) classified according to their reading frame shift relative
                        to the first orf (+0, +1, +2, -0, -1, -2)
    :ivar nt_sequence: list of Nucleotide objects
    """

    def __init__(self, original_seq, unsorted_orfs=None):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param original_seq: A list of Nucleotides
        :param unsorted_orfs: <option> A dictionary of ORFs, sorted by reading frame where:
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
        self.nt_sequence = DoubleLinkedList()
        for pos, nt in enumerate(self.original_seq):
            new_dict = self.create_nt_orf_dict(pos)
            self.nt_sequence.insert_nt(nt, pos, new_dict)

        self.nt_sequence.print_seq()

    def get_sequence(self):
        return self.nt_sequence

    def create_nt_orf_dict(self, pos):
        """
        Dictionary with keys as orfs (in which nt is involved) and items as a number (0,1 or 2) representing
        position of the nucleotide in the codon
        :param pos: position of the nucleotide in the sequence
        :return: position of the nucleotide in the codon (e.g. 0, 1 or 2), given an orf
        """
        nt_orf_dict = {'+0': None, '+1': None, '+2': None, '-0': None, '-1': None, '-2': None}

        for frame_shift, orfs_list in self.orfs.items():
            for orf_tuple in orfs_list:
                # if nucleotide is inside the orf
                if orf_tuple[0] <= pos <= orf_tuple[1] or orf_tuple[1] <= pos <= orf_tuple[0]:
                    position_in_codon = abs(orf_tuple[0] - pos) % 3
                    nt_orf_dict[frame_shift] = position_in_codon

        return nt_orf_dict

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
            if orf[1] > orf[0]:  # Forward strand
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
        :param: my_region <option> A sub-sequence of the original sequence
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
        :return sorted_orfs: List of ORFs classified according to their shift relative to
                    the plus zero reading frame  (+0, +1, +2, -0, -1, -2)
        """
        sorted_orfs = {'+0': [], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

        if unsorted_orfs:
            first_orf = unsorted_orfs[0]
            for orf in unsorted_orfs:
                difference = abs(orf[0] - first_orf[0]) % 3

                if first_orf[0] < first_orf[1]:
                    if orf[0] < orf[1]:  # positive strand
                        if difference == 0:
                            sorted_orfs['+0'].append(orf)
                        elif difference == 1:
                            sorted_orfs['+1'].append(orf)
                        elif difference == 2:
                            sorted_orfs['+2'].append(orf)

                    elif orf[0] > orf[1]:  # negative strand
                        if difference == 0:
                            sorted_orfs['-2'].append(orf)
                        elif difference == 1:
                            sorted_orfs['-1'].append(orf)
                        elif difference == 2:
                            sorted_orfs['-0'].append(orf)

                else:
                    if orf[0] < orf[1]:  # positive strand
                        if difference == 0:
                            sorted_orfs['+2'].append(orf)
                        elif difference == 1:  # plus one
                            sorted_orfs['+1'].append(orf)
                        elif difference == 2:  # plus two
                            sorted_orfs['+0'].append(orf)

                    elif orf[0] > orf[1]:  # negative strand
                        if difference == 0:
                            sorted_orfs['-0'].append(orf)
                        elif difference == 1:
                            sorted_orfs['-1'].append(orf)
                        elif difference == 2:
                            sorted_orfs['-2'].append(orf)

        return sorted_orfs


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    :ivar state: the nucleotide base (A, T, G, C)
    :ivar position: the position of the nucleotide relative to the start of the sequence
    :ivar pos_in_codon: a dictionary where:
                    - keys are the reading frames
                    - values are the nucleotide's position in the reading frame
    :ivar left_nt: reference to the nucleotide's left neighbour
    :ivar right_nt: reference to the nucleotide's right neighbour
    """

    def __init__(self, state, position, pos_in_codon={}, left_nt=None, right_nt=None):
        """
        :param state: Nucleotide A, C, G or T
        :param position: Position of the nucleotide in the sequence
        :param left_nt : reference to the adjacent nucleotide to the left (default to None)
        :param right_nt: reference to the adjacent nucleotide to the right (default to None)
        """
        self.state = state
        self.pos = position
        self.left_nt = left_nt
        self.right_nt = right_nt
        self.pos_in_codons = pos_in_codon
        # self.syn_nonsyn = self.get_nonsyn_subs()

    def get_state(self):
        return self.state

    def get_pos(self):
        return self.pos

    def get_left_nt(self):
        return self.left_nt

    def get_right_nt(self):
        return self.right_nt

    def get_pos_in_codons(self):
        return self.pos_in_codons

    def set_state(self, new_state):
        self.state = new_state

    def set_pos(self, new_pos):
        self.pos = new_pos

    def set_left_nt(self, new_left_nt):
        self.left_nt = new_left_nt

    def set_right_nt(self, new_right_nt):
        self.right_nt = new_right_nt

    def set_pos_in_codons(self, new_codons):
        self.pos_in_codons = new_codons
    #
    # def get_nonsyn_subs(self):
    #     """
    #     Get synonymous and non-synonymous substitutions given the open reading frame
    #     :return: a dictionary with possible mutations for each reading frame.
    #             The value is True if the substitution is non-synonymous, otherwise it is False
    #     """
    #
    #     nonsyn_syn = {}
    #     for frame in self.pos_in_codons:
    #         nonsyn_syn[frame] = {'A': None, 'T': None, 'G': None, 'C': None}
    #
    #         # Get the codon to which the Nucleotide belongs
    #         codon = []
    #         if self.pos_in_codons[frame] is not None:
    #
    #             if self.pos_in_codons[frame] == 0:
    #                 right = self.get_right_nt()
    #
    #                 # Check that right neighbours exits
    #                 if right is not None:
    #                     right_right = right.get_right_nt()
    #
    #                     if right_right is not None:
    #                         codon.append(self.get_state())
    #                         codon.append(right.get_state())
    #                         codon.append(right_right.get_state())
    #
    #             elif self.pos_in_codons[frame] == 1:
    #                 left = self.get_left_nt()
    #                 right = self.get_right_nt()
    #
    #                 # Check that left and right neighbours exist
    #                 if left and right is not None:
    #                     codon.append(left.get_state())
    #                     codon.append(self.get_state())
    #                     codon.append(right.get_state())
    #
    #             else:
    #                 left = self.get_left_nt()
    #                 # Check that left neighbour exists
    #                 if left is not None:
    #                     left_left = left.get_left_nt()
    #
    #                     if left_left is not None:
    #                         codon.append(left_left.get_state())
    #                         codon.append(left.get_state())
    #                         codon.append(self.get_state())
    #
    #             print(codon)
    #             if codon:
    #                 for nt in NUCLEOTIDES:
    #                     mutated_codon = codon
    #                     mutated_codon[self.pos_in_codons[frame]] = nt
    #
    #                     if CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]:
    #                         nonsyn_syn[frame][nt] = True    # Non-synonymous
    #                     else:
    #                         nonsyn_syn[frame][nt] = False   # Synonymous
    #
    #     return nonsyn_syn


class DoubleLinkedList:
    """
    Double linked list linking together objects of class Nucleotide
    default initialization with empty head node
    """

    def __init__(self):
        self.head = None  # head node (starting nucleotide)
        self.current_nt = None  # Pointer to current nt for insertion

    def insert_nt(self, state, position, nt_codon_dict):
        """
        Insert objects of class Nucleotide to the end of the DoubleLinkedList
        :param state: Nucleotide state in sequence
        :param position: Position of nt in sequence
        :param nt_codon_dict: Dictionary with nucleotide position inside the codon given a frame shift (0, 1 or 2)
        """
        new_nt = Nucleotide(state, position)  # create new Nucleotide object
        new_nt.set_pos_in_codons(nt_codon_dict)

        # Assign the first nucleotide as head
        if self.head is None:
            self.head = new_nt
            self.current_nt = new_nt

        else:
            new_nt.set_left_nt(self.current_nt)  # For the new nucleotide, create a left pointer towards the current one
            self.current_nt.set_right_nt(new_nt)  # Create the double link between current and new
            self.current_nt = new_nt

    def get_head(self):
        return self.head

    def print_seq(self):  # Print the string of nucleotides (check the class is working properly)
        temp = self.head
        while temp is not None:
            # print(temp.get_letter(), temp.get_pos())
            temp = temp.get_right_nt()
