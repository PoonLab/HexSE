# Store sequence information
import scipy
import scipy.stats as ss
import numpy as np
import re
import sys
from refact.event_tree import Event_tree


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

    def __init__(self, original_seq, mu, kappa, unsorted_orfs=None, pi=None, omega=None):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param original_seq: A list of Nucleotides
        :param mu: The global rate (substitutions/site/unit time)
        :param kappa: transition transversion rate ratio
        :param unsorted_orfs: <optional> A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        :param pi: <optional> a vector of stationary nucleotide frequencies.  If not specified, defaults
                   to the empirical frequencies in <seq>
        :param omega: <optional> a List of dN/dS (omega) ratios along sequence length as {dict}.
                      {dict} always contains key +0 (parent reading frame).  May contain omega for alternate
                      reading frame as (+1, +2, -0,sim_ -1 or -2). Codon position is determined by the nt's
                      location relative to start of <seq>
        :raises ValueError: If the sequence contains invalid characters or if it contains no ORFs
        """

        # Check if sequence is valid
        if not Sequence.valid_sequence(original_seq):
            print("Invalid sequence: {}".format(original_seq))
            sys.exit(0)

        self.original_seq = original_seq
        self.orfs = {}
        self.rcseq = self.reverse_and_complement(original_seq)
        self.mu = mu
        self.pi = pi
        self.kappa = kappa
        self.omega = omega

        # Calculate stationary frequency rates
        if self.pi is None:
            self.pi = self.get_frequency_rates(self.original_seq)

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

        # Create event tree for sequence (tree containing all possible mutation and the parameters that should be taken into account when calculating rates)
        self.tree = Event_tree(self.pi)
        self.event_tree = self.tree.create_event_tree()

        # Draw omega values from gamma distrubution
        omega_values = self.get_omega_values(2, 4)

        # Create Nucleotides
        self.nt_sequence = DoubleLinkedList()
        for pos_in_seq, nt in enumerate(self.original_seq):
            self.nt_sequence.insert_nt(nt, pos_in_seq)


        # Set Codons based on the reading frames
        for frame in self.orfs:
            orf_list = self.orfs[frame]
            for orf in orf_list:
                codons = self.find_codons(frame, orf)

                # Tell Nucleotide which Codon(s) it belongs to
                for codon in codons:
                    for i, nt in enumerate(codon.nts_in_codon):
                        nt.add_codon(codon)

        # TODO: Calculate mutation rates for each nucleotide according to event tree and codon info

    def get_sequence(self):
        return self.nt_sequence

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

    @staticmethod
    def codon_iterator(my_orf, start_pos, end_pos):
        """
        Generator to move every three nucleotides (codon)
        :param my_orf: A list of Nucleotides in the ORF
        :param start_pos: The start position of the ORF
        :param end_pos: The end position of the ORF
        :yield codon
        """
        if start_pos > end_pos:  # Negative strand
            my_orf.reverse()

        i = 0
        while i < len(my_orf):
            yield[my_orf[i:i + 3]]
            i += 3

    def find_codons(self, frame, orf):
        """
        Gets the Codon sequence
        :param frame: the frame of the ORF
        :param orf: tuple containing the coordinates of the ORF
        :return: a list of Codon objects for the specified ORF
        """
        codons = []
        start_pos = orf[0]
        end_pos = orf[1]
        my_orf = self.nt_sequence.slice_sequence(start_pos, end_pos)

        if start_pos < end_pos:
            # Iterate over list by threes and create Codons in the forward strand
            for cdn in self.codon_iterator(my_orf, start_pos, end_pos):
                codon = Codon(orf, frame, cdn)
                codons.append(codon)

        else:
            # TODO: Handle Nucleotides involved in reverse strand ORFs (Issue #41)
            pass

        return codons

    @staticmethod
    def get_frequency_rates(seq):
        """
        Frequency of nucleotides in the DNA sequence
        :param seq: the DNA sequence
        :return: a dictionary frequencies where the key is the nucleotide and the value is the frequency
        """
        frequencies = {
            'A': 0,
            'C': 0,
            'T': 0,
            'G': 0
        }

        for nucleotide in frequencies:
            frequencies[nucleotide] = round((float(seq.count(nucleotide)) / (len(seq))), 2)

        return frequencies

    def get_omega_values(self, alpha, ncat):
        """
        Draw ncat number of omega values from a discretized gamma distribution
        :param alpha: shape parameter
        :param ncat: Number of categories (expected omegas)
        :return: dictionary with number of categories as keys (e.i. {0: 0.29, 1: 0.65, 2: 1.06})
        """
        values = self.discretize_gamma(alpha = alpha, ncat = ncat)
        omega_values = list(values)
        return omega_values

    def discretize_gamma(self, alpha, ncat, dist=ss.gamma):
        """
        Divide the gamma distribution into a number of intervals with equal probability and get the mid point of those intervals
        From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
        :param alpha: shape parameter
        :param ncat: Number of categories
        :param dist: function from scipy stats
        :return: array with ncat number of values
        """
        if dist == ss.gamma:
            dist = dist(alpha, scale=1 / alpha)
        elif dist == ss.lognorm:
            dist = dist(s=alpha, scale=np.exp(0.5 * alpha ** 2))
        quantiles = dist.ppf(np.arange(0, ncat) / ncat)
        rates = np.zeros(ncat, dtype=np.double) # return a new array of shape ncat and type double
        for i in range(ncat - 1):
            rates[i] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                   quantiles[i], quantiles[i + 1])[0]
        rates[ncat - 1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                      quantiles[ncat - 1], np.inf)[0]
        return rates


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    :ivar state: the nucleotide base (A, T, G, C)
    :ivar pos_in_seq: the position of the nucleotide relative to the start of the sequence
    :ivar left_nt: reference to the nucleotide's left neighbour
    :ivar right_nt: reference to the nucleotide's right neighbour
    """

    def __init__(self, state, pos_in_seq, left_nt=None, right_nt=None):

        """
        :param state: Nucleotide A, C, G or T
        :param pos_in_seq: Position of the nucleotide in the sequence
        :param left_nt : reference to the adjacent nucleotide to the left (default to None)
        :param right_nt: reference to the adjacent nucleotide to the right (default to None)
        """
        self.state = state
        self.pos_in_seq = pos_in_seq
        self.left_nt = left_nt
        self.right_nt = right_nt
        self.codons = []
        self.complement_state = COMPLEMENT_DICT[self.state]

    def __repr__(self):
        return self.state

    def get_state(self):
        return self.state

    def get_pos_in_seq(self):
        return self.pos_in_seq

    def get_left_nt(self):
        return self.left_nt

    def get_right_nt(self):
        return self.right_nt

    def get_complement_state(self):
        return self.complement_state

    def set_state(self, new_state):
        self.state = new_state

    def set_pos_in_seq(self, new_pos):
        self.pos_in_seq = new_pos

    def set_left_nt(self, new_left_nt):
        self.left_nt = new_left_nt

    def set_right_nt(self, new_right_nt):
        self.right_nt = new_right_nt

    def add_codon(self, codon):
        self.codons.append(codon)


class DoubleLinkedList:
    """
    Double linked list linking together objects of class Nucleotide
    default initialization with empty head node
    """

    def __init__(self):
        self.head = None  # head node (starting nucleotide)
        self.current_nt = None  # Pointer to current nt for insertion

    def insert_nt(self, state, position):
        """
        Insert objects of class Nucleotide to the end of the DoubleLinkedList
        :param state: Nucleotide state in sequence
        :param position: Position of nt in sequence
        """
        new_nt = Nucleotide(state, position)  # create new Nucleotide object

        # Assign the first nucleotide as head
        if self.head is None:
            self.head = new_nt
            self.current_nt = new_nt

        else:
            new_nt.set_left_nt(self.current_nt)  # For the new nucleotide, create a left pointer towards the current one
            self.current_nt.set_right_nt(new_nt)  # Create the double link between current and new
            self.current_nt = new_nt

    def __iter__(self):
        return self.head()

    def __next__(self):
        current_nt = self.head()
        while current_nt is not None:
            return current_nt.get_right_nt()

    def get_head(self):
        return self.head

    def print_seq(self):  # Print the string of nucleotides (check the class is working properly)
        s = ''
        temp = self.head

        while temp is not None:
            s += temp.get_state()
            temp = temp.get_right_nt()
        print(s)

    def nucleotide_at_pos(self, position):
        """
        Traverse sequence to find the Nucleotide object at a specific position
        :param position: the position of the Nucleotide
        :return: the Nucleotide object in the specified position
        """
        current_nt = self.get_head()
        while current_nt is not None:
            if position == current_nt.get_pos_in_seq():
                return current_nt
            current_nt = current_nt.get_right_nt()
        return None

    def slice_sequence(self, start_pos, end_pos):
        """
        Slices the Nucleotide sequence
        :param start_pos: the start position
        :param end_pos: the end position
        :return sub_seq: a list of Nucleotides between the start and end positions
        """

        sub_seq = []
        curr_nt = self.nucleotide_at_pos(start_pos)

        if start_pos > end_pos:  # Positive strand
            # If there is a next Nucleotide and the position is within the range
            while curr_nt.get_right_nt() is not None and curr_nt.get_pos_in_seq() <= end_pos:
                sub_seq.append(curr_nt)
                curr_nt = curr_nt.get_right_nt()
        else:   # Negative strand
            # If there is a previous Nucleotide and the position is within the range
            while curr_nt.get_left_nt() is not None and curr_nt.get_pos_in_seq() >= end_pos:
                sub_seq.append(curr_nt)
                curr_nt = curr_nt.get_left_nt()

        return sub_seq


class Codon:
    """
    Stores information about the frameshift, ORF, omega, and pointers to 3 Nucleotide objects
    """

    def __init__(self, frame, orf, ptrs_to_nts):
        """
        Create a Codon
        :param frame: the reading frame (+0, +1, +2, -0, -1, -2)
        :param orf: a tuple containing the reading frame and the coordinates of the orf
        :param ptrs_to_nts: a list of pointers to the Nucleotides in the Codon
        """

        self.frame = frame
        self.orf = orf
        self.nts_in_codon = ptrs_to_nts

    def nt_in_pos(self, query_nt):
        """
        Finds the position of the Nucleotide in the Codon
        :param query_nt: the Nucleotide of interest
        :return: the position of the Nucleotide in the Codon
        """
        for idx, nt in enumerate(self.nts_in_codon):
            if query_nt == nt:
                return idx

    def is_nonsyn(self, pos, to_nt):
        """
        Finds if a substitution at the specified position results in a non-synonymous mutation
        :param pos: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to a non-synonymous mutation,
                 False if the substitution leads to a synonymous mutation
        """
        codon = [str(nt) for nt in self.nts_in_codon]    # Cast all Nucleotides in the Codon to strings

        mutated_codon = codon.copy()
        mutated_codon[pos] = to_nt

        if CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]:
            return True
        else:
            return False
