# Store sequence information

import random
import copy
import sys

import pprint
from tempfile import tempdir
import numpy as np

TRANSITIONS_DICT = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

NUCLEOTIDES = ['A', 'C', 'G', 'T']

AMBIGUOUS_NUCLEOTIDES = { 'R': ['G', 'A'], # Purine
                          'Y': ['C', 'T'], # Pyrimidine
                          'K': ['G', 'T'], # Ketone
                          'M': ['A', 'C'], # Amino
                          'S': ['C', 'G'], # Strong
                          'W': ['A', 'T'], # Weak
                          'B': ['C', 'G', 'T'], # Not A
                          'D': ['A', 'G', 'T'], # Not C
                          'H': ['A', 'C', 'T'], # Not G
                          'V': ['A', 'C', 'G'], # Not T
                          'N': ['A', 'C', 'G', 'T'], # Any one base
                        }

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
    """

    def __init__(self, str_sequence, orfs, kappa, global_rate, pi, cat_values, circular=False):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.

        :param str_sequence:  str, nucleotide sequence as a string object
        :param orfs:  dict, A dictionary of open reading frames (ORFs) in the sequence, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list containing the information for each ORF
                        - for example:
            {'+0': [{'coords': [[2849, 3182]], 'omega_shape': 1.5, 'omega_classes': 3,
                     'omega_values': [0.17, 0.48, 1.14]},
                    {'coords': [[3173, 3182]], 'omega_shape': 1.9, 'omega_classes': 5,
                     'omega_values': [0.18, 0.40, 0.63, 0.93, 1.63]}],
             '+1': [],
             '+2': [{'coords': [[0, 837]], 'omega_shape': 1.7, 'omega_classes': 4,
                     'omega_values': [0.17, 0.42, 0.72, 1.40]},
                    {'coords': [[156, 837]], 'omega_shape': 1.2, 'omega_classes': 6,
                     'omega_values': [0.05, 0.16, 0.28, 0.44, 0.67, 1.26]}],
             '-0': [],
             '-1': [],
             '-2': []}
        :param kappa:  float, transition/ transversion rate ratio
        :param global_rate:  float, the global substitution rate (/site/unit time)
        :param pi:  float, stationary frequencies of nucleotides, with nucleotide as keys
        :param cat_values:  dict, values drawn from a gamma distribution modeling rate variation among nucleotides, keyed by 'mu1', etc.
        :param circular:  bool, true if the genome is circular, false if the genome is linear (default: false)
        """
        self.orfs = orfs
        self.kappa = kappa
        self.global_rate = global_rate
        self.pi = pi
        self.cat_values = cat_values
        self.is_circular = circular

        self.nt_sequence = []
        self.__codons = []  # Store references to all codons
        self.total_omegas = {}  # Omega combos and their values. Populated when set_substitution_rates for nucleotides
        self.number_orfs = []  # List of "None"s with lenght equal to number of ORFs in sequence 

        pp = pprint.PrettyPrinter(indent=2)
        # Create Nucleotides
        for pos_in_seq, nt in enumerate(str_sequence):
            self.nt_sequence.append(Nucleotide(nt, pos_in_seq))

        # Set Codons based on the reading frames
        if self.orfs is not None:
            for frame, orf_list in self.orfs.items():
                for orf in orf_list:  # orf is a dictionary
                    codons = self.find_codons(frame, orf)  # generates Codon objects
                    self.number_orfs.append(None)  

                    # tell Nucleotide which Codon(s) it belongs to
                    for codon in codons:
                        for nt in codon.nts_in_codon:
                            nt.codons.append(codon)  # FIXME: shouldn't Codon __init__ do this?
                        self.__codons.append(codon)

        """
        # Dict to find ORF combination. 
        Keys are combos where 1 represent presence and 0 absence of the ORF assigned in that position of the tuple. 
        Values are coordinates and length of the fragment.
        E.g. { (0, 0, 1, 0): {'coords': [[837, 1374], [2454, 2848]], 'len': 933}, 
               (0, 0, 1, 1): {'coords': [[0, 836], [2849, 3181]], 'len': 1170}, 
               (0, 1, 0, 0): {'coords': [[1840, 2307]], 'len': 468}}
        """
        self.regions = {}  # Sequence regions, divided depending on the cooding context
        non_orf = True
 
        for nt in self.nt_sequence:
            if not nt.codons:  # If nt not in a coding region
                orf_map = [0]*len(self.number_orfs)
            
            else:  # Get the orf map for the nucleotide including overlapping regions (E.g. ((0,1,1,0))
                non_orf = False
                orf_map = sum([codon.orf['orf_map'] for codon in nt.codons])

            orf_map_key = tuple(orf_map)
            nt.set_orf_map_key(orf_map_key)  # Store the map key on the nt object. This value will never change

            # Save map key into all maps dict
            if orf_map_key not in self.regions:
                self.regions.update({orf_map_key: {'coords':[[nt.pos_in_seq, nt.pos_in_seq]], 'len': 0 }})
                if not nt.codons and non_orf is False:
                    non_orf = True
            
            else:  # Update orf coordinates and length by summing the position of new nucleotide
                if not nt.codons and non_orf is False:
                    self.regions[orf_map_key]['coords'].append([nt.pos_in_seq, nt.pos_in_seq])
                    non_orf = True
                last_pos = self.regions[orf_map_key]['coords'][-1][1]  # End position in 'coord'
                
                if nt.pos_in_seq - last_pos > 1:  # Splitted gene, start new coordinates
                    self.regions[orf_map_key]['coords'].append([nt.pos_in_seq, nt.pos_in_seq])
                else:  # Not splitted, update last nucleotide in coords
                    self.regions[orf_map_key]['coords'][-1][1] = nt.pos_in_seq
                
                self.regions[orf_map_key]['len'] = sum([abs(coord[0]-coord[1])+1 for coord in
                                                         self.regions[orf_map_key]['coords']])

        self.length = sum([self.regions[orf]['len'] for orf in self.regions.keys()])  # Store sequence length
        self.event_tree = self.create_event_tree()  # Nested dict containing info about all possible mutation events

        # Calculate mutation rates for each nucleotide in sequence, populate the event tree which each nucleotide
        for nt in self.nt_sequence:
            # Get substitution rates for the nucleotide, store omega combos in self.total_omegas
            self.set_substitution_rates(nt)
            self.nt_in_event_tree(nt)  # Locate nucleotide in the event tree

        self.count_events_per_layer()
        # pp.pprint(self.event_tree)
        # pp.pprint(self.regions)
        # sys.exit()

    def all_syn_values(self, nonsyn_values):
        for nonsyn_val in nonsyn_values:
            if any(nonsyn_val):
                return True
        return False

    def count_events_per_layer(self):
        """
        Modify Event Tree by traversing it calculating and storing the number of events in every branch, layer by layer
        Number of events are required to select a branch using weighted_random_choice
        """
        for to_nt in NUCLEOTIDES:
            to_events = 0  # count total number of events of this category
            for from_nt in NUCLEOTIDES:
                if to_nt != from_nt:
                    from_events = 0
                    branch = self.event_tree['to_nt'][to_nt]['from_nt'][from_nt]
                    for cat in self.cat_values.keys():
                        cat_events = 0
                        for orf_region in branch[cat].keys():
                            orf_region_events = 0
                            if type(orf_region) is not tuple:
                                continue

                            region_weight = 0
                            for omega_combo in branch[cat][orf_region].keys():
                                if type(omega_combo) is not tuple:
                                    continue
                                events = len(branch[cat][orf_region][omega_combo])
                                self.total_omegas[omega_combo]['nt_events'] = events
                                orf_region_events += events    
                                cat_events += events
                                from_events += events
                                to_events += events
                                # number of events multiplied by net effect of omegas
                                region_weight += self.total_omegas[omega_combo]['value'] * events

                            branch[cat][orf_region]['nt_events'] = orf_region_events
                            branch[cat][orf_region]['region_weight'] = region_weight

                        branch[cat]['nt_events'] = cat_events

                    self.event_tree['to_nt'][to_nt]['from_nt'][from_nt]['nt_events'] = from_events
            
            self.event_tree['to_nt'][to_nt]['nt_events'] = to_events         

    def __deepcopy__(self, memodict):
        """
        Creates a deepcopy of Sequence and sets the reference(s) for a Nucleotide's Codon(s)
        """

        # Creates a new Sequence
        cls = self.__class__
        new_sequence = cls.__new__(cls)
        memodict[id(self)] = new_sequence  # Avoid duplicate copying

        # Set attributes of new Sequence to the same as the original object
        for k, v in self.__dict__.items():
            setattr(new_sequence, k, copy.deepcopy(v, memodict))

        # Set references to Codons
        for codon in new_sequence.__codons:
            for i, nt in enumerate(codon.nts_in_codon):
                nt.codons.append(codon)

        return new_sequence

    def __str__(self):
        """
        Represents the Sequence as a string by casting nt_sequence to a string.
        """
        return ''.join(str(nt) for nt in self.nt_sequence)

    def get_sequence(self):
        return self.nt_sequence

    def get_right_nt(self, pos_in_seq):
        """
        Returns the next Nucleotide in the sequence
        :param pos_in_seq: the position of the Nucleotide in the sequence
        """
        if self.is_circular:
            if pos_in_seq == len(self.nt_sequence) - 1:
                return self.nt_sequence[0]
        else:
            return self.nt_sequence[pos_in_seq + 1]

    def get_left_nt(self, pos_in_seq):
        """
        Returns the previous Nucleotide in the sequence
        :param pos_in_seq: the position of the Nucleotide in the sequence
        """
        if pos_in_seq == 0:
            if self.is_circular:
                return self.nt_sequence[- 1]
            else:
                return None
        else:
            return self.nt_sequence[pos_in_seq - 1]

    def get_event_tree(self):
        return self.event_tree

    @staticmethod
    def complement(seq, rev=False):
        """
        Generates the complement of a DNA sequence
        :param seq: the input sequence
        :param <option> rev: option to find the reverse complement
        :return s: The complement (or reverse complement) of the sequence
        """
        if rev:
            s = reversed(seq.upper())
        else:
            s = seq

        result = ''
        for i in s:
            result += COMPLEMENT_DICT[i]

        return result

    @staticmethod
    def get_frequency_rates(seq):
        """
        Frequency of nucleotides in the DNA sequence
        :param seq: the DNA sequence
        :return frequencies: a dictionary frequencies where the key is the nucleotide and the value is the frequency
        """
        frequencies = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

        for nucleotide in frequencies:
            frequencies[nucleotide] = round((float(seq.count(nucleotide)) / (len(seq))), 2)

        return frequencies

    def create_event_tree(self):
        """
        Create an event tree (nested dictionaries) that stores pi, kappa, mu,
        and information about whether a mutation is a transition or transversion.
        :return event tree: a nested dictionary containing information about the mutation event
        """
        event_tree = {'to_nt': dict([(nuc, {'from_nt': {}}) for nuc in NUCLEOTIDES])}
        
        # dictionary keyed by binary tuples indicating combination of orfs (e.g., {(1, 0): {}, (0, 1): {} , (1, 1): {}})
        bin_orf_layer = {bin_code : {} for bin_code in [key for  key in self.regions.keys()]}
        
        # dictionary keyed by mutation rate categories
        # (e.g., {'mu1': {(1, 0): {}, (0, 1): {}}, 'mu2': {(1, 0): {}, (0, 1): {}}})
        cat_dict = {cat : copy.deepcopy(bin_orf_layer) for cat in self.cat_values.keys()}

        for to_nt in NUCLEOTIDES:
            # Update nucleotides with possible mutations
            for from_nt in NUCLEOTIDES:
                if from_nt == to_nt:
                    event_tree['to_nt'][to_nt]['from_nt'][from_nt] = None
                else: 
                    event_tree['to_nt'][to_nt]['from_nt'][from_nt] = copy.deepcopy(cat_dict)
           
        return event_tree

    def set_substitution_rates(self, nt):
        """
        Calculates substitution rates of a nucleotide
        Store newly created omega_combos in self.total_omegas as key with their value E.g {(None, 3,1) = (1*0.2*0.6)}
        Sets the sub-rates, omega values, rate classes, and total mutation rate of a Nucleotide
        :param nt: object of class Nucleotide
        """
        current_nt = nt.state
        sub_rates = {}
        selected_omegas = {}
        my_cat_keys = {}

        for to_nt in NUCLEOTIDES:
             
            if to_nt == current_nt:
                sub_rates[to_nt] = None
                selected_omegas[to_nt] = None
            
            else:
                # If mutation does not create STOP codons or affects START or STOP codons
                if not self.is_start_stop_codon(nt, to_nt):  

                    # Apply global substitution rate and stationary nucleotide frequency
                    sub_rates[to_nt] = self.global_rate * self.pi[current_nt]
                    if self.is_transv(current_nt, to_nt):
                        sub_rates[to_nt] *= self.kappa

                    chosen_omegas = list(self.number_orfs)  # Initialize chosen omegas as list of "None"s with len equal to number of ORFs
                    computed_omega = 1
                    # If nucleotide is part of at least one ORF, select omegas if mutation is non-syn or use "None" if it's syn
                    if nt.codons:
                    # For each codon, find the combination of omegas that will affect it according to the codons it is part of
                    # E.g.: [2, 3, None, None, None] represent a nucleotide in two reading frames where both mutations are non-syn. 
                    # Note: "None" means that nt is not part of that ORF

                        for codon in nt.codons:

                            codon_orf = codon.orf['orf_map']  # tuple of 1s and 0s Eg., (1, 0, 0)
                            orf_index = np.where(codon_orf==1)[0][0]  # The position with a 1 indicates the ORF for the codon
                            pos_in_codon = codon.nt_in_pos(nt)
                            
                            # If mutation is non-synonymous, apply codon omega
                            if codon.is_nonsyn(pos_in_codon, to_nt):
                                chosen_omegas[orf_index] = codon.omega
                                computed_omega *= codon.omega                          

                            # Use a -1 to indicate mutation is synonymous
                            else:
                                chosen_omegas[orf_index] = -1

                    # Store omega combination and calculated value in total_omega dict
                    if tuple(chosen_omegas) not in self.total_omegas.keys():
                        self.total_omegas[tuple(chosen_omegas)] = {'value' : computed_omega}

                    selected_omegas[to_nt] = tuple(chosen_omegas)  # Store omega keys used to describe the substitution
                    selected_cat = random.choice(list(self.cat_values))  # Randomly select one of the mu values (mutation rate) 
                    sub_rates[to_nt] *= self.cat_values[selected_cat]  # Apply my value over instant mutation rate
                    sub_rates[to_nt] *= self.total_omegas[tuple(chosen_omegas)]['value']  # Apply omega value over instant mutation rate
                    my_cat_keys[to_nt] = selected_cat
                
                else:  # Inform nucleotide that such subs cannot occur
                    selected_omegas[to_nt] = None
                    my_cat_keys[to_nt] = None
                    sub_rates[to_nt] = None

        # Set substitution rates and key values for the nucleotide object
        nt.set_rates(sub_rates)
        nt.set_categories(my_cat_keys)
        nt.set_omega(selected_omegas)
        nt.get_mutation_rate()  # Sum of mutation rates for all nucleotides is used to calculate rate at which mutations occurs on simulation.py

    @staticmethod
    def is_start_stop_codon(nt, to_nt):
        """"
        Check if mutation is a STOP codon or nucleotide belongs to a START codon
        :param nt: a Nucleotide object
        :param to_nt: the new state of the Nucleotide as a string
        :return: False if mutation does not create a STOP in any of the codons the nucleotide is part of
        """
        for codon in nt.codons:
            if codon.is_stop() or codon.is_start() or codon.creates_stop(codon.nt_in_pos(nt), to_nt):
                return True

        return False

    def nt_in_event_tree(self, nt):
        """
        Store nucleotide in each branch of the Event Tree where it belongs
        :param nt: a Nucleotide object
        :return: the new omega key
        """
        current_nt = nt.state

        for to_nt in NUCLEOTIDES:
            
            if to_nt != current_nt:
                
                if not self.is_start_stop_codon(nt, to_nt):  # Nucleotides that are part of start or stop should never mutate

                    # Find proper keys
                    category = nt.cat_keys[to_nt]
                    orf_map_key = nt.orf_map_key
                    omega_keys = nt.omega_keys[to_nt]

                    branch = self.event_tree['to_nt'][to_nt]['from_nt'][current_nt][category][orf_map_key]
                    
                    if omega_keys in branch.keys():  # Store nt on the branch
                        branch[omega_keys].append(nt)
                    
                    else:  # Create the omega layer when non existent
                        branch[omega_keys] = [nt]
                                    

    @staticmethod
    def is_transv(from_nt, to_nt):
        """
        Checks if a mutation is a transition or a transversion
        :param from_nt: the current nucleotide
        :param to_nt: the new nucleotide
        :return transv: True if the mutation is a transversion,
                        False if the mutation is a transition,
                        None if the current and new nucleotides are the same
        """
        if from_nt == to_nt:
            transv = None
        else:
            transv = True
            if TRANSITIONS_DICT[from_nt] == to_nt:
                transv = False
        return transv

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
            yield my_orf[i:i + 3]
            i += 3

    def find_codons(self, frame, orf):
        """
        Gets the Codon sequence
        :param frame:  str, the frame of the ORF, e.g., "+1"
        :param orf:  dict, containing the coordinates of the ORF, and the associated omega values
        :return: a list of Codon objects for the specified ORF
        """
        # extract coding sequence
        cds = []
        for start, stop in orf['coords']:
            cds.extend(self.nt_sequence[start:stop])  # concatenates spliced ORFs
        if frame.startswith('-'):
            cds = cds[::-1]  # negative strand ORF

        # Iterate over string by threes and create Codon objects
        codons = []
        for i in range(3, len(cds)+1, 3):
            codons.append(Codon(frame, orf, cds[(i-3):i]))

        return codons

    def check_event_tree(self):
        """
        When debugging, useful to check if nucleotides are being properly stored on the Event Tree
        """
        for key1, to_nt in self.event_tree['to_nt'].items():
            subset = to_nt['from_nt']

            for key2, from_nt in subset.items():
                if key2 != 'T' and from_nt and from_nt.get('nts_in_subs'):
                    nts_in_subs = list(from_nt['nts_in_subs'].keys())
                    if len([1 for tip in nts_in_subs if str(tip).lower() == 't0']) > 0:
                        meta2 = {'nts_in_subs': nts_in_subs}
                        print(f'>>>>>>>>>>>> meta2: from {key2}', meta2)
                        sys.exit(1)


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    """

    def __init__(self, state, pos_in_seq):
        """
        :param state:  str, nucleotide A, C, G or T
        :param pos_in_seq:  int, position of the nucleotide in the sequence
        """
        self.state = state
        self.pos_in_seq = pos_in_seq

        self.codons = []  # A list of codon objects the Nucleotide is part of
        self.complement_state = COMPLEMENT_DICT[self.state]  # The complement state
        self.rates = {}  # A dictionary of mutation rates
        self.omega_keys = {}  # For each substitution, tuple containing the omegas applied Eg. ({'A': None, 'C': [0, 3, -1, -1, -1], 'G': [1, 0, -1, -1, -1], 'T': [0, 0, -1, -1, -1]})
        self.cat_keys = {}  # category keys chosen when calculating rates
        self.mutation_rate = 0  # The total mutation rate
        self.orf_map_key = ()

    def __str__(self):
        return self.state

    def __repr__(self):
        return self.state.lower() + str(self.pos_in_seq)

    def __deepcopy__(self, memodict):
        """
        Creates a deepcopy of a Nucleotide.
        Note: A Nucleotide's reference(s) to its Codon(s) will be set in Sequence's deepcopy
        """

        # Creates a new Nucleotide
        cls = self.__class__
        new_nucleotide = cls.__new__(cls)
        memodict[id(self)] = new_nucleotide  # Avoid duplicate copying

        # Copy all attributes except the codons
        new_nucleotide.state = copy.deepcopy(self.state, memodict)
        new_nucleotide.pos_in_seq = copy.deepcopy(self.pos_in_seq, memodict)
        new_nucleotide.complement_state = copy.deepcopy(self.complement_state, memodict)
        new_nucleotide.rates = copy.deepcopy(self.rates, memodict)
        new_nucleotide.mutation_rate = copy.deepcopy(self.mutation_rate, memodict)
        new_nucleotide.omega_keys = copy.deepcopy(self.omega_keys, memodict)
        new_nucleotide.cat_keys = copy.deepcopy(self.cat_keys, memodict)
        new_nucleotide.codons = []  # References to Codons will be set when the Sequence is deep-copied
        new_nucleotide.orf_map_key = copy.deepcopy(self.orf_map_key, memodict)

        return new_nucleotide

    def set_orf_map_key(self, orf_map_key):
        self.orf_map_key = orf_map_key

    def set_state(self, new_state):
        self.state = new_state

    def get_complement_state(self):
        return self.complement_state

    def set_complement_state(self):
        self.complement_state = COMPLEMENT_DICT[self.state]

    def set_rates(self, rates):
        self.rates = rates

    def set_omega(self, omega_keys):
        self.omega_keys = omega_keys

    def set_categories(self, cat_keys):
        self.cat_keys = cat_keys

    def add_codon(self, codon):
        self.codons.append(codon)

    def set_mutation_rate(self, mutation_rate):
        self.mutation_rate = mutation_rate

    def get_mutation_rate(self):
        total_rate = 0
        for mutation_rate in self.rates.values():
            if mutation_rate:
                total_rate += mutation_rate
        self.mutation_rate = total_rate


class Codon:
    """
    Stores information about the frameshift, ORF, and pointers to 3 Nucleotide objects
    """

    def __init__(self, frame, orf, nts_in_codon):
        """
        Create a Codon
        :param frame:  str, the reading frame (+0, +1, +2, -0, -1, -2)
        :param orf:  tuple, containing the reading frame and the coordinates of the orf
        :param nts_in_codon:  a list of pointers to the Nucleotides in the Codon
        """
        self.frame = frame
        self.orf = orf
        self.nts_in_codon = nts_in_codon  # list of Nucleotide objects in the Codon
        self.omega = self.select_omega()  # Asign non-syn mutation rate 

    def __repr__(self):
        return ''.join(str(nt) for nt in self.nts_in_codon)

    def nt_in_pos(self, query_nt):
        """
        Finds the position of the Nucleotide in the Codon
        :param query_nt: the Nucleotide of interest
        :return: the position of the Nucleotide in the Codon
        """
        # FIXME: isn't this just self.nts_in_codon.index(query_nt) ?
        for idx, nt in enumerate(self.nts_in_codon):
            if query_nt is nt:
                return idx
    
    def select_omega(self):
        """
        Set rate at which non-synonymous substitutions occur in the codon
        """
        omega_values = self.orf['omega_values']
        codon_omega = random.choice(omega_values)
        return(codon_omega)

    def mutate_codon(self, pos_in_codon, to_nt):
        """
        Changes the state of the specified nucleotide in the codon
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide
        :return codon, mutated_codon: the codon and mutated codon represented as lists of strings
        """
        codon = [str(nt) for nt in self.nts_in_codon]  # Cast all Nucleotides in the Codon to strings
        mutated_codon = codon.copy()
        mutated_codon[pos_in_codon] = to_nt

        return codon, mutated_codon

    def is_nonsyn(self, pos_in_codon, to_nt):
        """
        Finds if a substitution at the specified position results in a non-synonymous mutation
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to a non-synonymous mutation,
                 False if the substitution leads to a synonymous mutation
        """
        codon, mutated_codon = self.mutate_codon(pos_in_codon, to_nt)
        return CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]

    def creates_stop(self, pos_in_codon, to_nt):
        """
        Finds if a substitution at the specified position results in a stop codon
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to stop codon,
                 False if the substitution doesn't lead to a stop codon
        """
        codon, mutated_codon = self.mutate_codon(pos_in_codon, to_nt)
        return CODON_DICT[''.join(mutated_codon)] == "*"

    def is_start(self):
        """
        Checks if the codon is a start codon
        :return True of the codon is a start codon, False otherwise
        """
        codon = ''.join(str(nt) for nt in self.nts_in_codon)  # Cast all Nucleotides in the Codon to strings
        if self.frame.startswith('+'):
            return codon == 'ATG' and self.nts_in_codon[0].pos_in_seq == self.orf['coords'][0][0]
            # Note: I don't think the codon needs to be ATG to be a start?
            #return self.nts_in_codon[0].pos_in_seq == self.orf['coords'][0][0]
        else:
            # +1 to account for non-inclusive indexing
            return codon == 'ATG' and self.nts_in_codon[0].pos_in_seq + 1 == self.orf['coords'][0][1]
            #return self.nts_in_codon[0].pos_in_seq + 1 == self.orf['coords'][0][1]
            
    def is_stop(self):
        """
        Checks if the codon is a STOP codon
        :return True of the codon is a STOP codon, False otherwise
        """
        codon = ''.join(str(nt) for nt in self.nts_in_codon)  # Cast all Nucleotides in the Codon to strings
        return codon == 'TAA' or codon == 'TGA' or codon == 'TAG'