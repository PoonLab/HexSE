# Simulate evolution in a sequence with overlapping reading frames according to the new data structure from new_sequence_info.py
import numpy as np
import random
from new_sequence_info import Sequence
from Bio import Phylo
import copy

class SimulateOnBranch:
    """
    Simulate evolution within a sequence throughout a branch in a phylogeny
    """

    def __init__(self, sequence , branch_length):
        """
        :param sequence: object of class Sequence (Double Linked list of Nucleotides)  in the parental node
        :param branch_length: length of the branch over which evolution is happening
        """
        self.sequence  = sequence
        self.branch_length = branch_length
        self.event_tree = self.sequence.event_tree # Tree of class EventTree (all possible mutations related to sequence)


    def get_substitution(self):
        """
        Select a substitution by moving over the event_tree according to the generation of random numbers
        """
        # Select: to nucleotide
        events = self.event_tree['total_events']
        to_mutation = self.select_weighted_values(self.event_tree['to_nt'], events, 'events_for_nt', 'stationary_frequency')

        # Select: possible from nucleotides
        from_dict = self.event_tree['to_nt'][to_mutation]['from_nt']
        events_for_nt = self.event_tree['to_nt'][to_mutation]['events_for_nt']
        from_mutation = self.select_weighted_values(from_dict, events_for_nt, 'number_of_events', 'kappa')
        final_mutation = self.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]

        # List of nucleotides that are candidates to mutate
        candidate_nts = final_mutation['nts_in_subs']
        rates_list = [nt.mutation_rate for nt in candidate_nts]
        nt_dict = dict(zip(candidate_nts, rates_list))

        # Select weighted nucleotide
        from_nucleotide = self.weighted_random_choice(nt_dict, sum(rates_list))

        return from_nucleotide, to_mutation

    def select_weighted_values(self, dictionary, number_of_total_events, key_for_local_events, key_to_weight):
        """
        Randomly selected a key from a dictionary of events with weighted values
        :param number_of_total_events: Number of total number of events on branch
        :param number_of_local_events: key to the number of events for each specific value
        :param key_to_weight: Depending on the level of the branch, this could be:
                - the key to stationary frequency
                - the key to transition/transversion rate ratio
        """

        total_events = number_of_total_events
        temp = {}
        sum_values = 0
        # Create a temp dictionary to store the weighted values
        for key, value in dictionary.items():
            if value:
                # weight values
                temp[key] = (value[key_for_local_events]/total_events)*value[key_to_weight]
                sum_values += temp[key]
            else:
                temp[key] = None

        # Randomly selected a key on the dictionary
        return self.weighted_random_choice(temp, sum_values)

    @staticmethod
    def weighted_random_choice(dictionary, sum_values):
        """
        Randomly select a key on a dictionary where values correspond to the weight for the key
        :param dictionary: Dictionary to select the value from
        :param sum_values: sum all values on dict to establish the limit for the mutation
        :return: random key from dict
        """
        iter_object = iter(dictionary.items())
        limit = random.uniform(0, sum_values)
        s = 0
        while s < limit:
            try:
                (key, value) = next(iter_object)
                if value:
                    s += value
                else:
                    pass
            except StopIteration:
                break

        return key

    def sum_rates(self):
        """
        Calculate the total mutation rate of sequence
        """
        total_rate = sum([nt.mutation_rate for nt in iter(self.sequence.nt_sequence)])
        return total_rate

    def draw_waiting_time(self, instant_rate):
        """
        Draw a time at which mutation occurs according to mutation rates.
        :return:
        """
        time = np.random.exponential(scale=instant_rate)
        return time

    def mutate_on_branch(self):
        """
        Simulate molecular evolution in sequence given a branch length
        """

        times_sum = 0
        instant_rate = self.sum_rates()

        while True:
            random_time = self.draw_waiting_time(instant_rate)
            times_sum += random_time
            if times_sum > self.branch_length:
                break

            # Draw a mutation
            mutation = self.get_substitution()
            my_nt = mutation[0]
            to_state = mutation[1]
            # Subtract the mutation rate of the nucleotide to mutate
            instant_rate = instant_rate - my_nt.mutation_rate

            # Remove the mutated nucleotide from the event tree and update information in Nucleotide
            self.remove_nt(my_nt)
            self.update_nucleotide(my_nt, to_state)
            # Add the mutation rate of the mutated nucletoide
            instant_rate = instant_rate + my_nt.mutation_rate

            # Update information about adjacent nucleotides
            if my_nt.codons:  # If the mutated nucleotide belongs to at least one codon
                adjacent_positions = [my_nt.pos_in_seq - 2, my_nt.pos_in_seq - 1, my_nt.pos_in_seq + 1, my_nt.pos_in_seq + 2]
                for i in adjacent_positions:
                    if 0 < i < len(self.sequence.original_seq):  # If position in sequence
                        adj_nt = self.sequence.nt_sequence.nucleotide_at_pos(i)
                        for codon in adj_nt.codons:
                            if codon in my_nt.codons:  # If adjacent nucleotide and mutated nucleotide share at least one codon
                                self.remove_nt(adj_nt)
                                instant_rate = instant_rate - adj_nt.mutation_rate
                                self.update_nucleotide(adj_nt, adj_nt.state)
                                # Update instant rate for adjacent nucleotide
                                instant_rate = instant_rate + adj_nt.mutation_rate
                                break

            # Update event_tree to include a list of nucleotides on the tips
            self.event_tree = self.sequence.get_nts_on_tips()

        return self.sequence

    def remove_nt(self, nt):
        """
        Find nucleotide selected to mutate in the event tree and remove it from every branch on the event tree
        """
        for key_to_nt, value_to_nt in self.event_tree['to_nt'].items():
            if key_to_nt != nt.state:
                # Find branches that contain my nucleotide
                my_branch = self.event_tree['to_nt'][key_to_nt]['from_nt'][nt.state]
                 # Remove nt from non-synonymous mutations
                for omega_key, nucleotide_list in my_branch['is_nonsyn'].items():
                    if nt in nucleotide_list:
                        nucleotide_list.remove(nt)
                # Remove nt from synonymous mutations and list of nucleotides in substitution
                if nt in my_branch['is_syn']: my_branch['is_syn'].remove(nt)
                if nt in my_branch['nts_in_subs']: my_branch['nts_in_subs'].remove(nt)

    def update_nucleotide(self, nt, to_state):
        """
        Update parameters on the mutated nucleotide
        """
        nt.set_state(to_state)  # Update the state of the nucleotide
        # Update rates, omega key and event tree with the nucleotide according to its new state
        nt.set_complement_state()  # Change complementary state given the mutation
        rates = self.sequence.get_substitution_rates(nt)  # Calculate new rates and update event tree
        nt.set_rates(rates[0])  # Update substitution rates
        nt.set_my_omegas(rates[1]) # Update omega keys
        nt.set_nt_rate()

class SimulateOnTree:
    """
    Simulate evolution within a sequence throughout an entire phylogeny
    """

    def __init__(self, root_sequence , phylo_tree):
        self.root_sequence  = root_sequence # Object of class Sequence (Double Linked list of Nucleotides)
        self.phylo_tree = phylo_tree # Phylogenetic tree over which sequence will evolve


    def get_parent_clade(self, child_clade):
        """
        :param child_clade: current clade
        :return: path to the parent clade
        """
        node_path = self.phylo_tree.get_path(child_clade)
        return node_path[-2] if len(node_path) > 1 else self.phylo_tree.root

    def traverse_tree(self):
        """
        Mutate a sequence along a phylogeny by traversing it in level-order
        @return Phylo tree with Clade objects annotated with sequences.
        """

        # Assign root_seq to root Clade
        self.phylo_tree.root.sequence = self.root_sequence
        root = self.phylo_tree.root

        for clade in self.phylo_tree.find_clades(order='level'):
            # skip the root
            if clade is root:
                continue

            parent = self.get_parent_clade(clade)
            #Create a deep copy of the parent sequence
            parent_sequence = copy.deepcopy(parent.sequence)
            print(parent_sequence)
            # Mutate sequence and store it on clade
            simulation = SimulateOnBranch(parent_sequence, clade.branch_length)
            clade.sequence = simulation.mutate_on_branch()

        return self.phylo_tree


    def get_alignment(self):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        """
        #aln = open('/home/lmunoz/Projects/ovrf/HBV/NEWtest_Output.txt', "w+")

        final_tree = self.traverse_tree()
        print(final_tree)
        for clade in final_tree.get_terminals():
            seq = clade.sequence.get_string_sequence()
            #aln.write(">Sequence_{} \n{}\n".format(clade, seq))
            print(seq)

        #aln.close()
