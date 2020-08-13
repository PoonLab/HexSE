# Simulate evolution in a sequence with overlapping reading frames

import copy
import random

import numpy as np


def is_stop(to_nt, from_nt):
    """
    Check if a STOP codon is being introduced in the sequence
    """
    for codon in from_nt.codons:
        codon, mutated_codon = codon.mutate_codon(codon.nt_in_pos(from_nt), to_nt)

        if str(mutated_codon) == 'TAA' or str(mutated_codon) == 'TGA' or str(mutated_codon) == 'TAG':
            return True

    return False


class SimulateOnBranch:
    """
    Simulate evolution within a sequence throughout a branch in a phylogeny
    """

    def __init__(self, sequence, branch_length):
        """
        :param sequence: object of class Sequence (Double Linked list of Nucleotides)  in the parental node
        :param branch_length: length of the branch over which evolution is happening
        """
        self.sequence = sequence  # Sequence object
        self.event_tree = sequence.get_event_tree()  # Tree of all possible mutations related to sequence
        self.probability_tree = sequence.get_probability_tree()
        self.branch_length = branch_length

    def get_substitution(self):
        """
        Select a substitution by moving over the event_tree according to the generation of random numbers
        """
        # Select: to nucleotide
        to_mutation = self.weighted_random_choice(self.sequence.pi, sum(self.sequence.pi.values()))

        # Select: from nt
        from_tree = self.probability_tree['to_nt'][to_mutation]['from_nt']
        trans_transv_dict = {}  # Dictionary contaning keys to from nucleotide weighted according to transition-transversion probability
        for from_nt, dict in from_tree.items():  # Create dictionary with probabilities for this branch
            if dict:  # Not none
                prob = from_tree[from_nt]['tr_p']
                trans_transv_dict[from_nt] = prob
        from_mutation = self.weighted_random_choice(trans_transv_dict, sum(trans_transv_dict.values()))

        # Select: mu class
        class_tree = self.probability_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['class_p']
        mu_dict = {}  # Dictionary contaning keys to mu classes weighted according to mu probability
        for mu_class, dict in class_tree.items():
            mu_dict[mu_class] = dict['mu_p']
        selected_class = self.weighted_random_choice(mu_dict, sum(mu_dict.values()))

        # Select omega branch
        omega_dict = self.probability_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['class_p'][selected_class]['omega_p']
        selected_omega = self.weighted_random_choice(omega_dict, sum(omega_dict.values()))

        # Select nucleotide
        nt_list = self.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['class'][selected_class][selected_omega]
        from_nucleotide = random.choice(nt_list)

        print(from_nucleotide, from_nucleotide.pos_in_seq)

        return from_nucleotide, to_mutation


    @staticmethod
    def weighted_random_choice(dictionary, sum_values):
        """
        Randomly select a key on a dictionary where values correspond to the weight for the key
        :param dictionary: Dictionary to select the value from
        :param sum_values: sum all values on dict to establish the limit for the mutation
        :return: random key from dict
        """
        key = None
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
        :return: the sum of the mutation rates
        """
        total_rate = sum([nt.mutation_rate for nt in iter(self.sequence.get_sequence())])
        return total_rate

    def mutate_on_branch(self):
        """
        Simulate molecular evolution in sequence given a branch length
        :return: the mutated sequence
        """
        times_sum = 0
        instant_rate = self.sum_rates()

        while True:
            # Draw a time at which mutation occurs according to mutation rates
            random_time = np.random.exponential(scale=(1 / instant_rate))
            times_sum += random_time
            if times_sum > self.branch_length:
                break

            # Draw a mutation
            mutation = self.get_substitution()
            my_nt = mutation[0]
            to_state = mutation[1]
            instant_rate = instant_rate - my_nt.mutation_rate

            # Remove the mutated nucleotide from the event tree and update information in Nucleotide
            self.remove_nt(my_nt)

            # Re-create substitution rates for the nucleotide
            self.update_nucleotide(my_nt, to_state)

            # Re-store the nucleotide on the Event Tree and update the total number of events
            self.update_nt_on_tree(my_nt, to_state)

            # Add the mutation rate of the mutated nucleotide
            instant_rate = instant_rate + my_nt.mutation_rate

            # Update information about adjacent nucleotides
            if my_nt.codons:  # If the mutated nucleotide belongs to at least one codon
                adjacent_positions = [my_nt.pos_in_seq - 2, my_nt.pos_in_seq - 1,
                                      my_nt.pos_in_seq + 1, my_nt.pos_in_seq + 2]

                for i in adjacent_positions:
                    if 0 < i < len(self.sequence.nt_sequence):  # If position in sequence
                        adj_nt = self.sequence.nt_sequence[i]

                        for codon in adj_nt.codons:
                            # If adjacent nucleotide and mutated nucleotide share at least one codon
                            if codon in my_nt.codons:
                                instant_rate = instant_rate - adj_nt.mutation_rate
                                self.remove_nt(adj_nt)
                                self.update_nucleotide(adj_nt, adj_nt.state)
                                self.update_nt_on_tree(adj_nt, adj_nt.state)
                                # Update instant rate for adjacent nucleotide
                                instant_rate = instant_rate + adj_nt.mutation_rate
                                break

            # Update event_tree to include a list of nucleotides on the tips
            self.event_tree = self.sequence.get_nts_on_tips()

        return self.sequence

    def remove_nt(self, nt):
        """
        Find nucleotide selected to mutate in the event tree and remove it from every branch on the event tree
        :param nt: The nucleotide to be removed from the event tree
        """
        for key_to_nt, value_to_nt in self.event_tree['to_nt'].items():
            if key_to_nt != nt.state:
                # Find branches that contain my nucleotide
                my_branch = self.event_tree['to_nt'][key_to_nt]['from_nt'][nt.state]

                # Remove nt from synonymous mutations and
                if nt in my_branch['is_syn']:
                    my_branch['is_syn'].remove(nt)

                    # Update the number of events
                    my_branch['number_of_events'] -= 1
                    self.event_tree['to_nt'][key_to_nt]['events_for_nt'] -= 1
                    self.event_tree['total_events'] -= 1

                # Remove nt from dictionary of nucleotides in substitution
                if nt in my_branch['nts_in_subs']:
                    my_branch['nts_in_subs'].pop(nt, None)

                # Remove nt from non-synonymous mutations
                for dN_key in list(my_branch['is_nonsyn']['dN'].keys()):
                    dN_nt_list = my_branch['is_nonsyn']['dN'][dN_key]
                    if nt in dN_nt_list:      # Remove nucleotide from dN list
                        dN_nt_list.remove(nt)
                    if not dN_nt_list:        # Remove dN keys with no associated nucleotides
                        del my_branch['is_nonsyn']['dN'][dN_key]

                for dS_key in list(my_branch['is_nonsyn']['dS'].keys()):
                    dS_nt_list = my_branch['is_nonsyn']['dS'][dS_key]
                    if nt in dS_nt_list:      # Remove nucleotide from dS list
                        dS_nt_list.remove(nt)
                    if not dS_nt_list:        # Remove dS keys with no associated nucleotides
                        del my_branch['is_nonsyn']['dS'][dS_key]

    def update_nucleotide(self, nt, to_state):
        """
        Update parameters on the mutated nucleotide
        """

        # Update the state of the nucleotide
        nt.set_state(to_state)
        # Change complementary state given the mutation
        nt.set_complement_state()
        # Update rates, omega key according to its new state
        self.sequence.get_substitution_rates(nt)

    def update_nt_on_tree(self, nt, to_state):
        """
        Update nucleotide on the Event Tree and their count in synonymous substitutions
        """
        self.sequence.nt_in_event_tree(nt) # Update the nucleotide on the branches of the Event Tree where it belongs

        # Update all occurrences of the nucleotide in the event tree
        for key_to_nt, val in self.event_tree['to_nt'].items():

            if key_to_nt != to_state:

                branch = self.event_tree['to_nt'][key_to_nt]['from_nt'][nt.state]['nts_in_subs']
                # Find branches that contain the nucleotide
                if nt in branch:

                    # Check if the mutation is synonymous
                    for codon in nt.codons:
                        pos_in_codon = codon.nt_in_pos(nt)

                        # Update the number of events
                        if codon.is_nonsyn(pos_in_codon, nt.state):
                            branch['number_of_events'] += 1
                            self.event_tree['to_nt'][key_to_nt]['events_for_nt'] += 1
                            self.event_tree['total_events'] += 1


class SimulateOnTree:
    """
    Simulate evolution within a sequence throughout an entire phylogeny
    """
    def __init__(self, root_sequence, phylo_tree, outfile=None):
        self.root_sequence = root_sequence  # Sequence object
        self.phylo_tree = phylo_tree  # Phylogenetic tree over which sequence will evolve
        self.outfile = outfile

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
        :return phylo_tree: A Phylo tree with Clade objects annotated with sequences.
        """
        # Assign root_seq to root Clade
        self.phylo_tree.root.sequence = self.root_sequence
        root = self.phylo_tree.root

        for clade in self.phylo_tree.find_clades(order='level'):
            # skip the root
            if clade is root:
                continue

            parent = self.get_parent_clade(clade)

            # Create a deep copy of the parent sequence
            parent_sequence = copy.deepcopy(parent.sequence)

            # Mutate sequence and store it on clade
            print("Simulating on Branch", clade)
            simulation = SimulateOnBranch(parent_sequence, clade.branch_length)
            clade.sequence = simulation.mutate_on_branch()

        return self.phylo_tree

    def get_alignment(self, outfile=None):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        """
        final_tree = self.traverse_tree()

        if outfile is not None:
            with open(outfile, 'w+') as out_handle:
                for clade in final_tree.get_terminals():
                    out_handle.write(">{} \n{}\n".format(clade, clade.sequence))
        else:
            for clade in final_tree.get_terminals():
                pass
                print(">{} \n{}".format(clade, clade.sequence))
