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
        self.branch_length = branch_length

    def get_substitution(self):
        """
        Select a substitution by moving over the event_tree according to the generation of random numbers
        """
        # Select: to nucleotide

        to_dict = {
                    to_nt: self.sequence.pi[to_nt]*self.sequence.probability_tree['to_nt'][to_nt]['number_of_events']
                    for to_nt in self.sequence.pi.keys()
                    }
        to_mutation = self.weighted_random_choice(to_dict, sum(to_dict.values()))

        # Select: from nt
        from_tree = self.sequence.probability_tree['to_nt'][to_mutation]['from_nt']
        from_mutation = self.select_key(from_tree)

        # Select: mu category
        def possible_cats():
            cat_tree = copy.copy(self.sequence.probability_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['cat'])
            for i in list(cat_tree.keys()):
                selected_cat = self.select_key(cat_tree)
                cat_tree.pop(selected_cat)  # Remove empty key

                yield selected_cat

        # Select omega branch
        def possible_nts(selected_cat):
            omega_dict = copy.copy(self.sequence.probability_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['cat'][selected_cat]['omega'])
            # Create a dictionary contaning only omega tuples and probability value
            omega_weights = {omega: omega_dict[omega]['prob']*omega_dict[omega]['number_of_events'] for omega in omega_dict.keys()}
            nt_dict = self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'][selected_cat]

            # Check if all tips of the tree for this class are empty
            if all(not value for value in nt_dict.values()):
                pass

            else:
                for i in range(1, len(omega_dict.keys())):
                    selected_omega = self.weighted_random_choice(omega_weights, sum(omega_weights.values()))
                    # Select nucleotide
                    nt_list = self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'][selected_cat][selected_omega]
                    omega_weights.pop(selected_omega)  # Remove empty key

                    yield nt_list

        for selected_cat in possible_cats():
            for nt_list in possible_nts(selected_cat): # Check if nt list is empty
                if nt_list:
                    from_nucleotide = random.choice(nt_list)
                    return from_nucleotide, to_mutation

        print("\n>>>> OH NO!!")
        print(self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation])

    def select_key(self, dictionary):
        """
        Select a random key using weighted_random_choice
        """
        temp_dict = {}  # Dictionary containing keys to 'from nucleotides' or 'categories' and their probability (weight)
        for key, value in dictionary.items():
            if value:
                temp_dict[key] = value['prob'] * value['number_of_events']
        selected_key = self.weighted_random_choice(temp_dict, sum(temp_dict.values()))
        return (selected_key)

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
                                # Update instant rate for adjacent nucleotide
                                instant_rate = instant_rate + adj_nt.mutation_rate
                                break

            # Update number of events in probability Tree
            self.sequence.populate_prob_tree_with_events()

        return self.sequence

    def remove_nt(self, nt):
        """
        Find nucleotide selected to mutate in the event tree and remove it from every branch on the event tree
        :param nt: The nucleotide to be removed from the event tree
        """

        for to_nt, value_to_nt in self.sequence.event_tree['to_nt'].items():

            if to_nt != nt.state:

                # Find branches that contain my nucleotide
                if to_nt in nt.omega_in_event_tree:  # Nucleotide may not have this key because it was not stored on the event tree (mutation would introduce a STOP or was part of a START)
                    cat_key = nt.cat_keys[to_nt]
                    omega_key = nt.omega_in_event_tree[to_nt]
                    my_branch = self.sequence.event_tree['to_nt'][to_nt]['from_nt'][nt.state]['category'][cat_key][omega_key]
                    my_branch.remove(nt)

    def update_nucleotide(self, nt, to_state):
        """
        Update parameters on the mutated nucleotide
        """

        # Update the state of the nucleotide
        nt.set_state(to_state)
        # Change complementary state given the mutation
        nt.set_complement_state()
        # RESET All vallues
        # Rates
        nt.set_rates({})
        # Omega keys
        nt.set_omega({})
        # Class Keys
        nt.set_categories({})
        # Omega in event tree
        nt.set_omega_in_event_tree({})
        # Mutation rate
        nt.set_mutation_rate(0)

        # Update rates, omega key according to its new state
        self.sequence.set_substitution_rates(nt)
        new_omega = self.sequence.nt_in_event_tree(nt)  # Update nucleotide on the Event Tree and return new key if created

        if new_omega: # If new omeka key is created in the event tree, update the probability tree
            self.sequence.probability_tree = self.sequence.create_probability_tree()

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
