# Simulate evolution in a sequence with overlapping reading frames

import copy
import random
import sys

import numpy as np

NUCLEOTIDES = ['A', 'C', 'G', 'T']

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
        :param sequence: object of class Sequence (Double Linked list of Nucleotides)
                         in the parental node
        :param branch_length: length of the branch over which evolution is happening
        """
        self.sequence = sequence  # Sequence object
        self.branch_length = branch_length

    @staticmethod
    def test_omega_tree(omega_tree):
        for keys in omega_tree.keys():
            if None in keys:
                print("KEYS FROM OMEGA TREE\n",keys)

    def get_substitution(self):
        """
        Select a substitution by moving over the event_tree according to the generation of random numbers
        """

        event_tree = self.sequence.event_tree

        # Select target nucleotide based on stationary base frequency (pi) and number of events for each nucleotide
        to_dict = {
                    to_nt: self.sequence.pi[to_nt]*event_tree['to_nt'][to_nt]['nt_events']
                    for to_nt in self.sequence.pi.keys()
                    }
        to_mutation = self.weighted_random_choice(to_dict, sum(to_dict.values()))

        # Select state of initial nucleotide with probability based on transition/transversion rate ratio
        from_tree = event_tree['to_nt'][to_mutation]['from_nt']
        from_dict = {}
        for from_nt in NUCLEOTIDES:
            if to_mutation != from_nt:
                if self.sequence.is_transv(from_nt, to_mutation):  # If its transversion apply kappa
                    from_dict[from_nt] = (self.sequence.kappa / (1 + 2 * self.sequence.kappa)) * from_tree[from_nt]['nt_events']

                else:
                    from_dict[from_nt] = (1 / (1 + 2 * self.sequence.kappa)) * from_tree[from_nt]['nt_events']

        from_mutation = self.weighted_random_choice(from_dict, sum(from_dict.values()))

        # Select category based on mu values
        cat_tree = from_tree[from_mutation]
        cat_sum = sum(self.sequence.cat_values.values())
        cat_dict = {
            cat: (self.sequence.cat_values[cat]/cat_sum)*cat_tree[cat]['nt_events']
            for cat in self.sequence.cat_values.keys()
            }

        selected_cat = self.weighted_random_choice(cat_dict, sum(cat_dict.values()))

        # Select sequence region to mutate based on the number of nucleotides
        orf_tree = cat_tree[selected_cat]
        all_maps = self.sequence.all_maps
        orf_dict = {}

        for orf_combo in orf_tree.keys():
            if type(orf_combo) == tuple:
                orf_dict[orf_combo] = (all_maps[orf_combo]['len']/self.sequence.length) * orf_tree[orf_combo]['nt_events']

        selected_orf_combo = self.weighted_random_choice(orf_dict, sum(orf_dict.values()))

        # Select omega combo when region with orfs. Else, select a random nucleotide
        omega_tree = orf_tree[selected_orf_combo]
        selected_nt = None
        selected_omega = None
        
        # If region has no ORFs, select any nucleotide on the list at random
        if all(v == 0 for v in selected_orf_combo):  
            tuple_key = tuple([None]*len(selected_orf_combo))
            selected_omega = tuple_key
            selected_nt = random.choice(omega_tree[tuple_key])
        
        else:  # If region has ORFs, select omega
            omega_dict = {}
            for omega_combo in omega_tree.keys():
                if type(omega_combo) == tuple:
                    omega_dict[omega_combo] = self.sequence.total_omegas[omega_combo]['value'] * len(omega_tree[omega_combo])
            
            selected_omega = self.weighted_random_choice(omega_dict, sum(omega_dict.values()))
            # Choose nucleotide
            selected_nt = random.choice(omega_tree[selected_omega])

        # print(">>",to_mutation, "\n", from_mutation, "\n", selected_cat, "\n", selected_orf_combo, "\n", selected_omega, "\n", selected_nt, "\n")
        # Return coordenates of the selection
        # sys.exit()
        return to_mutation, from_mutation, selected_cat, selected_orf_combo, selected_omega, selected_nt




        # # Select: mu category
        # def possible_cats():  # categories
        #     cat_tree = copy.copy(self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'])
        #     for _ in list(cat_tree.keys()):
        #         selected_cat = self.select_key(cat_tree)
        #         cat_tree.pop(selected_cat)  # Remove empty key
        #         yield selected_cat

        # # Select omega branch
        # def possible_nts(selected_cat):
        #     omega_dict = copy.copy(self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'][selected_cat]['omega'])
            
        #     # Create a dictionary containing only omega tuples and probability value
        #     omega_weights = {omega: omega_dict[omega]['prob']*omega_dict[omega]['number_of_events'] for omega in omega_dict.keys()}
        #     nt_dict = self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'][selected_cat]

        #     # Check if all tips of the tree for this class are empty
        #     if all(not value for value in nt_dict.values()):
        #         pass

        #     else:
        #         for _ in omega_dict.keys():
        #             selected_omega = self.weighted_random_choice(omega_weights, sum(omega_weights.values()))
        #             if selected_omega == None:
        #                 # element = list(omega_weights.keys())[0][0]
        #                 print("\n>>>>WHERE AM I\n",omega_weights)
        #             self.test_omega_tree(omega_dict) 
        #             # Select nucleotide
        #             try:
        #                 nt_list = self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]['category'][selected_cat]['omega'][selected_omega]
        #             except KeyError:
        #                 print(to_mutation, from_mutation, selected_cat, selected_omega)
        #                 print(self.sequence.event_tree)
        #                 raise
        #             omega_weights.pop(selected_omega)  # Remove empty key
        #             yield nt_list

        # for selected_cat in possible_cats():
        #     for nt_list in possible_nts(selected_cat):  # Check if nt list is empty
        #         if nt_list:
        #             from_nucleotide = random.choice(nt_list['nt'])
        #             return from_nucleotide, to_mutation

        # print("\n>>>> OH NO!!")
        #print(self.sequence.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation])

    # def select_key(self, dictionary):
    #     """
    #     Select a random key using weighted_random_choice
    #     """
    #     temp_dict = {}  # Dictionary containing keys to 'from nucleotides' or 'categories' and their probability (weight)
    #     for key, value in dictionary.items():
    #         if value:
    #             temp_dict[key] = value['prob'] * value['number_of_events']
    #     selected_key = self.weighted_random_choice(temp_dict, sum(temp_dict.values()))
    #     return selected_key

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
            to_mutation, from_mutation, selected_cat, selected_orf_combo, selected_omega, selected_nt = self.get_substitution()
            
            # Update instant rate by adding the rate at which this mutation occurs 
            instant_rate = instant_rate + selected_nt.rates[to_mutation]
            
            # Remove the mutated nucleotide from the event tree
            self.remove_nt(selected_nt, selected_orf_combo)

            # Calculate mutation rates for nucleotide under new states and place on the proper branches of the tree
            # Re-create substitution rates for the nucleotide
            self.update_nucleotide(selected_nt, to_mutation)
            # sys.exit()

            # Update information about adjacent nucleotides
            if selected_nt.codons:  # If the mutated nucleotide belongs to at least one codon
                adjacent_positions = [selected_nt.pos_in_seq - 2, selected_nt.pos_in_seq - 1,
                                      selected_nt.pos_in_seq + 1, selected_nt.pos_in_seq + 2]

                for i in adjacent_positions:
                    if 0 < i < len(self.sequence.nt_sequence):  # If position in sequence
                        adj_nt = self.sequence.nt_sequence[i]

                        for codon in adj_nt.codons:
                            # If adjacent nucleotide and mutated nucleotide share at least one codon
                            if codon in selected_nt.codons:
                                instant_rate = instant_rate - adj_nt.mutation_rate
                                self.remove_nt(adj_nt)
                                self.update_nucleotide(adj_nt, adj_nt.state)
                                # Update instant rate for adjacent nucleotide
                                instant_rate = instant_rate + adj_nt.mutation_rate
                                break

            # Update number of events in the Tree
            self.sequence.count_events_per_layer()

        return self.sequence

    def remove_nt(self,selected_nt, selected_orf_combo):
        """
        Find nucleotide selected to mutate in the event tree and remove it from every branch on the event tree
        Note: it should be in three branches (one nucleotide can change to three possible)   
        :param selected_nt: Nucleotide Object. 
        :param selected_ord_combo: Tuple.  1 for the orfs the nucleotide is part of. 0 for the ones it is not. 
        """

        for to_nt in NUCLEOTIDES:
            if to_nt != selected_nt.state:
                
                cat = selected_nt.cat_keys[to_nt]
                omega_combo = selected_nt.omega_keys[to_nt]
                my_branch = self.sequence.event_tree['to_nt'][to_nt]['from_nt'][selected_nt.state][cat][selected_orf_combo][omega_combo]
                my_branch.remove(selected_nt)
                

    def update_nucleotide(self, nt, to_state):
        """
        Update parameters on the mutated nucleotide
        """
        # Update the state of the nucleotide
        nt.set_state(to_state)
        # Change complementary state given the mutation
        nt.set_complement_state()
        # RESET All values
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

        if new_omega:  # If new omega key is created in the event tree, update the tree
            self.sequence.compute_probability()


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
