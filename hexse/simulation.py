# Simulate evolution in a sequence with overlapping reading frames

import copy
import random
import sys

import numpy as np
import pprint

from tqdm import tqdm

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
            to_nt: self.sequence.pi[to_nt] * event_tree['to_nt'][to_nt]['nt_events']
            for to_nt in self.sequence.pi.keys()  # A,C,G,T
        }
        to_mutation = self.weighted_random_choice(to_dict, sum(to_dict.values()))

        # Select state of initial nucleotide with probability based on transition/transversion rate ratio
        from_tree = event_tree['to_nt'][to_mutation]['from_nt']
        from_dict = {}
        for from_nt in NUCLEOTIDES:
            if to_mutation != from_nt:
                if self.sequence.is_transv(from_nt, to_mutation):  # If its transversion applies kappa
                    from_dict[from_nt] = (self.sequence.kappa / (1 + 2 * self.sequence.kappa)) * \
                                         from_tree[from_nt]['nt_events']
                else:
                    from_dict[from_nt] = (1 / (1 + 2 * self.sequence.kappa)) * from_tree[from_nt]['nt_events']

        from_mutation = self.weighted_random_choice(from_dict, sum(from_dict.values()))

        # Select category based on mu values
        cat_tree = from_tree[from_mutation]
        cat_sum = sum(self.sequence.cat_values.values())
        cat_dict = {
            cat: (self.sequence.cat_values[cat] / cat_sum) * cat_tree[cat]['nt_events']
            for cat in self.sequence.cat_values.keys()
        }

        selected_cat = self.weighted_random_choice(cat_dict, sum(cat_dict.values()))

        # Select sequence region to mutate based on the number of nucleotides (E.g, ((0, 1, 1, 0)))
        orf_tree = cat_tree[selected_cat]
        region_weights = {}

        for orf_combo in orf_tree.keys():
            if type(orf_combo) == tuple:
                region_weights[orf_combo] = (
                    orf_tree[orf_combo]['region_weight']  # weight calculated in count_events_per_layer()
                )
        selected_orf_combo = self.weighted_random_choice(region_weights, sum(region_weights.values()))

        # Select omega combo when region with orfs. Else, select a random nucleotide
        omega_tree = orf_tree[selected_orf_combo]
        selected_nt = None
        selected_omega = None
        
        # If region has no ORFs, select any nucleotide on the list at random
        if all(v == 0 for v in selected_orf_combo):  
            tuple_key = tuple([None]*len(selected_orf_combo))
            selected_omega = tuple_key
            selected_nt = random.choice(omega_tree[tuple_key])
        
        else:  # If region has ORFs, select omega combo
            omega_dict = {}
            for omega_combo in omega_tree.keys():
                if type(omega_combo) == tuple:  # Ignore the key 'nt_events' that contains the number on events on the branch
                    omega_dict[omega_combo] = self.sequence.total_omegas[omega_combo]['value'] * len(omega_tree[omega_combo])

            selected_omega = self.weighted_random_choice(omega_dict, sum(omega_dict.values()))
            # Choose nucleotide
            selected_nt = random.choice(omega_tree[selected_omega])

        # print(">>",to_mutation, "\n", from_mutation, "\n", selected_cat, "\n", selected_orf_combo, "\n", selected_omega, "\n", selected_nt, "\n")
        # Return coordenates of the selection
        return to_mutation, from_mutation, selected_cat, selected_orf_combo, selected_omega, selected_nt


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

    def mutate_on_branch(self, instant_rate):
        """
        Simulate molecular evolution in sequence given a branch length
        :return: the mutated sequence
        """
        times_sum = 0

        while True:
            # Draw a time at which mutation occurs according to mutation rates
            random_time = np.random.exponential(scale=(1 / instant_rate))
            times_sum += random_time
            
            if times_sum > self.branch_length:
                break

            # Draw a mutation
            new_state, from_mutation, selected_cat, selected_orf_combo, selected_omega, selected_nt = self.get_substitution()
            
            # Update instant rate by adding the rate at which this mutation occurs 
            instant_rate = instant_rate + selected_nt.rates[new_state]
            
            # Remove the mutated nucleotide from the Event Tree
            self.remove_nt(selected_nt, selected_orf_combo)  # ERROR OCCURS FROM HERE
            
            # For nt under new state, update it's attributes: state, complement, substitution_rate, omega_keys, cat_keys and mutation_rate
            self.update_nucleotide_info(selected_nt, new_state)
            updated_nts = [selected_nt] # Keep track of nucleotides that have been removed and updated on the Event Tree when a mutation occurs to not update them several times
            
            # Re-locate nucleotide on Event Tree according to it's new state
            self.sequence.nt_in_event_tree(selected_nt)

            # If nucleotide is part of a codon, the rest of nucleotides in the codon might be subject to different syn and non-syn mutations
            # Therefore, we will need to re calculate mutation rates for adjacent nucleotides and properly place them on the Event Tree
            
            if selected_nt.codons:
                for codon in selected_nt.codons:
                    # codon.omega = codon.select_omega()  # does mutations change the strength of selection against non-syn subs?
                    for adj_nt in codon.nts_in_codon:
                            # codons probably share nucleotides in common. Avoid to update them twice.
                            if adj_nt not in updated_nts:
                                # Remove, re calculate, and re populate Event Tree with adjacent nucleotides
                                self.remove_nt(adj_nt, adj_nt.orf_map_key)  # Adjacent nucleotides might have different orf maps since they can be not on the same overlapping region
                                self.update_nucleotide_info(adj_nt, adj_nt.state)
                                self.sequence.nt_in_event_tree(adj_nt)
                                updated_nts.append(adj_nt)
            
        
            # Update number of events in the Tree per branch, updates self.sequence.total_omegas and number of events per branch
            self.sequence.count_events_per_layer()

        return self.sequence

    def remove_nt(self,selected_nt, selected_orf_combo):
        """
        Find nucleotide on the Event Tree and remove it from the tips where it's stored
        Note: it should be on three branches (each nucleotide has three posible states that it can change to)   
        :param selected_nt: Nucleotide Object. 
        :param selected_orf_combo: Tuple.  1 for the orfs the nucleotide is part of. 0 for the ones it is not. (E.g (0,1,0))
        """

        for to_nt in NUCLEOTIDES:
            
            if to_nt != selected_nt.state:
                
                cat = selected_nt.cat_keys[to_nt]
                omega_combo = selected_nt.omega_keys[to_nt]

                if cat and omega_combo:  # If mutation did not introduced STOP codons on seq, therefore it IS stored on the tree:
                    # Note, when mutation causes STOP codons, both cat and omega_combo are Noneq
                    my_branch = self.sequence.event_tree['to_nt'][to_nt]['from_nt'][selected_nt.state][cat][selected_orf_combo][omega_combo]
                    my_branch.remove(selected_nt)
                

    def update_nucleotide_info(self, nt, new_state):
        """
        When a mutation occurs, the following attributes change for the nt object: state, complement_state, substitution_rate,
        omega_keys, cat_keys and mutation_rate
        """
        # Update the state of the nucleotide
        nt.set_state(new_state)
        # Change complementary state given the mutation
        nt.set_complement_state()
        # Sub_rates, omega_key, category_key, and mutation_rate are set on set_substition_rate
        self.sequence.set_substitution_rates(nt) 
        

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

    def traverse_tree(self, th):
        """
        Mutate a sequence along a phylogeny by traversing it in level-order
        :return phylo_tree: A Phylo tree with Clade objects annotated with sequences.
        """
        # Assign root_seq to root Clade
        self.phylo_tree.root.sequence = self.root_sequence
        root = self.phylo_tree.root
        clades = self.phylo_tree.find_clades(order='level')
        number_of_clades = sum(1 for _ in clades)
        
        with tqdm(total=number_of_clades - 1, unit='clades') as pbar:
            pbar.set_description("Traversing tree")
            
            for clade in self.phylo_tree.find_clades(order='level'):
                # skip the root
                if clade is root:
                    continue

                parent = self.get_parent_clade(clade)

                # Create a deep copy of the parent sequence
                parent_sequence = copy.deepcopy(parent.sequence)
                instant_rate = parent_sequence.get_instant_rate()
                
                if th and (clade.branch_length * instant_rate) > th:
                    raise TooManyEventsError(instant_rate, clade)

                # Mutate sequence and store it on clade
                simulation = SimulateOnBranch(parent_sequence, clade.branch_length)
                clade.sequence = simulation.mutate_on_branch(instant_rate)
            
                pbar.update(1)

        return self.phylo_tree

    def get_alignment(self, outfile=None, th=None):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        """
        final_tree = self.traverse_tree(th)

        if outfile is not None:
            with open(outfile, 'w+') as out_handle:
                for clade in final_tree.get_terminals():
                    out_handle.write(">{} \n{}\n".format(clade, clade.sequence))
        else:
            for clade in final_tree.get_terminals():
                pass
                print(">{} \n{}".format(clade, clade.sequence))

class TooManyEventsError(Exception):
    """Exception raised when there are too many events in a branch
    """

    def __init__(self, instant_rate, clade,
                message="Number of events is too high for branch "):
        
        self.instant_rate = instant_rate
        self.clade = clade
        branch_length = clade.branch_length
        n_events = clade.branch_length * instant_rate
        self.message = '\n'.join([
            message + str(self.clade),
            f"Instant rate: {self.instant_rate}",
            f"Branch length: {branch_length}",
            f"Number of events: {n_events}"
        ])

        super().__init__(self.message)