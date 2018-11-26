from Bio import Phylo


class Sim:
    def __init__(self, rate, matrix, alphabet='ACGT'):
        self.rate = rate
        self.matrix = matrix
        self.alphabet = alphabet
        # should include some checks on matrix here (e.g., square)

    def simulate_on_branch(self, seq0, time):
        """
        Simulate molecular evolution on the branch given starting sequence seq0
        :param seq0:
        :param time:
        :return:
        """
        seq1 = ''
        # TODO: evolution happens here
        return seq1

    def traverse_tree(self, tree, root_seq):
        """
        Post-order traversal of tree from the root.
        Call simulate_on_branch on each branch and feed the resulting sequence to initialize
        the next call.
        @return Phylo tree with Clade objects annotated with sequences.
        """
        # TODO: check if tree is rooted - if not, throw exception

        # assign root_seq to root Clade
        tree.root.sequence = root_seq

        # annotate Clades with parents
        for clade in tree.find_clades(order='level'):
            for child in clade:
                child.parent = clade

        for node in tree.find_clades(order='level'):
            # TODO: skip the root (sequence already assigned)
            node.sequence = self.simulate_on_branch(node.parent.sequence, node.branch_length)

        return tree

    def get_alignment(self, tree):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        attributes as list.
        :param tree:
        :return:
        """


sim = Sim()  # create an instance of object class Sim
tree = Phylo.read('<some file>', 'newick')
tree = sim.traverse_tree(tree, 'ACGT')
result = sim.get_alignment(tree)

