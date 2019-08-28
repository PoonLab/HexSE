import argparse

from Bio import Phylo

from src.evol_rates import Rates
from src.sequence_info import Sequence
from src.simulate import Simulate


def get_args(parser):
    parser.add_argument(
        'seq', type=argparse.FileType('r'),
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'tree',
        help='Path to file containing phylogenetic tree in Newick format.'
    )
    parser.add_argument(
        'outfile', type=argparse.FileType('w'),
        help='Path to the alignment file.'
    )
    parser.add_argument(
        '--orfs', type=argparse.FileType('r'), default=None,
        help='Path to a csv file containing the start and end coordinates of the open reading frames'
    )
    parser.add_argument(
        '--mu', type=float,
        help='Global substitution rate per site per unit time'
    )
    parser.add_argument(
        '--bias', nargs="6", type=float, default=[1, 1, 1, 1, 1, 1],
        help='List of transversion/ transition rates assuming time reversibility. '
             'Format: [AC, AG, AT, CG, CT, GT]'
    )
    parser.add_argument(
        '--pi', nargs="4", type=float, default=[None, None, None, None],
        help='Vector of stationary nucleotide frequencies. '
             'Format: [A, T, G, C]'
    )
    parser.add_argument(
        '--omega', nargs="6", type=float, default=[None, None, None, None, None, None],
        help='List of dN/dS ratios for each reading frame along the length of the sequence. '
             'Format: [+0, +1, +2, -0, -1, -2]'
    )

    return parser.parse_args()


def main():
    parser = argparse.ArgumentParser(
        description='Simulates and visualizes the evolution of a sequence through a phylogeny'
    )
    args = get_args(parser)

    # Read in the sequence
    s = ''
    for line in args.seq:
        # Skip header if the file is a FASTA file
        if not (line.startswith(">") or line().startswith("#")):
            s += line.strip('\n\r').upper()

    # Read in the tree
    tree = Phylo.read(args.tree, 'newick', rooted=True)

    # Read in ORFs as a list of tuples
    unsorted_orfs = []
    for line in args.orfs:
        line = line.split(',')
        orf = (line[0], line[1])
        unsorted_orfs.append(orf)

    # Make Sequence object
    seq = Sequence(s, unsorted_orfs)

    # Reformat pi into a dictionary
    # TODO: check that there is a check for all non-None values in pi
    keys = ['A', 'T', 'G', 'C']
    pi = dict(zip(keys, args.pi))

    # Reformat omega into a dictionary
    # TODO: check that there is a check for all non-None values in omega
    keys = ['+0', '+1', '+2', '-0', '-1', '-2']
    omega = dict(zip(keys, args.omega))

    # Get Rates
    rates = []
    for nt in seq.sequence:
        rate = Rates(args.mu, pi, args.bias, omega)
        rates.append(rate)

    # Run simulation
    sim = Simulate(rates, tree)

    sim.traverse_tree()
    sim.get_alignment(args.outfile)


if __name__ == '__main__':
    main()
