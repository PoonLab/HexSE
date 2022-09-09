# Get a specific reading frame from an alignment (simulation output)
import argparse
import os
import json
from urllib.parse import unquote
import ast
import sys


handle = "/home/laura/Projects/ovrf/temp/hiv_0.05.fa"
pos = (200,500)

def get_args(parser):
    parser.add_argument(
        'file',
        help = 'Path to the file containing the alignment'
    )
    parser.add_argument(
        'start', default=None, type=int,
        help = 'Position of the first nucleotide on the CDS'
    )

    parser.add_argument(
        'end', default=None, type=int,
        help = 'Position of the last nucleotide on the CDS'
    )

    return parser.parse_args()

def get_cds(alignment_file, start, end):

    file_name = os.path.basename(alignment_file)
    genome = file_name.split('_')[0]


    out_file = open("{}_{}_{}.fa".format(file_name, start, end),"w")

    with open(alignment_file) as alignment:
        for line in alignment:
            if line.startswith('>'):
                out_file.write(line)
            else:
                #cds = line[start : end-3]  # (-3) Remove the stop codon at the end of the reading frame
                cds = line[start : end]
                out_file.write(cds)
                out_file.write('\n\n')

def main():
    parser = argparse.ArgumentParser(
        description='Get a CDS region from a nucleotide alignment'
    )

    args = get_args(parser)
    alignment_file = args.file
    start, end = args.start, args.end
    
    print(start,end)
    get_cds(alignment_file, start, end)

if __name__ == '__main__':
    main()
