# Get a specific reading frame from an alignment (simulation output)
import argparse
import os
import json
from urllib.parse import unquote
import ast


handle = "/home/laura/Projects/ovrf/temp/hiv_0.05.fa"
pos = (200,500)

def get_args(parser):
    parser.add_argument(
        'file',
        help='Path to the file containing the alignment'
    )
    parser.add_argument(
        '-orfs', default=None,
        help='Dictionary with the valid reading frames in the sequence'
    )
    parser.add_argument(
        '--path_to_out', default=None, help='Path to the alignment file.'
    )

    return parser.parse_args()

def get_cds(alignment_file, start, end):

    file_name = os.path.basename(alignment_file)
    genome = file_name.split('_')[0]


    out_file = open("{}_{}_{}.fa".format(genome, start, end),"w")

    with open(alignment_file) as alignment:
        for line in alignment:
            if line.startswith('>'):
                out_file.write(line)
            else:
                cds = line[start : end-3]  # (-3) Remove the stop codon at the end of the reading frame
                out_file.write(cds)
                out_file.write('\n\n')

def main():
    parser = argparse.ArgumentParser(
        description='Get an ORF from an alignment'
    )

    args = get_args(parser)
    alignment_file = args.file
    orfs = ast.literal_eval(args.orfs)
    print(type(orfs))

    print("OPEN READING FRAMES", orfs)
    #path_to_out = args.path_to_out

    for key, value in orfs.items():
        #print(key, value)
        if value:
            #print(value)
            for orf in value:
                 start = orf[0]
                 end = orf[1]
                 get_cds(alignment_file, start, end)

if __name__ == '__main__':
    main()
