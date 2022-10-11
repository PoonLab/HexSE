# Get a specific reading frame from an alignment (simulation output)
import argparse
import os
from urllib.parse import unquote
import ast
import sys

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
    
    parser.add_argument(
        '--out', default=None,
        help = 'path to log file; defaults to stdout'
    )
    
    return parser.parse_args()

def get_cds(alignment_file, start, end, out_name):
    
    out_file = open(out_name, "w")
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
    out_file = args.out
    
    if not out_file:
        file_name = os.path.basename(alignment_file)
        out_file = f"{file_name}_{start}_{end}.fa"

    print(start,end)
    get_cds(alignment_file, start, end, out_file)

if __name__ == '__main__':
    main()
