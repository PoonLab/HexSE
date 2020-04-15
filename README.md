# Simulating evolution in Overlapping Reading Frames (OvRFs)

OvRFs is a bioinformatic pipeline to simulate evolution within Overlapping Reading Frames. The purpose of this pipeline is to account for effects of mutations in overlapping nucleotides when simulating evolution. 

1. [Usage](https://github.com/PoonLab/ovrf#usage)
2. [How it works]
3. [Outputs]

## Usage
To display the input options you can run:
```console
$ python3 run_simulation.py --help
usage: run_simulation.py [-h] [--outfile OUTFILE] [--orfs ORFS] [--mu MU]

[--kappa KAPPA] [--pi PI] [--omega OMEGA]

seq tree

  

Simulates and visualizes the evolution of a sequence through a phylogeny

  

positional arguments:

seq  Path to the file containing the query sequence.

tree Path to file containing phylogenetic tree in Newick

format.

  

optional arguments:

-h, --help show this help message and exit

--outfile OUTFILE  Path to the alignment file.

--orfs ORFS  Path to a csv file containing the start and end

coordinates of the open reading frames. Format:

start,endIf no ORFS are specified, the program will find

ORFs automatically

--mu MU  Global substitution rate per site per unit time

--kappa KAPPA  Transversion/ transition rate assuming time

reversibility.

--pi PI  Vector of stationary nucleotide frequencies. If no value

is specified, the program will use the empirical

frequencies in the sequence. Format: [A, T, G, C]

--omega OMEGA  List of dN/dS ratios along the length of the sequence.
 
```
### Input options
**Nucleotide sequence:** The sequence to evolve can be either in `fasta` of `genbank` format.

**Phylogenetic tree:** expect on `newick` format

**Open reading frames:** If the input sequence is a `fasta` file, you have two options:
1. Specify the ORFs as a list of tuples with the first and last nucleotide of the CDs `[(0, 8), (2,10)]` 
2. Let the program to recognize an start and stop codon as a CDs. 

If the input file is on `genbank` format, the program will retrieve the ORFs from according to the annotation of the CDs.

**Note** that in order to be valid, the ORFs mus be in the range of the genome length and be multiple of three (we are simulating evolution in nucleotides according to their position in the codon, so it is important for this information to be precise).  
 
## Unittest
To check that `ovrf` has been properly installed and that you have all the required dependencies, you can enter de `ovrf` folder and type:
```python
$python3 -m unittest
```



