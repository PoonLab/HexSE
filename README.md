# Simulating evolution in Overlapping Reading Frames (OvRFs)
Python module to simulate mutations on a nucleotide sequence along a phylogenetic tree taking into account the coding context of the nucleotide in the sequence.


## Usage
The minimum requirenments for the pipeline are:

**Nucleotide sequence** The sequence to evolve can be either in `fasta` of `genbank` format.

**Phylogenetic tree** expected on `newick` format

**Open reading frames**
If the input file is on `genbank` format, the program will retrieve the ORFs from the genome annotation. If the input sequence is a `fasta` file, you have two options:
1. Specify the ORFs with the first and last nucleotide of the CDs as an csv file
2. Let the program search for CDSs by recognizing start and stop codons

**Optional parameters**
Other substitution biases can be specified in a [YAML](https://en.wikipedia.org/wiki/YAML) file.

The pipeline can be run as a Python module from the terminal as following:

```console
$ python3 -m src.run_simulation sars-cov.fa phylo_tree.newick subs_rates.yaml
```

## Installation


## Unittest
To check that `ovrf` has been properly installed and that you have all the required dependencies, you can enter de `ovrf` folder and type:

```python
$python3 -m unittest
```
