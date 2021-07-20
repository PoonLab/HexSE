# HexSE: Simulating evolution in overlapping reading frames
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
$ python3 -m hexse.run_simulation sars-cov.fa phylo_tree.newick subs_rates.yaml
```

## Installation
`HexSE` is a Python package developed under the version 3.6. The following dependencies are required (as listed in the `setup.py` file):
- `scipy`
- `numpy`
- `biopython`
- `pyyaml`

To install `HexSE`, there are three options. To install using the cloned repository:
```console
$ git clone https://github.com/PoonLab/HexSE
$ cd HexSE
$ sudo python3 setup.py install
```

Ideally, consider using a virtual environment to ensure the right versions for the dependencies:
```console
$ git clone https://github.com/PoonLab/HexSE
$ cd HexSE 
$ python3 -m venv venv
$ source ./venv/bin/activate
$ python3 setup.py install
```

To install without cloning (Using `pip`):
```console
$ python3 -m pip install --upgrade git+https://github.com/PoonLab/HexSE
```

## Unittest
To check that `HexSE` has been properly installed and that you have all the required dependencies, use:

```python
$ python3 -m unittest
```
