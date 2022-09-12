# HexSE: Simulating evolution in overlapping reading frames
HexSE is a Python module designed to simulate sequence evolution along a phylogeny while considering the coding context the nucleotides. The ultimate porpuse of HexSE is to account for multiple selection preasures on Overlapping Reading Frames. 

![Pipeline_overview](./images/Pipeline_overview.png)

+ [System requirements](#system-requirements)  
+ [Installation](#installation)  
+ [Usage](#usage)
+ [Test](#test)  
+ [Output files](#output-files)  
+ [Additional features](#additional-features)
+ [Rationale](#rationale)

## System requirements
`HexSE` is a Python package developed under version 3.6.9.
The following dependencies are required (as listed in the `setup.py` file):

* `scipy`
* `numpy`
* `biopython`
* `pyyaml`

## Installation

To install `HexSE`, there are three options. 

* Clone repository:
```console
$ git clone https://github.com/PoonLab/HexSE
$ cd HexSE
$ sudo python3 setup.py install
```

* Create virtual environment (ideal to ensure the right versions for the dependencies):
```console
$ git clone https://github.com/PoonLab/HexSE
$ cd HexSE 
$ python3 -m venv venv
$ source ./venv/bin/activate
$ python3 setup.py install
```

* Install with `pip` (without clonning):
```console
$ python3 -m pip install --upgrade git+https://github.com/PoonLab/HexSE
```

## Usage
HexSe requires users to provide a nucleotide sequence (either on `fasta` or `genbank` format), a phylogenetic tree (on `newick` format), and a configuration file on [YAML](https://en.wikipedia.org/wiki/YAML) format specifying the substitution biases and the coordinates of the Open Reading Frames in the sequence.

HexSE runs as a Python module from the terminal with command line arguments as following:
```console
$ python3 -m hexse.run_simulation <path_to_sequence_file> <path_to_phylogenetic_tree> <path_to_configuration_file> --outfile <path_to_alignment_file> --logfile <path_to_log_file>
```

**Configuration file**

The configuration file contains the parameters that make the simulation possible. It must include a global substitution rate (`global_rate`), transition-transvertion rate ratio (`kappa`), mutation rates (`mu`), stationary nucleotide frequencies (`pi`), and locations of the coding sequences (`orfs`) with specifications of the distrbution from which dN/dS values (`omega`) are drawn (`omega_classes`, `omega_shape`, `omega_dist`):

```python
global_rate: 0.05
kappa: 0.3

pi: 
# Note that keys for pi values HAVE to be the exact nucleotides
  A: 0.25
  T: 0.25
  G: 0.25
  C: 0.25

orfs:   
  1902,2454:
    omega_classes: 5
    omega_shape: 2.5
    omega_dist: gamma

  2308, 3182; 0, 1625:  # From circular genome, this gene is specified as two separated fragments
    omega_classes: 4
    omega_shape: 1.7
    omega_dist: gamma

mu:
  classes: 2
  shape: 1.0
  dist: lognorm

circular: true

```
To create a configuration file from a `.gb` file, users can use `gb_to_yaml.py` available at `./accesory_scripts/`.

**Note:** To declare Coding Sequences resulting from spliced genes or circular genomes, coordinates must be specified as the same `orf` separated by a colon (`;`), in the exact order than the ribosome would read it. For example, a gene from a circular genome encoded from position `2308` to `1625`, is defined in the configuration file as:

```python
  2308, 3182; 0, 1625:
    omega_classes: 4
    omega_shape: 1.7
    omega_dist: gamma
```

## Test
To check that `HexSE` has been properly installed and that you have all the required dependencies, use:

```console
$ python3 -m unittest
```

Test files are provided in `.tests/fisxtures/`. For a run of *HexSE* on [HBV genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2) along a phylogeny with 100 tips, use:
```console
$ python3 -m hexse.run_simulation tests/fixtures/NC_003977.2_HBV.gb tests/fixtures/100_tree.newick tests/fixtures/conf_NC_00377.yaml --logfile test_HBV.log --outfile HBV_out.fasta
```

## Output Files
HeSE will output one alignment file in `fasta` format with as many mutated sequences as tips on the phylogeny. It will also create a log file specifying the run parameters. 

A log file for a test run on HBV looks includes the following information:
```python
INFO:root:
Simulation started at: 2022-09-09 12:30:53.815405

INFO:root:

FILES
	Sequence: tests/fixtures/NC_003977.2_HBV.gb
	Configuration: tests/fixtures/conf_NC_00377.yaml
	Phylo Tree: tests/fixtures/100_tree.newick
	Alignment: HBV_out.fasta

PARAMETERS: 
	Pi: {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
	Global rate: 0.05
	Kappa: 0.3
	Number of nucleotide classification classes: 2
	Nucleotide classification shape parameter: 1.0
	Rates classification values: {'mu1': 0.2615782918648644, 'mu2': 1.3871429788350027}
	
INFO:root:
	Valid ORFs: [[[0, 837]], [[156, 837]], [[1815, 2454]], [[1902, 2454]], [[1375, 1840]], [[1853, 1922]], [[2849, 3182]], [[3173, 3182]]]
	Total ORFs: 8

INFO:root:
	Simulation Ended at: 2022-09-09 15:37:49.275652
	Simulation lasted: 0:00:46.946687 seconds
```

## Aditional Features
*HexSE* includes a set of scripts located at `./accesory_scripts/` that might be usefult to prepare information for a run and to process the alignment afterwards. 

* `gb_to_yaml.py`: Creates a YAML configuration file by identifying the Coding Sequences (CDSs) on a `genbank` file and sets default values to `mu`, `kappa`, `global_rate`, and `pi`.

* `get_orfs.py`: Obtains the Open Reading Frames (ORF) of a genome based on its `genbank` annotations.

* `extract_cds`: From a sequence alignment in `fasta` file, extracts a specific region between two nucleotide locations.    

## Rationale

Gene overlap occurs when two or more genes are encoded by the same nucleotides.
This phenomenon is found in all taxonomic domains, but is particularly common in viruses, where it may provide a mechanism to increase the information content of compact genomes. The presence of overlapping reading frames (OvRFs) can skew estimates of selection based on the rates of non-synonymous and synonymous substitutions, since a substitution that is synonymous in one reading frame may be non-synonymous in another, and vice versa. 

To understand the impact of OvRFs on molecular evolution, HexSE implemented a versatile simulation model of nucleotide sequence evolution along a phylogeny with an arbitrary distribution of reading frames in linear or circular genomes. We use a custom data structure to track the substitution rates at every nucleotide site, which is determined by the stationary nucleotide frequencies, transition bias, and the distribution of selection biases (dN/dS) in the respective reading frames.

![Pipeline_overview](./images/Probability_tree.png)