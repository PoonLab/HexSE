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
+ [Workflow](#workflow)

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

## Test
To check that `HexSE` has been properly installed and that you have all the required dependencies, use:

```console
$ python3 -m unittest
```

Test files are provided in `.tests/fisxtures/`. For a test run of HexSE on [HBV genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2) along a phylogeny with 100 tips, use:
```console
$ python3 -m hexse.run_simulation tests/fixtures/NC_003977.2_HBV.gb tests/fixtures/100_tree.newick tests/fixtures/conf_NC_00377.yaml --logfile test_HBV.log --outfile HBV_out.fasta
```

## Output Files

## Aditional Features

## Rationale

## Workflow