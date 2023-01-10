from ..sequence_info import Sequence
from ..discretize import discretize
from Bio import SeqIO
from functools import reduce

import scipy.stats as ss
import numpy as np

import yaml
import sys

class Settings:
    """
    Store and contains all parameteres required from the run, 
    sorting them deppending from the file they come from
    TODO: detail the preference order between different sources (args, yaml or sequences)
    """
    def __init__(self, args):

        self.seq_path = args.seq
        self.seq = self.read_sequence(self.seq_path)  # Curated sequence
        self.file_paths = self.clean_file_paths([self.seq_path, args.config])
        self.yaml = self.read_settings_from_yaml(self.get_yaml_file_path(self.file_paths))
        self.tree = args.tree
        self.orfs = self.get_orfs()
        self.pi = self.get_pi()
        self.kappa = self.yaml.get('kappa', 0.3)
        self.global_rate = self.yaml.get('global_rate', 0.5)

        # Handling mu categories
        mu_keys = self.yaml['mu'].keys()
        if 'mu1' in mu_keys:  # mu values are specified as int
            self.mu_values = self.yaml['mu']
        else:  # mu values will be drawn from distribution
            mu_info = {
                            'classes': self.yaml['mu'].get('classes', 4), 
                            'dist': self.yaml['mu'].get('dist', 'gamma'), 
                            'shape': self.yaml['mu'].get('shape', 1),
                            'scale': self.yaml['mu'].get('scale', None)
                            }
        
            self.mu_values = self.draw_mu_values(
                                                mu_info['shape'], 
                                                mu_info['classes'], 
                                                "mu",
                                                mu_info['dist'], 
                                                mu_info['scale']
                                                )
            
    def get_pi(self):
        if 'pi' in self.yaml.keys():
            pi = self.yaml['pi']
     
        else:
            pi = self.calculate_pi(self.seq)

        return pi

    def get_kappa(self):
        return self.kappa

    def get_global_rate(self):
        return self.global_rate

    def get_orfs(self):
        orfs = {}
        if 'orfs' in self.yaml.keys():
            orfs = self.parse_orfs_from_yaml(self.yaml)
        else:
            orfs = self.parse_genbank_orfs(self.seq_path)
        
        return orfs

    def clean_file_paths(self, file_paths):
        return [i for i in file_paths if i]
        
    @staticmethod
    def parse_genbank_orfs(seq_path): #TODO: change name of input
        """
        Extract ORFs from the GenBank file
        """
        orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand

        # Loop through records
        for rec in SeqIO.parse(seq_path, format="genbank"):
            # Read ORFs from GenBank file
            cds = [feat for feat in rec.features if feat.type == "CDS"]
            # Record the first occurrence of the ORFs
            for cd in cds:
                orf = {'coords': []}
                strand = ''

                for loc in cd.location.parts:
                    # Get the strand
                    if loc.strand > 0:
                        strand = '+'
                    else:
                        strand = '-'
                    orf['coords'].append((int(loc.start), int(loc.end)))

                orf_locations[strand].append(orf)

        return orf_locations

    def get_gb_file_path(self, file_paths):
        for file_path in file_paths:
            if file_path.lower().endswith('.gb') or file_path.lower().endswith('genbank'):
                return file_path


    @staticmethod
    def read_settings_from_yaml(yaml_path):
        """
        Parse the config file if it exists
        :param config_path: path to the configuration file
        :return: a dictionary of settings from the configuration file (if no config file exists, the dictionary is empty)
        """
        settings = {}

        if yaml_path:
            with open(yaml_path, 'r') as stream:
                try:
                    settings = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print(e)

        return settings
    
    @staticmethod
    def get_yaml_file_path(file_paths):
        for file_path in file_paths:
            if file_path.lower().endswith('.yaml') or file_path.lower().endswith('.yml'):
                return file_path

    @staticmethod
    def read_sequence(seq_path):
        """
        Extract the sequence from the input file
        :param in_file: path to the file containing the sequence.
            Supported file types are GenBank (.gb, .genbank) and FASTA (.fa, .fasta)
        :return: the sequence as a string
  
        """
        seq = ''
        if seq_path.lower().endswith(".gb") or seq_path.lower().endswith("genbank"):
            # Loop through records
            for rec in SeqIO.parse(seq_path, format="genbank"):
                seq = rec.seq

        elif seq_path.lower().endswith(".fasta") or seq_path.lower().endswith(".fa"):
            # Read in the sequence
            with open(seq_path) as seq_file:
                seq = ''
                for line in seq_file:
                    # Skip header if the file is a FASTA file
                    if not (line.startswith(">") or line.startswith("#")):
                        seq += line.strip('\n\r').upper()

        else:
            pass
            #TODO: throw and exception, catch in main
            # print("Sequence files must end in '.fa', '.fasta', '.gb', 'genbank'")
            # logging.error("Invalid file type: files must end in '.fa', '.fasta', '.gb', 'genbank'")
            # sys.exit(1)

        return seq

    @staticmethod
    def calculate_pi(seq):
        """
        Extracts the value of pi from the command line or a configuration file
        :param pi: the value of pi input from the command line, can be 'None'
        :param s: the sequence as a string
        :return: the value of pi
        """

        return Sequence.get_frequency_rates(seq)

    def define_strand(self, start,end):
        """
        Define is a Reading frame is on + on - strand
        """
        strand = ''
        if start < end:
            strand = '+'
        else:
            strand = '-'

        return strand

    def parse_orfs_from_yaml(self, yaml):
        """
        Reads ORFs from a YAML file containing the ORF coordinates and the parameters of the dN/dS distribution for each ORF
        :param settings: dictionary representation of YAML file
        :return: orf_locations, a dictionary of ORFs sorted by strand (+ or -)
        """
        orf_locations = {'+': [], '-': []}

        raw_coords = list(yaml['orfs'].keys())

        for raw_coord in raw_coords:
            orf = {}

            # Spliced ORF
            if ';' in raw_coord:
                subset_raw = raw_coord.split(';')
                coords = []
                
                for set in subset_raw:
                    raw_set = list(map(int, set.split(',')))  # Convert string to integer
                    coords.append(raw_set)        
                
                # Use first set of coord to define strand
                strand = self.define_strand(coords[0][0], coords[0][1])
                orf['coords'] = coords


            else:
                orf = {}
                coords = raw_coord.split(',')
                coords = list(map(int, coords))  # Convert string to integer
                strand = self.define_strand(coords[0], coords[1])
                orf['coords'] = [coords]

            orf['omega_shape'] = yaml['orfs'][raw_coord]['omega_shape']
            orf['omega_classes'] = yaml['orfs'][raw_coord]['omega_classes']
            dist = yaml['orfs'][raw_coord]['omega_dist']
            orf['omega_scale'] = yaml['orfs'][raw_coord].get('scale', None)  # Default sets scale = 1/shape
            dist = '%s%s' % ('ss.', dist)
            orf['omega_values'] = list(discretize(yaml['orfs'][raw_coord]['omega_shape'],
                                                  yaml['orfs'][raw_coord]['omega_classes'], 
                                                  dist, 
                                                  orf['omega_scale']))

            orf_locations[strand].append(orf)
        
        #print(">> RUNING FROM SETTINGS")
        # pp = pprint.PrettyPrinter(indent=4)
        # pp.pprint(orf_locations)
    
        return orf_locations

    def draw_mu_values(self, alpha, ncat, key_name, dist, scale=None):
        """
        Creates dictionary with values (as rates) drawn from a discretized gamma or lognormal distribution
        :param alpha: int, In the gamma distribution, alpha is the shape. In lognormal distribution shape is the mean
        :param scale: int, Scale parameter
        :param ncat: int, Number of categories
        :param dist: str, the distribution (gamma or lognorm)
        :param key_name: str, i.e., "mu" to create the keys of the dictionary
        """

        nt_categories = discretize(alpha, ncat, dist, scale)
        nt_categories_dict = {}
        for i, item in enumerate(nt_categories):
            cat = f"{key_name}{i+1}"
            nt_categories_dict[cat] = item

        return nt_categories_dict

