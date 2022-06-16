# Convert input data into a yml object
import yaml
from Bio import SeqIO

def settings_from_files(file_paths):
    """
    :file_paths: list of data files (genbank, fasta, yml, filo_tree, csv)
    :return: simulation config object 
    """
    print(file_paths)
    settings = {}

    gb_file_path = get_gb_file_path(file_paths)
    if gb_file_path:
        settings = gb_settings(gb_file_path, get_yaml_file_path(file_paths))

    # if is_fasta(file_paths):
    #     return fasta_to_yaml(file_paths)
    
    # return parse_yaml(file)
    return settings
    
def get_gb_file_path(file_paths):
    for file_path in file_paths:
        if file_path.lower().endswith('.gb') or file_path.lower().endswith('genbank'):
            return file_path

def get_yaml_file_path(file_paths):
    for file_path in file_paths:
        if file_path.lower().endswith('.yaml') or file_path.lower().endswith('.yml'):
            return file_path

def gb_settings(gb_path, yaml_path=None):
    return {
        'orfs': parse_genbank_orfs(gb_path),
        **read_settings_from_yaml(yaml_path)
    }


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

    print (f">>>>>>>orfs Locations: {orf_locations}/n")
    return orf_locations

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

    print(">>> seetings from yaml\n")
    print(settings, "\n")

    return settings