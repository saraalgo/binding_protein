import json
import pandas as pd
import warnings
from Bio import SeqIO

# Omit warning messages from the following code
warnings.filterwarnings("ignore")

# Constant with the .json file path
CONFIG = '../../config/Brevundimonas-GCRA.json'

# Load .json file
def open_json(config_name):
    """
    Function to open json file and return data
    :params: config_name - json file name
    :return: json data
    """
    with open(config_name) as f:
        return json.load(f)

# Read getShapeR files
def read_DnaShapeR(algorithms,fasta_name):
    """
    Function to open fasta files and return DnaShapeR data
    :params: algorithms - functions used to calculate spacial rotations
    :params: fasta_name - fasta file name
    :return: shapes - pandas dict with the Sequence with the algorithm and the label
    """    
    shapes = {} # list that will hold them 
    for alg in algorithms:
        shapes[alg] = []
    for alg in algorithms:
        path = fasta_name[0] + alg
    # Open file with "with" statement to avoid problems with access 
    # to original file (in case computer hangs
    # or there will be any other problem)
        with open(path, mode='r') as handle:
            # Use Biopython's parse function to process individual
            # FASTA records (thus reducing memory footprint)
            for record in SeqIO.parse(handle, 'fasta'):
                # Extract individual parts of the FASTA record
                identifier = record.id
                sequence = record.seq
                shapes[alg].append([str(sequence),identifier])
    for alg in algorithms:
        shapes[alg] = pd.DataFrame(shapes[alg], columns=["Sequence", "Label"])
    return shapes

def extact_hist_descriptors(shapes):
    """
    Function to extract descriptors from desity in the histograms of each DnaShape function
    :params: shapes - pandas dict with the Sequence with the algorithm and the label
    :params: algorithms - functions used to calculate spacial rotations
    :return: hist_descriptors
    """ 
    for alg in shapes: # With this loop, 
        for i in range(len(shapes[alg]["Sequence"])):
            shapes[alg]["Sequence"][i] = (shapes[alg]["Sequence"][i].split(","))
            shapes[alg]["Sequence"][i] = [i for i in shapes[alg]["Sequence"][i] if not i.isalnum()]
            shapes[alg]["Sequence"][i] = [float(x) for x in shapes[alg]["Sequence"][i]]
    


# # Load .csv files after feature selection in R
# def load_csv_files(hist_features, alg_hist, descrip_features, n_descriptors):
# """
# Function to csv file and return data
# :params: config_name - json file name
# :return: json data
# """
# d_hist_features = {}  # dictionary that will hold them 
# d_descrip_featu = {}
# list_hist_features, list_descrip_features = [], []
# for alg in alg_hist:
#     list_hist_features.append(hist_features[0]+alg+".csv")
# for n in n_descriptors:
#     list_descrip_features.append(descrip_features[0]+str(n)+".csv")
# for file_name in list_hist_features:  # loop over files
# # read csv into a dataframe and add it to dict with file_name as it key
#     d_hist_features[file_name] = pd.read_csv(file_name)
# for file_name in list_descrip_features:  # loop over files
# # read csv into a dataframe and add it to dict with file_name as it key
#     d_descrip_featu[file_name] = pd.read_csv(file_name)
# return d_hist_features, d_descrip_featu


def data_preprocessing():
    # Read json
    config = open_json(CONFIG) 
    # Read csv files after feature selection in R
    # d_hist_features, d_descrip_featu = load_csv_files(config["Files with hist features"], config["Algorithm to extract features"], config["Files with descriptor features"], config["Number of descriptors selected"]) 
    # Read getShapeR files
    shapes = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR"])



if __name__ == '__main__':
    data_preprocessing()


# NECESSARY NOW TO USE model_test.py
# Read json
config = open_json(CONFIG) 
# Read csv files
#d_hist_features, d_descrip_featu = load_csv_files(config["Files with hist features"], config["Algorithm to extract features"], config["Files with descriptor features"], config["Number of descriptors selected"]) 
# Read getShapeR files
shapes = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR"])