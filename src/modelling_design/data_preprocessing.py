import json
import pandas as pd
import numpy as np
from itertools import combinations, permutations
import warnings
from Bio import SeqIO
from sklearn.decomposition import PCA

# Omit warning messages from the following code
warnings.filterwarnings("ignore")

# Constant with the .json file path
CONFIG = '../../config/Brevundimonas-GCRA.json'

# EXTRACT DESCRIPTORS
# ================================
## Load .json file
def open_json(config_name):
    """
    Function to open json file and return data
    :params: config_name - json file name
    :return: json data
    """
    with open(config_name) as f:
        return json.load(f)

## Read getShapeR files
def read_DnaShapeR(algorithms,fasta_name):
    """
    Function to open fasta files and return DnaShapeR data
    :params: algorithms - functions used to calculate spacial rotations
    :params: fasta_name - fasta file name
    :return: shapes_raw - pandas dict with the string of the Sequence of the algorithm and the label
    """    
    shapes_raw = {} # list that will hold them 
    for alg in algorithms:
        shapes_raw[alg] = []
    for alg in algorithms:
        path = fasta_name + alg
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
                shapes_raw[alg].append([str(sequence),identifier])
    for alg in algorithms:
        shapes_raw[alg] = pd.DataFrame(shapes_raw[alg], columns=["Sequence", "Label"])
    return shapes_raw

def extact_int_arrays(shapes_raw):
    """
    Function to transform sequences loaded from str to int arrays
    :params: shapes_raw - pandas dict with the string of the Sequence of the algorithm and the label
    :return: shapes, global_min, global_max - pandas dict with the Sequence of the algorithm and the label, global min of each algorithm in order, global max of each algorithm in order
    """ 
    global_min, global_max, alg_min, alg_max = [], [], [], []
    for alg in shapes_raw: # With this loop convert str array of the sequeces to array of int
        for i in range(len(shapes_raw[alg]["Sequence"])):
            splitted = (shapes_raw[alg]["Sequence"][i].split(",")) # split string
            drop_nas = [ x for x in splitted if "NA" not in x ] # drop NAs
            # deal with \n and numbers joined
            aux = []
            for value in drop_nas:
                if sum(c.isdigit() for c in value) > 5: #Because good ones have max 3 or 4 (30.68; -0.35...) 
                    beginning = value[0:(value.find(".")+3)]
                    remaining = value[len(beginning):]
                    aux.append(beginning) 
                    aux.append(remaining)
                else:
                    aux.append(value)
            shapes_raw[alg]["Sequence"][i] = [float(x) for x in aux]
            for value in shapes_raw[alg]["Sequence"][i]: # With this loop we save minimun and maximun value for each algorithm
                if alg_min == [] and alg_max == []:
                    alg_min = value
                    alg_max = value
                elif value < alg_min:
                    alg_min = value
                elif value > alg_max:
                    alg_max = value
        global_min.append(alg_min)
        global_max.append(alg_max)
    return shapes_raw, global_min, global_max

## Extract densities from histograms DNAshapeR
def extract_hist_descriptors(shapes, global_min, global_max):
    """
    Function to extract descriptors from desity in the histograms of each DnaShape function
    :params: shapes - pandas dict with the Sequence of the algorithm and the label
    :params: global_min - global min of each algorithm in order
    :params: global_max - global max of each algorithm in order
    :return: hist_descriptors - pandas dict with the density of the histogram descriptors for each algorithm
    """
    hist_descriptors = {}
    #Normalize algs Sequences values with globals min/max
    for idx, alg in enumerate(shapes):
        for i in range(len(shapes[alg]["Sequence"])):
            for idxx, value in enumerate(shapes[alg]["Sequence"][i]):
                shapes[alg]["Sequence"][i][idxx] = (value-global_min[idx])/(global_max[idx]-global_min[idx])
    #Extract densitites from histograms for each Sequence
    for alg in shapes:
        hist_descriptors[alg] = []
        for i in range(len(shapes[alg]["Sequence"])):
            histograms = np.histogram(shapes[alg]["Sequence"][i],density=True,bins=25)
            hist_descriptors[alg].append([histograms[0],shapes[alg]["Label"][i]]) 
    for alg in shapes:
        hist_descriptors[alg] = pd.DataFrame(hist_descriptors[alg], columns=["Hist_descriptor", "Label"])
    return hist_descriptors


### UNDER CONSTRUCTION, FOR THE MOMENT, USE R FUNCTION

# def extract_k_mers(data, namesmers, kmers):
#     """
#     Function to extract descriptors from desity in the histograms of each DnaShape function
#     :params: data - fasta with all Sequences and Labels
#     :params: kmers - number of nucleotids repeated we want to extract
#     :return: kmers_descriptors - pandas dict with the kmers descriptors for each sequence in each algorithm
#     """
#     data = []
#     # Create dict with keys for kmer descriptors
#     kmers_descriptors = {}
#     for names in namesmers:
#         kmers_descriptors[names] = None
#     # Read data.fa
#     with open(config["Data in fasta"], mode='r') as handle:
#         # Use Biopython's parse function to process individual
#         # FASTA records (thus reducing memory footprint)
#         for record in SeqIO.parse(handle, 'fasta'):
#             # Extract individual parts of the FASTA record
#             identifier = record.id
#             sequence = record.seq
#             data.append([str(sequence),identifier])
#     data = pd.DataFrame(data, columns=["Sequence", "Label"])
#     # Extract kmers combinations counts

# for idx, k in enumerate(config["Number of kmers selected"]):
#     # Set combinations for each k
#     comb = list(combinations("ACTG",k))
#     if k > 1:
#         comb = list(permutations(comb))
#     print(comb[0][0])

#     for seq in data["Sequence"]:


#         kmers_descriptors[namesmers[idx]].append()
                
#     for k in kmers:
#         kmers_descriptors[k] = pd.DataFrame(kmers_descriptors[alg], columns=[(str(k)+" mer_descriptor"), "Label"])
#     return hist_descriptors


# data =  extract_k_mers(config["Data in fasta"], config["Number of kmers selected"])
# kmers_descriptors = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"])

## MEANWHILE LOAD CSV DATA OF K-MERS FROM R OUTPUT
def load_csv_kmers(kmers, kmers_path):
    kmers_descriptors = {}
    for k in kmers:
        kmers_descriptors[k] = pd.read_csv(kmers_path+str(k)+"_descriptors.csv")
    return kmers_descriptors


# FEATURE SELECTION TECHNIQUES
# ================================
## PCA
def pca(ncomponents, descriptors):
    descriptors_PCA = {}
    pca = PCA(n_components=ncomponents)
    for var in descriptors:
        if len(descriptors[var].columns == 2): # for hist_descriptors
            x = descriptors[var][descriptors[var].columns[0]]
            x = pd.DataFrame(x.tolist(), index= x.index)
            y = descriptors[var][descriptors[var].columns[1]]
            principalComponents = pca.fit_transform(x)
        elif len(descriptors[var].columns) > ncomponents: # for kmers_descriptors   PROBLEMA AQUI!!!!!!!!!!!!!!!!!!!!!
            print(descriptors[var])
            x = descriptors[var].drop("Label",axis=1)
            print(x)
            y = descriptors[var]["Label"]
            print(y)
            principalComponents = pca.fit_transform(x)
        else: # for kmers_descriptors in case there are not enough
            principalComponents = descriptors[var].drop("Label",axis=1)
            y = descriptors[var]["Label"]
        descriptors_PCA[var] = [principalComponents,y]
    return descriptors_PCA

prueba = pca(15, kmers_descriptors)


def data_preprocessing():
    # Read json
    config = open_json(CONFIG) 
    # Read getShapeR files
    shapes_raw = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR"])
    # Convert to int arrays Sequences from getShapeR
    shapes, global_min, global_max = extact_int_arrays(shapes_raw)
    # Function to extract descriptors from desity in the histograms of each DnaShape function
    hist_descriptors = extract_hist_descriptors(shapes, global_min, global_max)
    # Function to extract descriptors from kmers
    # kmers_descriptors = extract_k_mers(config["Data in fasta"], config["Prefix of kmers selected"], config["Number of kmers selected"])
    # Load csv kmers descriptors from R
    kmers_descriptors = load_csv_kmers(config["Number of kmers selected"], config["Path to csv output k-mers"])
    # Apply FS techniques
    ## PCA
    hist_descriptors_PCA = 
    kmers_descriptors = 

if __name__ == '__main__':
    data_preprocessing()


# NECESSARY NOW TO USE model_test.py
# Read json
config = open_json(CONFIG) 
# Read csv files
#d_hist_features, d_descrip_featu = load_csv_files(config["Files with hist features"], config["Algorithm to extract features"], config["Files with descriptor features"], config["Number of descriptors selected"]) 
# Read getShapeR files
shapes_raw = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR"])
# Convert to int arrays Sequences from getShapeR
shapes, global_min, global_max = extact_int_arrays(shapes_raw)
# Function to extract descriptors from desity in the histograms of each DnaShape function
hist_descriptors = extract_hist_descriptors(shapes, global_min, global_max)
# Function to extract descriptors from kmers
# kmers_descriptors = extract_k_mers(config["Data in fasta"], config["Prefix of kmers selected"], config["Number of kmers selected"])
# Load csv kmers descriptors from R
kmers_descriptors = load_csv_kmers(config["Number of kmers selected"], config["Path to csv output k-mers"])
# Apply FS techniques
## PCA
hist_descriptors_PCA = pca(15, hist_descriptors)
kmers_descriptors_PCA = pca(15, kmers_descriptors)


pca = PCA(n_components=15)

x = hist_descriptors["HelT"]["Hist_descriptor"]
x = pd.DataFrame(x.tolist(), index= x.index)
