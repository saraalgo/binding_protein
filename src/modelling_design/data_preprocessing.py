import json
import pandas as pd
import numpy as np
import itertools
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif


# Omit warning messages from the following code
warnings.filterwarnings("ignore")

# Constant with the .json file path
CONFIG = '../../config/Brevundimonas-GCRA.json'

# READ CONFIG DATA
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

# EXTRACT DNASHAPER DESCRIPTORS
# ================================
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

# EXTRACT DNASHAPER DESCRIPTORS
# ================================
def charCombination(kmers):
    """
    Subfunciton to extract descriptors from desity in the histograms of each DnaShape function
    :params: kmers - number of nucleotids repeated we want to extract
    :return: all posible combinations for ACTG kmers
    """
    return ["".join(item) for item in itertools.product("ATCG", repeat=kmers)]

def extract_k_mers(datafa, kmers, complement):
    """
    Function to extract descriptors from desity in the histograms of each DnaShape function
    :params: datafa - fasta with all Sequences and Labels
    :params: kmers - number of nucleotids repeated we want to extract
    :params: complement - here you indicate if you want to count reverso complement too of each kmer
    :return: kmers_descriptors - pandas dict with the kmers descriptors for each sequence in each algorithm
    """
    data = []
    # Extract kmers combinations counts as key of the dict
    kmers_combinations = {}
    keys_combinations = charCombination(kmers)
    for i in keys_combinations:
        kmers_combinations[i] = int()
    # Read data.fa
    with open(datafa, mode='r') as handle:
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            data.append([str(sequence),identifier])
    data = pd.DataFrame(data, columns=["Sequence", "Label"])
    # Count how many strings of kmers each sequence has
    kmer_descriptors = []
    for seq in data["Sequence"]:
        for i in keys_combinations:
            kmers_combinations[i] = int()
        for idx, _ in enumerate(seq):
            if idx < len(seq)-kmers-1:
                kmers_combinations[seq[idx:idx+kmers]] += 1
                if complement == "Yes":
                    kmers_combinations[Seq(seq[idx:idx+kmers]).reverse_complement()] += 1 # Is counting the reverse complement of each kmer too
        kmer_descriptors.append([number / len(seq) for number in list(kmers_combinations.values())])
    
    kmer_descriptors=pd.DataFrame(kmer_descriptors,columns=keys_combinations)
    kmer_descriptors.index = data["Label"]
    kmer_descriptors = kmer_descriptors.loc[:, (kmer_descriptors != 0).any(axis=0)]
    return kmer_descriptors


# FEATURE SELECTION TECHNIQUES
# ================================
## PCA
def pca(ncomponents, descriptors):
    """
    Function to apply PCA FS technique
    :params: ncomponents - number os components to calculate by PCA
    :params: descriptors - descriptors we want to reduce features
    :return: descriptors_PCA - PCA features selected
    """
    pca = PCA(n_components=ncomponents)
    if type(descriptors) == dict:
        descriptors_PCA = {}
        for var in descriptors:
            if len(descriptors[var].columns == 2): # for hist_descriptors
                x = descriptors[var][descriptors[var].columns[0]]
                x = pd.DataFrame(x.tolist(), index= x.index)
                y = descriptors[var][descriptors[var].columns[1]]
                pca.fit_transform(x)
                columns = ['pca_%i' % i for i in range(ncomponents)]
            descriptors_PCA[var] = pd.DataFrame(pca.transform(x), columns=columns, index=y)
    else:
        pca.fit_transform(descriptors)
        columns = ['pca_%i' % i for i in range(ncomponents)]
        descriptors_PCA = pd.DataFrame(pca.transform(descriptors), columns=columns, index=descriptors.index)
    return descriptors_PCA

## Univariate method 
def univariate_FS(kbest,descriptors):
    if type(descriptors) == dict:
        descriptors_kbest = {}
        for var in descriptors:
            if len(descriptors[var].columns == 2): # for hist_descriptors
                x = descriptors[var][descriptors[var].columns[0]]
                x = pd.DataFrame(x.tolist(), index= x.index)
                y = descriptors[var][descriptors[var].columns[1]]
                y = y.str.replace('\d+', '')
                sel_anova = SelectKBest(f_classif, k=kbest).fit(x,y.ravel())
                kbest_fs = x.columns[sel_anova.get_support()]
                kbest_fs = x.iloc[:,kbest_fs]
                kbest_fs.index = y
            descriptors_kbest[var] = kbest_fs
    else:
        y = descriptors.index.str.replace('\d+', '')
        sel_anova = SelectKBest(f_classif, k=kbest).fit(descriptors,y)
        descriptors_kbest = descriptors.columns[sel_anova.get_support()]
        descriptors_kbest = descriptors.loc[:,descriptors_kbest]
    return descriptors_kbest


# SAVE DESCRIPTORS
# ================================
#Save descriptors
def save_csv(descriptor, filename, output_path):
    """
    :params: descriptor - descriptors we want to save
    :params: config - to call the path you want to save them
    :params: filename - name by which the file will be named
    :return: None - But a .csv file will be created in your output folder
    """
    #Separate data if they are in a dict to save them
    if len(descriptor.columns) == 2:
        descriptor = descriptor.set_index("Label")
        descriptor.iloc[:,0] = [list(row) for row in descriptor.iloc[:,0]]
        data = pd.DataFrame(descriptor.iloc[:,0].tolist(), columns=range(len(descriptor.iloc[:,0][0])), index=descriptor.index)
        #Save DataFrames as .csv
        data.to_csv(output_path+filename+".csv")
    else:
        descriptor.to_csv(output_path+filename+".csv")

     
# MAIN
# ================================
def data_preprocessing():
    # Read json
    config = open_json(CONFIG) 
    # Read getShapeR files
    print("Extracting DNAShapeR descriptors.....")
    shapes_raw = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR"])
    shapes_raw_b = read_DnaShapeR(config["Algorithm to extract features"],config["Files DnaShapeR_bootstraping"])
    # Convert to int arrays Sequences from getShapeR
    shapes, global_min, global_max = extact_int_arrays(shapes_raw)
    shapes_b, global_min_b, global_max_b = extact_int_arrays(shapes_raw_b)
    # Function to extract descriptors from desity in the histograms of each DnaShape function
    hist_descriptors = extract_hist_descriptors(shapes, global_min, global_max)
    hist_descriptors_b = extract_hist_descriptors(shapes_b, global_min_b, global_max_b)
    # Function to extract descriptors from kmers
    print("Extracting kmers descriptors.....")
    monomers = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][0], "Yes")
    dimers = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][1], "Yes")
    tetramers = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][2], "Yes")
    monomers_nocomp = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][0], "No")
    dimers_nocomp = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][1], "No")
    tetramers_nocomp = extract_k_mers(config["Data in fasta"], config["Number of kmers selected"][2], "No")
    monomers_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][0], "Yes")
    dimers_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][1], "Yes")
    tetramers_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][2], "Yes")
    monomers_nocomp_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][0], "No")
    dimers_nocomp_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][1], "No")
    tetramers_nocomp_b = extract_k_mers(config["Data in fasta_bootstraping"], config["Number of kmers selected"][2], "No")
    # Apply FS techniques
    ## PCA 15
    print("Applying FS techniques: PCA, kbest.....")
    hist_descriptors_PCA = pca(15, hist_descriptors)
    tetramers_descriptors_PCA = pca(15, tetramers)
    tetramers_nocomp_descriptors_PCA = pca(15, tetramers_nocomp)
    hist_descriptors_PCA_b = pca(15, hist_descriptors_b)
    tetramers_descriptors_PCA_b = pca(15, tetramers_b)
    tetramers_nocomp_descriptors_PCA_b = pca(15, tetramers_nocomp_b)
    ## Univariate kbest 15
    hist_descriptors_kbest = univariate_FS(15, hist_descriptors)
    tetramers_descriptors_kbest = univariate_FS(15, tetramers)
    tetramers_nocomp_descriptors_kbest = univariate_FS(15, tetramers_nocomp)
    hist_descriptors_kbest_b = univariate_FS(15, hist_descriptors_b)
    tetramers_descriptors_kbest_b = univariate_FS(15, tetramers_b)
    tetramers_nocomp_descriptors_kbest_b = univariate_FS(15, tetramers_nocomp_b)
    # Save results
    print("Saving descriptors features.....")
    save_csv(hist_descriptors["HelT"], "HelT_descriptors", config["Path to csv output descriptors"])
    save_csv(hist_descriptors["MGW"], "MGW_descriptors", config["Path to csv output descriptors"])
    save_csv(hist_descriptors["ProT"], "ProT_descriptors", config["Path to csv output descriptors"])
    save_csv(hist_descriptors["Roll"], "Roll_descriptors", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_PCA["HelT"], "HelT_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_PCA["MGW"], "MGW_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_PCA["ProT"], "ProT_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_PCA["Roll"], "Roll_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_kbest["HelT"], "HelT_descriptors_kbest", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_kbest["MGW"], "MGW_descriptors_kbest", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_kbest["ProT"], "ProT_descriptors_kbest", config["Path to csv output descriptors"])
    save_csv(hist_descriptors_kbest["Roll"], "Roll_descriptors_kbest", config["Path to csv output descriptors"])  
    save_csv(monomers, "monomers_descriptors", config["Path to csv output descriptors"])
    save_csv(dimers, "dimers_descriptors", config["Path to csv output descriptors"])
    save_csv(tetramers, "tetramers_descriptors", config["Path to csv output descriptors"])
    save_csv(monomers_nocomp, "monomers_nocomp_descriptors", config["Path to csv output descriptors"])
    save_csv(dimers_nocomp, "dimers_nocomp_descriptors", config["Path to csv output descriptors"])
    save_csv(tetramers_nocomp, "tetramers_nocomp_descriptors", config["Path to csv output descriptors"])
    save_csv(tetramers_descriptors_PCA, "tetramers_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(tetramers_nocomp_descriptors_PCA, "tetramers_nocomp_descriptors_PCA", config["Path to csv output descriptors"])
    save_csv(tetramers_descriptors_kbest, "tetramers_descriptors_kbest", config["Path to csv output descriptors"])
    save_csv(tetramers_nocomp_descriptors_kbest, "tetramers_nocomp_descriptors_kbest", config["Path to csv output descriptors"])
    print("Saving descriptors features with bootstraping.....")
    save_csv(hist_descriptors_b["HelT"], "HelT_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_b["MGW"], "MGW_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_b["ProT"], "ProT_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_b["Roll"], "Roll_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_PCA_b["HelT"], "HelT_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_PCA_b["MGW"], "MGW_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_PCA_b["ProT"], "ProT_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_PCA_b["Roll"], "Roll_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_kbest_b["HelT"], "HelT_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_kbest_b["MGW"], "MGW_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_kbest_b["ProT"], "ProT_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])
    save_csv(hist_descriptors_kbest_b["Roll"], "Roll_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])  
    save_csv(monomers_b, "monomers_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(dimers_b, "dimers_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_b, "tetramers_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(monomers_nocomp_b, "monomers_nocomp_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(dimers_nocomp_b, "dimers_nocomp_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_nocomp_b, "tetramers_nocomp_descriptors", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_descriptors_PCA_b, "tetramers_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_nocomp_descriptors_PCA_b, "tetramers_nocomp_descriptors_PCA", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_descriptors_kbest_b, "tetramers_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])
    save_csv(tetramers_nocomp_descriptors_kbest_b, "tetramers_nocomp_descriptors_kbest", config["Path to csv output descriptors bootstrapping"])

if __name__ == '__main__':
    data_preprocessing()

