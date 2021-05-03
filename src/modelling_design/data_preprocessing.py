import json
import pandas as pd
import warnings

# Omit warning messages from the following code
warnings.filterwarnings("ignore")

# Constant with the .json file path
CONFIG = '../../config/data.json'

# Load .json file
def open_json(config_name):
    """
    Function to open json file and return data
    :params: config_name - json file name
    :return: json data
    """
    with open(config_name) as f:
        return json.load(f)

# Load .csv files
def load_csv_files(hist_features, alg_hist, descrip_features, n_descriptors):
    """
    Function to csv file and return data
    :params: config_name - json file name
    :return: json data
    """
    d_hist_features = {}  # dictionary that will hold them 
    d_descrip_featu = {}
    list_hist_features, list_descrip_features = [], []
    for alg in alg_hist:
        list_hist_features.append(hist_features[0]+alg+".csv")
    for n in n_descriptors:
        list_descrip_features.append(descrip_features[0]+str(n)+".csv")
    for file_name in list_hist_features:  # loop over files
    # read csv into a dataframe and add it to dict with file_name as it key
        d_hist_features[file_name] = pd.read_csv(file_name)
    for file_name in list_descrip_features:  # loop over files
    # read csv into a dataframe and add it to dict with file_name as it key
        d_descrip_featu[file_name] = pd.read_csv(file_name)
    return d_hist_features, d_descrip_featu



def data_preprocessing():
    # Read json
    config = open_json(CONFIG) 
    # Read csv files
    d_hist_features, d_descrip_featu = load_csv_files(config["Files with hist features"], config["Algorithm to extract features"], config["Files with descriptor features"], config["Number of descriptors selected"]) 


if __name__ == '__main__':
    data_preprocessing()


# NECESSARY NOW TO USE model_test.py
# Read json
config = open_json(CONFIG) 
# Read csv files
d_hist_features, d_descrip_featu = load_csv_files(config["Files with hist features"], config["Algorithm to extract features"], config["Files with descriptor features"], config["Number of descriptors selected"]) 
