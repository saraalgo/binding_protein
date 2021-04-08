import json

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

# Load .csv and fasta files
def open_data(config_name):
    



def data_preprocessing():
    # Read json
    config = open_json(CONFIG) 
