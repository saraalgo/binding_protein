import json
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
import warnings

# Omit warning messages from the following code
warnings.filterwarnings("ignore")


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

# Read peaks
def read_peaks(csv_name):
    """
    Function to load peaks file and extract peaks
    :params: csv_name - path to your csv peaks file
    :return: (peaks, positions) - peaks info, where each peak starts and ends
    """
    peaks = []
    with open(csv_name)as f:
        for line in f:
            peaks.append(line.strip().split())
    return np.asarray(peaks), np.asarray(peaks)[:,1:3]

# Add no-peaks positions in peaks variable
def add_no_peaks(positions,genome_length):
    """
    Function add no peaks positions in positions variable
    :params: positions - positions variable where peaks positions are
    :params: genome_length - length of your genome
    :return: positions - positions variable uploaded with no peaks information added
    """
    positions = [[int(start),int(end)] for start,end in positions.tolist()] #transform all values to int 
    for i in range(len(positions)): #insert 'peak' label to peaks positions
        positions[i].append('peak')
    add = [0,positions[0][0]-1,'no-peak']
    positions.insert(0,add)
    for i in range(2,len(positions),2): #insert 'no peak' positions with their label in the gaps between peaks
        add = [positions[i-1][1]+1,positions[i][0]-1,'no-peak']
        positions.insert(i,add)
    add= [positions[-1][1],genome_length,'no-peak']
    positions.append(add)
    return positions


# Genome
## Download genome 
def download_genome(genebank_id, email):
    """
    Function to download the full genome of the organism
    :params: genebank_id - Genebank ID of the organism
    :params: email - email to inform to Entrez function (Entrez will send an email if we have exceeded the capacity allowed)
    :return: record - full genome
    """
    handle = Entrez.efetch(db="nucleotide", id=genebank_id, rettype="gb", retmode="text") #with "db" param you indicate that you want the full genome
    Entrez.email = email #let Entrez know who you are
    print('Downloading genome...')
    genome = SeqIO.read(handle, "genbank") # ID fixed to use genebank repository
    handle.close() #CAUTION: essential to close this connection before launching again
    return genome

## Select genic parts of the full genome
def select_genic_parts(genome):
    """
    :params: genome - full genome of the organism
    :return: genic_parts - genic positions
    """
    return [[feature.location.start.position,feature.location.end.position] for feature in genome.features if feature.type == "gene"] 

## Unify genic parts from genic parts contained in other genic parts 
def unify_genic_parts(genic_parts):
    new_genic_parts = [genic_parts[0]]
    for gen in genic_parts[1:]:
        if gen[0] < new_genic_parts[-1][1]:
            new_genic_parts[-1][1]=max(new_genic_parts[-1][1],gen[1])
        else:
            new_genic_parts.append(gen)
    return new_genic_parts

## Delete genic parts of the full genome
def delete_genic_parts(genic_parts, positions): 
    """
    :params: genic_parts - genic positions
    :params: positions - positions variable where peaks positions are
    :return: positions_aux - NEW positions variable where peaks positions are after filtering mode
    """
    positions_aux = positions.copy()
    # Delete positions of the array that are contained inside a genic part
    for start_gen, end_gen in genic_parts: 
        for pos in positions_aux:
            if start_gen <= pos[0] and end_gen >= pos[1]:
                positions_aux.remove(pos)
    # Update an end and the next start by removing a genic part that is placed between those two positions
    for start_gen, end_gen in genic_parts:
        for cnt_pos, pos in enumerate(positions_aux):
            if start_gen >= pos[0] and start_gen <= pos[1] and end_gen >= pos[1]:
                positions_aux[cnt_pos][1] = start_gen  
            if start_gen <= pos[0] and end_gen >= pos[0] and end_gen <= pos[1]:
                positions_aux[cnt_pos][0] = end_gen
                break          
    # Update and insert by splitting a genome part removing genic parts
    for start_gen, end_gen in genic_parts:
        for cnt_pos, pos in enumerate(positions_aux):
            if start_gen > pos[0] and start_gen < pos[1] and end_gen > pos[0] and end_gen < pos[1]:
                new_end = positions_aux[cnt_pos][1]
                positions_aux[cnt_pos][1] = start_gen
                positions_aux.insert(cnt_pos+1,[end_gen,new_end, positions_aux[cnt_pos][2]])
                break
    # Remove duplicated positions
    positions_aux = set(tuple(x) for x in positions_aux)
    # Sort and transform to list
    return sorted([list(x) for x in positions_aux])


## Calculate mean of peak length
def length_peaks(positions, fragment_fixed):
    """
    :params: positions - positions variable where peaks positions are
    :params: fragment_fixed - check if we want to use the mean of actual peaks, or we want to fix a number in our json file
    :return: fragment_size - size of range for extract negative cases from no-peaks positions
    """
    peaks_positions = [item for item in positions if item[2] == "peak"]
    fragment_size = fragment_fixed
    # Check if we have established the mean to grab, or an exact number
    if fragment_fixed == "mean":
        total_len = 0
        for i in range(len(peaks_positions)):
            total_len += peaks_positions[i][1]-peaks_positions[i][0]
        fragment_size = total_len/len(peaks_positions)
    
    return fragment_size

#Extract the sequences of the genome from peaks
def extract_pos_neg(genome, positions, fragment_size):  
    """
    :params: genome - full genome of the organism
    :params: positions - positions variable where peaks positions are
    :params: fragment_size - size of range for extract negative cases from no-peaks positions
    :return: (positives,negatives) - lists with positives and negatives sequences 
    """  
    positives, negatives = [], []
    #From positions where peaks are, extract the corresponding sequence from genome
    peaks_positions = [item for item in positions if item[2] == "peak"]
    for start, end, peak in peaks_positions:
        feature = SeqFeature(FeatureLocation(start,end), type="gene", strand=-1)
        feature_seq = genome.seq[feature.location.start:feature.location.end]
        positives.append([''.join(list(feature_seq)), peak])
    #From positions where no-peaks are, extract the corresponding sequence from genome
    no_peaks_positions = [item for item in positions if item[2] == "no-peak" and item[1]-item[0] >= fragment_size]
    for start, end, peak in no_peaks_positions:
        feature = SeqFeature(FeatureLocation(start,end), type="gene", strand=-1)
        feature_seq = genome.seq[feature.location.start:feature.location.end]
        negatives.append([''.join(list(feature_seq)), peak])

    return positives,negatives

#Save positive and negative data (peak/no-peak)
def save_csv(positives,negatives,config):
    """
    :params: positives - list of sequences for peaks regions
    :params: negatives - list of sequences for no-peaks regions
    :return: None - But a .csv file will be created in your output folder
    """
    #Write lists in a unique .csv file
    with open(config["Output folder"]+config["Output file name"], "w") as positives_file:
        csv_writer = csv.writer(positives_file)
        csv_writer.writerows(positives)
        csv_writer.writerows(negatives)
        

def data_extraction():
    # Read json
    config = open_json(CONFIG) 
    # Feedback of your .json config file
    print('You are working with the Specie: ', config["Species name"])
    print('The transcription factor you are analizing is: ', config["Transcription factor name"])
    print('The genome length of this organism is: ', config["Genome length"])
    # Read peaks and extract genome positions of the peaks
    peaks, positions = read_peaks(config["CSV file with peak center positions"])
    print('This is the first row in the peaks csv file:\n',peaks[0])
    # Add no-peaks positions to previous variable
    positions = add_no_peaks(positions,config["Genome length"])
    print('This is the head of the positions variable:\n',positions[0])
    # Download full genome
    genome = download_genome(config["Reference Genebank ID"],config["Email"])
    # Extract genic positions of the genome
    genic_parts = select_genic_parts(genome)
    # Unify genic positions of the genome, purging those overlapped
    new_genic_parts = unify_genic_parts(genic_parts)
    # IF FILTERING Delete genic parts from genome
    if config["Mode of operation"] == "filtering":
        print("Filtering mode has been used")
        positions = delete_genic_parts(new_genic_parts, positions)
    else: 
        print("Filtering mode has not been used")
    # Calculate fragment size to grab no-peak fragments 
    fragment_size = length_peaks(positions,config["Surrounding peak area to grab (i.e. fragment size)"])
    print("The fragment selected to grab negative data has a size of: ", fragment_size)
    # Extract positives and negatives sequences
    positives,negatives = extract_pos_neg(genome, positions,fragment_size)
    # Save final data to work with
    save_csv(positives,negatives,config)
    print("Your output file has been saved in the following path: ", config["Output folder"])
    print("You will find it with the name of: ", config["Output file name"])
    print("All functions has been correctly executed!")


if __name__ == '__main__':
    data_extraction()
