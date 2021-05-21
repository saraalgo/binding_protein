import json
import os
import numpy as np
import pandas as pd
import warnings
import math
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
import random
from pandas import DataFrame

# Omit warning messages from the following code
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings("ignore")

#0 = all messages are logged (default behavior)
#1 = INFO messages are not printed
#2 = INFO and WARNING messages are not printed
#3 = INFO, WARNING, and ERROR messages are not printed 

# Constant with the .json file path
CONFIG = '../../config/Brevundimonas-GcrA.json'

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
# ESTA FUNCIÓN TIENE QUE SER AUTOMATIZADA PARA LOS DATOS DE PEAKS QUE VAYAMOS A DESCARGAR!
def read_peaks(csv_name):
    """
    Function to load peaks file and extract peaks
    :params: csv_name - path to your csv peaks file
    :return: (peaks, peak_positions) - peaks info, where each peak starts and ends
    """
    peaks = []
    with open(csv_name)as f:
        for line in f:
            peaks.append(line.strip().split())
    return np.asarray(peaks), [item[1:3] for item in peaks]

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
    Function to select only genic regions from full genome from gene info
    :params: genome - full genome of the organism
    :return: genic_parts, preverse - genic positions, probability of a gene to be in the strand -1
    """
    genic_parts = [[feature.location.start.position,feature.location.end.position,feature.location.strand, "genic"] for feature in genome.features if feature.type == "gene"]
    cnt = 0
    for gen in genic_parts:
        if gen[2]==-1:
            cnt+=1
    preverse = cnt/len(genic_parts)
    return genic_parts, preverse

## Unify genic parts from genic parts contained in other genic parts 
def unify_genic_parts(genic_parts):
    """
    Function to unify genic positions which are overlapped or with a distance less than 100 bp and in the SAME strand
    :params: genic_parts - genic positions
    :return: new_genic_parts - genic positions unificated  
    """
    new_genic_parts = [genic_parts[0]]
    for gen in genic_parts[1:]: 
        if gen[0] > gen[1]:
            temp = gen[1]
            gen[1] = gen[0]
            gen[0] = temp
    for gen in genic_parts[1:]:
        if gen[0] < (new_genic_parts[-1][1]+100) and gen[2] == new_genic_parts[-1][2]:
            new_genic_parts[-1][1]=max(new_genic_parts[-1][1],gen[1])
        else:
            new_genic_parts.append(gen)
    return new_genic_parts

## Add intergenic positions in genic variable
def add_intergenic(new_genic_parts,genome_length):
    """
    Function add intergenic positions in genic_parts variable (start,end,genic/intergenic,strand)
    :params: new_genic_parts - positions variable where genic positions are
    :params: genome_length - length of your genome
    :return: genic_nogenic_parts - positions variable where genic and intergenic positions are
    """
    genic_nogenic_parts = new_genic_parts.copy()
    add = [0,genic_nogenic_parts[0][0]-1,1,'intergenic']
    genic_nogenic_parts.insert(0,add)
    for i in range(2,(len(new_genic_parts)*2),2): #insert 'intergenic' positions with their label in the gaps between genic parts
        add = [genic_nogenic_parts[i-1][1]+1,genic_nogenic_parts[i][0]-1,1,'intergenic']
        genic_nogenic_parts.insert(i,add)
    add= [genic_nogenic_parts[-1][1]+1,genome_length,1,'intergenic']
    genic_nogenic_parts.append(add)
    for gen in genic_nogenic_parts: 
        if gen[0] > gen[1]:
            temp = gen[1]
            gen[1] = gen[0]
            gen[0] = temp
    return genic_nogenic_parts

## Add prepromoter and postpromoter positions
### Treat no-genic parts to distinguish between upstream/downstream
def no_genic_treatment(positions, start, end, strand, label):
    """
    Subfunction to stablish how to divide upstream/downstream in intergenic
    :params: positions - new list to identify genic/intergenic/upstream/downstream parts in the genome
    :params: start - position where each an intergenic part starts
    :params: end - position where each an intergenic part ends
    :params: strand - strand of each intergenic part 
    :params: label - indication to be an intergenic part
    :return: positions - function to use to identify intergenic/upstream/downstream parts in the genome
    """
    if abs(end-start) > 250 and strand == 1:
        positions.append([start,start+100, 1, ["downstream"]])
        positions.append([start+101, end-149, 1, [label]])
        positions.append([end-150, end, 1, ["upstream"]])
    elif abs(end-start) == 250 and strand == 1:
        positions.append([start,start+100, 1, ["downstream"]])
        positions.append([end-150, end, 1, ["upstream"]])
    elif abs(end-start) < 250 and strand == 1:
        positions.append([start,(start+int(np.floor(abs((end-start)/2)))), 1, ["downstream"]])
        if int(abs(end-start)%2) == 0:
            positions.append([((start+int(np.floor(abs((end-start)/2))))+1),end, 1, ["upstream"]])
        else:
            positions.append([(start+int(np.ceil(abs((end-start)/2)))),end, 1, ["upstream"]])
    return positions

### Loop to create the final positions
def add_pre_post_positions(genic_nogenic_parts):
    """
    Function add upstream/downstream positions from intergenic regions
    :params: genic_nogenic_parts - positions variable where genic and intergenic positions are
    :return: positions - list to identify genic/intergenic/upstream/downstream parts in the genome
    """
    first_positions = genic_nogenic_parts[0][0],genic_nogenic_parts[0][1]-149, genic_nogenic_parts[0][2], [genic_nogenic_parts[0][3]]
    second_positions =  genic_nogenic_parts[0][1]-150,genic_nogenic_parts[0][1], genic_nogenic_parts[0][2], ["upstream"]
    positions = [first_positions,second_positions]
    for start,end,strand,label in genic_nogenic_parts[1:-1]:
        if label == "genic":
            positions.append([start,end,strand,[label]])
        else:
            no_genic_treatment(positions, start, end, strand, label)
    if len(genic_nogenic_parts[-1])<100:
        positions.append([genic_nogenic_parts[-1][0],genic_nogenic_parts[-1][1], strand, ["downstream"]])
    else:
        positions.append([genic_nogenic_parts[-1][0],genic_nogenic_parts[-1][0]+100, strand, ["downstream"]])
        positions.append([genic_nogenic_parts[-1][0]+101,genic_nogenic_parts[-1][1], strand, [genic_nogenic_parts[0][3]]])
    return positions


## Create variables of the sequences from each regions positions
def create_region_sequences(positions,genome):
    """
    Function to extract a unique sequence of nucleotids from correponding positions in each region 
    :params: positions - list to identify genic/intergenic/upstream/downstream parts in the genome
    :return: region_seq - dict with the whole sequence from correponding parts (genic, intergenic, upstream, downstream)
    """
    region_seq = {
    "genic": [],
    "intergenic": [],
    "upstream": [],
    "downstream": []    
    }
    for start, end, strand, label in positions:
        for l in label:
            if l != "genic" and start>end:
                break
            seqfeature = SeqFeature(FeatureLocation(start,end), strand=strand)
            seq = genome.seq[seqfeature.location.start:seqfeature.location.end]
            if seqfeature.strand == -1:
                seq = seq.reverse_complement()
            region_seq[l].append(seq)
    region_seq["genic"] = "".join([str(seq_rec) for seq_rec in region_seq["genic"]])
    region_seq["intergenic"] = "".join([str(seq_rec) for seq_rec in region_seq["intergenic"]])
    region_seq["upstream"] = "".join([str(seq_rec) for seq_rec in region_seq["upstream"]])
    region_seq["downstream"] = "".join([str(seq_rec) for seq_rec in region_seq["downstream"]])
    return region_seq

# Peak treatment
## Filter (if necessary) peaks to unify length according to the minimun peak length
def filter_peak_length(peak_positions, operation_mode):
    if operation_mode == "filtering":
        temp = [abs(int(end) - int(start)) for start, end in peak_positions]
        res = min(temp) # min peak length
        for idx,(start,end) in enumerate(peak_positions):
            diff = abs(int(end) - int(start))-res
            diff_side = diff/2
            peak_positions[idx][0] = str(int(start)+math.floor(diff_side))
            if (diff_side % 2) == 0 and diff != 0:
                peak_positions[idx][1] = str(int(end)-math.floor(diff_side))
            else:
                peak_positions[idx][1] = str(int(end)-math.ceil(diff_side))
    else: 
        pass
    return peak_positions

## Classify each peak
def classify_peak(peak_positions,positions):
    """
    Function to classify peaks according to the regions where they are found and counting bp in each region
    :params: peak_positions - start, end and strand of each peak
    :params: positions - list to identify genic/intergenic/upstream/downstream parts in the genome
    :return: peak_classification - info of each peak with the regions where they are found and how many bp they have in each one
    """
    peak_positions = [[int(start),int(end)] for start,end in peak_positions] #transform all values from peak list to int 
    peak_classification = []
    for start_p, end_p in peak_positions:
        bp_dict = {
        "genic": [],
        "intergenic": [],
        "upstream": [],
        "downstream": []
        }
        strand_aux,label_aux = None,[]
        for start, end, strand, label in positions:
            if start_p >= start and end_p <= end:
                for l in label:
                    bp_dict[l].append(end_p-start_p) 
                    label_aux.append(l)
                strand_aux = strand
                break
            elif start_p >= start and start_p <= end and end_p >= end:
                for l in label:
                    bp_dict[l].append(end-start_p)
                    label_aux.append(l)
            elif start_p <= start and end_p >= end:
                for l in label:
                    bp_dict[l].append(end-start)
                    label_aux.append(l)
            elif start_p <= start and end_p <= end and end_p >= start:
                for l in label:
                    bp_dict[l].append(end_p-start)
                    label_aux.append(l)
                strand_aux = strand
                break
        peak_classification.append([start_p,end_p,strand_aux,label_aux,bp_dict])
    return peak_classification

### Function to extract the nucleotids we need from the sequences
def pick_seq(peak_classification_aux, regions_seq,preverse):
    """
    Subfunction to indicate how to extract the sequences from the bp and the region of the peaks
    :params: peak_classification_aux - info of each peak with the regions where they are found and how many bp they have in each one
    :params: regions_seq - dict with the whole sequence from correponding parts (genic, intergenic, upstream, downstream)
    :params: preverse - probability of a gene to be in the strand -1
    :return: sublist - list of replicates of a list of sequences
    """
    sublist = []
    for peak in peak_classification_aux:
        replicate_sequence = [] 
        for par in peak:
            inicio_rdn = random.randint(0,len(regions_seq[par[0]])-par[1])
            if replicate_sequence == []:
                replicate_sequence = regions_seq[par[0]][inicio_rdn:inicio_rdn+par[1]]
            # Set that in the case of upstream/downstream, because we don't have strand info, the % of probability of the gene to be strand -1 and 50% of possibility 
            # of belong to each strand that part of the genome, will be flipped to decide if we let the sequence in the positive strand or we apply the reverse/complement
            else:
                prob = np.random.binomial(1, (0.5*preverse)) #50% of posibility to be on each strand * preverse% of genes = strand -1 
                if regions_seq[par[0]]=="downstream" or regions_seq[par[0]]=="upstream" and prob==1:
                    replicate_sequence = replicate_sequence + str(Seq(regions_seq[par[0]][inicio_rdn:inicio_rdn+par[1]]).reverse_complement())
                else:
                    replicate_sequence = replicate_sequence + (regions_seq[par[0]][inicio_rdn:inicio_rdn+par[1]])
        sublist.append(replicate_sequence)
    return sublist

## Replicate peaks (negatives)
def peak_replication(peak_classification, regions_seq,n_replicates,preverse):
    """
    Function to replicate peaks
    :params: peak_classification - info of each peak with the regions where they are found and how many bp they have in each one
    :params: regions_seq - dict with the whole sequence from correponding parts (genic, intergenic, upstream, downstream)
    :params: n_replicates - number of replicates we want to generate
    :params: preverse - probability of a gene to be in the strand -1
    :return: replicates - list of replicates of a list of sequences done n_replicates times (negative data)
    """
    replicates = []
    peak_classification_aux = []
    for idx, (_, _, _, label, bp) in enumerate(peak_classification):
        peak_classification_aux.append([(label,peak_classification[idx][4][label].pop(0)) for label in peak_classification[idx][3]]) # We are selecting label with 3 and bp with 4 
    for i in range(n_replicates):
        replicates += pick_seq(peak_classification_aux, regions_seq,preverse)
    return replicates

## Filter (if necessary) replicates to unify length according to the minimun peak length
def filter_replicate_length(peak_positions, operation_mode):
    if operation_mode == "filtering":
        temp = [abs(int(end) - int(start)) for start, end in peak_positions]
        res = min(temp) # min peak length
        for idx,(start,end) in enumerate(peak_positions):
            diff = abs(int(end) - int(start))-res
            diff_side = diff/2
            peak_positions[idx][0] = str(int(start)+math.floor(diff_side))
            if (diff_side % 2) == 0 and diff != 0:
                peak_positions[idx][1] = str(int(end)-math.floor(diff_side))
            else:
                peak_positions[idx][1] = str(int(end)-math.ceil(diff_side))
    else: 
        pass
    return peak_positions

## Extract nucleotics sequences from peaks data (positives)
def peak_seq(peak_positions, genome):
    """
    Function to extract sequences of the peaks from the genome
    :params: peak_positions - start, end and strand of each peak
    :params: genome - full genome of the organism
    :return: positives - list of sequences of the peak list
    """
    peak_positions = [[int(start),int(end)] for start,end in peak_positions] #transform all values from peak list to int
    positives = []
    for start, end in peak_positions:
        seqfeature = SeqFeature(FeatureLocation(start,end))
        seq = genome.seq[seqfeature.location.start:seqfeature.location.end]
        seq = str(seq)
        positives.append(seq)
    return positives

#Save positive and negative data (peaks/replicates)
def save_csv(positives,replicates,config):
    """
    :params: positives - list of sequences for peaks regions
    :params: negatives - list of sequences for no-peaks regions
    :return: None - But a .csv file will be created in your output folder
    """
    #Add label of peak or replicate
    positives = DataFrame(positives,columns=["Sequence"])
    positives['Label'] = 'Peak'
    negatives = DataFrame(replicates,columns=["Sequence"])
    negatives['Label'] = 'Replicate'
    data = positives.append(negatives)
    #Save DataFrames as .csv
    data.to_csv(config["Output folder"]+config["Output file name"])


def data_extraction_peak_replica():
    # Read json
    config = open_json(CONFIG) 
    # Feedback of your .json config file
    print('You are working with the Specie:', config["Species name"])
    print('The transcription factor you are analizing is:', config["Transcription factor name"])
    print('The genome length of this organism is:', config["Genome length"])
    # Read peaks and extract genome positions of the peaks
    peaks, peak_positions = read_peaks(config["CSV file with peak center positions"])
    print('This is the first row in the peaks csv file:\n',peaks[0])
    # Download full genome
    genome = download_genome(config["Reference Genebank ID"],config["Email"])
    # Extract genic positions of the genome
    genic_parts, preverse = select_genic_parts(genome)
    # Unify genic positions of the genome, purging those overlapped in the same strand
    new_genic_parts = unify_genic_parts(genic_parts)
    # Add no-genic positions to previous variable
    genic_nogenic_parts = add_intergenic(new_genic_parts,config["Genome length"])
    # Add upstream/downstream positions
    positions = add_pre_post_positions(genic_nogenic_parts)
    print('This is the head of the positions variable:\n',positions[0:4])
    # Create regions sequence dictionary
    regions_seq = create_region_sequences(positions,genome)
    # Filter (if necessary) peaks to unify length
    peak_positions = filter_peak_length(peak_positions, config["Mode of operation"])
    # Peak positions to replicate
    peak_classification = classify_peak(peak_positions,positions)
    # Replicates n times
    replicates = peak_replication(peak_classification,regions_seq,config["Number of replicates of each peak"],preverse)
    # Positive nucleotids, positives
    positives = peak_seq(peak_positions, genome)
    # Save final data to work with
    save_csv(positives,replicates,config)
    print("Your output file has been saved in the following path:", config["Output folder"])
    print("You will find it with the name of:", config["Output file name"])
    print("All functions has been correctly executed! :)")


if __name__ == '__main__':
    data_extraction_peak_replica()