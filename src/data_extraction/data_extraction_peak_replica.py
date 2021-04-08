import json
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
import warnings
import random
from pandas import DataFrame

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
    return np.asarray(peaks), np.asarray(peaks)[:,1:3]

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
    :return: genic_parts - genic positions
    """
    return [[feature.location.start.position,feature.location.end.position,feature.location.strand, "genic"] for feature in genome.features if feature.type == "gene"] 

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
    if end-start > 250 and strand == 1:
        positions.append([start,start+100, 1, ["downstream"]])
        positions.append([start+101, end-149, 1, [label]])
        positions.append([end-150, end, 1, ["upstream"]])
    elif end-start > 250:
        positions.append([start,start+150, -1, ["upstream"]])
        positions.append([start+151, end-99, -1, [label]])
        positions.append([end-100, end, -1, ["downstream"]])
    elif end-start == 250 and strand == 1:
        positions.append([start,start+100, 1, ["downstream"]])
        positions.append([end-150, end, 1, ["upstream"]])
    elif end-start == 250:
        positions.append([start,start+150, -1, ["upstream"]])
        positions.append([end-100, end, -1, ["downstream"]])
    else:
        positions.append([start,end, strand, ["upstream","downstream"]])
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
        positions.append([genic_nogenic_parts[-1][0],genic_nogenic_parts[-1][1], 1, ["downstream"]])
    else:
        positions.append([genic_nogenic_parts[-1][0],genic_nogenic_parts[-1][0]+100, 1, ["downstream"]])
        positions.append([genic_nogenic_parts[-1][0]+101,genic_nogenic_parts[-1][1], 1, [genic_nogenic_parts[0][3]]])
    return positions
                                                    
## Create variables of the sequences from each regions positions
def create_region_sequences(positions,genome):
    """
    Function to extract a unique sequence of nucleotids from correponding positions in each region 
    :params: positions - list to identify genic/intergenic/upstream/downstream parts in the genome
    :return: (genic_seq, intergenic_seq, upstream_seq, downstream_seq) - variables with the whole sequence from correponding parts unified from list of parts (genic, intergenic, upstream, downstream)
    """
    region_seq = {
    "genic": [],
    "intergenic": [],
    "upstream": [],
    "downstream": []    
    }
    for start, end, strand, label in positions:
        for l in label:
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
## Classify each peak
def classify_peak(peak_positions,positions):
    peak_positions = [[int(start),int(end)] for start,end in peak_positions.tolist()] #transform all values from peak list to int 
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
def pick_seq(peak_classification_aux, regions_seq):
    sublist = []
    for peak in peak_classification_aux:
        for par in peak:
            inicio_rdn = random.randint(0,len(regions_seq[par[0]])-par[1])
            sublist.append(regions_seq[par[0]][inicio_rdn:inicio_rdn+par[1]])
    return sublist


## Replicate peaks (negatives)
def peak_replication(peak_classification, regions_seq,config):
    n_replicates = config["Number of replicates of each peak"]
    replicates = []
    peak_classification_aux = []
    negatives = []
    for idx, (_, _, _, label, bp) in enumerate(peak_classification):
        peak_classification_aux.append([(label,peak_classification[idx][4][label].pop(0)) for label in peak_classification[idx][3]]) # We are selecting label with 3 and bp with 4 
    for i in range(n_replicates):
        replicates += pick_seq(peak_classification_aux, regions_seq)
    return replicates


## Extract nucleotics sequences from peaks data (positives)
def peak_seq(peak_positions, genome):
    peak_positions = [[int(start),int(end)] for start,end in peak_positions.tolist()] #transform all values from peak list to int
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
    genic_parts = select_genic_parts(genome)
    # Unify genic positions of the genome, purging those overlapped in the same strand
    new_genic_parts = unify_genic_parts(genic_parts)
    # Add no-genic positions to previous variable
    genic_nogenic_parts = add_intergenic(new_genic_parts,config["Genome length"])
    # Add upstream/downstream positions
    positions = add_pre_post_positions(genic_nogenic_parts)
    print('This is the head of the positions variable:\n',positions[0:4])
    # Create regions sequence dictionary
    regions_seq = create_region_sequences(positions,genome)
    # Peak positions to replicate
    peak_classification = classify_peak(peak_positions,positions)
    # Replicates n times
    replicates = peak_replication(peak_classification,regions_seq,config)
    # Positive nucleotids, positives
    positives = peak_seq(peak_positions, genome)
    # Save final data to work with
    save_csv(positives,replicates,config)
    print("Your output file has been saved in the following path:", config["Output folder"])
    print("You will find it with the name of:", config["Output file name"])
    print("All functions has been correctly executed! :)")


if __name__ == '__main__':
    data_extraction_peak_replica()