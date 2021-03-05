import numpy as np
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


#Read peaks
peaks = []
with open("peaks.bed")as f:
    for line in f:
        peaks.append(line.strip().split())

print(len(peaks))
peaks = np.asarray(peaks)
print(peaks.shape)
print(peaks[0])

#Extract genome positions of the peaks
positions = pd.DataFrame(data=peaks[:,1:3], columns=["start", "end"])
