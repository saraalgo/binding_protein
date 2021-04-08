# QSAR Models with DNA
require(DNAshapeR)
require(seqinr)
require(RBioinf)
library(rDNAse)

# Load Data
# =====================================
## Load positive data
dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/src/modelling_design"
setwd(dirData)

inFile = '../../data/trial1/peaks_data_generated.txt'
outPositives = '../../data/trial1/peaks_positive.fa'
seqs = read.delim2(inFile, sep = '\t', header = F)
seqs = as.list(seqs$V1)

## Save fasta with positive cases
write.fasta(seqs, names = paste0('positive', seq(1:length(seqs))), file.out = outPositives)

## Generate negative data randomly
outNegatives = '../../data/trial1/peaks_negatives.fa'
len = list()
for (i in seq_along(seqs)) {
  len[[i]] = nchar(seqs[[i]])
}
md = median(unlist(len))

set.seed(1995)
rdm = list()
for (j in seq_along(seqs)) {
  rdm[[j]] = randDNA(md)
}

## Save fasta with negative cases
write.fasta(rdm, names = paste0('negatives', seq(1:length(rdm))), file.out = outNegatives)

## Unify positive + negative cases in a fasta file
outUnified = '../../data/trial1/peaks_unified.fa'
list_names <- c(paste0('positives', seq(1:length(seqs))),paste0('negatives', seq(1:length(rdm))))
write.fasta(c(seqs,rdm), names = list_names, file.out = outUnified)

# Generate descriptors
# ===================================
## Generation of shape features and save them in fasta files  
pos = getShape(outPositives)
neg = getShape(outNegatives)
unified = getShape(outUnified)

plotShape(pos$MGW)
library(fields)
heatShape(pos$ProT,16)

for (i in seq_along(faPos)) {
  print(nchar(faPos[[i]]))
}

## Generation of sequence and shape features combined
library(Biostrings)
library(rjson)
json <- fromJSON(file = '../../config/data.json')

#featureType <- c(sprintf("%d-mer",json$`Nucleotids in a row`), "1-shape")
featureType <- c("2-mer", "1-shape")
featureVector <- encodeSeqShape(outUnified, unified, featureType)
head(featureVector)

write.csv(sprintf("featureVector_%d",json$`Nucleotids in a row`),sprintf("../../data/trial1/featureVector_%d.csv",json$`Nucleotids in a row`), row.names = FALSE)

# Generation of sequence features
source('adn_descriptors.r')

faPos = readFASTA(outPositives)
faNeg = readFASTA(outNegatives)

### Save octamers (or what has been set up in the json file) features
number_mers <- sprintf("prob_%d",json$`Nucleotids in a row`)
pos = as.data.frame(t(sapply(unlist(faPos), function(x) get(number_mers)(x))))
neg = as.data.frame(t(sapply(unlist(faNeg), function(x) get(number_mers)(x))))
pos = data.frame(pos, target = 'peak')
neg = data.frame(neg, target = 'no-peak')
dataSequence = rbind(pos, neg)
dataSequence = dataSequence[, colSums(dataSequence != 0) > 0]
write.csv(dataSequence,sprintf("../../data/trial1/dataSequence_%d.csv",json$`Nucleotids in a row`), row.names = FALSE) 