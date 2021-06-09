require(seqinr)
require(DNAshapeR)
require(dplyr)

# Load Data 
# =====================================
## Set directory
dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/data/trial2"
setwd(dirData)

## Load peaks/replicates data
inputPath = '../../data/trial2/data.csv'
data = read.delim2(inputPath, header = T, sep = ',')

# Save Data 
# =====================================
## Numerate peaks and replicates
reads = as.list(data$Sequence); names(reads) = data$Label
n = names(reads)
n1 = paste0(n[which(n == 'Peak')], c(1:length(which(n == 'Peak'))))
r = c(1:length(which(n == 'Peak')))
r = sort(rep(r, 10))
r = paste(r, c(1:10), sep = '.') #Because we have 10 replicates, change in function of what we want
rr = paste0(n[which(n == 'Replicate')], r)
names(reads) = c(n1, rr)

## Save fasta data
write.fasta(reads, names = names(reads), file.out = 'data.fa')

# DNAshapeR
# ===================================
shapes = getShape('data.fa')[1:4]
## Next step in Python script: data_preprocessing.py







# Load Data BOOTSTRAPING
# =====================================
## Load peaks/replicates data
inputPath = '../../data/trial2/data_bootstraping.csv'
data = read.delim2(inputPath, header = T, sep = ',')

# Save Data
# =====================================
## Numerate peaks and replicates
reads = as.list(data$Sequence); names(reads) = data$Label
n = names(reads)
n1 = paste0(n[which(n == 'Peak')], c(1:length(which(n == 'Peak'))))
r = c(1:length(which(n == 'Peak')))
r = sort(rep(r, 10))
r = paste(r, c(1:10), sep = '.') #Because we have 10 replicates, change in function of what we want
rr = paste0(n[which(n == 'Replicate')], r)
names(reads) = c(n1, rr)

## Save fasta data
write.fasta(reads, names = names(reads), file.out = 'data_bootstraping.fa')

# DNAshapeR
# ===================================
shapes = getShape('data_bootstrapping.fa')[1:4]
## Next step in Python script: data_preprocessing.py

