require(seqinr)

# Load Data
# =====================================
## Set directory
dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/data/trial2"
setwd(dirData)

## Load peaks/replicates data
inputPath = '../../data/trial2/data.csv'
data = read.delim2(inputPath, header = T, sep = ',')

## Filter replicates to ensure they are of the same lenght as peaks
# peaks = data[which(data$Label == 'Peak'), ]
# d = data[which(data$Label == 'Replicate'),]
# l = list()
# for (i in seq_along(d$Sequence)) {
#   l[[i]] = nchar(d$Sequence[i])
# }
# l = unlist(l)
# hist(l, 200)
# 
# pos = grep(290, l)
# d = d[pos, ]
# 
# data = rbind.data.frame(peaks, d)

# Save Data
# =====================================
## Save rds data
#saveRDS(data, file = 'data_annot_same_length.rds')
saveRDS(data, file = 'data_annot.rds')

## Save fasta data
reads = as.list(data$Sequence); names(reads) = data$Label
write.fasta(reads, names = names(reads), file.out = 'data.fa')



