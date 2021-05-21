require(DNAshapeR)
require(dplyr)
source('D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/src/modelling_design/utils_descriptors.r')

# Generate Descriptors
dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/data/trial2"
setwd(dirData)

# DNAshapeR
# ===================================
shapes = getShape('data.fa')[1:4]
## Next step in Python script: data_preprocessing.py


# Nucleotide descriptors
# ===================================
reads = rDNAse::readFASTA('data.fa')

nt1 = as.data.frame(t(sapply(unlist(reads), function(x) prob_1(x))))
nt1 <- tibble::rownames_to_column(nt1, "Label")
nt2 = as.data.frame(t(sapply(unlist(reads), function(x) prob_2(x))))
nt2 <- tibble::rownames_to_column(nt2, "Label")
nt4 = as.data.frame(t(sapply(unlist(reads), function(x) prob_4(x))))
nt4 <- tibble::rownames_to_column(nt4, "Label")

dir.create(file.path(dirData, "descriptors"), showWarnings = FALSE)
setwd(file.path(dirData, "descriptors"))

write.csv(nt1,"1_descriptors.csv" , row.names = FALSE)
write.csv(nt2,"2_descriptors.csv" , row.names = FALSE)
write.csv(nt4,"4_descriptors.csv" , row.names = FALSE)








