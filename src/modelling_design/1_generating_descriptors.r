require(DNAshapeR)
require(dplyr)
source('D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/src/modelling_design/utils_descriptors.r')

# Generate Descriptors
dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/data/trial2"
setwd(dirData)
#saveRDS(data, file = 'data_annot_same_length.rds')
saveRDS(data, file = 'data_annot.rds')
shapes = getShape('data.fa')[1:4]

# DNAShape R
# ===================================
# MGW
MGW.h = data.frame()
for (i in seq_along(1:nrow(shapes$MGW))) {
  histogram = hist(shapes$MGW[i, ], 50)
  MGW.h = rbind(MGW.h, histogram$density)
}

# HelT
helt.h = data.frame()
for (i in seq_along(1:nrow(shapes$HelT))) {
  histogram = hist(shapes$HelT[i,], 50)
  helt.h = rbind(helt.h, histogram$density)
}

# ProT
prot.h = data.frame()
for (i in seq_along(1:nrow(shapes$ProT))) {
  histogram = hist(shapes$ProT[i,], 50)
  prot.h = rbind(prot.h, histogram$density)
}

# Roll
roll.h = data.frame()
for (i in seq_along(1:nrow(shapes$Roll))) {
  histogram = hist(shapes$Roll[i,], 50)
  roll.h = rbind(roll.h, histogram$density)
}


# Nucleotide descriptors
# ===================================
reads = rDNAse::readFASTA('data.fa')

nt1 = as.data.frame(t(sapply(unlist(reads), function(x) prob_1(x))))
nt2 = as.data.frame(t(sapply(unlist(reads), function(x) prob_2(x))))
nt4 = as.data.frame(t(sapply(unlist(reads), function(x) prob_4(x))))



# Add Y to datasets
# ==================================
labeling = function(d, annot, seed = 1995){

  stopifnot('Label' %in% names(annot))
  print('Label correct!')
  stopifnot(length(annot$Label) == nrow(d))
  print('Same Length!')

  res = cbind.data.frame(d, target = annot$Label)
  print('done!')

  print('Balancing dataset ...!')

  pos = res[which(res$target == 'Peak'),]
  neg = res[which(res$target == 'Replicate'),]

  neg = sample_n(neg, nrow(pos))

  data = rbind.data.frame(pos, neg)

  print(table(data$target))

  return(data)


}

dir.create(file.path(dirData, "descriptors"), showWarnings = FALSE)
setwd(file.path(dirData, "descriptors"))

MGW.h = labeling(MGW.h, annot)
helt.h = labeling(helt.h, annot)
roll.h = labeling(roll.h, annot)
prot.h = labeling(prot.h, annot)

saveRDS(list(MGW.h, helt.h, roll.h, prot.h), file = 'DNAShapeR_descriptors.rds')

nt1 = labeling(nt1, annot)
nt2 = labeling(nt2, annot)
nt4 = labeling(nt4, annot)

saveRDS(list(nt1, nt2, nt4), file = 'Kmers_descriptors.rds')








