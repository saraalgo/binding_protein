source('D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/src/modelling_design/utils_fs.r')

dirData = "D:/Users/Sara/Documents/2019-21 Máster Bioinformática/2021 Febrero/TFM/Programacion/Git/data/trial2/descriptors"
setwd(dirData)
DNAshape = readRDS('DNAShapeR_descriptors.rds')
kmers = readRDS('Kmers_descriptors.rds')

dir.create(file.path(dirData, "../FS"), showWarnings = FALSE)
setwd(file.path(dirData, "../FS"))

# Making PCA - 15 variables
# ==========================================
MGW.pca = pca(DNAshape[[1]])
helt.pca = pca(DNAshape[[2]])
roll.pca = pca(DNAshape[[3]])
prot.pca = pca(DNAshape[[4]])

nt1 = kmers[[1]]
nt2 = kmers[[2]]
nt4.pca = kmers[[3]]

write.csv(MGW.pca, file = 'MGW_pca.csv')
write.csv(helt.pca, file = 'HelT_pca.csv')
write.csv(roll.pca, file = 'Roll_pca.csv')
write.csv(prot.pca, file = 'ProT_pca.csv')

write.csv(nt1, file = 'nt1.csv')
write.csv(nt2, file = 'nt2.csv')
write.csv(nt4.pca, file = 'nt4_pca.csv')


# Making kruskal
# ==========================================

mgw.kruskal = kruskal(DNAshape[[1]], perc = 0.1)
helt.kruskal = kruskal(DNAshape[[2]], perc = 0.1)
roll.kruskal = kruskal(DNAshape[[3]], perc = 0.1)
prot.kruskal = kruskal(DNAshape[[4]], perc = 0.1)

write.csv(mgw.kruskal, file = 'MGW_kruskal.csv')
write.csv(helt.kruskal, file = 'HelT_kruskal.csv')
write.csv(roll.kruskal, file = 'Roll_kruskal.csv')
write.csv(prot.kruskal, file = 'ProT_kruskal.csv')


nt4_kruskal = kruskal(kmers[[3]], fs.type = 'kruskal.test', perc = 0.1)
write.csv(nt4_kruskal, file = 'nt4_kruskal.csv')


# Making FCBF
# ==========================================
mgw.fcbf = fast.cor.FS(DNAshape[[1]], thresh = 0.00025, normalize = F)
helt.fcbf = fast.cor.FS(DNAshape[[2]], thresh = 0.00025, normalize = F)
roll.fcbf = fast.cor.FS(DNAshape[[3]], thresh = 0.00025, normalize = F)
prot.fcbf = fast.cor.FS(DNAshape[[4]], thresh = 0.00025, normalize = F)

dnashape.fcbf = cbind.data.frame(subset(mgw.fcbf, select = -c(target)),
                                 subset(helt.fcbf, select = -c(target)),
                                 subset(roll.fcbf, select = -c(target)),
                                 prot.fcbf)
n1 = paste0(c('mgw_'), c(1:(ncol(mgw.fcbf)-1)))
n2 = paste0(c('helt_'), c(1:(ncol(helt.fcbf)-1)))
n3 = paste0(c('roll_'), c(1:(ncol(roll.fcbf)-1)))
n4 = paste0(c('prot_'), c(1:(ncol(prot.fcbf)-1)))
names(dnashape.fcbf) = c(n1, n2, n3, n4, 'target')
write.csv(dnashape.fcbf, file = 'dnashape_fcbf.csv')

nt4_fcbf = fast.cor.FS(kmers[[3]], thresh = 0.0025, normalize = F)
write.csv(nt4_fcbf, file = 'nt4_fcbf.csv')



