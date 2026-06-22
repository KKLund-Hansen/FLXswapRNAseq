################################################################################################
#################### FLX SWAP MALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNESS ###################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/AssociationWithFitness")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("limma")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(limma)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
IMFemaleFit <- read.table(file = "InnocentiMorrowFemaleFitness.csv", h = T, sep = ",", stringsAsFactors = T)
IMMaleFit <- read.table(file = "InnocentiMorrowMaleFitness.csv", h = T, sep = ",", stringsAsFactors = T)
IMAnta <- read.table(file = "InnocentiMorrowAntagonistic.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# To be able to compare and subset the data sets from Innocenti & Morrow (2010), I first have to change the probe names to flybase names
# First add the probe names as row names
rownames(IMFemaleFit) <- IMFemaleFit[,1]
rownames(IMMaleFit) <- IMMaleFit[,1]
rownames(IMAnta) <- IMAnta[,1]

# Next add a new column with the flybase names
IMFemaleFit$gene.id <- mapIds(drosophila2.db, keys = row.names(IMFemaleFit), keytype = "PROBEID", column = "FLYBASE")
IMMaleFit$gene.id <- mapIds(drosophila2.db, keys = row.names(IMMaleFit), keytype = "PROBEID", column = "FLYBASE")
IMAnta$gene.id <- mapIds(drosophila2.db, keys = row.names(IMAnta), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the Innocenti & Morrow data set to match the male data set
IMFemaleFitM <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% RNAswapmale.data$gene.id)
IMMaleFitM <- subset(IMMaleFit, IMMaleFit$gene.id %in% RNAswapmale.data$gene.id)
IMAntaM <- subset(IMAnta, IMAnta$gene.id %in% RNAswapmale.data$gene.id)
# 14820 transcripts left



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAMupsigFF <- IMFemaleFitM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.1928399


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAMupsigMF <- IMMaleFitM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.4235699


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAMupsigA <- IMAntaM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigA, IMAntaM$t, alternative = "either")
# P = 0.4311733



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAMdownsigFF <- IMFemaleFitM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.0511769


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAMdownsigMF <- IMMaleFitM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.02137938


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAMdownsigA <- IMAntaM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.01478839




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBMupsigFF <- IMFemaleFitM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.06008038


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBMupsigMF <- IMMaleFitM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.4382402


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBMupsigA <- IMAntaM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigA, IMAntaM$t, alternative = "either")
# P = 0.2217262



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBMdownsigFF <- IMFemaleFitM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.9381228


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBMdownsigMF <- IMMaleFitM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.3003662


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBMdownsigA <- IMAntaM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.8105251




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCMupsigFF <- IMFemaleFitM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.1450573


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCMupsigMF <- IMMaleFitM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.2212441


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCMupsigA <- IMAntaM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigA, IMAntaM$t, alternative = "either")
# P = 0.7524221



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCMdownsigFF <- IMFemaleFitM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.09568322


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCMdownsigMF <- IMMaleFitM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.6233503


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCMdownsigA <- IMAntaM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.04884793




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDMupsigFF <- IMFemaleFitM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.1575967


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDMupsigMF <- IMMaleFitM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.9024935


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDMupsigA <- IMAntaM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigA, IMAntaM$t, alternative = "either")
# P = 0.1372531



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDMdownsigFF <- IMFemaleFitM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.0001679148

# Check if it is positive or negative associated
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "up")
# P = 8.395739e-05
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "down")
# P = 0.9999161
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDMdownsigMF <- IMMaleFitM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.03237659


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDMdownsigA <- IMAntaM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.0006116527

# Check if it is positive or negative associated
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "up")
# P = 0.9996942
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "down")
# P = 0.0003058263
## Negative




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEMupsigFF <- IMFemaleFitM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.7453963


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEMupsigMF <- IMMaleFitM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.7193883


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEMupsigA <- IMAntaM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigA, IMAntaM$t, alternative = "either")
# P = 0.5651469



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEMdownsigFF <- IMFemaleFitM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.138858


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEMdownsigMF <- IMMaleFitM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.1745708


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEMdownsigA <- IMAntaM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.1508073




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFMupsigFF <- IMFemaleFitM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.7276298


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFMupsigMF <- IMMaleFitM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.6639021


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFMupsigA <- IMAntaM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigA, IMAntaM$t, alternative = "either")
# P = 0.5502054



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFMdownsigFF <- IMFemaleFitM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.9078132


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFMdownsigMF <- IMMaleFitM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.5460534


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFMdownsigA <- IMAntaM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.8778293




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGMupsigFF <- IMFemaleFitM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.7872844


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGMupsigMF <- IMMaleFitM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.7690058


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGMupsigA <- IMAntaM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigA, IMAntaM$t, alternative = "either")
# P = 0.5378963



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGMdownsigFF <- IMFemaleFitM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.00703493


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGMdownsigMF <- IMMaleFitM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.6273459


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGMdownsigA <- IMAntaM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.4537439




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G2MupsigFF <- IMFemaleFitM$gene.id %in% G2Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.005372521


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G2MupsigMF <- IMMaleFitM$gene.id %in% G2Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.8151338


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G2MupsigA <- IMAntaM$gene.id %in% G2Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MupsigA, IMAntaM$t, alternative = "either")
# P = 0.1213271



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G2MdownsigFF <- IMFemaleFitM$gene.id %in% G2Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.6627506


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G2MdownsigMF <- IMMaleFitM$gene.id %in% G2Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.4128212


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G2MdownsigA <- IMAntaM$gene.id %in% G2Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G2MdownsigA, IMAntaM$t, alternative = "either")
# P = 0.7748828




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G11MupsigFF <- IMFemaleFitM$gene.id %in% G11Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.300943


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G11MupsigMF <- IMMaleFitM$gene.id %in% G11Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.3244463


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G11MupsigA <- IMAntaM$gene.id %in% G11Mupsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MupsigA, IMAntaM$t, alternative = "either")
# P = 0.879507



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G11MdownsigFF <- IMFemaleFitM$gene.id %in% G11Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.5744368


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G11MdownsigMF <- IMMaleFitM$gene.id %in% G11Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.1983263


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G11MdownsigA <- IMAntaM$gene.id %in% G11Mdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G11MdownsigA, IMAntaM$t, alternative = "either")
# P = 0.1590871

