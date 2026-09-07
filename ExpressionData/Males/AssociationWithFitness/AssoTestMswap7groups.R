################################################################################################
############## FLX SWAP MALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNESS - 7 GROUPS ##############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/AssociationWithFitness")

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
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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
IMFemaleFitM <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% FLXswap7gM.data$gene.id)
IMMaleFitM <- subset(IMMaleFit, IMMaleFit$gene.id %in% FLXswap7gM.data$gene.id)
IMAntaM <- subset(IMAnta, IMAnta$gene.id %in% FLXswap7gM.data$gene.id)
# 14527 transcripts left



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAMupsigFF <- IMFemaleFitM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.5549756


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAMupsigMF <- IMMaleFitM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.2366252


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAMupsigA <- IMAntaM$gene.id %in% GAMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMupsigA, IMAntaM$t, alternative = "either")
# P = 0.6276817



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAMdownsigFF <- IMFemaleFitM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.0488662


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAMdownsigMF <- IMMaleFitM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.02063879


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAMdownsigA <- IMAntaM$gene.id %in% GAMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.01417147




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBMupsigFF <- IMFemaleFitM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.0375347


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBMupsigMF <- IMMaleFitM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.3165787


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBMupsigA <- IMAntaM$gene.id %in% GBMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMupsigA, IMAntaM$t, alternative = "either")
# P = 0.1413982



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBMdownsigFF <- IMFemaleFitM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.9889031


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBMdownsigMF <- IMMaleFitM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.1568847


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBMdownsigA <- IMAntaM$gene.id %in% GBMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.4356018




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCMupsigFF <- IMFemaleFitM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.1381544


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCMupsigMF <- IMMaleFitM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.2313451


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCMupsigA <- IMAntaM$gene.id %in% GCMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMupsigA, IMAntaM$t, alternative = "either")
# P = 0.7303879



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCMdownsigFF <- IMFemaleFitM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.08847116


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCMdownsigMF <- IMMaleFitM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.5991248


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCMdownsigA <- IMAntaM$gene.id %in% GCMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.04503628




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDMupsigFF <- IMFemaleFitM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.07297216


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDMupsigMF <- IMMaleFitM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.9343955


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDMupsigA <- IMAntaM$gene.id %in% GDMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMupsigA, IMAntaM$t, alternative = "either")
# P = 0.1020662



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDMdownsigFF <- IMFemaleFitM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.0001567685

# Check if it is positive or negative associated
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "up")
# P = 7.838425e-05
wilcoxGST(GDMdownsigFF, IMFemaleFitM$t, alternative = "down")
# P = 0.9999216
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDMdownsigMF <- IMMaleFitM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.03155658


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDMdownsigA <- IMAntaM$gene.id %in% GDMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.0005819721

# Check if it is positive or negative associated
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "up")
# P = 0.9997091
wilcoxGST(GDMdownsigA, IMAntaM$t, alternative = "down")
# P = 0.0002909861
## Negative




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEMupsigFF <- IMFemaleFitM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.9271753


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEMupsigMF <- IMMaleFitM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.4556178


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEMupsigA <- IMAntaM$gene.id %in% GEMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMupsigA, IMAntaM$t, alternative = "either")
# P = 0.4423644



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEMdownsigFF <- IMFemaleFitM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.1335704


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEMdownsigMF <- IMMaleFitM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.2550671


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEMdownsigA <- IMAntaM$gene.id %in% GEMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.1904463




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFMupsigFF <- IMFemaleFitM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.7201297


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFMupsigMF <- IMMaleFitM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.6721884


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFMupsigA <- IMAntaM$gene.id %in% GFMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMupsigA, IMAntaM$t, alternative = "either")
# P = 0.5449601



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFMdownsigFF <- IMFemaleFitM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.9018726


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFMdownsigMF <- IMMaleFitM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.5386679


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFMdownsigA <- IMAntaM$gene.id %in% GFMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.8691493




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGMupsigFF <- IMFemaleFitM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.7999871


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGMupsigMF <- IMMaleFitM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.7801925


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGMupsigA <- IMAntaM$gene.id %in% GGMupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMupsigA, IMAntaM$t, alternative = "either")
# P = 0.5500717



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGMdownsigFF <- IMFemaleFitM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigFF, IMFemaleFitM$t, alternative = "either")
# P = 0.01062701


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGMdownsigMF <- IMMaleFitM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigMF, IMMaleFitM$t, alternative = "either")
# P = 0.4450067


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGMdownsigA <- IMAntaM$gene.id %in% GGMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGMdownsigA, IMAntaM$t, alternative = "either")
# P = 0.9357748

