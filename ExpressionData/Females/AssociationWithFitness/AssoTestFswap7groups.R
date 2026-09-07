################################################################################################
############# FLX SWAP FEMALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNESS - 7 GROUPS #############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/AssociationWithFitness")

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
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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

# Now I can subset the Innocenti & Morrow data set to match the female data set
IMFemaleFitF <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% FLXswap7gF.data$gene.id)
IMMaleFitF <- subset(IMMaleFit, IMMaleFit$gene.id %in% FLXswap7gF.data$gene.id)
IMAntaF <- subset(IMAnta, IMAnta$gene.id %in% FLXswap7gF.data$gene.id)
# 14527 transcripts left



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAFupsigFF <- IMFemaleFitF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.3159857


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAFupsigMF <- IMMaleFitF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.6029277


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAFupsigA <- IMAntaF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigA, IMAntaF$t, alternative = "either")
# P = 0.9303562



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAFdownsigFF <- IMFemaleFitF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.7462218


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAFdownsigMF <- IMMaleFitF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.003595479


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAFdownsigA <- IMAntaF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.2224542




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBFupsigFF <- IMFemaleFitF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 1.29158e-07

# Check if it is positive or negative associated
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9999999
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 6.457899e-08
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBFupsigMF <- IMMaleFitF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.001668793

# Check if it is positive or negative associated
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.0008343963
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.9991656
## Positive


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBFupsigA <- IMAntaF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "either")
# P = 4.076918e-06

# Check if it is positive or negative associated
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "up")
# P = 2.038459e-06
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "down")
# P = 0.999998
## Positive



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBFdownsigFF <- IMFemaleFitF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.000437776

# Check if it is positive or negative associated
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.000218888
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9997811
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBFdownsigMF <- IMMaleFitF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.01507385


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBFdownsigA <- IMAntaF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "either")
# P = 4.080176e-06

# Check if it is positive or negative associated
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "up")
# P = 0.999998
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "down")
# P = 2.040088e-06
## Negative




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCFupsigFF <- IMFemaleFitF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.0003899225

# Check if it is positive or negative associated
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9998051
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.0001949613
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCFupsigMF <- IMMaleFitF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.007822037


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCFupsigA <- IMAntaF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "either")
# P = 0.0001623454

# Check if it is positive or negative associated
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "up")
# P = 8.117269e-05
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "down")
# P = 0.9999188
## Positive



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCFdownsigFF <- IMFemaleFitF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.02952919


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCFdownsigMF <- IMMaleFitF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.748878


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCFdownsigA <- IMAntaF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.1301765




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDFupsigFF <- IMFemaleFitF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.0007564393

# Check if it is positive or negative associated
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9996218
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.0003782197
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDFupsigMF <- IMMaleFitF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.402238


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDFupsigA <- IMAntaF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigA, IMAntaF$t, alternative = "either")
# P = 0.05858759



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDFdownsigFF <- IMFemaleFitF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.001119372

# Check if it is positive or negative associated
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.0005596862
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9994403
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDFdownsigMF <- IMMaleFitF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.006971086


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDFdownsigA <- IMAntaF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "either")
# P = 5.041192e-05

# Check if it is positive or negative associated
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "up")
# P = 0.9999748
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "down")
# P = 2.520596e-05
## Negative




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEFupsigFF <- IMFemaleFitF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2471253


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEFupsigMF <- IMMaleFitF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.1100699


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEFupsigA <- IMAntaF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigA, IMAntaF$t, alternative = "either")
# P = 0.5397519



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEFdownsigFF <- IMFemaleFitF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.002446674

# Check if it is positive or negative associated
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9987767
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.001223337
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEFdownsigMF <- IMMaleFitF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.5220937


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEFdownsigA <- IMAntaF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.5417134




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFFupsigFF <- IMFemaleFitF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2578865


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFFupsigMF <- IMMaleFitF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.2461766


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFFupsigA <- IMAntaF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigA, IMAntaF$t, alternative = "either")
# P = 0.7082864



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFFdownsigFF <- IMFemaleFitF$gene.id %in% GFFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.3115428


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFFdownsigMF <- IMMaleFitF$gene.id %in% GFFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.8764864


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFFdownsigA <- IMAntaF$gene.id %in% GFFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.6359505




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGFupsigFF <- IMFemaleFitF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.03759636


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGFupsigMF <- IMMaleFitF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.03322205


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGFupsigA <- IMAntaF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigA, IMAntaF$t, alternative = "either")
# P = 0.04726396



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGFdownsigFF <- IMFemaleFitF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.8270173


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGFdownsigMF <- IMMaleFitF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.9975742


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGFdownsigA <- IMAntaF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.9846514

