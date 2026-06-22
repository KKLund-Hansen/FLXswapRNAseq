################################################################################################
################### FLX SWAP FEMALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNESS ##################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Females/AssociationWithFitness")

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
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
IMFemaleFitF <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% RNAswapfemale.data$gene.id)
IMMaleFitF <- subset(IMMaleFit, IMMaleFit$gene.id %in% RNAswapfemale.data$gene.id)
IMAntaF <- subset(IMAnta, IMAnta$gene.id %in% RNAswapfemale.data$gene.id)
# 14820 transcripts left



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAFupsigFF <- IMFemaleFitF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.5368313


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAFupsigMF <- IMMaleFitF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.4423676


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAFupsigA <- IMAntaF$gene.id %in% GAFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFupsigA, IMAntaF$t, alternative = "either")
# P = 0.5861587



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GAFdownsigFF <- IMFemaleFitF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.8478834


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GAFdownsigMF <- IMMaleFitF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.002192959

# Check if it is positive or negative associated
wilcoxGST(GAFdownsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.00109648
wilcoxGST(GAFdownsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.9989035
## Positive


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GAFdownsigA <- IMAntaF$gene.id %in% GAFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GAFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.1678788




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBFupsigFF <- IMFemaleFitF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 3.520648e-09

# Check if it is positive or negative associated
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 1
wilcoxGST(GBFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 1.760324e-09
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBFupsigMF <- IMMaleFitF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.002677438

# Check if it is positive or negative associated
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.001338719
wilcoxGST(GBFupsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.9986613
## Positive


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBFupsigA <- IMAntaF$gene.id %in% GBFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "either")
# P = 4.920873e-07

# Check if it is positive or negative associated
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "up")
# P = 2.460436e-07
wilcoxGST(GBFupsigA, IMAntaF$t, alternative = "down")
# P = 0.9999998
## Positive



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GBFdownsigFF <- IMFemaleFitF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 5.441237e-05

# Check if it is positive or negative associated
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 2.720619e-05
wilcoxGST(GBFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9999728
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GBFdownsigMF <- IMMaleFitF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.02919391


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GBFdownsigA <- IMAntaF$gene.id %in% GBFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "either")
# P = 3.186103e-06

# Check if it is positive or negative associated
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "up")
# P = 0.9999984
wilcoxGST(GBFdownsigA, IMAntaF$t, alternative = "down")
# P = 1.593051e-06
## Negative




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCFupsigFF <- IMFemaleFitF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.0004295017

# Check if it is positive or negative associated
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9997853
wilcoxGST(GCFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.0002147508
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCFupsigMF <- IMMaleFitF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.0007961291

# Check if it is positive or negative associated
wilcoxGST(GCFupsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.0003980645
wilcoxGST(GCFupsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.999602
## Positive


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCFupsigA <- IMAntaF$gene.id %in% GCFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "either")
# P = 3.061508e-05

# Check if it is positive or negative associated
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "up")
# P = 1.530754e-05
wilcoxGST(GCFupsigA, IMAntaF$t, alternative = "down")
# P = 0.9999847
## Positive



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GCFdownsigFF <- IMFemaleFitF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.0956634


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GCFdownsigMF <- IMMaleFitF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.7804812


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GCFdownsigA <- IMAntaF$gene.id %in% GCFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GCFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.1527451




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDFupsigFF <- IMFemaleFitF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.0008629556

# Check if it is positive or negative associated
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9995685
wilcoxGST(GDFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.0004314778
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDFupsigMF <- IMMaleFitF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.2366257


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDFupsigA <- IMAntaF$gene.id %in% GDFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFupsigA, IMAntaF$t, alternative = "either")
# P = 0.0763436



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GDFdownsigFF <- IMFemaleFitF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.002351781

# Check if it is positive or negative associated
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.00117589
wilcoxGST(GDFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9988242
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GDFdownsigMF <- IMMaleFitF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.02233869


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GDFdownsigA <- IMAntaF$gene.id %in% GDFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.0001992509

# Check if it is positive or negative associated
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "up")
# P = 0.9999004
wilcoxGST(GDFdownsigA, IMAntaF$t, alternative = "down")
# P = 9.962543e-05
## Negative




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEFupsigFF <- IMFemaleFitF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.1039677


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEFupsigMF <- IMMaleFitF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.195275


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEFupsigA <- IMAntaF$gene.id %in% GEFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFupsigA, IMAntaF$t, alternative = "either")
# P = 0.9021351



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GEFdownsigFF <- IMFemaleFitF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.001666808

# Check if it is positive or negative associated
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9991666
wilcoxGST(GEFdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.0008334041
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GEFdownsigMF <- IMMaleFitF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.5874009


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GEFdownsigA <- IMAntaF$gene.id %in% GEFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GEFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.4608472




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFFupsigFF <- IMFemaleFitF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.3438662


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GFFupsigMF <- IMMaleFitF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.05938513


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GFFupsigA <- IMAntaF$gene.id %in% GFFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFupsigA, IMAntaF$t, alternative = "either")
# P = 0.8235733



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GFFdownsigFF <- IMFemaleFitF$gene.id %in% GFFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GFFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.1245142


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
# P = 0.5373169




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGFupsigFF <- IMFemaleFitF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2419909


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGFupsigMF <- IMMaleFitF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.07568146


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGFupsigA <- IMAntaF$gene.id %in% GGFupsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFupsigA, IMAntaF$t, alternative = "either")
# P = 0.1711562



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group F
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
GGFdownsigFF <- IMFemaleFitF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.7893756


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
GGFdownsigMF <- IMMaleFitF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.9443706


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
GGFdownsigA <- IMAntaF$gene.id %in% GGFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(GGFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.9461502




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G23FupsigFF <- IMFemaleFitF$gene.id %in% G23Fupsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.7295095


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G23FupsigMF <- IMMaleFitF$gene.id %in% G23Fupsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.3908414


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G23FupsigA <- IMAntaF$gene.id %in% G23Fupsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FupsigA, IMAntaF$t, alternative = "either")
# P = 0.6030832



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
G23FdownsigFF <- IMFemaleFitF$gene.id %in% G23Fdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2640873


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
G23FdownsigMF <- IMMaleFitF$gene.id %in% G23Fdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.3907626


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
G23FdownsigA <- IMAntaF$gene.id %in% G23Fdownsig$gene.id

# Next do a wilcox test
wilcoxGST(G23FdownsigA, IMAntaF$t, alternative = "either")
# P = 0.2208549

