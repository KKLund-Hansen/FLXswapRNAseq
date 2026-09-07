################################################################################################
################### FLX SWAP FEMALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNES ###################
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
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
IMFemaleFitF <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% FLXswapF.data$gene.id)
IMMaleFitF <- subset(IMMaleFit, IMMaleFit$gene.id %in% FLXswapF.data$gene.id)
IMAntaF <- subset(IMAnta, IMAnta$gene.id %in% FLXswapF.data$gene.id)
# 14527 transcripts left



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FvCfupsigFF <- IMFemaleFitF$gene.id %in% FvCfupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.07098094


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FvCfupsigMF <- IMMaleFitF$gene.id %in% FvCfupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.5739838


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FvCfupsigA <- IMAntaF$gene.id %in% FvCfupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfupsigA, IMAntaF$t, alternative = "either")
# P = 0.2326507



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FvCfdownsigFF <- IMFemaleFitF$gene.id %in% FvCfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.4550993


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FvCfdownsigMF <- IMMaleFitF$gene.id %in% FvCfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.4860495


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FvCfdownsigA <- IMAntaF$gene.id %in% FvCfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCfdownsigA, IMAntaF$t, alternative = "either")
# P = 0.9929588




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
AvSfupsigFF <- IMFemaleFitF$gene.id %in% AvSfupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.07216242


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
AvSfupsigMF <- IMMaleFitF$gene.id %in% AvSfupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.4861625


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
AvSfupsigA <- IMAntaF$gene.id %in% AvSfupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfupsigA, IMAntaF$t, alternative = "either")
# P = 0.07655991




##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
AvSfdownsigFF <- IMFemaleFitF$gene.id %in% AvSfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2021551


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
AvSfdownsigMF <- IMMaleFitF$gene.id %in% AvSfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.7881229


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
AvSfdownsigA <- IMAntaF$gene.id %in% AvSfdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSfdownsigA, IMAntaF$t, alternative = "either")
# P = 0.273919




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FsvCsFupsigFF <- IMFemaleFitF$gene.id %in% FsvCsFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.001359505

# Check if it is positive or negative associated
wilcoxGST(FsvCsFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.0006797526
wilcoxGST(FsvCsFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9993203
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FsvCsFupsigMF <- IMMaleFitF$gene.id %in% FsvCsFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.0006089157

# Check if it is positive or negative associated
wilcoxGST(FsvCsFupsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.9996956
wilcoxGST(FsvCsFupsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.0003044578
## Negative


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FsvCsFupsigA <- IMAntaF$gene.id %in% FsvCsFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFupsigA, IMAntaF$t, alternative = "either")
# P = 9.390236e-05

# Check if it is positive or negative associated
wilcoxGST(FsvCsFupsigA, IMAntaF$t, alternative = "up")
# P = 0.9999531
wilcoxGST(FsvCsFupsigA, IMAntaF$t, alternative = "down")
# P = 4.695118e-05
## Negative



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FsvCsFdownsigFF <- IMFemaleFitF$gene.id %in% FsvCsFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2589071


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FsvCsFdownsigMF <- IMMaleFitF$gene.id %in% FsvCsFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.3026093


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FsvCsFdownsigA <- IMAntaF$gene.id %in% FsvCsFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.1531668




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FavCaFupsigFF <- IMFemaleFitF$gene.id %in% FavCaFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 4.430738e-05

# Check if it is positive or negative associated
wilcoxGST(FavCaFupsigFF, IMFemaleFitF$t, alternative = "up")
# P = 0.9999778
wilcoxGST(FavCaFupsigFF, IMFemaleFitF$t, alternative = "down")
# P = 2.215369e-05
## Negative


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FavCaFupsigMF <- IMMaleFitF$gene.id %in% FavCaFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.06660911


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FavCaFupsigA <- IMAntaF$gene.id %in% FavCaFupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFupsigA, IMAntaF$t, alternative = "either")
# P = 0.000291904

# Check if it is positive or negative associated
wilcoxGST(FavCaFupsigA, IMAntaF$t, alternative = "up")
# P = 0.000145952
wilcoxGST(FavCaFupsigA, IMAntaF$t, alternative = "down")
# P = 0.9998541
## Positive



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FavCaFdownsigFF <- IMFemaleFitF$gene.id %in% FavCaFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.03106073


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FavCaFdownsigMF <- IMMaleFitF$gene.id %in% FavCaFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.2229372


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FavCaFdownsigA <- IMAntaF$gene.id %in% FavCaFdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaFdownsigA, IMAntaF$t, alternative = "either")
# P = 0.3007947

