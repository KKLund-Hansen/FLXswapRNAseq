################################################################################################
#################### FLX SWAP MALE ASSOCIATION TEST FOR SEX-SPECIFIC FITNESS ###################
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
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
IMFemaleFitF <- subset(IMFemaleFit, IMFemaleFit$gene.id %in% FLXswapM.data$gene.id)
IMMaleFitF <- subset(IMMaleFit, IMMaleFit$gene.id %in% FLXswapM.data$gene.id)
IMAntaF <- subset(IMAnta, IMAnta$gene.id %in% FLXswapM.data$gene.id)
# 14527 transcripts left



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FvCmupsigFF <- IMFemaleFitF$gene.id %in% FvCmupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.7810446


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FvCmupsigMF <- IMMaleFitF$gene.id %in% FvCmupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.3505292


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FvCmupsigA <- IMAntaF$gene.id %in% FvCmupsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmupsigA, IMAntaF$t, alternative = "either")
# P = 0.2094304



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FvCmdownsigFF <- IMFemaleFitF$gene.id %in% FvCmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.04476385


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FvCmdownsigMF <- IMMaleFitF$gene.id %in% FvCmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.8195455


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FvCmdownsigA <- IMAntaF$gene.id %in% FvCmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FvCmdownsigA, IMAntaF$t, alternative = "either")
# P = 0.3165164



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
AvSmupsigFF <- IMFemaleFitF$gene.id %in% AvSmupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.07376327


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
AvSmupsigMF <- IMMaleFitF$gene.id %in% AvSmupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.8246837


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
AvSmupsigA <- IMAntaF$gene.id %in% AvSmupsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmupsigA, IMAntaF$t, alternative = "either")
# P = 0.4284283



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
AvSmdownsigFF <- IMFemaleFitF$gene.id %in% AvSmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.00570763


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
AvSmdownsigMF <- IMMaleFitF$gene.id %in% AvSmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.003192757

# Check if it is positive or negative associated
wilcoxGST(AvSmdownsigMF, IMMaleFitF$t, alternative = "up")
# P = 0.9984037
wilcoxGST(AvSmdownsigMF, IMMaleFitF$t, alternative = "down")
# P = 0.001596379
## Negative


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
AvSmdownsigA <- IMAntaF$gene.id %in% AvSmdownsig$gene.id

# Next do a wilcox test
wilcoxGST(AvSmdownsigA, IMAntaF$t, alternative = "either")
# P = 0.0006799345

# Check if it is positive or negative associated
wilcoxGST(AvSmdownsigA, IMAntaF$t, alternative = "up")
# P = 0.9996601
wilcoxGST(AvSmdownsigA, IMAntaF$t, alternative = "down")
# P = 0.0003399673
## Negative




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FsvCsMupsigFF <- IMFemaleFitF$gene.id %in% FsvCsMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.2623722


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FsvCsMupsigMF <- IMMaleFitF$gene.id %in% FsvCsMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.915865


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FsvCsMupsigA <- IMAntaF$gene.id %in% FsvCsMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMupsigA, IMAntaF$t, alternative = "either")
# P = 0.5992161



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FsvCsMdownsigFF <- IMFemaleFitF$gene.id %in% FsvCsMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.7538138


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FsvCsMdownsigMF <- IMMaleFitF$gene.id %in% FsvCsMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.660664


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FsvCsMdownsigA <- IMAntaF$gene.id %in% FsvCsMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FsvCsMdownsigA, IMAntaF$t, alternative = "either")
# P = 0.6708729




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FavCaMupsigFF <- IMFemaleFitF$gene.id %in% FavCaMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMupsigFF, IMFemaleFitF$t, alternative = "either")
# P = 0.1522607


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FavCaMupsigMF <- IMMaleFitF$gene.id %in% FavCaMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMupsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.8139491


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FavCaMupsigA <- IMAntaF$gene.id %in% FavCaMupsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMupsigA, IMAntaF$t, alternative = "either")
# P = 0.076287



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to do the association test for down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121


## Test for association with female fitness

# Make an logical vector for the transcripts which are significant for female fitness
FavCaMdownsigFF <- IMFemaleFitF$gene.id %in% FavCaMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMdownsigFF, IMFemaleFitF$t, alternative = "either")
# P = 2.455171e-06

# Check if it is positive or negative associated
wilcoxGST(FavCaMdownsigFF, IMFemaleFitF$t, alternative = "up")
# P = 1.227586e-06
wilcoxGST(FavCaMdownsigFF, IMFemaleFitF$t, alternative = "down")
# P = 0.9999988
## Positive


## Test for association with male fitness

# Make an logical vector for the transcripts which are significant for male fitness
FavCaMdownsigMF <- IMMaleFitF$gene.id %in% FavCaMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMdownsigMF, IMMaleFitF$t, alternative = "either")
# P = 0.01971592


## Test for association with antagonistic fitness

# Make an logical vector for the transcripts which are significant for antagonistic fitness
FavCaMdownsigA <- IMAntaF$gene.id %in% FavCaMdownsig$gene.id

# Next do a wilcox test
wilcoxGST(FavCaMdownsigA, IMAntaF$t, alternative = "either")
# P = 6.327688e-06

# Check if it is positive or negative associated
wilcoxGST(FavCaMdownsigA, IMAntaF$t, alternative = "up")
# P = 0.9999968
wilcoxGST(FavCaMdownsigA, IMAntaF$t, alternative = "down")
# P = 3.163844e-06
## Negative

