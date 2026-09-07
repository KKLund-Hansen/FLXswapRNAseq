################################################################################################
########## CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES AGAINST PSEUDOOBSCURA #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexualConflict")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
ImmonenSC.data <- read.table(file = "ImmonenetalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)
DpGenes.data <- read.table(file = "DpseudoobscuraGenes.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add CG gene ID to the Immonene et al. data set, as adding the flybase name doesn't work.
ImmonenSC.data$gene.id <- DpGenes.data$CG.ID[match(ImmonenSC.data$ID, DpGenes.data$Glean.Accession)]

# Next remove the NAs
ImmonenSC.data <- ImmonenSC.data[!(is.na(ImmonenSC.data$gene.id)),]
# 1937 genes left

# To be able to match the two transcript list I have to add CG gene IDs to my data set
# First add the flybase gene names as row names
rownames(FLXswapM.data) <- FLXswapM.data[,1]
# The add a new column with the CG gene names
FLXswapM.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(FLXswapM.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonen et al. data set to match the male data set
ImmonenSCM <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% FLXswapM.data$CGgene.id)
# 1797 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropM <- length(ImmonenSCM$gene.id)/length(FLXswapM.data$CGgene.id)
# 0.12
# And for the chisq.test I also have to calculate
NotSCpropM <- 1-SCpropM
# 0.88



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many transcripts are related to sexual conflict?
FvCmupsigSC <- sum(FvCmupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 3
# How many are not?
FvCmupsigNotSC <- length(FvCmupsig$CGgene.id)-FvCmupsigSC
# 54

#The number of expected transcripts which are related to sexual conflict
length(FvCmupsig$CGgene.id)*SCpropM
# 6.606192

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCmupsigSC, FvCmupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.1357



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many transcripts are related to sexual conflict?
FvCmdownsigSC <- sum(FvCmdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
FvCmdownsigNotSC <- length(FvCmdownsig$CGgene.id)-FvCmdownsigSC
# 54

# The number of expected transcripts which are related to sexual conflict
length(FvCmdownsig$CGgene.id)*SCpropM
# 6.606192

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCmdownsigSC, FvCmdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.05665




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many transcripts are related to sexual conflict?
AvSmupsigSC <- sum(AvSmupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 22
# How many are not?
AvSmupsigNotSC <- length(AvSmupsig$CGgene.id)-AvSmupsigSC
# 193

#The number of expected transcripts which are related to sexual conflict
length(AvSmupsig$CGgene.id)*SCpropM
# 24.91809

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSmupsigSC, AvSmupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5341



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many transcripts are related to sexual conflict?
AvSmdownsigSC <- sum(AvSmdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 25
# How many are not?
AvSmdownsigNotSC <- length(AvSmdownsig$CGgene.id)-AvSmdownsigSC
# 237

# The number of expected transcripts which are related to sexual conflict
length(AvSmdownsig$CGgene.id)*SCpropM
# 30.3653

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSmdownsigSC, AvSmdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.3004




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many transcripts are related to sexual conflict?
FsvCsMupsigSC <- sum(FsvCsMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 7
# How many are not?
FsvCsMupsigNotSC <- length(FsvCsMupsig$CGgene.id)-FsvCsMupsigSC
# 49

#The number of expected transcripts which are related to sexual conflict
length(FsvCsMupsig$CGgene.id)*SCpropM
# 6.490293

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsMupsigSC, FsvCsMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.8315



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many transcripts are related to sexual conflict?
FsvCsMdownsigSC <- sum(FsvCsMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
FsvCsMdownsigNotSC <- length(FsvCsMdownsig$CGgene.id)-FsvCsMdownsigSC
# 47

# The number of expected transcripts which are related to sexual conflict
length(FsvCsMdownsig$CGgene.id)*SCpropM
# 5.679007

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsMdownsigSC, FsvCsMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.1006




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many transcripts are related to sexual conflict?
FavCaMupsigSC <- sum(FavCaMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 10
# How many are not?
FavCaMupsigNotSC <- length(FavCaMupsig$CGgene.id)-FavCaMupsigSC
# 55

#The number of expected transcripts which are related to sexual conflict
length(FavCaMupsig$CGgene.id)*SCpropM
# 7.533376

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaMupsigSC, FavCaMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.3392



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many transcripts are related to sexual conflict?
FavCaMdownsigSC <- sum(FavCaMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 23
# How many are not?
FavCaMdownsigNotSC <- length(FavCaMdownsig$CGgene.id)-FavCaMdownsigSC
# 98

# The number of expected transcripts which are related to sexual conflict
length(FavCaMdownsig$CGgene.id)*SCpropM
# 14.02367

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaMdownsigSC, FavCaMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.01079

