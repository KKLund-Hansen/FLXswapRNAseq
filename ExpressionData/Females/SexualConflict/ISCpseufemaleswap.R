################################################################################################
######### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES AGAINST PSEUDOOBSCURA ########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexualConflict")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
rownames(FLXswapF.data) <- FLXswapF.data[,1]
# The add a new column with the CG gene names
FLXswapF.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(FLXswapF.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonene et al. data set to match the female data set
ImmonenSCF <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% FLXswapF.data$CGgene.id)
# 1797 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropF <- length(ImmonenSCF$gene.id)/length(FLXswapF.data$CGgene.id)
# 0.12
# And for the chisq.test I also have to calculate
NotSCpropF <- 1-SCpropF
# 0.88



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many transcripts are related to sexual conflict?
FvCfupsigSC <- sum(FvCfupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 1
# How many are not?
FvCfupsigNotSC <- length(FvCfupsig$CGgene.id)-FvCfupsigSC
# 34

#The number of expected transcripts which are related to sexual conflict
length(FvCfupsig$CGgene.id)*SCpropF
# 4.056433

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCfupsigSC, FvCfupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1065

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.FvCfupF <- rbind(c(FvCfupsigSC, FvCfupsigNotSC),
                    c(length(ImmonenSCF$gene.id)-FvCfupsigSC, (length(FLXswapF.data$CGgene.id)-length(ImmonenSCF$gene.id))-FvCfupsigNotSC))
#Add names to the columns and rows
colnames(SC.FvCfupF) <- c("SC", "NotSC")
rownames(SC.FvCfupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.FvCfupF)
# P = 0.178

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many transcripts are related to sexual conflict?
FvCfdownsigSC <- sum(FvCfdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 5
# How many are not?
FvCfdownsigNotSC <- length(FvCfdownsig$CGgene.id)-FvCfdownsigSC
# 35

# The number of expected transcripts which are related to sexual conflict
length(FvCfdownsig$CGgene.id)*SCpropF
# 4.635924

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCfdownsigSC, FvCfdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.8573

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.FvCfdownF <- rbind(c(FvCfdownsigSC, FvCfdownsigNotSC),
                      c(length(ImmonenSCF$gene.id)-FvCfdownsigSC, (length(FLXswapF.data$CGgene.id)-length(ImmonenSCF$gene.id))-FvCfdownsigNotSC))
#Add names to the columns and rows
colnames(SC.FvCfdownF) <- c("SC", "NotSC")
rownames(SC.FvCfdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.FvCfdownF)
# P = 0.8042

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many transcripts are related to sexual conflict?
AvSfupsigSC <- sum(AvSfupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 25
# How many are not?
AvSfupsigNotSC <- length(AvSfupsig$CGgene.id)-AvSfupsigSC
# 219

#The number of expected transcripts which are related to sexual conflict
length(AvSfupsig$CGgene.id)*SCpropF
# 17.88073

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSfupsigSC, AvSfupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.5119



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many transcripts are related to sexual conflict?
AvSfdownsigSC <- sum(AvSfdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 17
# How many are not?
AvSfdownsigNotSC <- length(AvSfdownsig$CGgene.id)-AvSfdownsigSC
# 144

# The number of expected transcripts which are related to sexual conflict
length(AvSfdownsig$CGgene.id)*SCpropF
# 18.65959

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSfdownsigSC, AvSfdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.6828




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many transcripts are related to sexual conflict?
FsvCsFupsigSC <- sum(FsvCsFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 3
# How many are not?
FsvCsFupsigNotSC <- length(FsvCsFupsig$CGgene.id)-FsvCsFupsigSC
# 134

#The number of expected transcripts which are related to sexual conflict
length(FsvCsFupsig$CGgene.id)*SCpropF
# 15.87804

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsFupsigSC, FsvCsFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.0005878



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many transcripts are related to sexual conflict?
FsvCsFdownsigSC <- sum(FsvCsFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 8
# How many are not?
FsvCsFdownsigNotSC <- length(FsvCsFdownsig$CGgene.id)-FsvCsFdownsigSC
# 47

# The number of expected transcripts which are related to sexual conflict
length(FsvCsFdownsig$CGgene.id)*SCpropF
# 6.374395

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsFdownsigSC, FsvCsFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.4935




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many transcripts are related to sexual conflict?
FavCaFupsigSC <- sum(FavCaFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 86
# How many are not?
FavCaFupsigNotSC <- length(FavCaFupsig$CGgene.id)-FavCaFupsigSC
# 439

#The number of expected transcripts which are related to sexual conflict
length(FavCaFupsig$CGgene.id)*SCpropF
# 60.8465

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaFupsigSC, FavCaFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.0006047



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many transcripts are related to sexual conflict?
FavCaFdownsigSC <- sum(FavCaFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 408
# How many are not?
FavCaFdownsigNotSC <- length(FavCaFdownsig$CGgene.id)-FavCaFdownsigSC
# 2586

# The number of expected transcripts which are related to sexual conflict
length(FavCaFdownsig$CGgene.id)*SCpropF
# 346.9989

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaFdownsigSC, FavCaFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.0004963

