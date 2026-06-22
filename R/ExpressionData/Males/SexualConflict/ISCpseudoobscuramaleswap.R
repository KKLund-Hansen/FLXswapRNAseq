################################################################################################
########## CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES AGAINST PSEUDOOBSCURA #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/SexualConflict")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
rownames(RNAswapmale.data) <- RNAswapmale.data[,1]
# The add a new column with the CG gene names
RNAswapmale.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(RNAswapmale.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonene et al. data set to match the male data set
ImmonenSCM <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% RNAswapmale.data$CGgene.id)
# 1808 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropM <- length(ImmonenSCM$gene.id)/length(RNAswapmale.data$CGgene.id)
# 0.11
# And for the chisq.test I also have to calculate
NotSCpropM <- 1-SCpropM
# 0.89



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many transcripts are related to sexual conflict?
GAMupsigSC <- sum(GAMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GAMupsigNotSC <- length(GAMupsig$CGgene.id)-GAMupsigSC
# 4

# The number of expected transcripts which are related to sexual conflict
length(GAMupsig$CGgene.id)*SCpropM
# 0.4308609

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMupsigSC, GAMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.4871

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GAMupF <- rbind(c(GAMupsigSC, GAMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GAMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GAMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GAMupF) <- c("SA", "NotSA")
rownames(SC.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GAMdownsigSC <- sum(GAMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 1
# How many are not?
GAMdownsigNotSC <- length(GAMdownsig$CGgene.id)-GAMdownsigSC
# 15

# The number of expected transcripts which are related to sexual conflict
length(GAMdownsig$CGgene.id)*SCpropM
# 1.723444

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigSC, GAMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5596

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GAMdownF <- rbind(c(GAMdownsigSC, GAMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GAMdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GAMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GAMdownF) <- c("SA", "NotSA")
rownames(SC.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GAMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many transcripts are related to sexual conflict?
GBMupsigSC <- sum(GBMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 4
# How many are not?
GBMupsigNotSC <- length(GBMupsig$CGgene.id)-GBMupsigSC
# 50

# The number of expected transcripts which are related to sexual conflict
length(GBMupsig$CGgene.id)*SCpropM
# 5.816622

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMupsigSC, GBMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.4252



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many transcripts are related to sexual conflict?
GBMdownsigSC <- sum(GBMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 10
# How many are not?
GBMdownsigNotSC <- length(GBMdownsig$CGgene.id)-GBMdownsigSC
# 48

# The number of expected transcripts which are related to sexual conflict
length(GBMdownsig$CGgene.id)*SCpropM
# 6.247483

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigSC, GBMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.112




####################### GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many transcripts are related to sexual conflict?
GCMupsigSC <- sum(GCMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 5
# How many are not?
GCMupsigNotSC <- length(GCMupsig$CGgene.id)-GCMupsigSC
# 22

# The number of expected transcripts which are related to sexual conflict
length(GCMupsig$CGgene.id)*SCpropM
# 2.908311

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMupsigSC, GCMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.1941

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GCMupF <- rbind(c(GCMupsigSC, GCMupsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GCMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GCMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GCMupF) <- c("SA", "NotSA")
rownames(SC.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GCMupF)
# P = 0.2048

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many transcripts are related to sexual conflict?
GCMdownsigSC <- sum(GCMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 6
# How many are not?
GCMdownsigNotSC <- length(GCMdownsig$CGgene.id)-GCMdownsigSC
# 31

# The number of expected transcripts which are related to sexual conflict
length(GCMdownsig$CGgene.id)*SCpropM
# 3.985463

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigSC, GCMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.2854

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GCMdownF <- rbind(c(GCMdownsigSC, GCMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GCMdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GCMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GCMdownF) <- c("SA", "NotSA")
rownames(SC.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GCMdownF)
# P = 0.2833

# We get the same results, so that is good




####################### GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many transcripts are related to sexual conflict?
GDMupsigSC <- sum(GDMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 1
# How many are not?
GDMupsigNotSC <- length(GDMupsig$CGgene.id)-GDMupsigSC
# 13

# The number of expected transcripts which are related to sexual conflict
length(GDMupsig$CGgene.id)*SCpropM
# 1.508013

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMupsigSC, GDMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.6614

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GDMupF <- rbind(c(GDMupsigSC, GDMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GDMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GDMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GDMupF) <- c("SA", "NotSA")
rownames(SC.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GDMdownsigSC <- sum(GDMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
GDMdownsigNotSC <- length(GDMdownsig$CGgene.id)-GDMdownsigSC
# 14

# The number of expected transcripts which are related to sexual conflict
length(GDMdownsig$CGgene.id)*SCpropM
# 1.723444

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigSC, GDMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.8235

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GDMdownF <- rbind(c(GDMdownsigSC, GDMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GDMdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GDMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GDMdownF) <- c("SA", "NotSA")
rownames(SC.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GDMdownF)
# P = 0.6881

# We get the same results, so that is good




####################### GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many transcripts are related to sexual conflict?
GEMupsigSC <- sum(GEMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 1
# How many are not?
GEMupsigNotSC <- length(GEMupsig$CGgene.id)-GEMupsigSC
# 16

# The number of expected transcripts which are related to sexual conflict
length(GEMupsig$CGgene.id)*SCpropM
# 1.831159

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMupsigSC, GEMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5155

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GEMupF <- rbind(c(GEMupsigSC, GEMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GEMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GEMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GEMupF) <- c("SA", "NotSA")
rownames(SC.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GEMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many transcripts are related to sexual conflict?
GEMdownsigSC <- sum(GEMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 5
# How many are not?
GEMdownsigNotSC <- length(GEMdownsig$CGgene.id)-GEMdownsigSC
# 21

# The number of expected transcripts which are related to sexual conflict
length(GEMdownsig$CGgene.id)*SCpropM
# 2.800596

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigSC, GEMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.1641

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GEMdownF <- rbind(c(GEMdownsigSC, GEMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GEMdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GEMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GEMdownF) <- c("SA", "NotSA")
rownames(SC.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GEMdownF)
# P = 0.1927

# We get the same results, so that is good




####################### GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many transcripts are related to sexual conflict?
GFMupsigSC <- sum(GFMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GFMupsigNotSC <- length(GFMupsig$CGgene.id)-GFMupsigSC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMupsig$CGgene.id)*SCpropM
# 0.2154304

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMupsigSC, GFMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.6232

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFMupF <- rbind(c(GFMupsigSC, GFMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GFMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GFMupsigNotSC))
#Add names to the columns and rows
colnames(SC.GFMupF) <- c("SC", "NotSC")
rownames(SC.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many transcripts are related to sexual conflict?
GFMdownsigSC <- sum(GFMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GFMdownsigNotSC <- length(GFMdownsig$CGgene.id)-GFMdownsigSC
# 3

# The number of expected transcripts which are related to sexual conflict
length(GFMdownsig$CGgene.id)*SCpropM
# 0.3231457

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigSC, GFMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5473

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFMupF <- rbind(c(GFMupsigSC, GFMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GFMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GFMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GFMupF) <- c("SC", "NotSC")
rownames(SC.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GFMupF)
# P = 1

# We get the same results, so that is goo




####################### GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many transcripts are related to sexual conflict?
GGMupsigSC <- sum(GGMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GGMupsigNotSC <- length(GGMupsig$CGgene.id)-GGMupsigSC
# 10

# The number of expected transcripts which are related to sexual conflict
length(GGMupsig$CGgene.id)*SCpropM
# 1.077152

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMupsigSC, GGMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.2719

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GGMupF <- rbind(c(GGMupsigSC, GGMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GGMupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GGMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GGMupF) <- c("SC", "NotSC")
rownames(SC.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GGMupF)
# P = 0.6137

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many transcripts are related to sexual conflict?
GGMdownsigSC <- sum(GGMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
GGMdownsigNotSC <- length(GGMdownsig$CGgene.id)-GGMdownsigSC
# 17

# The number of expected transcripts which are related to sexual conflict
length(GGMdownsig$CGgene.id)*SCpropM
# 2.046589

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigSC, GGMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.9207

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GGMdownF <- rbind(c(GGMdownsigSC, GGMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GGMdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-GGMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GGMdownF) <- c("SC", "NotSC")
rownames(SC.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GGMdownF)
# P = 1

# We get the same results, so that is good




####################### GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many transcripts are related to sexual conflict?
G2MupsigSC <- sum(G2Mupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 4
# How many are not?
G2MupsigNotSC <- length(G2Mupsig$CGgene.id)-G2MupsigSC
# 25

# The number of expected transcripts which are related to sexual conflict
length(G2Mupsig$CGgene.id)*SCpropM
# 3.123741

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G2MupsigSC, G2MupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5997

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G2MupF <- rbind(c(G2MupsigSC, G2MupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-G2MupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-G2MupsigNotSC))
# Add names to the columns and rows
colnames(SC.G2MupF) <- c("SC", "NotSC")
rownames(SC.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G2MupF)
# P = 0.5468

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many transcripts are related to sexual conflict?
G2MdownsigSC <- sum(G2Mdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 5
# How many are not?
G2MdownsigNotSC <- length(G2Mdownsig$CGgene.id)-G2MdownsigSC
# 31

# The number of expected transcripts which are related to sexual conflict
length(G2Mdownsig$CGgene.id)*SCpropM
# 3.877748

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G2MdownsigSC, G2MdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5463

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G2MdownF <- rbind(c(G2MdownsigSC, G2MdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-G2MdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-G2MdownsigNotSC))
# Add names to the columns and rows
colnames(SC.G2MdownF) <- c("SC", "NotSC")
rownames(SC.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G2MdownF)
# P = 0.5854

# We get the same results, so that is good




####################### GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many transcripts are related to sexual conflict?
G11MupsigSC <- sum(G11Mupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 7
# How many are not?
G11MupsigNotSC <- length(G11Mupsig$CGgene.id)-G11MupsigSC
# 14

# The number of expected transcripts which are related to sexual conflict
length(G11Mupsig$CGgene.id)*SCpropM
# 2.26202

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G11MupsigSC, G11MupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.0008531

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G11MupF <- rbind(c(G11MupsigSC, G11MupsigNotSC),
                    c(length(ImmonenSCM$gene.id)-G11MupsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-G11MupsigNotSC))
# Add names to the columns and rows
colnames(SC.G11MupF) <- c("SC", "NotSC")
rownames(SC.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G11MupF)
# P = 0.004946

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many transcripts are related to sexual conflict?
G11MdownsigSC <- sum(G11Mdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 4
# How many are not?
G11MdownsigNotSC <- length(G11Mdownsig$CGgene.id)-G11MdownsigSC
# 3

# The number of expected transcripts which are related to sexual conflict
length(G11Mdownsig$CGgene.id)*SCpropM
# 0.7540066

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G11MdownsigSC, G11MdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 7.577e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G11MdownF <- rbind(c(G11MdownsigSC, G11MdownsigNotSC),
                      c(length(ImmonenSCM$gene.id)-G11MdownsigSC, (length(RNAswapmale.data$CGgene.id)-length(ImmonenSCM$gene.id))-G11MdownsigNotSC))
# Add names to the columns and rows
colnames(SC.G11MdownF) <- c("SC", "NotSC")
rownames(SC.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G11MdownF)
# P = 0.003591

# It goes from significant to non-significant
