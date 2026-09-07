################################################################################################
#### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES AGAINST PSEUDOOBSCURA - 7 GROUPS ####
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
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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
rownames(FLXswap7gM.data) <- FLXswap7gM.data[,1]
# The add a new column with the CG gene names
FLXswap7gM.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(FLXswap7gM.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonene et al. data set to match the male data set
ImmonenSCM <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% FLXswap7gM.data$CGgene.id)
# 1797 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropM <- length(ImmonenSCM$gene.id)/length(FLXswap7gM.data$CGgene.id)
# 0.12
# And for the chisq.test I also have to calculate
NotSCpropM <- 1-SCpropM
# 0.88



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many transcripts are related to sexual conflict?
GAMupsigSC <- sum(GAMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GAMupsigNotSC <- length(GAMupsig$CGgene.id)-GAMupsigSC
# 5

# The number of expected transcripts which are related to sexual conflict
length(GAMupsig$CGgene.id)*SCpropM
# 0.5794905

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMupsigSC, GAMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.4182

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GAMupF <- rbind(c(GAMupsigSC, GAMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GAMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GAMupsigNotSC))
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
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GAMdownsigSC <- sum(GAMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 1
# How many are not?
GAMdownsigNotSC <- length(GAMdownsig$CGgene.id)-GAMdownsigSC
# 15

# The number of expected transcripts which are related to sexual conflict
length(GAMdownsig$CGgene.id)*SCpropM
# 1.85437

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigSC, GAMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5046

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GAMdownF <- rbind(c(GAMdownsigSC, GAMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GAMdownsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GAMdownsigNotSC))
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
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many transcripts are related to sexual conflict?
GBMupsigSC <- sum(GBMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 4
# How many are not?
GBMupsigNotSC <- length(GBMupsig$CGgene.id)-GBMupsigSC
# 47

# The number of expected transcripts which are related to sexual conflict
length(GBMupsig$CGgene.id)*SCpropM
# 5.910803

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMupsigSC, GBMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.4032



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many transcripts are related to sexual conflict?
GBMdownsigSC <- sum(GBMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 10
# How many are not?
GBMdownsigNotSC <- length(GBMdownsig$CGgene.id)-GBMdownsigSC
# 46

# The number of expected transcripts which are related to sexual conflict
length(GBMdownsig$CGgene.id)*SCpropM
# 6.490293

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigSC, GBMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.1429




####################### GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many transcripts are related to sexual conflict?
GCMupsigSC <- sum(GCMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 5
# How many are not?
GCMupsigNotSC <- length(GCMupsig$CGgene.id)-GCMupsigSC
# 22

# The number of expected transcripts which are related to sexual conflict
length(GCMupsig$CGgene.id)*SCpropM
# 3.129249

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMupsigSC, GCMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.2607

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GCMupF <- rbind(c(GCMupsigSC, GCMupsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GCMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GCMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GCMupF) <- c("SA", "NotSA")
rownames(SC.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GCMupF)
# P = 0.2327

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many transcripts are related to sexual conflict?
GCMdownsigSC <- sum(GCMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 6
# How many are not?
GCMdownsigNotSC <- length(GCMdownsig$CGgene.id)-GCMdownsigSC
# 32

# The number of expected transcripts which are related to sexual conflict
length(GCMdownsig$CGgene.id)*SCpropM
# 4.404128

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigSC, GCMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.4187

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GCMdownF <- rbind(c(GCMdownsigSC, GCMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GCMdownsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GCMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GCMdownF) <- c("SA", "NotSA")
rownames(SC.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GCMdownF)
# P = 0.4411

# We get the same results, so that is good




####################### GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many transcripts are related to sexual conflict?
GDMupsigSC <- sum(GDMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 1
# How many are not?
GDMupsigNotSC <- length(GDMupsig$CGgene.id)-GDMupsigSC
# 11

# The number of expected transcripts which are related to sexual conflict
length(GDMupsig$CGgene.id)*SCpropM
# 1.390777

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMupsigSC, GDMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.7245

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GDMupF <- rbind(c(GDMupsigSC, GDMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GDMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GDMupsigNotSC))
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
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GDMdownsigSC <- sum(GDMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
GDMdownsigNotSC <- length(GDMdownsig$CGgene.id)-GDMdownsigSC
# 14

# The number of expected transcripts which are related to sexual conflict
length(GDMdownsig$CGgene.id)*SCpropM
# 1.85437

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigSC, GDMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.9094

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GDMdownF <- rbind(c(GDMdownsigSC, GDMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GDMdownsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GDMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GDMdownF) <- c("SA", "NotSA")
rownames(SC.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GDMdownF)
# P = 0.7077

# We get the same results, so that is good




####################### GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many transcripts are related to sexual conflict?
GEMupsigSC <- sum(GEMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
GEMupsigNotSC <- length(GEMupsig$CGgene.id)-GEMupsigSC
# 18

# The number of expected transcripts which are related to sexual conflict
length(GEMupsig$CGgene.id)*SCpropM
# 2.317962

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMupsigSC, GEMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.8242

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GEMupF <- rbind(c(GEMupsigSC, GEMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GEMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GEMupsigNotSC))
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
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many transcripts are related to sexual conflict?
GEMdownsigSC <- sum(GEMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 5
# How many are not?
GEMdownsigNotSC <- length(GEMdownsig$CGgene.id)-GEMdownsigSC
# 21

# The number of expected transcripts which are related to sexual conflict
length(GEMdownsig$CGgene.id)*SCpropM
# 3.013351

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigSC, GEMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.2235

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GEMdownF <- rbind(c(GEMdownsigSC, GEMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GEMdownsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GEMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GEMdownF) <- c("SA", "NotSA")
rownames(SC.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GEMdownF)
# P = 0.2166

# We get the same results, so that is good




####################### GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many transcripts are related to sexual conflict?
GFMupsigSC <- sum(GFMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GFMupsigNotSC <- length(GFMupsig$CGgene.id)-GFMupsigSC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMupsig$CGgene.id)*SCpropM
# 0.2317962

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMupsigSC, GFMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.6086

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFMupF <- rbind(c(GFMupsigSC, GFMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GFMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GFMupsigNotSC))
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
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many transcripts are related to sexual conflict?
GFMdownsigSC <- sum(GFMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GFMdownsigNotSC <- length(GFMdownsig$CGgene.id)-GFMdownsigSC
# 3

# The number of expected transcripts which are related to sexual conflict
length(GFMdownsig$CGgene.id)*SCpropM
# 0.3476943

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigSC, GFMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.5306

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFMupF <- rbind(c(GFMupsigSC, GFMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GFMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GFMupsigNotSC))
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
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many transcripts are related to sexual conflict?
GGMupsigSC <- sum(GGMupsig$CGgene.id %in% ImmonenSCM$gene.id)
# 0
# How many are not?
GGMupsigNotSC <- length(GGMupsig$CGgene.id)-GGMupsigSC
# 10

# The number of expected transcripts which are related to sexual conflict
length(GGMupsig$CGgene.id)*SCpropM
# 1.158981

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMupsigSC, GGMupsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.2522

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GGMupF <- rbind(c(GGMupsigSC, GGMupsigNotSC),
                   c(length(ImmonenSCM$gene.id)-GGMupsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GGMupsigNotSC))
# Add names to the columns and rows
colnames(SC.GGMupF) <- c("SC", "NotSC")
rownames(SC.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GGMupF)
# P = 0.6174

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many transcripts are related to sexual conflict?
GGMdownsigSC <- sum(GGMdownsig$CGgene.id %in% ImmonenSCM$gene.id)
# 2
# How many are not?
GGMdownsigNotSC <- length(GGMdownsig$CGgene.id)-GGMdownsigSC
# 20

# The number of expected transcripts which are related to sexual conflict
length(GGMdownsig$CGgene.id)*SCpropM
# 2.549758

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigSC, GGMdownsigNotSC), p = c(SCpropM, NotSCpropM))
# P = 0.7142

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GGMdownF <- rbind(c(GGMdownsigSC, GGMdownsigNotSC),
                     c(length(ImmonenSCM$gene.id)-GGMdownsigSC, (length(FLXswap7gM.data$CGgene.id)-length(ImmonenSCM$gene.id))-GGMdownsigNotSC))
# Add names to the columns and rows
colnames(SC.GGMdownF) <- c("SC", "NotSC")
rownames(SC.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GGMdownF)
# P = 1

# We get the same results, so that is good

