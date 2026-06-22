################################################################################################
##################### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES ####################
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
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the male data set
InnocentiSCM <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% RNAswapmale.data$gene.id)
# 1114 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropM <- length(InnocentiSCM$gene.id)/length(RNAswapmale.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropM <- 1-ISCpropM
# 0.93



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many transcripts are related to sexual conflict?
GAMupsigISC <- sum(GAMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GAMupsigNotISC <- length(GAMupsig$gene.id)-GAMupsigISC
# 4

# The number of expected transcripts which are related to sexual conflict
length(GAMupsig$gene.id)*ISCpropM
# 0.2654751

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMupsigISC, GAMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.5939

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GAMupF <- rbind(c(GAMupsigISC, GAMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GAMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GAMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GAMupF) <- c("SC", "NotSC")
rownames(ISC.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GAMdownsigISC <- sum(GAMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 4
# How many are not?
GAMdownsigNotISC <- length(GAMdownsig$gene.id)-GAMdownsigISC
# 12

# The number of expected transcripts which are related to sexual conflict
length(GAMdownsig$gene.id)*ISCpropM
# 1.061901

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigISC, GAMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.00317

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GAMdownF <- rbind(c(GAMdownsigISC, GAMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GAMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GAMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GAMdownF) <- c("SC", "NotSC")
rownames(ISC.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GAMdownF)
# P = 0.01846

# We went from significant to non-significant




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many transcripts are related to sexual conflict?
GBMupsigISC <- sum(GBMupsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GBMupsigNotISC <- length(GBMupsig$gene.id)-GBMupsigISC
# 51

# The number of expected transcripts which are related to sexual conflict
length(GBMupsig$gene.id)*ISCpropM
# 3.583914

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMupsigISC, GBMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.7496

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GBMupF <- rbind(c(GBMupsigISC, GBMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GBMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GBMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GBMupF) <- c("SC", "NotSC")
rownames(ISC.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GBMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many transcripts are related to sexual conflict?
GBMdownsigISC <- sum(GBMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 11
# How many are not?
GBMdownsigNotISC <- length(GBMdownsig$gene.id)-GBMdownsigISC
# 47

# The number of expected transcripts which are related to sexual conflict
length(GBMdownsig$gene.id)*ISCpropM
# 3.849389

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigISC, GBMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.000162

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GBMdownF <- rbind(c(GBMdownsigISC, GBMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GBMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GBMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GBMdownF) <- c("SC", "NotSC")
rownames(ISC.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GBMdownF)
# P = 0.001329

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many transcripts are related to sexual conflict?
GCMupsigISC <- sum(GCMupsig$gene.id %in% InnocentiSCM$gene.id)
# 6
# How many are not?
GCMupsigNotISC <- length(GCMupsig$gene.id)-GCMupsigISC
# 21

# The number of expected transcripts which are related to sexual conflict
length(GCMupsig$gene.id)*ISCpropM
# 1.791957

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMupsigISC, GCMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.001141

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GCMupF <- rbind(c(GCMupsigISC, GCMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GCMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GCMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GCMupF) <- c("SC", "NotSC")
rownames(ISC.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GCMupF)
# P = 0.007469

# We went from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many transcripts are related to sexual conflict?
GCMdownsigISC <- sum(GCMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 12
# How many are not?
GCMdownsigNotISC <- length(GCMdownsig$gene.id)-GCMdownsigISC
# 25

# The number of expected transcripts which are related to sexual conflict
length(GCMdownsig$gene.id)*ISCpropM
# 2.455645

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigISC, GCMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 2.912e-10

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GCMdownF <- rbind(c(GCMdownsigISC, GCMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GCMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GCMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GCMdownF) <- c("SC", "NotSC")
rownames(ISC.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GCMdownF)
# P = 2.697e-06

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many transcripts are related to sexual conflict?
GDMupsigISC <- sum(GDMupsig$gene.id %in% InnocentiSCM$gene.id)
# 1
# How many are not?
GDMupsigNotISC <- length(GDMupsig$gene.id)-GDMupsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(GDMupsig$gene.id)*ISCpropM
# 0.9291629

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMupsigISC, GDMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.9394

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GDMupF <- rbind(c(GDMupsigISC, GDMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GDMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GDMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GDMupF) <- c("SC", "NotSC")
rownames(ISC.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GDMupF)
# P = 0.6178

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GDMdownsigISC <- sum(GDMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GDMdownsigNotISC <- length(GDMdownsig$gene.id)-GDMdownsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(GDMdownsig$gene.id)*ISCpropM
# 1.061901

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigISC, GDMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.0516

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GDMdownF <- rbind(c(GDMdownsigISC, GDMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GDMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GDMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GDMdownF) <- c("SC", "NotSC")
rownames(ISC.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GDMdownF)
# P = 0.08547

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many transcripts are related to sexual conflict?
GEMupsigISC <- sum(GEMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GEMupsigNotISC <- length(GEMupsig$gene.id)-GEMupsigISC
# 17

# The number of expected transcripts which are related to sexual conflict
length(GEMupsig$gene.id)*ISCpropM
# 1.128269

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMupsigISC, GEMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.2716

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GEMupF <- rbind(c(GEMupsigISC, GEMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GEMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GEMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GEMupF) <- c("SC", "NotSC")
rownames(ISC.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GEMupF)
# P = 0.6238

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many transcripts are related to sexual conflict?
GEMdownsigISC <- sum(GEMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 7
# How many are not?
GEMdownsigNotISC <- length(GEMdownsig$gene.id)-GEMdownsigISC
# 19

# The number of expected transcripts which are related to sexual conflict
length(GEMdownsig$gene.id)*ISCpropM
# 1.725588

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigISC, GEMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 3.247e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GEMdownF <- rbind(c(GEMdownsigISC, GEMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GEMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GEMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GEMdownF) <- c("SC", "NotSC")
rownames(ISC.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GEMdownF)
# P = 0.001197

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many transcripts are related to sexual conflict?
GFMupsigISC <- sum(GFMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GFMupsigNotISC <- length(GFMupsig$gene.id)-GFMupsigISC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMupsig$gene.id)*ISCpropM
# 0.1327376

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMupsigISC, GFMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.7061

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFMupF <- rbind(c(GFMupsigISC, GFMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GFMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GFMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GFMupF) <- c("SC", "NotSC")
rownames(ISC.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many transcripts are related to sexual conflict?
GFMdownsigISC <- sum(GFMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 1
# How many are not?
GFMdownsigNotISC <- length(GFMdownsig$gene.id)-GFMdownsigISC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMdownsig$gene.id)*ISCpropM
# 0.1991063

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigISC, GFMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.06323

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFMdownF <- rbind(c(GFMdownsigISC, GFMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GFMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GFMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GFMdownF) <- c("SC", "NotSC")
rownames(ISC.GFMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GFMdownF)
# P = 0.1862

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many transcripts are related to sexual conflict?
GGMupsigISC <- sum(GGMupsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GGMupsigNotISC <- length(GGMupsig$gene.id)-GGMupsigISC
# 7

# The number of expected transcripts which are related to sexual conflict
length(GGMupsig$gene.id)*ISCpropM
# 0.6636878

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMupsigISC, GGMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.002998

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GGMupF <- rbind(c(GGMupsigISC, GGMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GGMupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GGMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GGMupF) <- c("SC", "NotSC")
rownames(ISC.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GGMupF)
# P = 0.02459

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many transcripts are related to sexual conflict?
GGMdownsigISC <- sum(GGMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 9
# How many are not?
GGMdownsigNotISC <- length(GGMdownsig$gene.id)-GGMdownsigISC
# 10

# The number of expected transcripts which are related to sexual conflict
length(GGMdownsig$gene.id)*ISCpropM
# 1.261007

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigISC, GGMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 9.861e-13

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GGMdownF <- rbind(c(GGMdownsigISC, GGMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GGMdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-GGMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GGMdownF) <- c("SC", "NotSC")
rownames(ISC.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GGMdownF)
# P = 1.217e-06

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many transcripts are related to sexual conflict?
G2MupsigISC <- sum(G2Mupsig$gene.id %in% InnocentiSCM$gene.id)
# 10
# How many are not?
G2MupsigNotISC <- length(G2Mupsig$gene.id)-G2MupsigISC
# 19

# The number of expected transcripts which are related to sexual conflict
length(G2Mupsig$gene.id)*ISCpropM
# 1.924695

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G2MupsigISC, G2MupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 1.701e-09

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G2MupF <- rbind(c(G2MupsigISC, G2MupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-G2MupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-G2MupsigNotISC))
# Add names to the columns and rows
colnames(ISC.G2MupF) <- c("SC", "NotSC")
rownames(ISC.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G2MupF)
# P = 9.959e-06

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many transcripts are related to sexual conflict?
G2MdownsigISC <- sum(G2Mdownsig$gene.id %in% InnocentiSCM$gene.id)
# 5
# How many are not?
G2MdownsigNotISC <- length(G2Mdownsig$gene.id)-G2MdownsigISC
# 31

# The number of expected transcripts which are related to sexual conflict
length(G2Mdownsig$gene.id)*ISCpropM
# 2.389276

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G2MdownsigISC, G2MdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.08046

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G2MdownF <- rbind(c(G2MdownsigISC, G2MdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-G2MdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-G2MdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.G2MdownF) <- c("SC", "NotSC")
rownames(ISC.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G2MdownF)
# P = 0.08733

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many transcripts are related to sexual conflict?
G11MupsigISC <- sum(G11Mupsig$gene.id %in% InnocentiSCM$gene.id)
# 8
# How many are not?
G11MupsigNotISC <- length(G11Mupsig$gene.id)-G11MupsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(G11Mupsig$gene.id)*ISCpropM
# 1.393744

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G11MupsigISC, G11MupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 6.984e-09

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G11MupF <- rbind(c(G11MupsigISC, G11MupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-G11MupsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-G11MupsigNotISC))
# Add names to the columns and rows
colnames(ISC.G11MupF) <- c("SC", "NotSC")
rownames(ISC.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G11MupF)
# P = 3.425e-05

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many transcripts are related to sexual conflict?
G11MdownsigISC <- sum(G11Mdownsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
G11MdownsigNotISC <- length(G11Mdownsig$gene.id)-G11MdownsigISC
# 7

# The number of expected transcripts which are related to sexual conflict
length(G11Mdownsig$gene.id)*ISCpropM
# 0.4645815

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G11MdownsigISC, G11MdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.4806

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G11MdownF <- rbind(c(G11MdownsigISC, G11MdownsigNotISC),
                       c(length(InnocentiSCM$gene.id)-G11MdownsigISC, (length(RNAswapmale.data$gene.id)-length(InnocentiSCM$gene.id))-G11MdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.G11MdownF) <- c("SC", "NotSC")
rownames(ISC.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G11MdownF)
# P = 1

# We get the same results, so that is good

