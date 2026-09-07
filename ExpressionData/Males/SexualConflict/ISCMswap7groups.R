################################################################################################
############### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES - 7 GROUPS ###############
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
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the male data set
InnocentiSCM <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% FLXswap7gM.data$gene.id)
# 1106 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropM <- length(InnocentiSCM$gene.id)/length(FLXswap7gM.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropM <- 1-ISCpropM
# 0.93



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many transcripts are related to sexual conflict?
GAMupsigISC <- sum(GAMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GAMupsigNotISC <- length(GAMupsig$gene.id)-GAMupsigISC
# 5

# The number of expected transcripts which are related to sexual conflict
length(GAMupsig$gene.id)*ISCpropM
# 0.3566591

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMupsigISC, GAMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.5354

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GAMupF <- rbind(c(GAMupsigISC, GAMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GAMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GAMupsigNotISC))
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
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GAMdownsigISC <- sum(GAMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 4
# How many are not?
GAMdownsigNotISC <- length(GAMdownsig$gene.id)-GAMdownsigISC
# 12

# The number of expected transcripts which are related to sexual conflict
length(GAMdownsig$gene.id)*ISCpropM
# 1.141309

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigISC, GAMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.005491

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GAMdownF <- rbind(c(GAMdownsigISC, GAMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GAMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GAMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GAMdownF) <- c("SC", "NotSC")
rownames(ISC.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GAMdownF)
# P = 0.02346

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many transcripts are related to sexual conflict?
GBMupsigISC <- sum(GBMupsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GBMupsigNotISC <- length(GBMupsig$gene.id)-GBMupsigISC
# 48

# The number of expected transcripts which are related to sexual conflict
length(GBMupsig$gene.id)*ISCpropM
# 3.637923

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMupsigISC, GBMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.7285

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GBMupF <- rbind(c(GBMupsigISC, GBMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GBMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GBMupsigNotISC))
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
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many transcripts are related to sexual conflict?
GBMdownsigISC <- sum(GBMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 12
# How many are not?
GBMdownsigNotISC <- length(GBMdownsig$gene.id)-GBMdownsigISC
# 44

# The number of expected transcripts which are related to sexual conflict
length(GBMdownsig$gene.id)*ISCpropM
# 3.994582

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigISC, GBMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 3.233e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GBMdownF <- rbind(c(GBMdownsigISC, GBMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GBMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GBMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GBMdownF) <- c("SC", "NotSC")
rownames(ISC.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GBMdownF)
# P = 0.0004853

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many transcripts are related to sexual conflict?
GCMupsigISC <- sum(GCMupsig$gene.id %in% InnocentiSCM$gene.id)
# 6
# How many are not?
GCMupsigNotISC <- length(GCMupsig$gene.id)-GCMupsigISC
# 21

# The number of expected transcripts which are related to sexual conflict
length(GCMupsig$gene.id)*ISCpropM
# 1.925959

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMupsigISC, GCMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.002317

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GCMupF <- rbind(c(GCMupsigISC, GCMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GCMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GCMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GCMupF) <- c("SC", "NotSC")
rownames(ISC.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GCMupF)
# P = 0.01051

# We went from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many transcripts are related to sexual conflict?
GCMdownsigISC <- sum(GCMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 12
# How many are not?
GCMdownsigNotISC <- length(GCMdownsig$gene.id)-GCMdownsigISC
# 26

# The number of expected transcripts which are related to sexual conflict
length(GCMdownsig$gene.id)*ISCpropM
# 2.710609

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigISC, GCMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 4.771e-09

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GCMdownF <- rbind(c(GCMdownsigISC, GCMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GCMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GCMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GCMdownF) <- c("SC", "NotSC")
rownames(ISC.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GCMdownF)
# P = 7.772e-066

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many transcripts are related to sexual conflict?
GDMupsigISC <- sum(GDMupsig$gene.id %in% InnocentiSCM$gene.id)
# 1
# How many are not?
GDMupsigNotISC <- length(GDMupsig$gene.id)-GDMupsigISC
# 11

# The number of expected transcripts which are related to sexual conflict
length(GDMupsig$gene.id)*ISCpropM
# 0.8559819

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMupsigISC, GDMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.8717

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GDMupF <- rbind(c(GDMupsigISC, GDMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GDMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GDMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GDMupF) <- c("SC", "NotSC")
rownames(ISC.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GDMupF)
# P = 0.5887

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many transcripts are related to sexual conflict?
GDMdownsigISC <- sum(GDMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GDMdownsigNotISC <- length(GDMdownsig$gene.id)-GDMdownsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(GDMdownsig$gene.id)*ISCpropM
# 1.141309

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigISC, GDMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.07101

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GDMdownF <- rbind(c(GDMdownsigISC, GDMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GDMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GDMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GDMdownF) <- c("SC", "NotSC")
rownames(ISC.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GDMdownF)
# P = 0.1011

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many transcripts are related to sexual conflict?
GEMupsigISC <- sum(GEMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GEMupsigNotISC <- length(GEMupsig$gene.id)-GEMupsigISC
# 20

# The number of expected transcripts which are related to sexual conflict
length(GEMupsig$gene.id)*ISCpropM
# 1.426637

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMupsigISC, GEMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.2152

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GEMupF <- rbind(c(GEMupsigISC, GEMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GEMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GEMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GEMupF) <- c("SC", "NotSC")
rownames(ISC.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GEMupF)
# P = 0.3948

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many transcripts are related to sexual conflict?
GEMdownsigISC <- sum(GEMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 7
# How many are not?
GEMdownsigNotISC <- length(GEMdownsig$gene.id)-GEMdownsigISC
# 19

# The number of expected transcripts which are related to sexual conflict
length(GEMdownsig$gene.id)*ISCpropM
# 1.854628

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigISC, GEMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 8.831e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GEMdownF <- rbind(c(GEMdownsigISC, GEMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GEMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GEMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GEMdownF) <- c("SC", "NotSC")
rownames(ISC.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GEMdownF)
# P = 0.00182

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many transcripts are related to sexual conflict?
GFMupsigISC <- sum(GFMupsig$gene.id %in% InnocentiSCM$gene.id)
# 0
# How many are not?
GFMupsigNotISC <- length(GFMupsig$gene.id)-GFMupsigISC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMupsig$gene.id)*ISCpropM
# 0.1426637

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMupsigISC, GFMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.6951

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFMupF <- rbind(c(GFMupsigISC, GFMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GFMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GFMupsigNotISC))
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
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many transcripts are related to sexual conflict?
GFMdownsigISC <- sum(GFMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 1
# How many are not?
GFMdownsigNotISC <- length(GFMdownsig$gene.id)-GFMdownsigISC
# 2

# The number of expected transcripts which are related to sexual conflict
length(GFMdownsig$gene.id)*ISCpropM
# 0.2139955

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigISC, GFMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.07787

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFMdownF <- rbind(c(GFMdownsigISC, GFMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GFMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GFMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GFMdownF) <- c("SC", "NotSC")
rownames(ISC.GFMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GFMdownF)
# P = 0.1991

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many transcripts are related to sexual conflict?
GGMupsigISC <- sum(GGMupsig$gene.id %in% InnocentiSCM$gene.id)
# 3
# How many are not?
GGMupsigNotISC <- length(GGMupsig$gene.id)-GGMupsigISC
# 7

# The number of expected transcripts which are related to sexual conflict
length(GGMupsig$gene.id)*ISCpropM
# 0.7133183

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMupsigISC, GGMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.004961

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GGMupF <- rbind(c(GGMupsigISC, GGMupsigNotISC),
                    c(length(InnocentiSCM$gene.id)-GGMupsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GGMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.GGMupF) <- c("SC", "NotSC")
rownames(ISC.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GGMupF)
# P = 0.02973

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many transcripts are related to sexual conflict?
GGMdownsigISC <- sum(GGMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 9
# How many are not?
GGMdownsigNotISC <- length(GGMdownsig$gene.id)-GGMdownsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(GGMdownsig$gene.id)*ISCpropM
# 1.5693

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigISC, GGMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 7.495e-10

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GGMdownF <- rbind(c(GGMdownsigISC, GGMdownsigNotISC),
                      c(length(InnocentiSCM$gene.id)-GGMdownsigISC, (length(FLXswap7gM.data$gene.id)-length(InnocentiSCM$gene.id))-GGMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.GGMdownF) <- c("SC", "NotSC")
rownames(ISC.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GGMdownF)
# P = 9.839e-06

# We get the same results, so that is good

