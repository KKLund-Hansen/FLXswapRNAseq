################################################################################################
##################### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP MALES ####################
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
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the female data set
InnocentiSCM <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% FLXswapM.data$gene.id)
# 1106 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropM <- length(InnocentiSCM$gene.id)/length(FLXswapM.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropM <- 1-ISCpropM
# 0.93



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many transcripts are related to sexual conflict?
FvCmupsigISC <- sum(FvCmupsig$gene.id %in% InnocentiSCM$gene.id)
# 7
# How many are not?
FvCmupsigNotISC <- length(FvCmupsig$gene.id)-FvCmupsigISC
# 50

# The number of expected transcripts which are related to sexual conflict
length(FvCmupsig$gene.id)*ISCpropM
# 4.065914

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCmupsigISC, FvCmupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.1311

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FvCmupF <- rbind(c(FvCmupsigISC, FvCmupsigNotISC),
                     c(length(InnocentiSCM$gene.id)-FvCmupsigISC, (length(FLXswapM.data$gene.id)-length(InnocentiSCM$gene.id))-FvCmupsigNotISC))
# Add names to the columns and rows
colnames(ISC.FvCmupF) <- c("SA", "NotSA")
rownames(ISC.FvCmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FvCmupF)
# P = 0.1241

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many transcripts are related to sexual conflict?
FvCmdownsigISC <- sum(FvCmdownsig$gene.id %in% InnocentiSCM$gene.id)
# 11
# How many are not?
FvCmdownsigNotISC <- length(FvCmdownsig$gene.id)-FvCmdownsigISC
# 46

# The number of expected transcripts which are related to sexual conflict
length(FvCmdownsig$gene.id)*ISCpropM
# 4.065914

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCmdownsigISC, FvCmdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.0003591

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FvCmdownF <- rbind(c(FvCmdownsigISC, FvCmdownsigNotISC),
                       c(length(InnocentiSCM$gene.id)-FvCmdownsigISC, (length(FLXswapM.data$gene.id)-length(InnocentiSCM$gene.id))-FvCmdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.FvCmdownF) <- c("SA", "NotSA")
rownames(ISC.FvCmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FvCmdownF)
# P = 0.00204

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many transcripts are related to sexual conflict?
AvSmupsigISC <- sum(AvSmupsig$gene.id %in% InnocentiSCM$gene.id)
# 34
# How many are not?
AvSmupsigNotISC <- length(AvSmupsig$gene.id)-AvSmupsigISC
# 181

# The number of expected transcripts which are related to sexual conflict
length(AvSmupsig$gene.id)*ISCpropM
# 15.33634

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSmupsigISC, AvSmupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 7.597e-07



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many transcripts are related to sexual conflict?
AvSmdownsigISC <- sum(AvSmdownsig$gene.id %in% InnocentiSCM$gene.id)
# 28
# How many are not?
AvSmdownsigNotISC <- length(AvSmdownsig$gene.id)-AvSmdownsigISC
# 234

# The number of expected transcripts which are related to sexual conflict
length(AvSmdownsig$gene.id)*ISCpropM
# 18.68894

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSmdownsigISC, AvSmdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.02542




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many transcripts are related to sexual conflict?
FsvCsMupsigISC <- sum(FsvCsMupsig$gene.id %in% InnocentiSCM$gene.id)
# 8
# How many are not?
FsvCsMupsigNotISC <- length(FsvCsMupsig$gene.id)-FsvCsMupsigISC
# 48

# The number of expected transcripts which are related to sexual conflict
length(FsvCsMupsig$gene.id)*ISCpropM
# 3.994582

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsMupsigISC, FsvCsMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.03756

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FsvCsMupF <- rbind(c(FsvCsMupsigISC, FsvCsMupsigNotISC),
                       c(length(InnocentiSCM$gene.id)-FsvCsMupsigISC, (length(FLXswapM.data$gene.id)-length(InnocentiSCM$gene.id))-FsvCsMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.FsvCsMupF) <- c("SA", "NotSA")
rownames(ISC.FsvCsMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FsvCsMupF)
# P = 0.05973

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many transcripts are related to sexual conflict?
FsvCsMdownsigISC <- sum(FsvCsMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 4
# How many are not?
FsvCsMdownsigNotISC <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigISC
# 45

# The number of expected transcripts which are related to sexual conflict
length(FsvCsMdownsig$gene.id)*ISCpropM
# 3.49526

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsMdownsigISC, FsvCsMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.7794

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FsvCsMdownF <- rbind(c(FsvCsMdownsigISC, FsvCsMdownsigNotISC),
                         c(length(InnocentiSCM$gene.id)-FsvCsMdownsigISC, (length(FLXswapM.data$gene.id)-length(InnocentiSCM$gene.id))-FsvCsMdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.FsvCsMdownF) <- c("SA", "NotSA")
rownames(ISC.FsvCsMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FsvCsMdownF)
# P = 0.7775

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many transcripts are related to sexual conflict?
FavCaMupsigISC <- sum(FavCaMupsig$gene.id %in% InnocentiSCM$gene.id)
# 8
# How many are not?
FavCaMupsigNotISC <- length(FavCaMupsig$gene.id)-FavCaMupsigISC
# 57

# The number of expected transcripts which are related to sexual conflict
length(FavCaMupsig$gene.id)*ISCpropM
# 4.636569

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaMupsigISC, FavCaMupsigNotISC), p = c(ISCpropM, NotISCpropM))
# P = 0.105

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FavCaMupF <- rbind(c(FavCaMupsigISC, FavCaMupsigNotISC),
                       c(length(InnocentiSCM$gene.id)-FavCaMupsigISC, (length(FLXswapM.data$gene.id)-length(InnocentiSCM$gene.id))-FavCaMupsigNotISC))
# Add names to the columns and rows
colnames(ISC.FavCaMupF) <- c("SA", "NotSA")
rownames(ISC.FavCaMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FavCaMupF)
# P = 0.1385

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many transcripts are related to sexual conflict?
FavCaMdownsigISC <- sum(FavCaMdownsig$gene.id %in% InnocentiSCM$gene.id)
# 34
# How many are not?
FavCaMdownsigNotISC <- length(FavCaMdownsig$gene.id)-FavCaMdownsigISC
# 87

# The number of expected transcripts which are related to sexual conflict
length(FavCaMdownsig$gene.id)*ISCpropM
# 8.631151

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaMdownsigISC, FavCaMdownsigNotISC), p = c(ISCpropM, NotISCpropM))
# P < 2.2e-16

