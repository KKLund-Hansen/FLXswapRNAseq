################################################################################################
################### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES  ###################
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
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the female data set
InnocentiSCF <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% FLXswapF.data$gene.id)
# 1106 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropF <- length(InnocentiSCF$gene.id)/length(FLXswapF.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropF <- 1-ISCpropF
# 0.93



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many transcripts are related to sexual conflict?
FvCfupsigISC <- sum(FvCfupsig$gene.id %in% InnocentiSCF$gene.id)
# 9
# How many are not?
FvCfupsigNotISC <- length(FvCfupsig$gene.id)-FvCfupsigISC
# 26

# The number of expected transcripts which are related to sexual conflict
length(FvCfupsig$gene.id)*ISCpropF
# 2.496614

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCfupsigISC, FvCfupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 1.946e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FvCfupF <- rbind(c(FvCfupsigISC, FvCfupsigNotISC),
                     c(length(InnocentiSCF$gene.id)-FvCfupsigISC, (length(FLXswapF.data$gene.id)-length(InnocentiSCF$gene.id))-FvCfupsigNotISC))
# Add names to the columns and rows
colnames(ISC.FvCfupF) <- c("SA", "NotSA")
rownames(ISC.FvCfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FvCfupF)
# P = 0.0006005

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many transcripts are related to sexual conflict?
FvCfdownsigISC <- sum(FvCfdownsig$gene.id %in% InnocentiSCF$gene.id)
# 5
# How many are not?
FvCfdownsigNotISC <- length(FvCfdownsig$gene.id)-FvCfdownsigISC
# 35

# The number of expected transcripts which are related to sexual conflict
length(FvCfdownsig$gene.id)*ISCpropF
# 2.853273

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FvCfdownsigISC, FvCfdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.1872

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FvCfdownF <- rbind(c(FvCfdownsigISC, FvCfdownsigNotISC),
                       c(length(InnocentiSCF$gene.id)-FvCfdownsigISC, (length(FLXswapF.data$gene.id)-length(InnocentiSCF$gene.id))-FvCfdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.FvCfdownF) <- c("SA", "NotSA")
rownames(ISC.FvCfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FvCfdownF)
# P = 0.2051

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many transcripts are related to sexual conflict?
AvSfupsigISC <- sum(AvSfupsig$gene.id %in% InnocentiSCF$gene.id)
# 36
# How many are not?
AvSfupsigNotISC <- length(AvSfupsig$gene.id)-AvSfupsigISC
# 208

# The number of expected transcripts which are related to sexual conflict
length(AvSfupsig$gene.id)*ISCpropF
# 17.40497

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSfupsigISC, AvSfupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 3.742e-06



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig") 
# 161

# How many transcripts are related to sexual conflict?
AvSfdownsigISC <- sum(AvSfdownsig$gene.id %in% InnocentiSCF$gene.id)
# 37
# How many are not?
AvSfdownsigNotISC <- length(AvSfdownsig$gene.id)-AvSfdownsigISC
# 124

# The number of expected transcripts which are related to sexual conflict
length(AvSfdownsig$gene.id)*ISCpropF
# 11.48442

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(AvSfdownsigISC, AvSfdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 5.582e-15




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many transcripts are related to sexual conflict?
FsvCsFupsigISC <- sum(FsvCsFupsig$gene.id %in% InnocentiSCF$gene.id)
# 15
# How many are not?
FsvCsFupsigNotISC <- length(FsvCsFupsig$gene.id)-FsvCsFupsigISC
# 122

# The number of expected transcripts which are related to sexual conflict
length(FsvCsFupsig$gene.id)*ISCpropF
# 9.77246

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsFupsigISC, FsvCsFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.08269



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many transcripts are related to sexual conflict?
FsvCsFdownsigISC <- sum(FsvCsFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 4
# How many are not?
FsvCsFdownsigNotISC <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigISC
# 51

# The number of expected transcripts which are related to sexual conflict
length(FsvCsFdownsig$gene.id)*ISCpropF
# 3.923251

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FsvCsFdownsigISC, FsvCsFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.9679

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.FsvCsFdownF <- rbind(c(FsvCsFdownsigISC, FsvCsFdownsigNotISC),
                         c(length(InnocentiSCF$gene.id)-FsvCsFdownsigISC, (length(FLXswapF.data$gene.id)-length(InnocentiSCF$gene.id))-FsvCsFdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.FsvCsFdownF) <- c("SA", "NotSA")
rownames(ISC.FsvCsFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.FsvCsFdownF)
# P = 0.7968

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many transcripts are related to sexual conflict?
FavCaFupsigISC <- sum(FavCaFupsig$gene.id %in% InnocentiSCF$gene.id)
# 47
# How many are not?
FavCaFupsigNotISC <- length(FavCaFupsig$gene.id)-FavCaFupsigISC
# 478

# The number of expected transcripts which are related to sexual conflict
length(FavCaFupsig$gene.id)*ISCpropF
# 37.44921

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaFupsigISC, FavCaFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.1053



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts  in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many transcripts are related to sexual conflict?
FavCaFdownsigISC <- sum(FavCaFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 295
# How many are not?
FavCaFdownsigNotISC <- length(FavCaFdownsig$gene.id)-FavCaFdownsigISC
# 2699

# The number of expected transcripts which are related to sexual conflict
length(FavCaFdownsig$gene.id)*ISCpropF
# 213.5675

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(FavCaFdownsigISC, FavCaFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 7.369e-09

