################################################################################################
#################### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES ###################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/SexualConflict")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the female data set
InnocentiSCF <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% RNAswapfemale.data$gene.id)
# 1114 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropF <- length(InnocentiSCF$gene.id)/length(RNAswapfemale.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropF <- 1-ISCpropF
# 0.93



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many transcripts are related to sexual conflict?
GAFupsigISC <- sum(GAFupsig$gene.id %in% InnocentiSCF$gene.id)
# 19
# How many are not?
GAFupsigNotISC <- length(GAFupsig$gene.id)-GAFupsigISC
# 110

# The number of expected transcripts which are related to sexual conflict
length(GAFupsig$gene.id)*ISCpropF
# 8.561573

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFupsigISC, GAFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.0002224



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many transcripts are related to sexual conflict?
GAFdownsigISC <- sum(GAFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 178
# How many are not?
GAFdownsigNotISC <- length(GAFdownsig$gene.id)-GAFdownsigISC
# 1843

# The number of expected transcripts which are related to sexual conflict
length(GAFdownsig$gene.id)*ISCpropF
# 134.1313

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigISC, GAFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 8.85e-05


######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many transcripts are related to sexual conflict?
GBFupsigISC <- sum(GBFupsig$gene.id %in% InnocentiSCF$gene.id)
# 51
# How many are not?
GBFupsigNotISC <- length(GBFupsig$gene.id)-GBFupsigISC
# 519

# The number of expected transcripts which are related to sexual conflict
length(GBFupsig$gene.id)*ISCpropF
# 37.83021

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFupsigISC, GBFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.02669



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many transcripts are related to sexual conflict?
GBFdownsigISC <- sum(GBFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 49
# How many are not?
GBFdownsigNotISC <- length(GBFdownsig$gene.id)-GBFdownsigISC
# 333

# The number of expected transcripts which are related to sexual conflict
length(GBFdownsig$gene.id)*ISCpropF
# 25.35287

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigISC, GBFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 1.171e-06




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many transcripts are related to sexual conflict?
GCFupsigISC <- sum(GCFupsig$gene.id %in% InnocentiSCF$gene.id)
# 12
# How many are not?
GCFupsigNotISC <- length(GCFupsig$gene.id)-GCFupsigISC
# 165

# The number of expected transcripts which are related to sexual conflict
length(GCFupsig$gene.id)*ISCpropF
# 11.74727

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFupsigISC, GCFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.9392



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many transcripts are related to sexual conflict?
GCFdownsigISC <- sum(GCFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 42
# How many are not?
GCFdownsigNotISC <- length(GCFdownsig$gene.id)-GCFdownsigISC
# 330

# The number of expected transcripts which are related to sexual conflict
length(GCFdownsig$gene.id)*ISCpropF
# 24.68919

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigISC, GCFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.0003114




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many transcripts are related to sexual conflict?
GDFupsigISC <- sum(GDFupsig$gene.id %in% InnocentiSCF$gene.id)
# 10
# How many are not?
GDFupsigNotISC <- length(GDFupsig$gene.id)-GDFupsigISC
# 181

# The number of expected transcripts which are related to sexual conflict
length(GDFupsig$gene.id)*ISCpropF
# 12.67644

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFupsigISC, GDFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.4366



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many transcripts are related to sexual conflict?
GDFdownsigISC <- sum(GDFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 28
# How many are not?
GDFdownsigNotISC <- length(GDFdownsig$gene.id)-GDFdownsigISC
# 256

# The number of expected transcripts which are related to sexual conflict
length(GDFdownsig$gene.id)*ISCpropF
# 18.84873

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigISC, GDFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.02915




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many transcripts are related to sexual conflict?
GEFupsigISC <- sum(GEFupsig$gene.id %in% InnocentiSCF$gene.id)
# 14
# How many are not?
GEFupsigNotISC <- length(GEFupsig$gene.id)-GEFupsigISC
# 84

# The number of expected transcripts which are related to sexual conflict
length(GEFupsig$gene.id)*ISCpropF
# 6.504141

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFupsigISC, GEFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.002351



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many transcripts are related to sexual conflict?
GEFdownsigISC <- sum(GEFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 19
# How many are not?
GEFdownsigNotISC <- length(GEFdownsig$gene.id)-GEFdownsigISC
# 245

# The number of expected transcripts which are related to sexual conflict
length(GEFdownsig$gene.id)*ISCpropF
# 17.52136

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigISC, GEFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.7147




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many transcripts are related to sexual conflict?
GFFupsigISC <- sum(GFFupsig$gene.id %in% InnocentiSCF$gene.id)
# 5
# How many are not?
GFFupsigNotISC <- length(GFFupsig$gene.id)-GFFupsigISC
# 35

# The number of expected transcripts which are related to sexual conflict
length(GFFupsig$gene.id)*ISCpropF
# 2.654751

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFupsigISC, GFFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.1363

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFFupF <- rbind(c(GFFupsigISC, GFFupsigNotISC),
                    c(length(InnocentiSCF$gene.id)-GFFupsigISC, (length(RNAswapfemale.data$gene.id)-length(InnocentiSCF$gene.id))-GFFupsigNotISC))
#Add names to the columns and rows
colnames(ISC.GFFupF) <- c("SA", "NotSA")
rownames(ISC.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GFFupF)
# P = 0.1872

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many transcripts are related to sexual conflict?
GFFdownsigISC <- sum(GFFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 30
# How many are not?
GFFdownsigNotISC <- length(GFFdownsig$gene.id)-GFFdownsigISC
# 234

# The number of expected transcripts which are related to sexual conflict
length(GFFdownsig$gene.id)*ISCpropF
# 17.52136

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigISC, GFFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.002034




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many transcripts are related to sexual conflict?
GGFupsigISC <- sum(GGFupsig$gene.id %in% InnocentiSCF$gene.id)
# 10
# How many are not?
GGFupsigNotISC <- length(GGFupsig$gene.id)-GGFupsigISC
# 82

# The number of expected transcripts which are related to sexual conflict
length(GGFupsig$gene.id)*ISCpropF
# 6.105928

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFupsigISC, GGFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.1029



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many transcripts are related to sexual conflict?
GGFdownsigISC <- sum(GGFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 20
# How many are not?
GGFdownsigNotISC <- length(GGFdownsig$gene.id)-GGFdownsigISC
# 179

# The number of expected transcripts which are related to sexual conflict
length(GGFdownsig$gene.id)*ISCpropF
# 13.20739

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigISC, GGFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.05307




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many transcripts are related to sexual conflict?
G23FupsigISC <- sum(G23Fupsig$gene.id %in% InnocentiSCF$gene.id)
# 3
# How many are not?
G23FupsigNotISC <- length(G23Fupsig$gene.id)-G23FupsigISC
# 13

# The number of expected transcripts which are related to sexual conflict
length(G23Fupsig$gene.id)*ISCpropF
# 1.061901

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G23FupsigISC, G23FupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.0516

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G23FupF <- rbind(c(G23FupsigISC, G23FupsigNotISC),
                     c(length(InnocentiSCF$gene.id)-G23FupsigISC, (length(RNAswapfemale.data$gene.id)-length(InnocentiSCF$gene.id))-G23FupsigNotISC))
#Add names to the columns and rows
colnames(ISC.G23FupF) <- c("SA", "NotSA")
rownames(ISC.G23FupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G23FupF)
# P = 0.08547

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many transcripts are related to sexual conflict?
G23FdownsigISC <- sum(G23Fdownsig$gene.id %in% InnocentiSCF$gene.id)
# 0
# How many are not?
G23FdownsigNotISC <- length(G23Fdownsig$gene.id)-G23FdownsigISC
# 6

# The number of expected transcripts which are related to sexual conflict
length(G23Fdownsig$gene.id)*ISCpropF
# 0.3982127

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G23FdownsigISC, G23FdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.5137

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.G23FdownF <- rbind(c(G23FdownsigISC, G23FdownsigNotISC),
                       c(length(InnocentiSCF$gene.id)-G23FdownsigISC, (length(RNAswapfemale.data$gene.id)-length(InnocentiSCF$gene.id))-G23FdownsigNotISC))
# Add names to the columns and rows
colnames(ISC.G23FdownF) <- c("SA", "NotSA")
rownames(ISC.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.G23FdownF)
# P = 1

# We get the same results, so that is good

