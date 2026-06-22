################################################################################################
######### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES AGAINST PSEUDOOBSCURA ########
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
rownames(RNAswapfemale.data) <- RNAswapfemale.data[,1]
# The add a new column with the CG gene names
RNAswapfemale.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(RNAswapfemale.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonene et al. data set to match the female data set
ImmonenSCF <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% RNAswapfemale.data$CGgene.id)
# 1808 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropF <- length(ImmonenSCF$gene.id)/length(RNAswapfemale.data$CGgene.id)
# 0.11
# And for the chisq.test I also have to calculate
NotSCpropF <- 1-SCpropF
# 0.89



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many transcripts are related to sexual conflict?
GAFupsigSC <- sum(GAFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 16
# How many are not?
GAFupsigNotSC <- length(GAFupsig$CGgene.id)-GAFupsigSC
# 113

#The number of expected transcripts which are related to sexual conflict
length(GAFupsig$CGgene.id)*SCpropF
# 13.89526

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFupsigSC, GAFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.55



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many transcripts are related to sexual conflict?
GAFdownsigSC <- sum(GAFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 292
# How many are not?
GAFdownsigNotSC <- length(GAFdownsig$CGgene.id)-GAFdownsigSC
# 1792

# The number of expected transcripts which are related to sexual conflict
length(GAFdownsig$CGgene.id)*SCpropF
# 217.6925

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigSC, GAFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 9.734e-08




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many transcripts are related to sexual conflict?
GBFupsigSC <- sum(GBFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 109
# How many are not?
GBFupsigNotSC <- length(GBFupsig$CGgene.id)-GBFupsigSC
# 461

#The number of expected transcripts which are related to sexual conflict
length(GBFupsig$CGgene.id)*SCpropF
# 61.39768

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFupsigSC, GBFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 1.265e-10



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many transcripts are related to sexual conflict?
GBFdownsigSC <- sum(GBFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 34
# How many are not?
GBFdownsigNotSC <- length(GBFdownsig$CGgene.id)-GBFdownsigSC
# 348

# The number of expected transcripts which are related to sexual conflict
length(GBFdownsig$CGgene.id)*SCpropF
# 41.14721

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigSC, GBFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.2382




####################### GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many transcripts are related to sexual conflict?
GCFupsigSC <- sum(GCFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 29
# How many are not?
GCFupsigNotSC <- length(GCFupsig$CGgene.id)-GCFupsigSC
# 148

#The number of expected transcripts which are related to sexual conflict
length(GCFupsig$CGgene.id)*SCpropF
# 19.06559

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFupsigSC, GCFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.01601



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many transcripts are related to sexual conflict?
GCFdownsigSC <- sum(GCFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 44
# How many are not?
GCFdownsigNotSC <- length(GCFdownsig$CGgene.id)-GCFdownsigSC
# 328

# The number of expected transcripts which are related to sexual conflict
length(GCFdownsig$CGgene.id)*SCpropF
# 40.07006

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigSC, GCFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.511




####################### GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many transcripts are related to sexual conflict?
GDFupsigSC <- sum(GDFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 40
# How many are not?
GDFupsigNotSC <- length(GDFupsig$CGgene.id)-GDFupsigSC
# 151

#The number of expected transcripts which are related to sexual conflict
length(GDFupsig$CGgene.id)*SCpropF
# 20.57361

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFupsigSC, GDFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 5.787e-06



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many transcripts are related to sexual conflict?
GDFdownsigSC <- sum(GDFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 39
# How many are not?
GDFdownsigNotSC <- length(GDFdownsig$CGgene.id)-GDFdownsigSC
# 245

# The number of expected transcripts which are related to sexual conflict
length(GDFdownsig$CGgene.id)*SCpropF
# 30.59112

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigSC, GDFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1075




####################### GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many transcripts are related to sexual conflict?
GEFupsigSC <- sum(GEFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 7
# How many are not?
GEFupsigNotSC <- length(GEFupsig$CGgene.id)-GEFupsigSC
# 91

#The number of expected transcripts which are related to sexual conflict
length(GEFupsig$CGgene.id)*SCpropF
# 10.55609

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFupsigSC, GEFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.2466



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many transcripts are related to sexual conflict?
GEFdownsigSC <- sum(GEFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 35
# How many are not?
GEFdownsigNotSC <- length(GEFdownsig$CGgene.id)-GEFdownsigSC
# 229

# The number of expected transcripts which are related to sexual conflict
length(GEFdownsig$CGgene.id)*SCpropF
# 28.43682

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigSC, GEFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1926




####################### GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many transcripts are related to sexual conflict?
GFFupsigSC <- sum(GFFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 5
# How many are not?
GFFupsigNotSC <- length(GFFupsig$CGgene.id)-GFFupsigSC
# 35

#The number of expected transcripts which are related to sexual conflict
length(GFFupsig$CGgene.id)*SCpropF
# 4.308609

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFupsigSC, GFFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.7244

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFFupF <- rbind(c(GFFupsigSC, GFFupsigNotSC),
                   c(length(ImmonenSCF$gene.id)-GFFupsigSC, (length(RNAswapfemale.data$CGgene.id)-length(ImmonenSCF$gene.id))-GFFupsigNotSC))
#Add names to the columns and rows
colnames(SC.GFFupF) <- c("SC", "NotSC")
rownames(SC.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GFFupF)
# P = 0.6144

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many transcripts are related to sexual conflict?
GFFdownsigSC <- sum(GFFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 33
# How many are not?
GFFdownsigNotSC <- length(GFFdownsig$CGgene.id)-GFFdownsigSC
# 231

# The number of expected transcripts which are related to sexual conflict
length(GFFdownsig$CGgene.id)*SCpropF
# 28.43682

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigSC, GFFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.365




####################### GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many transcripts are related to sexual conflict?
GGFupsigSC <- sum(GGFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 8
# How many are not?
GGFupsigNotSC <- length(GGFupsig$CGgene.id)-GGFupsigSC
# 84

#The number of expected transcripts which are related to sexual conflict
length(GGFupsig$CGgene.id)*SCpropF
# 9.9098

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFupsigSC, GGFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.5207



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many transcripts are related to sexual conflict?
GGFdownsigSC <- sum(GGFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 21
# How many are not?
GGFdownsigNotSC <- length(GGFdownsig$CGgene.id)-GGFdownsigSC
# 178

# The number of expected transcripts which are related to sexual conflict
length(GGFdownsig$CGgene.id)*SCpropF
# 21.43533

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigSC, GGFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.9207




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many transcripts are related to sexual conflict?
G23FupsigSC <- sum(G23Fupsig$gene.id %in% ImmonenSCF$gene.id)
# 0
# How many are not?
G23FupsigNotSC <- length(G23Fupsig$gene.id)-G23FupsigSC
# 16

# The number of expected transcripts which are related to sexual conflict
length(G23Fupsig$gene.id)*SCpropF
# 1.723444

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G23FupsigSC, G23FupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1646

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G23FupF <- rbind(c(G23FupsigSC, G23FupsigNotSC),
                    c(length(ImmonenSCF$gene.id)-G23FupsigSC, (length(RNAswapfemale.data$gene.id)-length(ImmonenSCF$gene.id))-G23FupsigNotSC))
# Add names to the columns and rows
colnames(SC.G23FupF) <- c("SC", "NotSC")
rownames(SC.G23FupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G23FupF)
# P = 0.4056

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many transcripts are related to sexual conflict?
G23FdownsigSC <- sum(G23Fdownsig$gene.id %in% ImmonenSCF$gene.id)
# 0
# How many are not?
G23FdownsigNotSC <- length(G23Fdownsig$gene.id)-G23FdownsigSC
# 6

# The number of expected transcripts which are related to sexual conflict
length(G23Fdownsig$gene.id)*SCpropF
# 0.6462913

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(G23FdownsigSC, G23FdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.3947

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.G23FdownF <- rbind(c(G23FdownsigSC, G23FdownsigNotSC),
                      c(length(ImmonenSCF$gene.id)-G23FdownsigSC, (length(RNAswapfemale.data$gene.id)-length(ImmonenSCF$gene.id))-G23FdownsigNotSC))
# Add names to the columns and rows
colnames(SC.G23FdownF) <- c("SC", "NotSC")
rownames(SC.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.G23FdownF)
# P = 1

# We get the same results, so that is good

