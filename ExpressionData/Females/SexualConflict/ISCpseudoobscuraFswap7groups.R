################################################################################################
###### CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES v PSEUDOOBSCURA - 7 GROUPS ######
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
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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
rownames(FLXswap7gF.data) <- FLXswap7gF.data[,1]
# The add a new column with the CG gene names
FLXswap7gF.data$CGgene.id <- mapIds(drosophila2.db, keys = row.names(FLXswap7gF.data), keytype = "FLYBASE", column = "FLYBASECG")


# Now I can finally subset the Immonene et al. data set to match the female data set
ImmonenSCF <- subset(ImmonenSC.data, ImmonenSC.data$gene.id %in% FLXswap7gF.data$CGgene.id)
# 1797 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
SCpropF <- length(ImmonenSCF$gene.id)/length(FLXswap7gF.data$CGgene.id)
# 0.12
# And for the chisq.test I also have to calculate
NotSCpropF <- 1-SCpropF
# 0.88



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many transcripts are related to sexual conflict?
GAFupsigSC <- sum(GAFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 13
# How many are not?
GAFupsigNotSC <- length(GAFupsig$CGgene.id)-GAFupsigSC
# 109

#The number of expected transcripts which are related to sexual conflict
length(GAFupsig$CGgene.id)*SCpropF
# 14.13957

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFupsigSC, GAFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.7472



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many transcripts are related to sexual conflict?
GAFdownsigSC <- sum(GAFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 288
# How many are not?
GAFdownsigNotSC <- length(GAFdownsig$CGgene.id)-GAFdownsigSC
# 1683

# The number of expected transcripts which are related to sexual conflict
length(GAFdownsig$CGgene.id)*SCpropF
# 228.4351

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigSC, GAFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 2.773e-05




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many transcripts are related to sexual conflict?
GBFupsigSC <- sum(GBFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 112
# How many are not?
GBFupsigNotSC <- length(GBFupsig$CGgene.id)-GBFupsigSC
# 451

#The number of expected transcripts which are related to sexual conflict
length(GBFupsig$CGgene.id)*SCpropF
# 65.25063

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFupsigSC, GBFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 7.505e-10



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many transcripts are related to sexual conflict?
GBFdownsigSC <- sum(GBFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 35
# How many are not?
GBFdownsigNotSC <- length(GBFdownsig$CGgene.id)-GBFdownsigSC
# 347

# The number of expected transcripts which are related to sexual conflict
length(GBFdownsig$CGgene.id)*SCpropF
# 44.27307

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigSC, GBFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1383




####################### GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many transcripts are related to sexual conflict?
GCFupsigSC <- sum(GCFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 25
# How many are not?
GCFupsigNotSC <- length(GCFupsig$CGgene.id)-GCFupsigSC
# 137

#The number of expected transcripts which are related to sexual conflict
length(GCFupsig$CGgene.id)*SCpropF
# 18.77549

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFupsigSC, GCFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.1266



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many transcripts are related to sexual conflict?
GCFdownsigSC <- sum(GCFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 41
# How many are not?
GCFdownsigNotSC <- length(GCFdownsig$CGgene.id)-GCFdownsigSC
# 324

# The number of expected transcripts which are related to sexual conflict
length(GCFdownsig$CGgene.id)*SCpropF
# 42.30281

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigSC, GCFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.8313




####################### GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many transcripts are related to sexual conflict?
GDFupsigSC <- sum(GDFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 36
# How many are not?
GDFupsigNotSC <- length(GDFupsig$CGgene.id)-GDFupsigSC
# 140

#The number of expected transcripts which are related to sexual conflict
length(GDFupsig$CGgene.id)*SCpropF
# 20.39807

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFupsigSC, GDFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.0002388



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many transcripts are related to sexual conflict?
GDFdownsigSC <- sum(GDFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 39
# How many are not?
GDFdownsigNotSC <- length(GDFdownsig$CGgene.id)-GDFdownsigSC
# 244

# The number of expected transcripts which are related to sexual conflict
length(GDFdownsig$CGgene.id)*SCpropF
# 32.79916

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigSC, GDFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.2495




####################### GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many transcripts are related to sexual conflict?
GEFupsigSC <- sum(GEFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 8
# How many are not?
GEFupsigNotSC <- length(GEFupsig$CGgene.id)-GEFupsigSC
# 83

#The number of expected transcripts which are related to sexual conflict
length(GEFupsig$CGgene.id)*SCpropF
# 10.54673

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFupsigSC, GEFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.4043



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many transcripts are related to sexual conflict?
GEFdownsigSC <- sum(GEFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 37
# How many are not?
GEFdownsigNotSC <- length(GEFdownsig$CGgene.id)-GEFdownsigSC
# 226

# The number of expected transcripts which are related to sexual conflict
length(GEFdownsig$CGgene.id)*SCpropF
# 30.4812

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigSC, GEFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.2092




####################### GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many transcripts are related to sexual conflict?
GFFupsigSC <- sum(GFFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 5
# How many are not?
GFFupsigNotSC <- length(GFFupsig$CGgene.id)-GFFupsigSC
# 32

#The number of expected transcripts which are related to sexual conflict
length(GFFupsig$CGgene.id)*SCpropF
# 4.28823

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFupsigSC, GFFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.7147

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SC.GFFupF <- rbind(c(GFFupsigSC, GFFupsigNotSC),
                   c(length(ImmonenSCF$gene.id)-GFFupsigSC, (length(FLXswap7gF.data$CGgene.id)-length(ImmonenSCF$gene.id))-GFFupsigNotSC))
#Add names to the columns and rows
colnames(SC.GFFupF) <- c("SC", "NotSC")
rownames(SC.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SC.GFFupF)
# P = 0.6117

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many transcripts are related to sexual conflict?
GFFdownsigSC <- sum(GFFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 31
# How many are not?
GFFdownsigNotSC <- length(GFFdownsig$CGgene.id)-GFFdownsigSC
# 232

# The number of expected transcripts which are related to sexual conflict
length(GFFdownsig$CGgene.id)*SCpropF
# 30.4812

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigSC, GFFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.9204




####################### GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many transcripts are related to sexual conflict?
GGFupsigSC <- sum(GGFupsig$CGgene.id %in% ImmonenSCF$gene.id)
# 10
# How many are not?
GGFupsigNotSC <- length(GGFupsig$CGgene.id)-GGFupsigSC
# 85

#The number of expected transcripts which are related to sexual conflict
length(GGFupsig$CGgene.id)*SCpropF
# 11.01032

#Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFupsigSC, GGFupsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.7461



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

#First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many transcripts are related to sexual conflict?
GGFdownsigSC <- sum(GGFdownsig$CGgene.id %in% ImmonenSCF$gene.id)
# 24
# How many are not?
GGFdownsigNotSC <- length(GGFdownsig$CGgene.id)-GGFdownsigSC
# 170

# The number of expected transcripts which are related to sexual conflict
length(GGFdownsig$CGgene.id)*SCpropF
# 22.48423

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigSC, GGFdownsigNotSC), p = c(SCpropF, NotSCpropF))
# P = 0.7339


