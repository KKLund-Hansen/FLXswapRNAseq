################################################################################################
########################## MITO-SENSITIVE TRANSCRIPTS FLX SWAP FEMALES ######################### 
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the female data set
MSfemaleswap <- subset(MS.data, MS.data$gene.id %in% RNAswapfemale.data$gene.id)
# 1221 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropF <- length(MSfemaleswap$gene.id)/length(RNAswapfemale.data$gene.id)
# 0.07

# And for the chisq.test I also have to calculate
NotMSpropF <- 1-MSpropF
# 0.93




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many mito-sensitive transcripts are there?
GAFupsigMS <- sum(GAFupsig$gene.id %in% MSfemaleswap$gene.id)
# 6
# How many are not?
GAFupsigNotMS <- length(GAFupsig$gene.id)-GAFupsigMS
# 123

# The number of expected mito-sensitive transcripts
length(GAFupsig$gene.id)*MSpropF
# 9.383914

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMS, GAFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2513



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many mito-sensitive transcripts are there?
GAFdownsigMS <- sum(GAFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 217
# How many are not?
GAFdownsigNotMS <- length(GAFdownsig$gene.id)-GAFdownsigMS
# 1804

# The number of expected mito-sensitive transcripts
length(GAFdownsig$gene.id)*MSpropF
# 147.0147

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMS, GAFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 2.046e-09




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many mito-sensitive transcripts are there?
GBFupsigMS <- sum(GBFupsig$gene.id %in% MSfemaleswap$gene.id)
# 55
# How many are not?
GBFupsigNotMS <- length(GBFupsig$gene.id)-GBFupsigMS
# 515

# The number of expected mito-sensitive transcripts
length(GBFupsig$gene.id)*MSpropF
# 41.46381

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMS, GBFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.02903



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many mito-sensitive transcripts are there?
GBFdownsigMS <- sum(GBFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 20
# How many are not?
GBFdownsigNotMS <- length(GBFdownsig$gene.id)-GBFdownsigMS
# 362

# The number of expected mito-sensitive transcripts
length(GBFdownsig$gene.id)*MSpropF
# 27.78803

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMS, GBFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.125




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many mito-sensitive transcripts are there?
GCFupsigMS <- sum(GCFupsig$gene.id %in% MSfemaleswap$gene.id)
# 16
# How many are not?
GCFupsigNotMS <- length(GCFupsig$gene.id)-GCFupsigMS
# 161

# The number of expected mito-sensitive transcripts
length(GCFupsig$gene.id)*MSpropF
# 12.8756

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMS, GCFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.3659



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many mito-sensitive transcripts are there?
GCFdownsigMS <- sum(GCFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 14
# How many are not?
GCFdownsigNotMS <- length(GCFdownsig$gene.id)-GCFdownsigMS
# 358

# The number of expected mito-sensitive transcripts
length(GCFdownsig$gene.id)*MSpropF
# 27.06059

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMS, GCFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.009125




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many mito-sensitive transcripts are there?
GDFupsigMS <- sum(GDFupsig$gene.id %in% MSfemaleswap$gene.id)
# 15
# How many are not?
GDFupsigNotMS <- length(GDFupsig$gene.id)-GDFupsigMS
# 176

# The number of expected mito-sensitive transcripts
length(GDFupsig$gene.id)*MSpropF
# 13.89401

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMS, GDFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.758



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many mito-sensitive transcripts are there?
GDFdownsigMS <- sum(GDFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 9
# How many are not?
GDFdownsigNotMS <- length(GDFdownsig$gene.id)-GDFdownsigMS
# 275

# The number of expected mito-sensitive transcripts
length(GDFdownsig$gene.id)*MSpropF
# 20.65916

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMS, GDFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.007725




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many mito-sensitive transcripts are there?
GEFupsigMS <- sum(GEFupsig$gene.id %in% MSfemaleswap$gene.id)
# 4
# How many are not?
GEFupsigNotMS <- length(GEFupsig$gene.id)-GEFupsigMS
# 94

# The number of expected mito-sensitive transcripts
length(GEFupsig$gene.id)*MSpropF
# 7.128865

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMS, GEFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2236



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many mito-sensitive transcripts are there?
GEFdownsigMS <- sum(GEFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 16
# How many are not?
GEFdownsigNotMS <- length(GEFdownsig$gene.id)-GEFdownsigMS
# 248

# The number of expected mito-sensitive transcripts
length(GEFdownsig$gene.id)*MSpropF
# 19.20429

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMS, GEFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.4477




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many mito-sensitive transcripts are there?
GFFupsigMS <- sum(GFFupsig$gene.id %in% MSfemaleswap$gene.id)
# 3
# How many are not?
GFFupsigNotMS <- length(GFFupsig$gene.id)-GFFupsigMS
# 37

# The number of expected mito-sensitive transcripts
length(GFFupsig$gene.id)*MSpropF
# 2.909741

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMS, GFFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.9562

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFFupF <- rbind(c(GFFupsigMS, GFFupsigNotMS), 
                   c(length(MSfemaleswap$gene.id)-GFFupsigMS, (length(RNAswapfemale.data$gene.id)-length(MSfemaleswap$gene.id))-GFFupsigNotMS))
#Add names to the columns and rows
colnames(MS.GFFupF) <- c("MS", "NotMS")
rownames(MS.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFFupF)
# P = 0.7659

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many mito-sensitive transcripts are there?
GFFdownsigMS <- sum(GFFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 21
# How many are not?
GFFdownsigNotMS <- length(GFFdownsig$gene.id)-GFFdownsigMS
# 243

# The number of expected mito-sensitive transcripts
length(GFFdownsig$gene.id)*MSpropF
# 19.20429

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMS, GFFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.6704




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many mito-sensitive transcripts are there?
GGFupsigMS <- sum(GGFupsig$gene.id %in% MSfemaleswap$gene.id)
# 4
# How many are not?
GGFupsigNotMS <- length(GGFupsig$gene.id)-GGFupsigMS
# 88

# The number of expected mito-sensitive transcripts
length(GGFupsig$gene.id)*MSpropF
# 6.692404

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMS, GGFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2798



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many mito-sensitive transcripts are there?
GGFdownsigMS <- sum(GGFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 18
# How many are not?
GGFdownsigNotMS <- length(GGFdownsig$gene.id)-GGFdownsigMS
# 181

# The number of expected mito-sensitive transcripts
length(GGFdownsig$gene.id)*MSpropF
# 14.47596

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMS, GGFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.3361




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many mito-sensitive transcripts are there?
G23FupsigMS <- sum(G23Fupsig$gene.id %in% MSfemaleswap$gene.id)
# 2
# How many are not?
G23FupsigNotMS <- length(G23Fupsig$gene.id)-G23FupsigMS
# 14

# The number of expected mito-sensitive transcripts
length(G23Fupsig$gene.id)*MSpropF
# 1.163896

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FupsigMS, G23FupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.4209

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G23FupF <- rbind(c(G23FupsigMS, G23FupsigNotMS), 
                    c(length(MSfemaleswap$gene.id)-G23FupsigMS, (length(RNAswapfemale.data$gene.id)-length(MSfemaleswap$gene.id))-G23FupsigNotMS))
#Add names to the columns and rows
colnames(MS.G23FupF) <- c("MS", "NotMS")
rownames(MS.G23FupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G23FupF)
# P = 0.3264

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many mito-sensitive transcripts are there?
G23FdownsigMS <- sum(G23Fdownsig$gene.id %in% MSfemaleswap$gene.id)
# 0
# How many are not?
G23FdownsigNotMS <- length(G23Fdownsig$gene.id)-G23FdownsigMS
# 6

# The number of expected mito-sensitive transcripts
length(G23Fdownsig$gene.id)*MSpropF
# 0.4364611

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FdownsigMS, G23FdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.4927

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G23FdownF <- rbind(c(G23FdownsigMS, G23FdownsigNotMS),
                      c(length(MSfemaleswap$gene.id)-G23FdownsigMS, (length(RNAswapfemale.data$gene.id)-length(MSfemaleswap$gene.id))-G23FdownsigNotMS))
#Add names to the columns and rows
colnames(MS.G23FdownF) <- c("MS", "NotMS")
rownames(MS.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G23FdownF)
# P = 1

# We get the same results, so that is good


