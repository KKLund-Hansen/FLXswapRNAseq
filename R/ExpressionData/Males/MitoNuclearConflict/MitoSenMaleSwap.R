################################################################################################
########################### MITO-SENSITIVE TRANSCRIPTS FLX SWAP MALES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the male data set
MSmaleswap <- subset(MS.data, MS.data$gene.id %in% RNAswapmale.data$gene.id)
# 1221 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropM <- length(MSmaleswap$gene.id)/length(RNAswapmale.data$gene.id)
# 0.07

# And for the chisq.test I also have to calculate
NotMSpropM <- 1-MSpropM
# 0.93




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many mito-sensitive transcripts are there?
GAMupsigMS <- sum(GAMupsig$gene.id %in% MSmaleswap$gene.id)
# 0
# How many are not?
GAMupsigNotMS <- length(GAMupsig$gene.id)-GAMupsigMS
# 4

# The number of expected mito-sensitive transcripts
length(GAMupsig$gene.id)*MSpropM
# 0.2909741

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMupsigMS, GAMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5754

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GAMupF <- rbind(c(GAMupsigMS, GAMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GAMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GAMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GAMupF) <- c("MS", "NotMS")
rownames(MS.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many mito-sensitive transcripts are there?
GAMdownsigMS <- sum(GAMdownsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
GAMdownsigNotMS <- length(GAMdownsig$gene.id)-GAMdownsigMS
# 14

# The number of expected mito-sensitive transcripts
length(GAMdownsig$gene.id)*MSpropM
# 1.163896

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigMS, GAMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4209

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GAMdownF <- rbind(c(GAMdownsigMS, GAMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GAMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GAMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GAMdownF) <- c("MS", "NotMS")
rownames(MS.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GAMdownF)
# P = 0.3264

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many mito-sensitive transcripts are there?
GBMupsigMS <- sum(GBMupsig$gene.id %in% MSmaleswap$gene.id)
# 4
# How many are not?
GBMupsigNotMS <- length(GBMupsig$gene.id)-GBMupsigMS
# 50

# The number of expected mito-sensitive transcripts
length(GBMupsig$gene.id)*MSpropM
# 3.92815

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMupsigMS, GBMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.97

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GBMupF <- rbind(c(GBMupsigMS, GBMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GBMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GBMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GBMupF) <- c("MS", "NotMS")
rownames(MS.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GBMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many mito-sensitive transcripts are there?
GBMdownsigMS <- sum(GBMdownsig$gene.id %in% MSmaleswap$gene.id)
# 3
# How many are not?
GBMdownsigNotMS <- length(GBMdownsig$gene.id)-GBMdownsigMS
# 55

# The number of expected mito-sensitive transcripts
length(GBMdownsig$gene.id)*MSpropM
# 4.219124

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigMS, GBMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5377

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GBMdownF <- rbind(c(GBMdownsigMS, GBMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GBMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GBMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GBMdownF) <- c("MS", "NotMS")
rownames(MS.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GBMdownF)
# P = 0.7985

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many mito-sensitive transcripts are there?
GCMupsigMS <- sum(GCMupsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
GCMupsigNotMS <- length(GCMupsig$gene.id)-GCMupsigMS
# 25

# The number of expected mito-sensitive transcripts
length(GCMupsig$gene.id)*MSpropM
# 1.964075

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMupsigMS, GCMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9788

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GCMupF <- rbind(c(GCMupsigMS, GCMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GCMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GCMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GCMupF) <- c("MS", "NotMS")
rownames(MS.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GCMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many mito-sensitive transcripts are there?
GCMdownsigMS <- sum(GCMdownsig$gene.id %in% MSmaleswap$gene.id)
# 3
# How many are not?
GCMdownsigNotMS <- length(GCMdownsig$gene.id)-GCMdownsigMS
# 34

# The number of expected mito-sensitive transcripts
length(GCMdownsig$gene.id)*MSpropM
# 2.69151

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigMS, GCMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.8452

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GCMdownF <- rbind(c(GCMdownsigMS, GCMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GCMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GCMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GCMdownF) <- c("MS", "NotMS")
rownames(MS.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GCMdownF)
# P = 0.7492

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many mito-sensitive transcripts are there?
GDMupsigMS <- sum(GDMupsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
GDMupsigNotMS <- length(GDMupsig$gene.id)-GDMupsigMS
# 13

# The number of expected mito-sensitive transcripts
length(GDMupsig$gene.id)*MSpropM
# 1.018409

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMupsigMS, GDMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9849

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GDMupF <- rbind(c(GDMupsigMS, GDMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GDMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GDMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GDMupF) <- c("MS", "NotMS")
rownames(MS.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many mito-sensitive transcripts are there?
GDMdownsigMS <- sum(GDMdownsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
GDMdownsigNotMS <- length(GDMdownsig$gene.id)-GDMdownsigMS
# 14

# The number of expected mito-sensitive transcripts
length(GDMdownsig$gene.id)*MSpropM
# 1.163896

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigMS, GDMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4209

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GDMdownF <- rbind(c(GDMdownsigMS, GDMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GDMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GDMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GDMdownF) <- c("MS", "NotMS")
rownames(MS.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GDMdownF)
# P = 0.3264

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many mito-sensitive transcripts are there?
GEMupsigMS <- sum(GEMupsig$gene.id %in% MSmaleswap$gene.id)
# 0
# How many are not?
GEMupsigNotMS <- length(GEMupsig$gene.id)-GEMupsigMS
# 17

# The number of expected mito-sensitive transcripts
length(GEMupsig$gene.id)*MSpropM
# 1.23664

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMupsigMS, GEMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.2482

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GEMupF <- rbind(c(GEMupsigMS, GEMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GEMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GEMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GEMupF) <- c("MS", "NotMS")
rownames(MS.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GEMupF)
# P = 0.6305

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many mito-sensitive transcripts are there?
GEMdownsigMS <- sum(GEMdownsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
GEMdownsigNotMS <- length(GEMdownsig$gene.id)-GEMdownsigMS
# 25

# The number of expected mito-sensitive transcripts
length(GEMdownsig$gene.id)*MSpropM
# 1.891332

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigMS, GEMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5009

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GEMdownF <- rbind(c(GEMdownsigMS, GEMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GEMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GEMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GEMdownF) <- c("MS", "NotMS")
rownames(MS.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GEMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many mito-sensitive transcripts are there?
GFMupsigMS <- sum(GFMupsig$gene.id %in% MSmaleswap$gene.id)
# 0
# How many are not?
GFMupsigNotMS <- length(GFMupsig$gene.id)-GFMupsigMS
# 2

# The number of expected mito-sensitive transcripts
length(GFMupsig$gene.id)*MSpropM
# 0.145487

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMupsigMS, GFMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.692

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFMupF <- rbind(c(GFMupsigMS, GFMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GFMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GFMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GFMupF) <- c("MS", "NotMS")
rownames(MS.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many mito-sensitive transcripts are there?
GFMdownsigMS <- sum(GFMdownsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
GFMdownsigNotMS <- length(GFMdownsig$gene.id)-GFMdownsigMS
# 1

# The number of expected mito-sensitive transcripts
length(GFMdownsig$gene.id)*MSpropM
# 0.2182306

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigMS, GFMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 7.467e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFMdownF <- rbind(c(GFMdownsigMS, GFMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GFMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GFMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GFMdownF) <- c("MS", "NotMS")
rownames(MS.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFMdownF)
# P = 0.01509

# It goes from significant to non-significant




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many mito-sensitive transcripts are there?
GGMupsigMS <- sum(GGMupsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
GGMupsigNotMS <- length(GGMupsig$gene.id)-GGMupsigMS
# 9

# The number of expected mito-sensitive transcripts
length(GGMupsig$gene.id)*MSpropM
# 0.7274352

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMupsigMS, GGMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.74

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GGMupF <- rbind(c(GGMupsigMS, GGMupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-GGMupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GGMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GGMupF) <- c("MS", "NotMS")
rownames(MS.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GGMupF)
# P = 0.5302

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many mito-sensitive transcripts are there?
GGMdownsigMS <- sum(GGMdownsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
GGMdownsigNotMS <- length(GGMdownsig$gene.id)-GGMdownsigMS
# 18

# The number of expected mito-sensitive transcripts
length(GGMdownsig$gene.id)*MSpropM
# 1.382127

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigMS, GGMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.7357

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GGMdownF <- rbind(c(GGMdownsigMS, GGMdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-GGMdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-GGMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GGMdownF) <- c("MS", "NotMS")
rownames(MS.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GGMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many mito-sensitive transcripts are there?
G2MupsigMS <- sum(G2Mupsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
G2MupsigNotMS <- length(G2Mupsig$gene.id)-G2MupsigMS
# 27

# The number of expected mito-sensitive transcripts
length(G2Mupsig$gene.id)*MSpropM
# 2.109562

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MupsigMS, G2MupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9376

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G2MupF <- rbind(c(G2MupsigMS, G2MupsigNotMS), 
                   c(length(MSmaleswap$gene.id)-G2MupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-G2MupsigNotMS))
#Add names to the columns and rows
colnames(MS.G2MupF) <- c("MS", "NotMS")
rownames(MS.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G2MupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many mito-sensitive transcripts are there?
G2MdownsigMS <- sum(G2Mdownsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
G2MdownsigNotMS <- length(G2Mdownsig$gene.id)-G2MdownsigMS
# 34

# The number of expected mito-sensitive transcripts
length(G2Mdownsig$gene.id)*MSpropM
# 2.618767

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MdownsigMS, G2MdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.6913

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G2MdownF <- rbind(c(G2MdownsigMS, G2MdownsigNotMS), 
                     c(length(MSmaleswap$gene.id)-G2MdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-G2MdownsigNotMS))
#Add names to the columns and rows
colnames(MS.G2MdownF) <- c("MS", "NotMS")
rownames(MS.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G2MdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many mito-sensitive transcripts are there?
G11MupsigMS <- sum(G11Mupsig$gene.id %in% MSmaleswap$gene.id)
# 0
# How many are not?
G11MupsigNotMS <- length(G11Mupsig$gene.id)-G11MupsigMS
# 21

# The number of expected mito-sensitive transcripts
length(G11Mupsig$gene.id)*MSpropM
# 1.527614

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MupsigMS, G11MupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.1993

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G11MupF <- rbind(c(G11MupsigMS, G11MupsigNotMS), 
                    c(length(MSmaleswap$gene.id)-G11MupsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-G11MupsigNotMS))
#Add names to the columns and rows
colnames(MS.G11MupF) <- c("MS", "NotMS")
rownames(MS.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G11MupF)
# P = 0.3978

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many mito-sensitive transcripts are there?
G11MdownsigMS <- sum(G11Mdownsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
G11MdownsigNotMS <- length(G11Mdownsig$gene.id)-G11MdownsigMS
# 6

# The number of expected mito-sensitive transcripts
length(G11Mdownsig$gene.id)*MSpropM
# 0.5092046

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MdownsigMS, G11MdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4751

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.G11MdownF <- rbind(c(G11MdownsigMS, G11MdownsigNotMS), 
                      c(length(MSmaleswap$gene.id)-G11MdownsigMS, (length(RNAswapmale.data$gene.id)-length(MSmaleswap$gene.id))-G11MdownsigNotMS))
#Add names to the columns and rows
colnames(MS.G11MdownF) <- c("MS", "NotMS")
rownames(MS.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.G11MdownF)
# P = 0.4107

# We get the same results, so that is good

