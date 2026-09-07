################################################################################################
##### SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP MALES - 7 GROUPS #####
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexualAntagonisticTranscripts")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the male data set
IMsaM.data <- subset(IMsa.data, IMsa.data$gene.id %in% FLXswap7gM.data$gene.id)
# 1337

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMM <- length(IMsaM.data$gene.id)/length(FLXswap7gM.data$gene.id)
# 0.09
# And for the chisq.test I also have to calculate
NotSApropIMM <- 1-SApropIMM
# 0.91



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many transcripts are sexual antagonistic?
GAMupsigSAIM <- sum(GAMupsig$gene.id %in% IMsaM.data$gene.id)
# 1
# How many are not?
GAMupsigNotSAIM <- length(GAMupsig$gene.id)-GAMupsigSAIM
# 4

# The number of expected transcripts which are SA
length(GAMupsig$gene.id)*SApropIMM
# 0.4311512

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMupsigSAIM, GAMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.3648

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GAMupF <- rbind(c(GAMupsigSAIM, GAMupsigNotSAIM),
                      c(length(IMsaM.data$gene.id)-GAMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GAMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GAMupF) <- c("SA", "NotSA")
rownames(SAIM.GAMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GAMupF)
# P = 0.363

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GAMdownsigSAIM <- sum(GAMdownsig$gene.id %in% IMsaM.data$gene.id)
# 5
# How many are not?
GAMdownsigNotSAIM <- length(GAMdownsig$gene.id)-GAMdownsigSAIM
# 11

# The number of expected transcripts which are SA
length(GAMdownsig$gene.id)*SApropIMM
# 1.379684

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMdownsigSAIM, GAMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.001263

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GAMdownF <- rbind(c(GAMdownsigSAIM, GAMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GAMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GAMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GAMdownF) <- c("SA", "NotSA")
rownames(SAIM.GAMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GAMdownF)
# P = 0.009221

# It went from significant to not




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many transcripts are sexual antagonistic?
GBMupsigSAIM <- sum(GBMupsig$gene.id %in% IMsaM.data$gene.id)
# 8
# How many are not?
GBMupsigNotSAIM <- length(GBMupsig$gene.id)-GBMupsigSAIM
# 43

# The number of expected transcripts which are SA
length(GBMupsig$gene.id)*SApropIMM
# 4.397743

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMupsigSAIM, GBMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.07234

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GBMupF <- rbind(c(GBMupsigSAIM, GBMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GBMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GBMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GBMupF) <- c("SA", "NotSA")
rownames(SAIM.GBMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GBMupF)
# P = 0.07895

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many transcripts are sexual antagonistic?
GBMdownsigSAIM <- sum(GBMdownsig$gene.id %in% IMsaM.data$gene.id)
# 13
# How many are not?
GBMdownsigNotSAIM <- length(GBMdownsig$gene.id)-GBMdownsigSAIM
# 43

# The number of expected transcripts which are SA
length(GBMdownsig$gene.id)*SApropIMM
# 4.828894

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMdownsigSAIM, GBMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.0001003

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GBMdownF <- rbind(c(GBMdownsigSAIM, GBMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GBMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GBMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GBMdownF) <- c("SA", "NotSA")
rownames(SAIM.GBMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GBMdownF)
# P = 0.0007726

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many transcripts are sexual antagonistic?
GCMupsigSAIM <- sum(GCMupsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GCMupsigNotSAIM <- length(GCMupsig$gene.id)-GCMupsigSAIM
# 23

# The number of expected transcripts which are SA
length(GCMupsig$gene.id)*SApropIMM
# 2.328217

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMupsigSAIM, GCMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2517

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GCMupF <- rbind(c(GCMupsigSAIM, GCMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GCMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GCMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GCMupF) <- c("SA", "NotSA")
rownames(SAIM.GCMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GCMupF)
# P = 0.2872

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many transcripts are sexual antagonistic?
GCMdownsigSAIM <- sum(GCMdownsig$gene.id %in% IMsaM.data$gene.id)
# 2
# How many are not?
GCMdownsigNotSAIM <- length(GCMdownsig$gene.id)-GCMdownsigSAIM
# 36

# The number of expected transcripts which are SA
length(GCMdownsig$gene.id)*SApropIMM
# 3.276749

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMdownsigSAIM, GCMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.4606

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GCMdownF <- rbind(c(GCMdownsigSAIM, GCMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GCMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GCMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GCMdownF) <- c("SA", "NotSA")
rownames(SAIM.GCMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GCMdownF)
# P = 0.7694

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many transcripts are sexual antagonistic?
GDMupsigSAIM <- sum(GDMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GDMupsigNotSAIM <- length(GDMupsig$gene.id)-GDMupsigSAIM
# 12

# The number of expected transcripts which are SA
length(GDMupsig$gene.id)*SApropIMM
# 1.034763

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMupsigSAIM, GDMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2873

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GDMupF <- rbind(c(GDMupsigSAIM, GDMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GDMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GDMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GDMupF) <- c("SA", "NotSA")
rownames(SAIM.GDMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GDMupF)
# P = 0.6161

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GDMdownsigSAIM <- sum(GDMdownsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GDMdownsigNotSAIM <- length(GDMdownsig$gene.id)-GDMdownsigSAIM
# 12

# The number of expected transcripts which are SA
length(GDMdownsig$gene.id)*SApropIMM
# 1.379684

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMdownsigSAIM, GDMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.01961

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GDMdownF <- rbind(c(GDMdownsigSAIM, GDMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GDMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GDMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GDMdownF) <- c("SA", "NotSA")
rownames(SAIM.GDMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GDMdownF)
# P = 0.04327

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many transcripts are sexual antagonistic?
GEMupsigSAIM <- sum(GEMupsig$gene.id %in% IMsaM.data$gene.id)
# 1
# How many are not?
GEMupsigNotSAIM <- length(GEMupsig$gene.id)-GEMupsigSAIM
# 19

# The number of expected transcripts which are SA
length(GEMupsig$gene.id)*SApropIMM
# 1.724605

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMupsigSAIM, GEMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.5638

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GEMupF <- rbind(c(GEMupsigSAIM, GEMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GEMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GEMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GEMupF) <- c("SA", "NotSA")
rownames(SAIM.GEMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GEMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many transcripts are sexual antagonistic?
GEMdownsigSAIM <- sum(GEMdownsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GEMdownsigNotSAIM <- length(GEMdownsig$gene.id)-GEMdownsigSAIM
# 22

# The number of expected transcripts which are SA
length(GEMdownsig$gene.id)*SApropIMM
# 2.241986

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMdownsigSAIM, GEMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2194

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GEMdownF <- rbind(c(GEMdownsigSAIM, GEMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GEMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GEMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GEMdownF) <- c("SA", "NotSA")
rownames(SAIM.GEMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GEMdownF)
# P = 0.2774

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many transcripts are sexual antagonistic?
GFMupsigSAIM <- sum(GFMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GFMupsigNotSAIM <- length(GFMupsig$gene.id)-GFMupsigSAIM
# 2

# The number of expected transcripts which are SA
length(GFMupsig$gene.id)*SApropIMM
# 0.1724605

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMupsigSAIM, GFMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.664

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFMupF <- rbind(c(GFMupsigSAIM, GFMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GFMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GFMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GFMupF) <- c("SA", "NotSA")
rownames(SAIM.GFMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many transcripts are sexual antagonistic?
GFMdownsigSAIM <- sum(GFMdownsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GFMdownsigNotSAIM <- length(GFMdownsig$gene.id)-GFMdownsigSAIM
# 3

# The number of expected transcripts which are SA
length(GFMdownsig$gene.id)*SApropIMM
# 0.2586907

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMdownsigSAIM, GFMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.5947

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFMdownF <- rbind(c(GFMdownsigSAIM, GFMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GFMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GFMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GFMdownF) <- c("SA", "NotSA")
rownames(SAIM.GFMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GFMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many transcripts are sexual antagonistic?
GGMupsigSAIM <- sum(GGMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GGMupsigNotSAIM <- length(GGMupsig$gene.id)-GGMupsigSAIM
# 10

# The number of expected transcripts which are SA
length(GGMupsig$gene.id)*SApropIMM
# 0.8623025

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMupsigSAIM, GGMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.3313

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GGMupF <- rbind(c(GGMupsigSAIM, GGMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GGMupsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GGMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GGMupF) <- c("SA", "NotSA")
rownames(SAIM.GGMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GGMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many transcripts are sexual antagonistic?
GGMdownsigSAIM <- sum(GGMdownsig$gene.id %in% IMsaM.data$gene.id)
# 3
# How many are not?
GGMdownsigNotSAIM <- length(GGMdownsig$gene.id)-GGMdownsigSAIM
# 19

# The number of expected transcripts which are SA
length(GGMdownsig$gene.id)*SApropIMM
# 1.897065

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMdownsigSAIM, GGMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.4022

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GGMdownF <- rbind(c(GGMdownsigSAIM, GGMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GGMdownsigSAIM, (length(FLXswap7gM.data$gene.id)-length(IMsaM.data$gene.id))-GGMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GGMdownF) <- c("SA", "NotSA")
rownames(SAIM.GGMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GGMdownF)
# P = 0.4313

# We get the same results, so that is good

