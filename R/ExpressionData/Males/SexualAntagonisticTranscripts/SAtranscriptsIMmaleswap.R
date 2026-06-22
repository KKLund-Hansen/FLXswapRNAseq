################################################################################################
########### SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP MALES ##########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/SexualAntagonisticTranscripts")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the male data set
IMsaM.data <- subset(IMsa.data, IMsa.data$gene.id %in% RNAswapmale.data$gene.id)
# 1353

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMM <- length(IMsaM.data$gene.id)/length(RNAswapmale.data$gene.id)
# 0.08
# And for the chisq.test I also have to calculate
NotSApropIMM <- 1-SApropIMM
# 0.92



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many transcripts are sexual antagonistic?
GAMupsigSAIM <- sum(GAMupsig$gene.id %in% IMsaM.data$gene.id)
# 1
# How many are not?
GAMupsigNotSAIM <- length(GAMupsig$gene.id)-GAMupsigSAIM
# 3

# The number of expected transcripts which are SA
length(GAMupsig$gene.id)*SApropIMM
# 0.3224307

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMupsigSAIM, GAMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2133

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GAMupF <- rbind(c(GAMupsigSAIM, GAMupsigNotSAIM),
                      c(length(IMsaM.data$gene.id)-GAMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GAMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GAMupF) <- c("SA", "NotSA")
rownames(SAIM.GAMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GAMupF)
# P = 0.2855

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GAMdownsigSAIM <- sum(GAMdownsig$gene.id %in% IMsaM.data$gene.id)
# 5
# How many are not?
GAMdownsigNotSAIM <- length(GAMdownsig$gene.id)-GAMdownsigSAIM
# 11

# The number of expected transcripts which are SA
length(GAMdownsig$gene.id)*SApropIMM
# 1.289723

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMdownsigSAIM, GAMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.0006561

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GAMdownF <- rbind(c(GAMdownsigSAIM, GAMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GAMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GAMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GAMdownF) <- c("SA", "NotSA")
rownames(SAIM.GAMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GAMdownF)
# P = 0.006947

# It went from significant to not




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many transcripts are sexual antagonistic?
GBMupsigSAIM <- sum(GBMupsig$gene.id %in% IMsaM.data$gene.id)
# 9
# How many are not?
GBMupsigNotSAIM <- length(GBMupsig$gene.id)-GBMupsigSAIM
# 45

# The number of expected transcripts which are SA
length(GBMupsig$gene.id)*SApropIMM
# 4.352815

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMupsigSAIM, GBMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.02018

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GBMupF <- rbind(c(GBMupsigSAIM, GBMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GBMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GBMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GBMupF) <- c("SA", "NotSA")
rownames(SAIM.GBMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GBMupF)
# P = 0.03812

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many transcripts are sexual antagonistic?
GBMdownsigSAIM <- sum(GBMdownsig$gene.id %in% IMsaM.data$gene.id)
# 12
# How many are not?
GBMdownsigNotSAIM <- length(GBMdownsig$gene.id)-GBMdownsigSAIM
# 46

# The number of expected transcripts which are SA
length(GBMdownsig$gene.id)*SApropIMM
# 4.675246

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMdownsigSAIM, GBMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.0004109

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GBMdownF <- rbind(c(GBMdownsigSAIM, GBMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GBMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GBMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GBMdownF) <- c("SA", "NotSA")
rownames(SAIM.GBMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GBMdownF)
# P = 0.001966

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many transcripts are sexual antagonistic?
GCMupsigSAIM <- sum(GCMupsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GCMupsigNotSAIM <- length(GCMupsig$gene.id)-GCMupsigSAIM
# 23

# The number of expected transcripts which are SA
length(GCMupsig$gene.id)*SApropIMM
# 2.176408

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMupsigSAIM, GCMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.1973

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GCMupF <- rbind(c(GCMupsigSAIM, GCMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GCMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GCMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GCMupF) <- c("SA", "NotSA")
rownames(SAIM.GCMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GCMupF)
# P = 0.2721

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many transcripts are sexual antagonistic?
GCMdownsigSAIM <- sum(GCMdownsig$gene.id %in% IMsaM.data$gene.id)
# 2
# How many are not?
GCMdownsigNotSAIM <- length(GCMdownsig$gene.id)-GCMdownsigSAIM
# 35

# The number of expected transcripts which are SA
length(GCMdownsig$gene.id)*SApropIMM
# 2.982484

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMdownsigSAIM, GCMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.553

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GCMdownF <- rbind(c(GCMdownsigSAIM, GCMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GCMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GCMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GCMdownF) <- c("SA", "NotSA")
rownames(SAIM.GCMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GCMdownF)
# P = 0.7661

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many transcripts are sexual antagonistic?
GDMupsigSAIM <- sum(GDMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GDMupsigNotSAIM <- length(GDMupsig$gene.id)-GDMupsigSAIM
# 14

# The number of expected transcripts which are SA
length(GDMupsig$gene.id)*SApropIMM
# 1.128508

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMupsigSAIM, GDMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2679

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GDMupF <- rbind(c(GDMupsigSAIM, GDMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GDMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GDMupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GDMupF) <- c("SA", "NotSA")
rownames(SAIM.GDMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GDMupF)
# P = 0.6214

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GDMdownsigSAIM <- sum(GDMdownsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GDMdownsigNotSAIM <- length(GDMdownsig$gene.id)-GDMdownsigSAIM
# 12

# The number of expected transcripts which are SA
length(GDMdownsig$gene.id)*SApropIMM
# 1.289723

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMdownsigSAIM, GDMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.01281

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GDMdownF <- rbind(c(GDMdownsigSAIM, GDMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GDMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GDMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GDMdownF) <- c("SA", "NotSA")
rownames(SAIM.GDMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GDMdownF)
# P = 0.03493

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many transcripts are sexual antagonistic?
GEMupsigSAIM <- sum(GEMupsig$gene.id %in% IMsaM.data$gene.id)
# 1
# How many are not?
GEMupsigNotSAIM <- length(GEMupsig$gene.id)-GEMupsigSAIM
# 16

# The number of expected transcripts which are SA
length(GEMupsig$gene.id)*SApropIMM
# 1.370331

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMupsigSAIM, GEMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.7415

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GEMupF <- rbind(c(GEMupsigSAIM, GEMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GEMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GEMupsigNotSAIM))
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
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many transcripts are sexual antagonistic?
GEMdownsigSAIM <- sum(GEMdownsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
GEMdownsigNotSAIM <- length(GEMdownsig$gene.id)-GEMdownsigSAIM
# 22

# The number of expected transcripts which are SA
length(GEMdownsig$gene.id)*SApropIMM
# 2.0958

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMdownsigSAIM, GEMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.1701

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GEMdownF <- rbind(c(GEMdownsigSAIM, GEMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GEMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GEMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GEMdownF) <- c("SA", "NotSA")
rownames(SAIM.GEMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GEMdownF)
# P = 0.153

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many transcripts are sexual antagonistic?
GFMupsigSAIM <- sum(GFMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GFMupsigNotSAIM <- length(GFMupsig$gene.id)-GFMupsigSAIM
# 2

# The number of expected transcripts which are SA
length(GFMupsig$gene.id)*SApropIMM
# 0.1612154

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMupsigSAIM, GFMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.6754

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFMupF <- rbind(c(GFMupsigSAIM, GFMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GFMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GFMupsigNotSAIM))
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
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many transcripts are sexual antagonistic?
GFMdownsigSAIM <- sum(GFMdownsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GFMdownsigNotSAIM <- length(GFMdownsig$gene.id)-GFMdownsigSAIM
# 3

# The number of expected transcripts which are SA
length(GFMdownsig$gene.id)*SApropIMM
# 0.2418231

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMdownsigSAIM, GFMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.608

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFMdownF <- rbind(c(GFMdownsigSAIM, GFMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GFMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GFMdownsigNotSAIM))
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
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many transcripts are sexual antagonistic?
GGMupsigSAIM <- sum(GGMupsig$gene.id %in% IMsaM.data$gene.id)
# 0
# How many are not?
GGMupsigNotSAIM <- length(GGMupsig$gene.id)-GGMupsigSAIM
# 10

# The number of expected transcripts which are SA
length(GGMupsig$gene.id)*SApropIMM
# 0.8060769

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMupsigSAIM, GGMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.3491

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GGMupF <- rbind(c(GGMupsigSAIM, GGMupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-GGMupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GGMupsigNotSAIM))
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
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many transcripts are sexual antagonistic?
GGMdownsigSAIM <- sum(GGMdownsig$gene.id %in% IMsaM.data$gene.id)
# 3
# How many are not?
GGMdownsigNotSAIM <- length(GGMdownsig$gene.id)-GGMdownsigSAIM
# 16

# The number of expected transcripts which are SA
length(GGMdownsig$gene.id)*SApropIMM
# 1.531546

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMdownsigSAIM, GGMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.2159

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GGMdownF <- rbind(c(GGMdownsigSAIM, GGMdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-GGMdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-GGMdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GGMdownF) <- c("SA", "NotSA")
rownames(SAIM.GGMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GGMdownF)
# P = 0.1937

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many transcripts are sexual antagonistic?
G2MupsigSAIM <- sum(G2Mupsig$gene.id %in% IMsaM.data$gene.id)
# 6
# How many are not?
G2MupsigNotSAIM <- length(G2Mupsig$gene.id)-G2MupsigSAIM
# 23

# The number of expected transcripts which are SA
length(G2Mupsig$gene.id)*SApropIMM
# 2.337623

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G2MupsigSAIM, G2MupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.01248

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G2MupF <- rbind(c(G2MupsigSAIM, G2MupsigNotSAIM),
                     c(length(IMsaM.data$gene.id)-G2MupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-G2MupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G2MupF) <- c("SA", "NotSA")
rownames(SAIM.G2MupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G2MupF)
# P = 0.02581

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many transcripts are sexual antagonistic?
G2MdownsigSAIM <- sum(G2Mdownsig$gene.id %in% IMsaM.data$gene.id)
# 5
# How many are not?
G2MdownsigNotSAIM <- length(G2Mdownsig$gene.id)-G2MdownsigSAIM
# 31

# The number of expected transcripts which are SA
length(G2Mdownsig$gene.id)*SApropIMM
# 2.901877

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G2MdownsigSAIM, G2MdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.199

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G2MdownF <- rbind(c(G2MdownsigSAIM, G2MdownsigNotSAIM),
                       c(length(IMsaM.data$gene.id)-G2MdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-G2MdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G2MdownF) <- c("SA", "NotSA")
rownames(SAIM.G2MdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G2MdownF)
# P = 0.209

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many transcripts are sexual antagonistic?
G11MupsigSAIM <- sum(G11Mupsig$gene.id %in% IMsaM.data$gene.id)
# 4
# How many are not?
G11MupsigNotSAIM <- length(G11Mupsig$gene.id)-G11MupsigSAIM
# 17

# The number of expected transcripts which are SA
length(G11Mupsig$gene.id)*SApropIMM
# 1.692761

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G11MupsigSAIM, G11MupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.06439

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G11MupF <- rbind(c(G11MupsigSAIM, G11MupsigNotSAIM),
                      c(length(IMsaM.data$gene.id)-G11MupsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-G11MupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G11MupF) <- c("SA", "NotSA")
rownames(SAIM.G11MupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G11MupF)
# P = 0.08362

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many transcripts are sexual antagonistic?
G11MdownsigSAIM <- sum(G11Mdownsig$gene.id %in% IMsaM.data$gene.id)
# 2
# How many are not?
G11MdownsigNotSAIM <- length(G11Mdownsig$gene.id)-G11MdownsigSAIM
# 5

# The number of expected transcripts which are SA
length(G11Mdownsig$gene.id)*SApropIMM
# 0.5642538

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G11MdownsigSAIM, G11MdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.04622

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G11MdownF <- rbind(c(G11MdownsigSAIM, G11MdownsigNotSAIM),
                        c(length(IMsaM.data$gene.id)-G11MdownsigSAIM, (length(RNAswapmale.data$gene.id)-length(IMsaM.data$gene.id))-G11MdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G11MdownF) <- c("SA", "NotSA")
rownames(SAIM.G11MdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G11MdownF)
# P = 0.1039

# We get the same results, so that is good

