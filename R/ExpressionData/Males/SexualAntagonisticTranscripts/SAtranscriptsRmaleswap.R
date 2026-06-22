################################################################################################
############## SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP MALES #############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/SexualAntagonisticTranscripts")

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the male data set
RsaM.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% RNAswapmale.data$gene.id)
# 490

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRM <- length(RsaM.data$FBgene.id)/length(RNAswapmale.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRM <- 1-SApropRM
# 0.97


######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many transcripts are sexual antagonistic?
GAMupsigSAR <- sum(GAMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GAMupsigNotSAR <- length(GAMupsig$gene.id)-GAMupsigSAR
# 4

# The number of expected transcripts which are SA
length(GAMupsig$gene.id)*SApropRM
# 0.1167709

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMupsigSAR, GAMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.7287

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAMupF <- rbind(c(GAMupsigSAR, GAMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GAMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GAMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GAMupF) <- c("SA", "NotSA")
rownames(SAR.GAMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GAMdownsigSAR <- sum(GAMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GAMdownsigNotSAR <- length(GAMdownsig$gene.id)-GAMdownsigSAR
# 15

# The number of expected transcripts which are SA
length(GAMdownsig$gene.id)*SApropRM
# 0.4670837

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMdownsigSAR, GAMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4287

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAMdownF <- rbind(c(GAMdownsigSAR, GAMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GAMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GAMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GAMdownF) <- c("SA", "NotSA")
rownames(SAR.GAMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GAMdownF)
# P = 0.3776

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many transcripts are sexual antagonistic?
GBMupsigSAR <- sum(GBMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GBMupsigNotSAR <- length(GBMupsig$gene.id)-GBMupsigSAR
# 54

# The number of expected transcripts which are SA
length(GBMupsig$gene.id)*SApropRM
# 1.576408

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMupsigSAR, GBMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.2026

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GBMupF <- rbind(c(GBMupsigSAR, GBMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GBMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GBMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GBMupF) <- c("SA", "NotSA")
rownames(SAR.GBMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GBMupF)
# P = 0.4101

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many transcripts are sexual antagonistic?
GBMdownsigSAR <- sum(GBMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 4
# How many are not?
GBMdownsigNotSAR <- length(GBMdownsig$gene.id)-GBMdownsigSAR
# 54

# The number of expected transcripts which are SA
length(GBMdownsig$gene.id)*SApropRM
# 1.693178

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMdownsigSAR, GBMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.07198

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GBMdownF <- rbind(c(GBMdownsigSAR, GBMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GBMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GBMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GBMdownF) <- c("SA", "NotSA")
rownames(SAR.GBMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GBMdownF)
# P = 0.08893

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many transcripts are sexual antagonistic?
GCMupsigSAR <- sum(GCMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GCMupsigNotSAR <- length(GCMupsig$gene.id)-GCMupsigSAR
# 26

# The number of expected transcripts which are SA
length(GCMupsig$gene.id)*SApropRM
# 0.7882038

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMupsigSAR, GCMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.8087

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GCMupF <- rbind(c(GCMupsigSAR, GCMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GCMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GCMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GCMupF) <- c("SA", "NotSA")
rownames(SAR.GCMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GCMupF)
# P = 0.5509

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many transcripts are sexual antagonistic?
GCMdownsigSAR <- sum(GCMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 3
# How many are not?
GCMdownsigNotSAR <- length(GCMdownsig$gene.id)-GCMdownsigSAR
# 34

# The number of expected transcripts which are SA
length(GCMdownsig$gene.id)*SApropRM
# 1.080131

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMdownsigSAR, GCMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.06081

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GCMdownF <- rbind(c(GCMdownsigSAR, GCMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GCMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GCMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GCMdownF) <- c("SA", "NotSA")
rownames(SAR.GCMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GCMdownF)
# P = 0.09267

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many transcripts are sexual antagonistic?
GDMupsigSAR <- sum(GDMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GDMupsigNotSAR <- length(GDMupsig$gene.id)-GDMupsigSAR
# 13

# The number of expected transcripts which are SA
length(GDMupsig$gene.id)*SApropRM
# 0.4086982

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMupsigSAR, GDMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.3479

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GDMupF <- rbind(c(GDMupsigSAR, GDMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GDMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GDMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GDMupF) <- c("SA", "NotSA")
rownames(SAR.GDMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GDMupF)
# P = 0.3396

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GDMdownsigSAR <- sum(GDMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GDMdownsigNotSAR <- length(GDMdownsig$gene.id)-GDMdownsigSAR
# 16

# The number of expected transcripts which are SA
length(GDMdownsig$gene.id)*SApropRM
# 0.4670837

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMdownsigSAR, GDMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4879

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GDMdownF <- rbind(c(GDMdownsigSAR, GDMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GDMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GDMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GDMdownF) <- c("SA", "NotSA")
rownames(SAR.GDMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GDMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many transcripts are sexual antagonistic?
GEMupsigSAR <- sum(GEMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GEMupsigNotSAR <- length(GEMupsig$gene.id)-GEMupsigSAR
# 17

# The number of expected transcripts which are SA
length(GEMupsig$gene.id)*SApropRM
# 0.4962764

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMupsigSAR, GEMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4746

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEMupF <- rbind(c(GEMupsigSAR, GEMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GEMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GEMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GEMupF) <- c("SA", "NotSA")
rownames(SAR.GEMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GEMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many transcripts are sexual antagonistic?
GEMdownsigSAR <- sum(GEMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GEMdownsigNotSAR <- length(GEMdownsig$gene.id)-GEMdownsigSAR
# 26

# The number of expected transcripts which are SA
length(GEMdownsig$gene.id)*SApropRM
# 0.759011

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMdownsigSAR, GEMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.3766

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEMdownF <- rbind(c(GEMdownsigSAR, GEMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GEMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GEMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GEMdownF) <- c("SA", "NotSA")
rownames(SAR.GEMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GEMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many transcripts are sexual antagonistic?
GFMupsigSAR <- sum(GFMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GFMupsigNotSAR <- length(GFMupsig$gene.id)-GFMupsigSAR
# 2

# The number of expected transcripts which are SA
length(GFMupsig$gene.id)*SApropRM
# 0.05838546

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMupsigSAR, GFMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.8063

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFMupF <- rbind(c(GFMupsigSAR, GFMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GFMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GFMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GFMupF) <- c("SA", "NotSA")
rownames(SAR.GFMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many transcripts are sexual antagonistic?
GFMdownsigSAR <- sum(GFMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GFMdownsigNotSAR <- length(GFMdownsig$gene.id)-GFMdownsigSAR
# 3

# The number of expected transcripts which are SA
length(GFMdownsig$gene.id)*SApropRM
# 0.08757819

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMdownsigSAR, GFMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.7639

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFMdownF <- rbind(c(GFMdownsigSAR, GFMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GFMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GFMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GFMdownF) <- c("SA", "NotSA")
rownames(SAR.GFMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GFMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many transcripts are sexual antagonistic?
GGMupsigSAR <- sum(GGMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GGMupsigNotSAR <- length(GGMupsig$gene.id)-GGMupsigSAR
# 9

# The number of expected transcripts which are SA
length(GGMupsig$gene.id)*SApropRM
# 0.2919273

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMupsigSAR, GGMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.1835

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGMupF <- rbind(c(GGMupsigSAR, GGMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GGMupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GGMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGMupF) <- c("SA", "NotSA")
rownames(SAR.GGMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGMupF)
# P = 0.2565

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many transcripts are sexual antagonistic?
GGMdownsigSAR <- sum(GGMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GGMdownsigNotSAR <- length(GGMdownsig$gene.id)-GGMdownsigSAR
# 18

# The number of expected transcripts which are SA
length(GGMdownsig$gene.id)*SApropRM
# 0.5546619

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMdownsigSAR, GGMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.5439

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGMdownF <- rbind(c(GGMdownsigSAR, GGMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GGMdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-GGMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGMdownF) <- c("SA", "NotSA")
rownames(SAR.GGMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGMdownF)
# P = 0.4306

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many transcripts are sexual antagonistic?
G2MupsigSAR <- sum(G2Mupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
G2MupsigNotSAR <- length(G2Mupsig$gene.id)-G2MupsigSAR
# 29

# The number of expected transcripts which are SA
length(G2Mupsig$gene.id)*SApropRM
# 0.8465892

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G2MupsigSAR, G2MupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.3504

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G2MupF <- rbind(c(G2MupsigSAR, G2MupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-G2MupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-G2MupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G2MupF) <- c("SA", "NotSA")
rownames(SAR.G2MupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G2MupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many transcripts are sexual antagonistic?
G2MdownsigSAR <- sum(G2Mdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
G2MdownsigNotSAR <- length(G2Mdownsig$gene.id)-G2MdownsigSAR
# 36

# The number of expected transcripts which are SA
length(G2Mdownsig$gene.id)*SApropRM
# 1.050938

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G2MdownsigSAR, G2MdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.2981

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G2MdownF <- rbind(c(G2MdownsigSAR, G2MdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-G2MdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-G2MdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G2MdownF) <- c("SA", "NotSA")
rownames(SAR.G2MdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G2MdownF)
# P = 0.627

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many transcripts are sexual antagonistic?
G11MupsigSAR <- sum(G11Mupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
G11MupsigNotSAR <- length(G11Mupsig$gene.id)-G11MupsigSAR
# 21

# The number of expected transcripts which are SA
length(G11Mupsig$gene.id)*SApropRM
# 0.6130474

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G11MupsigSAR, G11MupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4268

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G11MupF <- rbind(c(G11MupsigSAR, G11MupsigNotSAR),
                     c(length(RsaM.data$FBgene.id)-G11MupsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-G11MupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G11MupF) <- c("SA", "NotSA")
rownames(SAR.G11MupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G11MupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many transcripts are sexual antagonistic?
G11MdownsigSAR <- sum(G11Mdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
G11MdownsigNotSAR <- length(G11Mdownsig$gene.id)-G11MdownsigSAR
# 7

# The number of expected transcripts which are SA
length(G11Mdownsig$gene.id)*SApropRM
# 0.2043491

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G11MdownsigSAR, G11MdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.6464

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G11MdownF <- rbind(c(G11MdownsigSAR, G11MdownsigNotSAR),
                       c(length(RsaM.data$FBgene.id)-G11MdownsigSAR, (length(RNAswapmale.data$gene.id)-length(RsaM.data$FBgene.id))-G11MdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G11MdownF) <- c("SA", "NotSA")
rownames(SAR.G11MdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G11MdownF)
# P = 1

# We get the same results, so that is good

