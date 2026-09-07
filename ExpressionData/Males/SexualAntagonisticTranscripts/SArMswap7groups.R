################################################################################################
######## SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP MALES - 7 GROUPS ########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexualAntagonisticTranscripts")

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the male data set
RsaM.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% FLXswap7gM.data$gene.id)
# 468

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRM <- length(RsaM.data$FBgene.id)/length(FLXswap7gM.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRM <- 1-SApropRM
# 0.97


######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many transcripts are sexual antagonistic?
GAMupsigSAR <- sum(GAMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GAMupsigNotSAR <- length(GAMupsig$gene.id)-GAMupsigSAR
# 5

# The number of expected transcripts which are SA
length(GAMupsig$gene.id)*SApropRM
# 0.1509191

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMupsigSAR, GAMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.6932

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAMupF <- rbind(c(GAMupsigSAR, GAMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GAMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GAMupsigNotSAR))
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
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GAMdownsigSAR <- sum(GAMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GAMdownsigNotSAR <- length(GAMdownsig$gene.id)-GAMdownsigSAR
# 15

# The number of expected transcripts which are SA
length(GAMdownsig$gene.id)*SApropRM
# 0.482941

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAMdownsigSAR, GAMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4499

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAMdownF <- rbind(c(GAMdownsigSAR, GAMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GAMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GAMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GAMdownF) <- c("SA", "NotSA")
rownames(SAR.GAMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GAMdownF)
# P = 0.3878

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many transcripts are sexual antagonistic?
GBMupsigSAR <- sum(GBMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GBMupsigNotSAR <- length(GBMupsig$gene.id)-GBMupsigSAR
# 51

# The number of expected transcripts which are SA
length(GBMupsig$gene.id)*SApropRM
# 1.539374

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMupsigSAR, GBMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.2077

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GBMupF <- rbind(c(GBMupsigSAR, GBMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GBMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GBMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GBMupF) <- c("SA", "NotSA")
rownames(SAR.GBMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GBMupF)
# P = 0.408

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many transcripts are sexual antagonistic?
GBMdownsigSAR <- sum(GBMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 4
# How many are not?
GBMdownsigNotSAR <- length(GBMdownsig$gene.id)-GBMdownsigSAR
# 52

# The number of expected transcripts which are SA
length(GBMdownsig$gene.id)*SApropRM
# 1.690293

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBMdownsigSAR, GBMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.07123

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GBMdownF <- rbind(c(GBMdownsigSAR, GBMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GBMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GBMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GBMdownF) <- c("SA", "NotSA")
rownames(SAR.GBMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GBMdownF)
# P = 0.08838

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many transcripts are sexual antagonistic?
GCMupsigSAR <- sum(GCMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GCMupsigNotSAR <- length(GCMupsig$gene.id)-GCMupsigSAR
# 26

# The number of expected transcripts which are SA
length(GCMupsig$gene.id)*SApropRM
# 0.8149629

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMupsigSAR, GCMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.8351

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GCMupF <- rbind(c(GCMupsigSAR, GCMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GCMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GCMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GCMupF) <- c("SA", "NotSA")
rownames(SAR.GCMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GCMupF)
# P = 0.5632

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many transcripts are sexual antagonistic?
GCMdownsigSAR <- sum(GCMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 3
# How many are not?
GCMdownsigNotSAR <- length(GCMdownsig$gene.id)-GCMdownsigSAR
# 35

# The number of expected transcripts which are SA
length(GCMdownsig$gene.id)*SApropRM
# 1.146985

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCMdownsigSAR, GCMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.07893

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GCMdownF <- rbind(c(GCMdownsigSAR, GCMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GCMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GCMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GCMdownF) <- c("SA", "NotSA")
rownames(SAR.GCMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GCMdownF)
# P = 0.1062

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many transcripts are sexual antagonistic?
GDMupsigSAR <- sum(GDMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GDMupsigNotSAR <- length(GDMupsig$gene.id)-GDMupsigSAR
# 11

# The number of expected transcripts which are SA
length(GDMupsig$gene.id)*SApropRM
# 0.3622057

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMupsigSAR, GDMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.2819

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GDMupF <- rbind(c(GDMupsigSAR, GDMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GDMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GDMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GDMupF) <- c("SA", "NotSA")
rownames(SAR.GDMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GDMupF)
# P = 0.3078

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many transcripts are sexual antagonistic?
GDMdownsigSAR <- sum(GDMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GDMdownsigNotSAR <- length(GDMdownsig$gene.id)-GDMdownsigSAR
# 16

# The number of expected transcripts which are SA
length(GDMdownsig$gene.id)*SApropRM
# 0.482941

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDMdownsigSAR, GDMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4804

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GDMdownF <- rbind(c(GDMdownsigSAR, GDMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GDMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GDMdownsigNotSAR))
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
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many transcripts are sexual antagonistic?
GEMupsigSAR <- sum(GEMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GEMupsigNotSAR <- length(GEMupsig$gene.id)-GEMupsigSAR
# 20

# The number of expected transcripts which are SA
length(GEMupsig$gene.id)*SApropRM
# 0.6036762

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMupsigSAR, GEMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4301

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEMupF <- rbind(c(GEMupsigSAR, GEMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GEMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GEMupsigNotSAR))
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
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many transcripts are sexual antagonistic?
GEMdownsigSAR <- sum(GEMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GEMdownsigNotSAR <- length(GEMdownsig$gene.id)-GEMdownsigSAR
# 26

# The number of expected transcripts which are SA
length(GEMdownsig$gene.id)*SApropRM
# 0.7847791

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEMdownsigSAR, GEMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.3684

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEMdownF <- rbind(c(GEMdownsigSAR, GEMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GEMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GEMdownsigNotSAR))
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
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many transcripts are sexual antagonistic?
GFMupsigSAR <- sum(GFMupsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GFMupsigNotSAR <- length(GFMupsig$gene.id)-GFMupsigSAR
# 2

# The number of expected transcripts which are SA
length(GFMupsig$gene.id)*SApropRM
# 0.06036762

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMupsigSAR, GFMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.803

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFMupF <- rbind(c(GFMupsigSAR, GFMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GFMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GFMupsigNotSAR))
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
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many transcripts are sexual antagonistic?
GFMdownsigSAR <- sum(GFMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 0
# How many are not?
GFMdownsigNotSAR <- length(GFMdownsig$gene.id)-GFMdownsigSAR
# 3

# The number of expected transcripts which are SA
length(GFMdownsig$gene.id)*SApropRM
# 0.09055144

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFMdownsigSAR, GFMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.7599

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFMdownF <- rbind(c(GFMdownsigSAR, GFMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GFMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GFMdownsigNotSAR))
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
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many transcripts are sexual antagonistic?
GGMupsigSAR <- sum(GGMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GGMupsigNotSAR <- length(GGMupsig$gene.id)-GGMupsigSAR
# 9

# The number of expected transcripts which are SA
length(GGMupsig$gene.id)*SApropRM
# 0.3018381

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMupsigSAR, GGMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.1969

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGMupF <- rbind(c(GGMupsigSAR, GGMupsigNotSAR),
                    c(length(RsaM.data$FBgene.id)-GGMupsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GGMupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGMupF) <- c("SA", "NotSA")
rownames(SAR.GGMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGMupF)
# P = 0.264

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many transcripts are sexual antagonistic?
GGMdownsigSAR <- sum(GGMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
GGMdownsigNotSAR <- length(GGMdownsig$gene.id)-GGMdownsigSAR
# 21

# The number of expected transcripts which are SA
length(GGMdownsig$gene.id)*SApropRM
# 0.6640439

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGMdownsigSAR, GGMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.6755

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGMdownF <- rbind(c(GGMdownsigSAR, GGMdownsigNotSAR),
                      c(length(RsaM.data$FBgene.id)-GGMdownsigSAR, (length(FLXswap7gM.data$gene.id)-length(RsaM.data$FBgene.id))-GGMdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGMdownF) <- c("SA", "NotSA")
rownames(SAR.GGMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGMdownF)
# P = 0.4907

# We get the same results, so that is good

