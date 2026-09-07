################################################################################################
############## SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP MALES #############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexualAntagonisticTranscripts")

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the male data set
RsaM.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% FLXswapM.data$gene.id)
# 468

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRM <- length(RsaM.data$FBgene.id)/length(FLXswapM.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRM <- 1-SApropRM
# 0.97



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many transcripts are sexual antagonistic?
FvCmupsigSAR <- sum(FvCmupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
FvCmupsigNotSAR <- length(FvCmupsig$gene.id)-FvCmupsigSAR
# 56

# The number of expected transcripts which are SA
length(FvCmupsig$gene.id)*SApropRM
# 1.720477

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCmupsigSAR, FvCmupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.577

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FvCmupF <- rbind(c(FvCmupsigSAR, FvCmupsigNotSAR),
                     c(length(RsaM.data$FBgene.id)-FvCmupsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FvCmupsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FvCmupF) <- c("SA", "NotSA")
rownames(SAR.FvCmupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FvCmupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many transcripts are sexual antagonistic?
FvCmdownsigSAR <- sum(FvCmdownsig$gene.id %in% RsaM.data$FBgene.id)
# 4
# How many are not?
FvCmdownsigNotSAR <- length(FvCmdownsig$gene.id)-FvCmdownsigSAR
# 53

# The number of expected transcripts which are SA
length(FvCmdownsig$gene.id)*SApropRM
# 1.720477

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCmdownsigSAR, FvCmdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.07761

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FvCmdownF <- rbind(c(FvCmdownsigSAR, FvCmdownsigNotSAR),
                       c(length(RsaM.data$FBgene.id)-FvCmdownsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FvCmdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FvCmdownF) <- c("SA", "NotSA")
rownames(SAR.FvCmdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FvCmdownF)
# P = 0.09291

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many transcripts are sexual antagonistic?
AvSmupsigSAR <- sum(AvSmupsig$gene.id %in% RsaM.data$FBgene.id)
# 4
# How many are not?
AvSmupsigNotSAR <- length(AvSmupsig$gene.id)-AvSmupsigSAR
# 211

# The number of expected transcripts which are SA
length(AvSmupsig$gene.id)*SApropRM
# 6.48952

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSmupsigSAR, AvSmupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.321



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many transcripts are sexual antagonistic?
AvSmdownsigSAR <- sum(AvSmdownsig$gene.id %in% RsaM.data$FBgene.id)
# 8
# How many are not?
AvSmdownsigNotSAR <- length(AvSmdownsig$gene.id)-AvSmdownsigSAR
# 254

# The number of expected transcripts which are SA
length(AvSmdownsig$gene.id)*SApropRM
# 7.908159

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSmdownsigSAR, AvSmdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.3564




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many transcripts are sexual antagonistic?
FsvCsMupsigSAR <- sum(FsvCsMupsig$gene.id %in% RsaM.data$FBgene.id)
# 1
# How many are not?
FsvCsMupsigNotSAR <- length(FsvCsMupsig$gene.id)-FsvCsMupsigSAR
# 55

# The number of expected transcripts which are SA
length(FsvCsMupsig$gene.id)*SApropRM
# 1.690293

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsMupsigSAR, FsvCsMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.5898

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FsvCsMupF <- rbind(c(FsvCsMupsigSAR, FsvCsMupsigNotSAR),
                       c(length(RsaM.data$FBgene.id)-FsvCsMupsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FsvCsMupsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FsvCsMupF) <- c("SA", "NotSA")
rownames(SAR.FsvCsMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FsvCsMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many transcripts are sexual antagonistic?
FsvCsMdownsigSAR <- sum(FsvCsMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 4
# How many are not?
FsvCsMdownsigNotSAR <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigSAR
# 45

# The number of expected transcripts which are SA
length(FsvCsMdownsig$gene.id)*SApropRM
# 1.479007

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsMdownsigSAR, FsvCsMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.0353

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FsvCsMdownF <- rbind(c(FsvCsMdownsigSAR, FsvCsMdownsigNotSAR),
                         c(length(RsaM.data$FBgene.id)-FsvCsMdownsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FsvCsMdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FsvCsMdownF) <- c("SA", "NotSA")
rownames(SAR.FsvCsMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FsvCsMdownF)
# P = 0.05986

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many transcripts are sexual antagonistic?
FavCaMupsigSAR <- sum(FavCaMupsig$gene.id %in% RsaM.data$FBgene.id)
# 3
# How many are not?
FavCaMupsigNotSAR <- length(FavCaMupsig$gene.id)-FavCaMupsigSAR
# 62

# The number of expected transcripts which are SA
length(FavCaMupsig$gene.id)*SApropRM
# 1.961948

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaMupsigSAR, FavCaMupsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4517

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FavCaMupF <- rbind(c(FavCaMupsigSAR, FavCaMupsigNotSAR),
                       c(length(RsaM.data$FBgene.id)-FavCaMupsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FavCaMupsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FavCaMupF) <- c("SA", "NotSA")
rownames(SAR.FavCaMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FavCaMupF)
# P = 0.4486

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many transcripts are sexual antagonistic?
FavCaMdownsigSAR <- sum(FavCaMdownsig$gene.id %in% RsaM.data$FBgene.id)
# 5
# How many are not?
FavCaMdownsigNotSAR <- length(FavCaMdownsig$gene.id)-FavCaMdownsigSAR
# 116

# The number of expected transcripts which are SA
length(FavCaMdownsig$gene.id)*SApropRM
# 3.652241

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaMdownsigSAR, FavCaMdownsigNotSAR), p = c(SApropRM, NotSApropRM))
# P = 0.4739

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FavCaMdownF <- rbind(c(FavCaMdownsigSAR, FavCaMdownsigNotSAR),
                         c(length(RsaM.data$FBgene.id)-FavCaMdownsigSAR, (length(FLXswapM.data$gene.id)-length(RsaM.data$FBgene.id))-FavCaMdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FavCaMdownF) <- c("SA", "NotSA")
rownames(SAR.FavCaMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FavCaMdownF)
# P = 0.4181

# We get the same results, so that is good

