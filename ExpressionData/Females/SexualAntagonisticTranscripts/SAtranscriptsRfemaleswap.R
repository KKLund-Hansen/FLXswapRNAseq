################################################################################################
############# SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP FEMALES ############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexualAntagonisticTranscripts")

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the female data set
RsaF.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% FLXswapF.data$gene.id)
# 468

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRF <- length(RsaF.data$FBgene.id)/length(FLXswapF.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRF <- 1-SApropRF
# 0.97



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many transcripts are sexual antagonistic?
FvCfupsigSAR <- sum(FvCfupsig$gene.id %in% RsaF.data$FBgene.id)
# 4
# How many are not?
FvCfupsigNotSAR <- length(FvCfupsig$gene.id)-FvCfupsigSAR
# 31

# The number of expected transcripts which are SA
length(FvCfupsig$gene.id)*SApropRF
# 1.056433

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCfupsigSAR, FvCfupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.003636

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FvCfupF <- rbind(c(FvCfupsigSAR, FvCfupsigNotSAR),
                      c(length(RsaF.data$FBgene.id)-FvCfupsigSAR, (length(FLXswapF.data$gene.id)-length(RsaF.data$FBgene.id))-FvCfupsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FvCfupF) <- c("SA", "NotSA")
rownames(SAR.FvCfupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FvCfupF)
# P = 0.02048

# We go from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many transcripts are sexual antagonistic?
FvCfdownsigSAR <- sum(FvCfdownsig$gene.id %in% RsaF.data$FBgene.id)
# 1
# How many are not?
FvCfdownsigNotSAR <- length(FvCfdownsig$gene.id)-FvCfdownsigSAR
# 39

# The number of expected transcripts which are SA
length(FvCfdownsig$gene.id)*SApropRF
# 1.207352

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCfdownsigSAR, FvCfdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.848

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FvCfdownF <- rbind(c(FvCfdownsigSAR, FvCfdownsigNotSAR),
                       c(length(RsaF.data$FBgene.id)-FvCfdownsigSAR, (length(FLXswapF.data$gene.id)-length(RsaF.data$FBgene.id))-FvCfdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FvCfdownF) <- c("SA", "NotSA")
rownames(SAR.FvCfdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FvCfdownF)
# P = 1

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many transcripts are sexual antagonistic?
AvSfupsigSAR <- sum(AvSfupsig$gene.id %in% RsaF.data$FBgene.id)
# 11
# How many are not?
AvSfupsigNotSAR <- length(AvSfupsig$gene.id)-AvSfupsigSAR
# 233

# The number of expected transcripts which are SA
length(AvSfupsig$gene.id)*SApropRF
# 7.36485

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSfupsigSAR, AvSfupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.1738



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many transcripts are sexual antagonistic?
AvSfdownsigSAR <- sum(AvSfdownsig$gene.id %in% RsaF.data$FBgene.id)
# 4
# How many are not?
AvSfdownsigNotSAR <- length(AvSfdownsig$gene.id)-AvSfdownsigSAR
# 157

# The number of expected transcripts which are SA
length(AvSfdownsig$gene.id)*SApropRF
# 4.859594

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSfdownsigSAR, AvSfdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.6921

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.AvSfdownF <- rbind(c(AvSfdownsigSAR, AvSfdownsigNotSAR),
                       c(length(RsaF.data$FBgene.id)-AvSfdownsigSAR, (length(FLXswapF.data$gene.id)-length(RsaF.data$FBgene.id))-AvSfdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.AvSfdownF) <- c("SA", "NotSA")
rownames(SAR.AvSfdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.AvSfdownF)
# P = 1

# We get the same results, so that is good




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many transcripts are sexual antagonistic?
FsvCsFupsigSAR <- sum(FsvCsFupsig$gene.id %in% RsaF.data$FBgene.id)
# 8
# How many are not?
FsvCsFupsigNotSAR <- length(FsvCsFupsig$gene.id)-FsvCsFupsigSAR
# 129

# The number of expected transcripts which are SA
length(FsvCsFupsig$gene.id)*SApropRF
# 4.135182

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsFupsigSAR, FsvCsFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.05362

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FsvCsFupF <- rbind(c(FsvCsFupsigSAR, FsvCsFupsigNotSAR),
                       c(length(RsaF.data$FBgene.id)-FsvCsFupsigSAR, (length(FLXswapF.data$gene.id)-length(RsaF.data$FBgene.id))-FsvCsFupsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FsvCsFupF) <- c("SA", "NotSA")
rownames(SAR.FsvCsFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FsvCsFupF)
# P = 0.07063

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many transcripts are sexual antagonistic?
FsvCsFdownsigSAR <- sum(FsvCsFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 0
# How many are not?
FsvCsFdownsigNotSAR <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigSAR
# 55

# The number of expected transcripts which are SA
length(FsvCsFdownsig$gene.id)*SApropRF
# 1.66011

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsFdownsigSAR, FsvCsFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.1908

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.FsvCsFdownF <- rbind(c(FsvCsFdownsigSAR, FsvCsFdownsigNotSAR),
                         c(length(RsaF.data$FBgene.id)-FsvCsFdownsigSAR, (length(FLXswapF.data$gene.id)-length(RsaF.data$FBgene.id))-FsvCsFdownsigNotSAR))
# Add names to the columns and rows
colnames(SAR.FsvCsFdownF) <- c("SA", "NotSA")
rownames(SAR.FsvCsFdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.FsvCsFdownF)
# P = 0.4155

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many transcripts are sexual antagonistic?
FavCaFupsigSAR <- sum(FavCaFupsig$gene.id %in% RsaF.data$FBgene.id)
# 20
# How many are not?
FavCaFupsigNotSAR <- length(FavCaFupsig$gene.id)-FavCaFupsigSAR
# 505

# The number of expected transcripts which are SA
length(FavCaFupsig$gene.id)*SApropRF
# 15.8465

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaFupsigSAR, FavCaFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.2894



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many transcripts are sexual antagonistic?
FavCaFdownsigSAR <- sum(FavCaFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 115
# How many are not?
FavCaFdownsigNotSAR <- length(FavCaFdownsig$gene.id)-FavCaFdownsigSAR
# 2879

# The number of expected transcripts which are SA
length(FavCaFdownsig$gene.id)*SApropRF
# 90.37033

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaFdownsigSAR, FavCaFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.008516


