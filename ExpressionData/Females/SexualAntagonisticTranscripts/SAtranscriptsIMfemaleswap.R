################################################################################################
########## SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP FEMALES #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexualAntagonisticTranscripts")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the female data set
IMsaF.data <- subset(IMsa.data, IMsa.data$gene.id %in% FLXswapF.data$gene.id)
# 1337

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMF <- length(IMsaF.data$gene.id)/length(FLXswapF.data$gene.id)
# 0.09
# And for the chisq.test I also have to calculate
NotSApropIMF <- 1-SApropIMF
# 0.91



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many transcripts are sexual antagonistic?
FvCfupsigSAIM <- sum(FvCfupsig$gene.id %in% IMsaF.data$gene.id)
# 3
# How many are not?
FvCfupsigNotSAIM <- length(FvCfupsig$gene.id)-FvCfupsigSAIM
# 32

# The number of expected transcripts which are SA
length(FvCfupsig$gene.id)*SApropIMF
# 3.018059

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCfupsigSAIM, FvCfupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.9913

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FvCfupF <- rbind(c(FvCfupsigSAIM, FvCfupsigNotSAIM),
                      c(length(IMsaF.data$gene.id)-FvCfupsigSAIM, (length(FLXswapF.data$gene.id)-length(IMsaF.data$gene.id))-FvCfupsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FvCfupF) <- c("SA", "NotSA")
rownames(SAIM.FvCfupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FvCfupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many transcripts are sexual antagonistic?
FvCfdownsigSAIM <- sum(FvCfdownsig$gene.id %in% IMsaF.data$gene.id)
# 5
# How many are not?
FvCfdownsigNotSAIM <- length(FvCfdownsig$gene.id)-FvCfdownsigSAIM
# 35

# The number of expected transcripts which are SA
length(FvCfdownsig$gene.id)*SApropIMF
# 3.44921

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCfdownsigSAIM, FvCfdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.3824

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FvCfdownF <- rbind(c(FvCfdownsigSAIM, FvCfdownsigNotSAIM),
                        c(length(IMsaF.data$gene.id)-FvCfdownsigSAIM, (length(FLXswapF.data$gene.id)-length(IMsaF.data$gene.id))-FvCfdownsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FvCfdownF) <- c("SA", "NotSA")
rownames(SAIM.FvCfdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FvCfdownF)
# P = 0.3891

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many transcripts are sexual antagonistic?
AvSfupsigSAIM <- sum(AvSfupsig$gene.id %in% IMsaF.data$gene.id)
# 33
# How many are not?
AvSfupsigNotSAIM <- length(AvSfupsig$gene.id)-AvSfupsigSAIM
# 211

# The number of expected transcripts which are SA
length(AvSfupsig$gene.id)*SApropIMF
# 21.04018

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSfupsigSAIM, AvSfupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.00638



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many transcripts are sexual antagonistic?
AvSfdownsigSAIM <- sum(AvSfdownsig$gene.id %in% IMsaF.data$gene.id)
# 30
# How many are not?
AvSfdownsigNotSAIM <- length(AvSfdownsig$gene.id)-AvSfdownsigSAIM
# 131

# The number of expected transcripts which are SA
length(AvSfdownsig$gene.id)*SApropIMF
# 13.88307

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSfdownsigSAIM, AvSfdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 6.039e-06




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many transcripts are sexual antagonistic?
FsvCsFupsigSAIM <- sum(FsvCsFupsig$gene.id %in% IMsaF.data$gene.id)
# 19
# How many are not?
FsvCsFupsigNotSAIM <- length(FsvCsFupsig$gene.id)-FsvCsFupsigSAIM
# 118

# The number of expected transcripts which are SA
length(FsvCsFupsig$gene.id)*SApropIMF
# 11.81354

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsFupsigSAIM, FsvCsFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.02872



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many transcripts are sexual antagonistic?
FsvCsFdownsigSAIM <- sum(FsvCsFdownsig$gene.id %in% IMsaF.data$gene.id)
# 3
# How many are not?
FsvCsFdownsigNotSAIM <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigSAIM
# 52

# The number of expected transcripts which are SA
length(FsvCsFdownsig$gene.id)*SApropIMF
# 4.742664

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsFdownsigSAIM, FsvCsFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.4025

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FsvCsFdownF <- rbind(c(FsvCsFdownsigSAIM, FsvCsFdownsigNotSAIM),
                          c(length(IMsaF.data$gene.id)-FsvCsFdownsigSAIM, (length(FLXswapF.data$gene.id)-length(IMsaF.data$gene.id))-FsvCsFdownsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FsvCsFdownF) <- c("SA", "NotSA")
rownames(SAIM.FsvCsFdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FsvCsFdownF)
# P = 0.6271

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many transcripts are sexual antagonistic?
FavCaFupsigSAIM <- sum(FavCaFupsig$gene.id %in% IMsaF.data$gene.id)
# 56
# How many are not?
FavCaFupsigNotSAIM <- length(FavCaFupsig$gene.id)-FavCaFupsigSAIM
# 469

# The number of expected transcripts which are SA
length(FavCaFupsig$gene.id)*SApropIMF
# 45.27088

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaFupsigSAIM, FavCaFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.09529



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many transcripts are sexual antagonistic?
FavCaFdownsigSAIM <- sum(FavCaFdownsig$gene.id %in% IMsaF.data$gene.id)
# 346
# How many are not?
FavCaFdownsigNotSAIM <- length(FavCaFdownsig$gene.id)-FavCaFdownsigSAIM
# 2648

# The number of expected transcripts which are SA
length(FavCaFdownsig$gene.id)*SApropIMF
# 258.1734

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaFdownsigSAIM, FavCaFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 1.077e-08

