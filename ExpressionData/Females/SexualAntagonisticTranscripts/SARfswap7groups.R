################################################################################################
####### SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP FEMALES - 7 GROUPS #######
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexualAntagonisticTranscripts")

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the female data set
RsaF.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% FLXswap7gF.data$gene.id)
# 468

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRF <- length(RsaF.data$FBgene.id)/length(FLXswap7gF.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRF <- 1-SApropRF
# 0.97


######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many transcripts are sexual antagonistic?
GAFupsigSAR <- sum(GAFupsig$gene.id %in% RsaF.data$FBgene.id)
# 5
# How many are not?
GAFupsigNotSAR <- length(GAFupsig$gene.id)-GAFupsigSAR
# 117

# The number of expected transcripts which are SA
length(GAFupsig$gene.id)*SApropRF
# 3.682425

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFupsigSAR, GAFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.4857

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAFupF <- rbind(c(GAFupsigSAR, GAFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GAFupsigSAR, (length(FLXswap7gF.data$gene.id)-length(RsaF.data$FBgene.id))-GAFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GAFupF) <- c("SA", "NotSA")
rownames(SAR.GAFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GAFupF)
# P = 0.4212

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many transcripts are sexual antagonistic?
GAFdownsigSAR <- sum(GAFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 72
# How many are not?
GAFdownsigNotSAR <- length(GAFdownsig$gene.id)-GAFdownsigSAR
# 1899

# The number of expected transcripts which are SA
length(GAFdownsig$gene.id)*SApropRF
# 59.49229

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFdownsigSAR, GAFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.09963




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many transcripts are sexual antagonistic?
GBFupsigSAR <- sum(GBFupsig$gene.id %in% RsaF.data$FBgene.id)
# 13
# How many are not?
GBFupsigNotSAR <- length(GBFupsig$gene.id)-GBFupsigSAR
# 550

# The number of expected transcripts which are SA
length(GBFupsig$gene.id)*SApropRF
# 16.99349

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFupsigSAR, GBFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.3253



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many transcripts are sexual antagonistic?
GBFdownsigSAR <- sum(GBFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 10
# How many are not?
GBFdownsigNotSAR <- length(GBFdownsig$gene.id)-GBFdownsigSAR
# 372

# The number of expected transcripts which are SA
length(GBFdownsig$gene.id)*SApropRF
# 11.53022

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFdownsigSAR, GBFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.6472




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many transcripts are sexual antagonistic?
GCFupsigSAR <- sum(GCFupsig$gene.id %in% RsaF.data$FBgene.id)
# 6
# How many are not?
GCFupsigNotSAR <- length(GCFupsig$gene.id)-GCFupsigSAR
# 156

# The number of expected transcripts which are SA
length(GCFupsig$gene.id)*SApropRF
# 4.889777

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFupsigSAR, GCFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.6102

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GCFupF <- rbind(c(GCFupsigSAR, GCFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GCFupsigSAR, (length(FLXswap7gF.data$gene.id)-length(RsaF.data$FBgene.id))-GCFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GCFupF) <- c("SA", "NotSA")
rownames(SAR.GCFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GCFupF)
# P = 0.4932

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many transcripts are sexual antagonistic?
GCFdownsigSAR <- sum(GCFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 11
# How many are not?
GCFdownsigNotSAR <- length(GCFdownsig$gene.id)-GCFdownsigSAR
# 354

# The number of expected transcripts which are SA
length(GCFdownsig$gene.id)*SApropRF
# 11.01709

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFdownsigSAR, GCFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.9958




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many transcripts are sexual antagonistic?
GDFupsigSAR <- sum(GDFupsig$gene.id %in% RsaF.data$FBgene.id)
# 7
# How many are not?
GDFupsigNotSAR <- length(GDFupsig$gene.id)-GDFupsigSAR
# 169

# The number of expected transcripts which are SA
length(GDFupsig$gene.id)*SApropRF
# 5.312351

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFupsigSAR, GDFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.4572



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many transcripts are sexual antagonistic?
GDFdownsigSAR <- sum(GDFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 14
# How many are not?
GDFdownsigNotSAR <- length(GDFdownsig$gene.id)-GDFdownsigSAR
# 269

# The number of expected transcripts which are SA
length(GDFdownsig$gene.id)*SApropRF
# 8.542019

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFdownsigSAR, GDFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.05792




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many transcripts are sexual antagonistic?
GEFupsigSAR <- sum(GEFupsig$gene.id %in% RsaF.data$FBgene.id)
# 5
# How many are not?
GEFupsigNotSAR <- length(GEFupsig$gene.id)-GEFupsigSAR
# 86

# The number of expected transcripts which are SA
length(GEFupsig$gene.id)*SApropRF
# 2.746727

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFupsigSAR, GEFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.1674

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEFupF <- rbind(c(GEFupsigSAR, GEFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GEFupsigSAR, (length(FLXswap7gF.data$gene.id)-length(RsaF.data$FBgene.id))-GEFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GEFupF) <- c("SA", "NotSA")
rownames(SAR.GEFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GEFupF)
# P = 0.2016

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many transcripts are sexual antagonistic?
GEFdownsigSAR <- sum(GEFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 7
# How many are not?
GEFdownsigNotSAR <- length(GEFdownsig$gene.id)-GEFdownsigSAR
# 256

# The number of expected transcripts which are SA
length(GEFdownsig$gene.id)*SApropRF
# 7.938342

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFdownsigSAR, GEFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.7352




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many transcripts are sexual antagonistic?
GFFupsigSAR <- sum(GFFupsig$gene.id %in% RsaF.data$FBgene.id)
# 2
# How many are not?
GFFupsigNotSAR <- length(GFFupsig$gene.id)-GFFupsigSAR
# 35

# The number of expected transcripts which are SA
length(GFFupsig$gene.id)*SApropRF
# 1.116801

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFupsigSAR, GFFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.3961

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFFupF <- rbind(c(GFFupsigSAR, GFFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GFFupsigSAR, (length(FLXswap7gF.data$gene.id)-length(RsaF.data$FBgene.id))-GFFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GFFupF) <- c("SA", "NotSA")
rownames(SAR.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GFFupF)
# P = 0.3078

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many transcripts are sexual antagonistic?
GFFdownsigSAR <- sum(GFFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 16
# How many are not?
GFFdownsigNotSAR <- length(GFFdownsig$gene.id)-GFFdownsigSAR
# 247

# The number of expected transcripts which are SA
length(GFFdownsig$gene.id)*SApropRF
# 7.938342

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFdownsigSAR, GFFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.003667




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many transcripts are sexual antagonistic?
GGFupsigSAR <- sum(GGFupsig$gene.id %in% RsaF.data$FBgene.id)
# 3
# How many are not?
GGFupsigNotSAR <- length(GGFupsig$gene.id)-GGFupsigSAR
# 92

# The number of expected transcripts which are SA
length(GGFupsig$gene.id)*SApropRF
# 2.867462

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFupsigSAR, GGFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.9367

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGFupF <- rbind(c(GGFupsigSAR, GGFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GGFupsigSAR, (length(FLXswap7gF.data$gene.id)-length(RsaF.data$FBgene.id))-GGFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGFupF) <- c("SA", "NotSA")
rownames(SAR.GGFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGFupF)
# P = 0.7645

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many transcripts are sexual antagonistic?
GGFdownsigSAR <- sum(GGFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 3
# How many are not?
GGFdownsigNotSAR <- length(GGFdownsig$gene.id)-GGFdownsigSAR
# 191

# The number of expected transcripts which are SA
length(GGFdownsig$gene.id)*SApropRF
# 5.855659

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFdownsigSAR, GGFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.2308

