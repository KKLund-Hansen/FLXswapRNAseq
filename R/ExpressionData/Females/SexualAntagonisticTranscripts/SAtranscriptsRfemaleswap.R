################################################################################################
############# SEXUAL ANTAGONISTIC TRANSCRIPTS USING RUZICAK et al. FLX SWAP FEMALES ############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/SexualAntagonisticTranscripts")

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
Rsa.data <- read.table(file = "RuzickaetalAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

#First subset the Ruzicka et al. gene set to match the female data set
RsaF.data <- subset(Rsa.data, Rsa.data$FBgene.id %in% RNAswapfemale.data$gene.id)
# 490

# Next I'll calculate the proportion of the number of SA genes out of all genes
SApropRF <- length(RsaF.data$FBgene.id)/length(RNAswapfemale.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotSApropRF <- 1-SApropRF
# 0.97


######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many transcripts are sexual antagonistic?
GAFupsigSAR <- sum(GAFupsig$gene.id %in% RsaF.data$FBgene.id)
# 5
# How many are not?
GAFupsigNotSAR <- length(GAFupsig$gene.id)-GAFupsigSAR
# 124

# The number of expected transcripts which are SA
length(GAFupsig$gene.id)*SApropRF
# 3.765862

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFupsigSAR, GAFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.5186

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GAFupF <- rbind(c(GAFupsigSAR, GAFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GAFupsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-GAFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GAFupF) <- c("SA", "NotSA")
rownames(SAR.GAFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GAFupF)
# P = 0.4305

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many transcripts are sexual antagonistic?
GAFdownsigSAR <- sum(GAFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 73
# How many are not?
GAFdownsigNotSAR <- length(GAFdownsig$gene.id)-GAFdownsigSAR
# 1948

# The number of expected transcripts which are SA
length(GAFdownsig$gene.id)*SApropRF
# 58.99851

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFdownsigSAR, GAFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.0643




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many transcripts are sexual antagonistic?
GBFupsigSAR <- sum(GBFupsig$gene.id %in% RsaF.data$FBgene.id)
# 15
# How many are not?
GBFupsigNotSAR <- length(GBFupsig$gene.id)-GBFupsigSAR
# 555

# The number of expected transcripts which are SA
length(GBFupsig$gene.id)*SApropRF
# 16.63986

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFupsigSAR, GBFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.6833



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many transcripts are sexual antagonistic?
GBFdownsigSAR <- sum(GBFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 12
# How many are not?
GBFdownsigNotSAR <- length(GBFdownsig$gene.id)-GBFdownsigSAR
# 370

# The number of expected transcripts which are SA
length(GBFdownsig$gene.id)*SApropRF
# 11.15162

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFdownsigSAR, GBFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.7965




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many transcripts are sexual antagonistic?
GCFupsigSAR <- sum(GCFupsig$gene.id %in% RsaF.data$FBgene.id)
# 5
# How many are not?
GCFupsigNotSAR <- length(GCFupsig$gene.id)-GCFupsigSAR
# 172

# The number of expected transcripts which are SA
length(GCFupsig$gene.id)*SApropRF
# 5.167113

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFupsigSAR, GCFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.9405



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many transcripts are sexual antagonistic?
GCFdownsigSAR <- sum(GCFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 11
# How many are not?
GCFdownsigNotSAR <- length(GCFdownsig$gene.id)-GCFdownsigSAR
# 361

# The number of expected transcripts which are SA
length(GCFdownsig$gene.id)*SApropRF
# 10.8597

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFdownsigSAR, GCFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.9655




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many transcripts are sexual antagonistic?
GDFupsigSAR <- sum(GDFupsig$gene.id %in% RsaF.data$FBgene.id)
# 7
# How many are not?
GDFupsigNotSAR <- length(GDFupsig$gene.id)-GDFupsigSAR
# 184

# The number of expected transcripts which are SA
length(GDFupsig$gene.id)*SApropRF
# 5.575812

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFupsigSAR, GDFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.5404



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many transcripts are sexual antagonistic?
GDFdownsigSAR <- sum(GDFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 16
# How many are not?
GDFdownsigNotSAR <- length(GDFdownsig$gene.id)-GDFdownsigSAR
# 268

# The number of expected transcripts which are SA
length(GDFdownsig$gene.id)*SApropRF
# 8.290736

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFdownsigSAR, GDFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.00658




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many transcripts are sexual antagonistic?
GEFupsigSAR <- sum(GEFupsig$gene.id %in% RsaF.data$FBgene.id)
# 6
# How many are not?
GEFupsigNotSAR <- length(GEFupsig$gene.id)-GEFupsigSAR
# 92

# The number of expected transcripts which are SA
length(GEFupsig$gene.id)*SApropRF
# 2.860888

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFupsigSAR, GEFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.05962

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GEFupF <- rbind(c(GEFupsigSAR, GEFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GEFupsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-GEFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GEFupF) <- c("SA", "NotSA")
rownames(SAR.GEFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GEFupF)
# P = 0.06708

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many transcripts are sexual antagonistic?
GEFdownsigSAR <- sum(GEFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 8
# How many are not?
GEFdownsigNotSAR <- length(GEFdownsig$gene.id)-GEFdownsigSAR
# 256

# The number of expected transcripts which are SA
length(GEFdownsig$gene.id)*SApropRF
# 7.706881

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFdownsigSAR, GEFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.9147




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many transcripts are sexual antagonistic?
GFFupsigSAR <- sum(GFFupsig$gene.id %in% RsaF.data$FBgene.id)
# 3
# How many are not?
GFFupsigNotSAR <- length(GFFupsig$gene.id)-GFFupsigSAR
# 37

# The number of expected transcripts which are SA
length(GFFupsig$gene.id)*SApropRF
# 1.167709

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFupsigSAR, GFFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.08527

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GFFupF <- rbind(c(GFFupsigSAR, GFFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GFFupsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-GFFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GFFupF) <- c("SA", "NotSA")
rownames(SAR.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GFFupF)
# P = 0.1107

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many transcripts are sexual antagonistic?
GFFdownsigSAR <- sum(GFFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 16
# How many are not?
GFFdownsigNotSAR <- length(GFFdownsig$gene.id)-GFFdownsigSAR
# 248

# The number of expected transcripts which are SA
length(GFFdownsig$gene.id)*SApropRF
# 7.706881

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFdownsigSAR, GFFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.00243




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many transcripts are sexual antagonistic?
GGFupsigSAR <- sum(GGFupsig$gene.id %in% RsaF.data$FBgene.id)
# 3
# How many are not?
GGFupsigNotSAR <- length(GGFupsig$gene.id)-GGFupsigSAR
# 89

# The number of expected transcripts which are SA
length(GGFupsig$gene.id)*SApropRF
# 2.685731

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFupsigSAR, GGFupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.8457

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.GGFupF <- rbind(c(GGFupsigSAR, GGFupsigNotSAR),
                    c(length(RsaF.data$FBgene.id)-GGFupsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-GGFupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.GGFupF) <- c("SA", "NotSA")
rownames(SAR.GGFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.GGFupF)
# P = 0.7517

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many transcripts are sexual antagonistic?
GGFdownsigSAR <- sum(GGFdownsig$gene.id %in% RsaF.data$FBgene.id)
# 3
# How many are not?
GGFdownsigNotSAR <- length(GGFdownsig$gene.id)-GGFdownsigSAR
# 196

# The number of expected transcripts which are SA
length(GGFdownsig$gene.id)*SApropRF
# 5.809354

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFdownsigSAR, GGFdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.2368




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many transcripts are sexual antagonistic?
G23FupsigSAR <- sum(G23Fupsig$gene.id %in% RsaF.data$FBgene.id)
# 1
# How many are not?
G23FupsigNotSAR <- length(G23Fupsig$gene.id)-G23FupsigSAR
# 15

# The number of expected transcripts which are SA
length(G23Fupsig$gene.id)*SApropRF
# 0.4670837

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G23FupsigSAR, G23FupsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.4287

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G23FupF <- rbind(c(G23FupsigSAR, G23FupsigNotSAR),
                     c(length(RsaF.data$FBgene.id)-G23FupsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-G23FupsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G23FupF) <- c("SA", "NotSA")
rownames(SAR.G23FupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G23FupF)
# P = 0.3776

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many transcripts are sexual antagonistic?
G23FdownsigSAR <- sum(G23Fdownsig$gene.id %in% RsaF.data$FBgene.id)
# 0
# How many are not?
G23FdownsigNotSAR <- length(G23Fdownsig$gene.id)-G23FdownsigSAR
# 6

# The number of expected transcripts which are SA
length(G23Fdownsig$gene.id)*SApropRF
# 0.1751564

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G23FdownsigSAR, G23FdownsigNotSAR), p = c(SApropRF, NotSApropRF))
# P = 0.671

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAR.G23FdownF <- rbind(c(G23FdownsigSAR, G23FdownsigNotSAR),
                       c(length(RsaF.data$FBgene.id)-G23FdownsigSAR, (length(RNAswapfemale.data$gene.id)-length(RsaF.data$FBgene.id))-G23FdownsigNotSAR))
#Add names to the columns and rows
colnames(SAR.G23FdownF) <- c("SA", "NotSA")
rownames(SAR.G23FdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAR.G23FdownF)
# P = 1

# We get the same results, so that is good

