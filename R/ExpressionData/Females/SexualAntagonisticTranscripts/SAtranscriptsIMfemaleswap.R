################################################################################################
########## SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP FEMALES #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/SexualAntagonisticTranscripts")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the female data set
IMsaF.data <- subset(IMsa.data, IMsa.data$gene.id %in% RNAswapfemale.data$gene.id)
# 1353

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMF <- length(IMsaF.data$gene.id)/length(RNAswapfemale.data$gene.id)
# 0.08
# And for the chisq.test I also have to calculate
NotSApropIMF <- 1-SApropIMF
# 0.92



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many transcripts are sexual antagonistic?
GAFupsigSAIM <- sum(GAFupsig$gene.id %in% IMsaF.data$gene.id)
# 15
# How many are not?
GAFupsigNotSAIM <- length(GAFupsig$gene.id)-GAFupsigSAIM
# 114

# The number of expected transcripts which are SA
length(GAFupsig$gene.id)*SApropIMF
# 10.39839

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFupsigSAIM, GAFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.1367



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many transcripts are sexual antagonistic?
GAFdownsigSAIM <- sum(GAFdownsig$gene.id %in% IMsaF.data$gene.id)
# 233
# How many are not?
GAFdownsigNotSAIM <- length(GAFdownsig$gene.id)-GAFdownsigSAIM
# 1788

# The number of expected transcripts which are SA
length(GAFdownsig$gene.id)*SApropIMF
# 162.9081

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFdownsigSAIM, GAFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 1.021e-08




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many transcripts are sexual antagonistic?
GBFupsigSAIM <- sum(GBFupsig$gene.id %in% IMsaF.data$gene.id)
# 48
# How many are not?
GBFupsigNotSAIM <- length(GBFupsig$gene.id)-GBFupsigSAIM
# 522

# The number of expected transcripts which are SA
length(GBFupsig$gene.id)*SApropIMF
# 45.94638

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFupsigSAIM, GBFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.752



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many transcripts are sexual antagonistic?
GBFdownsigSAIM <- sum(GBFdownsig$gene.id %in% IMsaF.data$gene.id)
# 58
# How many are not?
GBFdownsigNotSAIM <- length(GBFdownsig$gene.id)-GBFdownsigSAIM
# 324

# The number of expected transcripts which are SA
length(GBFdownsig$gene.id)*SApropIMF
# 30.79214

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFdownsigSAIM, GBFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 3.161e-07




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many transcripts are sexual antagonistic?
GCFupsigSAIM <- sum(GCFupsig$gene.id %in% IMsaF.data$gene.id)
# 21
# How many are not?
GCFupsigNotSAIM <- length(GCFupsig$gene.id)-GCFupsigSAIM
# 156

# The number of expected transcripts which are SA
length(GCFupsig$gene.id)*SApropIMF
# 14.26756

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFupsigSAIM, GCFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.06305



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many transcripts are sexual antagonistic?
GCFdownsigSAIM <- sum(GCFdownsig$gene.id %in% IMsaF.data$gene.id)
# 41
# How many are not?
GCFdownsigNotSAIM <- length(GCFdownsig$gene.id)-GCFdownsigSAIM
# 331

# The number of expected transcripts which are SA
length(GCFdownsig$gene.id)*SApropIMF
# 29.98606

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFdownsigSAIM, GCFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.03594




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many transcripts are sexual antagonistic?
GDFupsigSAIM <- sum(GDFupsig$gene.id %in% IMsaF.data$gene.id)
# 14
# How many are not?
GDFupsigNotSAIM <- length(GDFupsig$gene.id)-GDFupsigSAIM
# 177

# The number of expected transcripts which are SA
length(GDFupsig$gene.id)*SApropIMF
# 15.39607

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFupsigSAIM, GDFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.7106



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many transcripts are sexual antagonistic?
GDFdownsigSAIM <- sum(GDFdownsig$gene.id %in% IMsaF.data$gene.id)
# 37
# How many are not?
GDFdownsigNotSAIM <- length(GDFdownsig$gene.id)-GDFdownsigSAIM
# 247

# The number of expected transcripts which are SA
length(GDFdownsig$gene.id)*SApropIMF
# 22.89258

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFdownsigSAIM, GDFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.002105




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many transcripts are sexual antagonistic?
GEFupsigSAIM <- sum(GEFupsig$gene.id %in% IMsaF.data$gene.id)
# 15
# How many are not?
GEFupsigNotSAIM <- length(GEFupsig$gene.id)-GEFupsigSAIM
# 83

# The number of expected transcripts which are SA
length(GEFupsig$gene.id)*SApropIMF
# 7.899553

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFupsigSAIM, GEFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.008421



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many transcripts are sexual antagonistic?
GEFdownsigSAIM <- sum(GEFdownsig$gene.id %in% IMsaF.data$gene.id)
# 18
# How many are not?
GEFdownsigNotSAIM <- length(GEFdownsig$gene.id)-GEFdownsigSAIM
# 246

# The number of expected transcripts which are SA
length(GEFdownsig$gene.id)*SApropIMF
# 21.28043

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFdownsigSAIM, GEFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.4583




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many transcripts are sexual antagonistic?
GFFupsigSAIM <- sum(GFFupsig$gene.id %in% IMsaF.data$gene.id)
# 3
# How many are not?
GFFupsigNotSAIM <- length(GFFupsig$gene.id)-GFFupsigSAIM
# 37

# The number of expected transcripts which are SA
length(GFFupsig$gene.id)*SApropIMF
# 3.224307

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFupsigSAIM, GFFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.8963

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFFupF <- rbind(c(GFFupsigSAIM, GFFupsigNotSAIM),
                     c(length(IMsaF.data$gene.id)-GFFupsigSAIM, (length(RNAswapfemale.data$gene.id)-length(IMsaF.data$gene.id))-GFFupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.GFFupF) <- c("SA", "NotSA")
rownames(SAIM.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.GFFupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many transcripts are sexual antagonistic?
GFFdownsigSAIM <- sum(GFFdownsig$gene.id %in% IMsaF.data$gene.id)
# 30
# How many are not?
GFFdownsigNotSAIM <- length(GFFdownsig$gene.id)-GFFdownsigSAIM
# 234

# The number of expected transcripts which are SA
length(GFFdownsig$gene.id)*SApropIMF
# 21.28043

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFdownsigSAIM, GFFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.04869




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many transcripts are sexual antagonistic?
GGFupsigSAIM <- sum(GGFupsig$gene.id %in% IMsaF.data$gene.id)
# 9
# How many are not?
GGFupsigNotSAIM <- length(GGFupsig$gene.id)-GGFupsigSAIM
# 83

# The number of expected transcripts which are SA
length(GGFupsig$gene.id)*SApropIMF
# 7.415907

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFupsigSAIM, GGFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.5441



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many transcripts are sexual antagonistic?
GGFdownsigSAIM <- sum(GGFdownsig$gene.id %in% IMsaF.data$gene.id)
# 25
# How many are not?
GGFdownsigNotSAIM <- length(GGFdownsig$gene.id)-GGFdownsigSAIM
# 174

# The number of expected transcripts which are SA
length(GGFdownsig$gene.id)*SApropIMF
# 16.04093

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFdownsigSAIM, GGFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.01965




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many transcripts are sexual antagonistic?
G23FupsigSAIM <- sum(G23Fupsig$gene.id %in% IMsaF.data$gene.id)
# 2
# How many are not?
G23FupsigNotSAIM <- length(G23Fupsig$gene.id)-G23FupsigSAIM
# 14

# The number of expected transcripts which are SA
length(G23Fupsig$gene.id)*SApropIMF
# 1.289723

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G23FupsigSAIM, G23FupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.5142

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G23FupF <- rbind(c(G23FupsigSAIM, G23FupsigNotSAIM),
                      c(length(IMsaF.data$gene.id)-G23FupsigSAIM, (length(RNAswapfemale.data$gene.id)-length(IMsaF.data$gene.id))-G23FupsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G23FupF) <- c("SA", "NotSA")
rownames(SAIM.G23FupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G23FupF)
# P = 0.3738

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many transcripts are sexual antagonistic?
G23FdownsigSAIM <- sum(G23Fdownsig$gene.id %in% IMsaF.data$gene.id)
# 1
# How many are not?
G23FdownsigNotSAIM <- length(G23Fdownsig$gene.id)-G23FdownsigSAIM
# 5

# The number of expected transcripts which are SA
length(G23Fdownsig$gene.id)*SApropIMF
# 0.4836461

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(G23FdownsigSAIM, G23FdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.4387

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.G23FdownF <- rbind(c(G23FdownsigSAIM, G23FdownsigNotSAIM),
                        c(length(IMsaF.data$gene.id)-G23FdownsigSAIM, (length(RNAswapfemale.data$gene.id)-length(IMsaF.data$gene.id))-G23FdownsigNotSAIM))
#Add names to the columns and rows
colnames(SAIM.G23FdownF) <- c("SA", "NotSA")
rownames(SAIM.G23FdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.G23FdownF)
# P = 0.3961

# We get the same results, so that is good

