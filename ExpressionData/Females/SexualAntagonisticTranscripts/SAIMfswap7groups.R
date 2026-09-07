################################################################################################
#### SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP FEMALES - 7 GROUPS ####
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
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the female data set
IMsaF.data <- subset(IMsa.data, IMsa.data$gene.id %in% FLXswap7gF.data$gene.id)
# 1337

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMF <- length(IMsaF.data$gene.id)/length(FLXswap7gF.data$gene.id)
# 0.09
# And for the chisq.test I also have to calculate
NotSApropIMF <- 1-SApropIMF
# 0.91



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many transcripts are sexual antagonistic?
GAFupsigSAIM <- sum(GAFupsig$gene.id %in% IMsaF.data$gene.id)
# 14
# How many are not?
GAFupsigNotSAIM <- length(GAFupsig$gene.id)-GAFupsigSAIM
# 108

# The number of expected transcripts which are SA
length(GAFupsig$gene.id)*SApropIMF
# 10.52009

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFupsigSAIM, GAFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.2617



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many transcripts are sexual antagonistic?
GAFdownsigSAIM <- sum(GAFdownsig$gene.id %in% IMsaF.data$gene.id)
# 230
# How many are not?
GAFdownsigNotSAIM <- length(GAFdownsig$gene.id)-GAFdownsigSAIM
# 1741

# The number of expected transcripts which are SA
length(GAFdownsig$gene.id)*SApropIMF
# 169.9598

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GAFdownsigSAIM, GAFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 1.451e-06




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many transcripts are sexual antagonistic?
GBFupsigSAIM <- sum(GBFupsig$gene.id %in% IMsaF.data$gene.id)
# 50
# How many are not?
GBFupsigNotSAIM <- length(GBFupsig$gene.id)-GBFupsigSAIM
# 513

# The number of expected transcripts which are SA
length(GBFupsig$gene.id)*SApropIMF
# 48.54763

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFupsigSAIM, GBFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.8274



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many transcripts are sexual antagonistic?
GBFdownsigSAIM <- sum(GBFdownsig$gene.id %in% IMsaF.data$gene.id)
# 59
# How many are not?
GBFdownsigNotSAIM <- length(GBFdownsig$gene.id)-GBFdownsigSAIM
# 323

# The number of expected transcripts which are SA
length(GBFdownsig$gene.id)*SApropIMF
# 32.93995

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GBFdownsigSAIM, GBFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 2.034e-06




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many transcripts are sexual antagonistic?
GCFupsigSAIM <- sum(GCFupsig$gene.id %in% IMsaF.data$gene.id)
# 21
# How many are not?
GCFupsigNotSAIM <- length(GCFupsig$gene.id)-GCFupsigSAIM
# 141

# The number of expected transcripts which are SA
length(GCFupsig$gene.id)*SApropIMF
# 13.9693

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFupsigSAIM, GCFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.04908



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many transcripts are sexual antagonistic?
GCFdownsigSAIM <- sum(GCFdownsig$gene.id %in% IMsaF.data$gene.id)
# 41
# How many are not?
GCFdownsigNotSAIM <- length(GCFdownsig$gene.id)-GCFdownsigSAIM
# 324

# The number of expected transcripts which are SA
length(GCFdownsig$gene.id)*SApropIMF
# 31.47404

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GCFdownsigSAIM, GCFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.07568




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many transcripts are sexual antagonistic?
GDFupsigSAIM <- sum(GDFupsig$gene.id %in% IMsaF.data$gene.id)
# 12
# How many are not?
GDFupsigNotSAIM <- length(GDFupsig$gene.id)-GDFupsigSAIM
# 164

# The number of expected transcripts which are SA
length(GDFupsig$gene.id)*SApropIMF
# 15.17652

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFupsigSAIM, GDFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.3937



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many transcripts are sexual antagonistic?
GDFdownsigSAIM <- sum(GDFdownsig$gene.id %in% IMsaF.data$gene.id)
# 36
# How many are not?
GDFdownsigNotSAIM <- length(GDFdownsig$gene.id)-GDFdownsigSAIM
# 247

# The number of expected transcripts which are SA
length(GDFdownsig$gene.id)*SApropIMF
# 24.40316

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GDFdownsigSAIM, GDFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.01406




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many transcripts are sexual antagonistic?
GEFupsigSAIM <- sum(GEFupsig$gene.id %in% IMsaF.data$gene.id)
# 15
# How many are not?
GEFupsigNotSAIM <- length(GEFupsig$gene.id)-GEFupsigSAIM
# 76

# The number of expected transcripts which are SA
length(GEFupsig$gene.id)*SApropIMF
# 7.846953

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFupsigSAIM, GEFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.007556



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many transcripts are sexual antagonistic?
GEFdownsigSAIM <- sum(GEFdownsig$gene.id %in% IMsaF.data$gene.id)
# 17
# How many are not?
GEFdownsigNotSAIM <- length(GEFdownsig$gene.id)-GEFdownsigSAIM
# 246

# The number of expected transcripts which are SA
length(GEFdownsig$gene.id)*SApropIMF
# 22.67856

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GEFdownsigSAIM, GEFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.2122




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many transcripts are sexual antagonistic?
GFFupsigSAIM <- sum(GFFupsig$gene.id %in% IMsaF.data$gene.id)
# 3
# How many are not?
GFFupsigNotSAIM <- length(GFFupsig$gene.id)-GFFupsigSAIM
# 34

# The number of expected transcripts which are SA
length(GFFupsig$gene.id)*SApropIMF
# 3.190519

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFupsigSAIM, GFFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.9112

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.GFFupF <- rbind(c(GFFupsigSAIM, GFFupsigNotSAIM),
                     c(length(IMsaF.data$gene.id)-GFFupsigSAIM, (length(FLXswap7gF.data$gene.id)-length(IMsaF.data$gene.id))-GFFupsigNotSAIM))
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
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many transcripts are sexual antagonistic?
GFFdownsigSAIM <- sum(GFFdownsig$gene.id %in% IMsaF.data$gene.id)
# 32
# How many are not?
GFFdownsigNotSAIM <- length(GFFdownsig$gene.id)-GFFdownsigSAIM
# 231

# The number of expected transcripts which are SA
length(GFFdownsig$gene.id)*SApropIMF
# 22.67856

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GFFdownsigSAIM, GFFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.04059




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many transcripts are sexual antagonistic?
GGFupsigSAIM <- sum(GGFupsig$gene.id %in% IMsaF.data$gene.id)
# 11
# How many are not?
GGFupsigNotSAIM <- length(GGFupsig$gene.id)-GGFupsigSAIM
# 84

# The number of expected transcripts which are SA
length(GGFupsig$gene.id)*SApropIMF
# 8.191874

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFupsigSAIM, GGFupsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.3047



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many transcripts are sexual antagonistic?
GGFdownsigSAIM <- sum(GGFdownsig$gene.id %in% IMsaF.data$gene.id)
# 26
# How many are not?
GGFdownsigNotSAIM <- length(GGFdownsig$gene.id)-GGFdownsigSAIM
# 168

# The number of expected transcripts which are SA
length(GGFdownsig$gene.id)*SApropIMF
# 16.72867

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(GGFdownsigSAIM, GGFdownsigNotSAIM), p = c(SApropIMF, NotSApropIMF))
# P = 0.01772

