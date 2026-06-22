################################################################################################
########################## MITO-ANNOTATED TRANSCRIPTS FLX SWAP FEMALES ######################### 
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/MitoNuclearConflict")

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
MA.data <- read.table(file = "mitoannotated.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First I'll subset the mito-annotated data set to match the female data set
MAfemaleswap <- subset(MA.data, MA.data$From %in% RNAswapfemale.data$gene.id)
# 468 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-annotated transcripts out of all my transcripts
MApropF <- length(MAfemaleswap$From)/length(RNAswapfemale.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotMApropF <- 1-MApropF
# 0.97



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many mito-annotated transcripts are there?
GAFupsigMA <- sum(GAFupsig$gene.id %in% MAfemaleswap$From)
# 6
# How many are not?
GAFupsigNotMA <- length(GAFupsig$gene.id)-GAFupsigMA
# 123

# The number of expected mito-annotated transcripts
length(GAFupsig$gene.id)*MApropF
# 3.596783 

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMA, GAFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.1987

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GAFupF <- rbind(c(GAFupsigMA, GAFupsigNotMA), 
                       c(length(MAfemaleswap$From)-GAFupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GAFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GAFupF) <- c("MA", "NotMA")
rownames(MA.GAFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GAFupF)
# P = 0.1776

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many mito-annotated transcripts are there?
GAFdownsigMA <- sum(GAFdownsig$gene.id %in% MAfemaleswap$From)
# 41
# How many are not?
GAFdownsigNotMA <- length(GAFdownsig$gene.id)-GAFdownsigMA
# 1980

# The number of expected mito-annotated transcripts
length(GAFdownsig$gene.id)*MApropF
# 56.3496

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMA, GAFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.03809




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many mito-annotated transcripts are there?
GBFupsigMA <- sum(GBFupsig$gene.id %in% MAfemaleswap$From)
# 24
# How many are not?
GBFupsigNotMA <- length(GBFupsig$gene.id)-GBFupsigMA
# 546

# The number of expected mito-annotated transcripts
length(GBFupsig$gene.id)*MApropF
# 15.89276 

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMA, GBFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.03915



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many mito-annotated transcripts are there?
GBFdownsigMA <- sum(GBFdownsig$gene.id %in% MAfemaleswap$From)
# 21
# How many are not?
GBFdownsigNotMA <- length(GBFdownsig$gene.id)-GBFdownsigMA
# 361

# The number of expected mito-annotated transcripts
length(GBFdownsig$gene.id)*MApropF
# 10.65094

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMA, GBFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.001299




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many mito-annotated transcripts are there?
GCFupsigMA <- sum(GCFupsig$gene.id %in% MAfemaleswap$From)
# 10
# How many are not?
GCFupsigNotMA <- length(GCFupsig$gene.id)-GCFupsigMA
# 167

# The number of expected mito-annotated transcripts
length(GCFupsig$gene.id)*MApropF
# 4.935121 

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMA, GCFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.02076

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCFupF <- rbind(c(GCFupsigMA, GCFupsigNotMA), 
                   c(length(MAfemaleswap$From)-GCFupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GCFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GCFupF) <- c("MA", "NotMA")
rownames(MA.GCFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCFupF)
# P = 0.0334

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many mito-annotated transcripts are there?
GCFdownsigMA <- sum(GCFdownsig$gene.id %in% MAfemaleswap$From)
# 3
# How many are not?
GCFdownsigNotMA <- length(GCFdownsig$gene.id)-GCFdownsigMA
# 369

# The number of expected mito-annotated transcripts
length(GCFdownsig$gene.id)*MApropF
# 10.37212

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMA, GCFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.02025

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCFdownF <- rbind(c(GCFdownsigMA, GCFdownsigNotMA), 
                     c(length(MAfemaleswap$From)-GCFdownsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GCFdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GCFdownF) <- c("MA", "NotMA")
rownames(MA.GCFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCFdownF)
# P = 0.01536

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many mito-annotated transcripts are there?
GDFupsigMA <- sum(GDFupsig$gene.id %in% MAfemaleswap$From)
# 18
# How many are not?
GDFupsigNotMA <- length(GDFupsig$gene.id)-GDFupsigMA
# 173

# The number of expected mito-annotated transcripts
length(GDFupsig$gene.id)*MApropF
# 5.325469

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMA, GDFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 2.54e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many mito-annotated transcripts are there?
GDFdownsigMA <- sum(GDFdownsig$gene.id %in% MAfemaleswap$From)
# 4
# How many are not?
GDFdownsigNotMA <- length(GDFdownsig$gene.id)-GDFdownsigMA
# 280

# The number of expected mito-annotated transcripts
length(GDFdownsig$gene.id)*MApropF
# 7.918499

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMA, GDFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.1579




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many mito-annotated transcripts are there?
GEFupsigMA <- sum(GEFupsig$gene.id %in% MAfemaleswap$From)
# 1
# How many are not?
GEFupsigNotMA <- length(GEFupsig$gene.id)-GEFupsigMA
# 97

# The number of expected mito-annotated transcripts
length(GEFupsig$gene.id)*MApropF
# 2.73244

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMA, GEFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.2878

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GEFupF <- rbind(c(GEFupsigMA, GEFupsigNotMA), 
                   c(length(MAfemaleswap$From)-GEFupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GEFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GEFupF) <- c("MA", "NotMA")
rownames(MA.GEFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GEFupF)
# P = 0.5296

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many mito-annotated transcripts are there?
GEFdownsigMA <- sum(GEFdownsig$gene.id %in% MAfemaleswap$From)
# 14
# How many are not?
GEFdownsigNotMA <- length(GEFdownsig$gene.id)-GEFdownsigMA
# 250

# The number of expected mito-annotated transcripts
length(GEFdownsig$gene.id)*MApropF
# 7.360858

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMA, GEFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.01307




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many mito-annotated transcripts are there?
GFFupsigMA <- sum(GFFupsig$gene.id %in% MAfemaleswap$From)
# 2
# How many are not?
GFFupsigNotMA <- length(GFFupsig$gene.id)-GFFupsigMA
# 38

# The number of expected mito-annotated transcripts
length(GFFupsig$gene.id)*MApropF
# 1.115282

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMA, GFFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.3955

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GFFupF <- rbind(c(GFFupsigMA, GFFupsigNotMA), 
                   c(length(MAfemaleswap$From)-GFFupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GFFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GFFupF) <- c("MA", "NotMA")
rownames(MA.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GFFupF)
# P = 0.3072

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many mito-annotated transcripts are there?
GFFdownsigMA <- sum(GFFdownsig$gene.id %in% MAfemaleswap$From)
# 5
# How many are not?
GFFdownsigNotMA <- length(GFFdownsig$gene.id)-GFFdownsigMA
# 259

# The number of expected mito-annotated transcripts
length(GFFdownsig$gene.id)*MApropF
# 7.360858

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMA, GFFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.3775




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many mito-annotated transcripts are there?
GGFupsigMA <- sum(GGFupsig$gene.id %in% MAfemaleswap$From)
# 4
# How many are not?
GGFupsigNotMA <- length(GGFupsig$gene.id)-GGFupsigMA
# 88

# The number of expected mito-annotated transcripts
length(GGFupsig$gene.id)*MApropF
# 2.565147

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMA, GGFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.3635

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GGFupF <- rbind(c(GGFupsigMA, GGFupsigNotMA), 
                   c(length(MAfemaleswap$From)-GGFupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-GGFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GGFupF) <- c("MA", "NotMA")
rownames(MA.GGFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GGFupF)
# P = 0.3284

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many mito-annotated transcripts are there?
GGFdownsigMA <- sum(GGFdownsig$gene.id %in% MAfemaleswap$From)
# 13
# How many are not?
GGFdownsigNotMA <- length(GGFdownsig$gene.id)-GGFdownsigMA
# 186

# The number of expected mito-annotated transcripts
length(GGFdownsig$gene.id)*MApropF
# 5.548525

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMA, GGFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.001335




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many mito-annotated transcripts are there?
G23FupsigMA <- sum(G23Fupsig$gene.id %in% MAfemaleswap$From)
# 1
# How many are not?
G23FupsigNotMA <- length(G23Fupsig$gene.id)-G23FupsigMA
# 15

# The number of expected mito-annotated transcripts
length(G23Fupsig$gene.id)*MApropF
# 0.4461126

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FupsigMA, G23FupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.4003

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G23FupF <- rbind(c(G23FupsigMA, G23FupsigNotMA), 
                    c(length(MAfemaleswap$From)-G23FupsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-G23FupsigNotMA))
#Add names to the columns and rows
colnames(MA.G23FupF) <- c("MA", "NotMA")
rownames(MA.G23FupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G23FupF)
# P = 0.3641

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many mito-annotated transcripts are there?
G23FdownsigMA <- sum(G23Fdownsig$gene.id %in% MAfemaleswap$From)
# 0
# How many are not?
G23FdownsigNotMA <- length(G23Fdownsig$gene.id)-G23FdownsigMA
# 6

# The number of expected mito-annotated transcripts
length(G23Fdownsig$gene.id)*MApropF
# 0.1672922

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FdownsigMA, G23FdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.6783

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G23FdownF <- rbind(c(G23FdownsigMA, G23FdownsigNotMA),
                      c(length(MAfemaleswap$From)-G23FdownsigMA, (length(RNAswapfemale.data$gene.id)-length(MAfemaleswap$From))-G23FdownsigNotMA))
#Add names to the columns and rows
colnames(MA.G23FdownF) <- c("MA", "NotMA")
rownames(MA.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G23FdownF)
# P = 1

# We get the same results, so that is good

