################################################################################################
#################### MITO-ANNOTATED TRANSCRIPTS FLX SWAP FEMALES - 7 GROUPS ####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/MitoNuclearConflict")

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
MA.data <- read.table(file = "mitoannotated.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First I'll subset the mito-annotated data set to match the female data set
MAfswap7g <- subset(MA.data, MA.data$From %in% FLXswap7gF.data$gene.id)
# 463 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-annotated transcripts out of all my transcripts
MApropF <- length(MAfswap7g$From)/length(FLXswap7gF.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotMApropF <- 1-MApropF
# 0.97



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many mito-annotated transcripts are there?
GAFupsigMA <- sum(GAFupsig$gene.id %in% MAfswap7g$From)
# 6
# How many are not?
GAFupsigNotMA <- length(GAFupsig$gene.id)-GAFupsigMA
# 116

# The number of expected mito-annotated transcripts
length(GAFupsig$gene.id)*MApropF
# 3.643083


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMA, GAFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.21

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GAFupF <- rbind(c(GAFupsigMA, GAFupsigNotMA), 
                       c(length(MAfswap7g$From)-GAFupsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GAFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GAFupF) <- c("MA", "NotMA")
rownames(MA.GAFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GAFupF)
# P = 0.1826

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many mito-annotated transcripts are there?
GAFdownsigMA <- sum(GAFdownsig$gene.id %in% MAfswap7g$From)
# 40
# How many are not?
GAFdownsigNotMA <- length(GAFdownsig$gene.id)-GAFdownsigMA
# 1931

# The number of expected mito-annotated transcripts
length(GAFdownsig$gene.id)*MApropF
# 56.3496


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMA, GAFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.01258




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many mito-annotated transcripts are there?
GBFupsigMA <- sum(GBFupsig$gene.id %in% MAfswap7g$From)
# 25
# How many are not?
GBFupsigNotMA <- length(GBFupsig$gene.id)-GBFupsigMA
# 538

# The number of expected mito-annotated transcripts
length(GBFupsig$gene.id)*MApropF
# 16.81193


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMA, GBFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.04261



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many mito-annotated transcripts are there?
GBFdownsigMA <- sum(GBFdownsig$gene.id %in% MAfswap7g$From)
# 22
# How many are not?
GBFdownsigNotMA <- length(GBFdownsig$gene.id)-GBFdownsigMA
# 360

# The number of expected mito-annotated transcripts
length(GBFdownsig$gene.id)*MApropF
# 11.40703


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMA, GBFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.001451




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many mito-annotated transcripts are there?
GCFupsigMA <- sum(GCFupsig$gene.id %in% MAfswap7g$From)
# 9
# How many are not?
GCFupsigNotMA <- length(GCFupsig$gene.id)-GCFupsigMA
# 153

# The number of expected mito-annotated transcripts
length(GCFupsig$gene.id)*MApropF
# 4.837536


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMA, GCFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.05468

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCFupF <- rbind(c(GCFupsigMA, GCFupsigNotMA), 
                   c(length(MAfswap7g$From)-GCFupsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GCFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GCFupF) <- c("MA", "NotMA")
rownames(MA.GCFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCFupF)
# P = 0.0614

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many mito-annotated transcripts are there?
GCFdownsigMA <- sum(GCFdownsig$gene.id %in% MAfswap7g$From)
# 3
# How many are not?
GCFdownsigNotMA <- length(GCFdownsig$gene.id)-GCFdownsigMA
# 362

# The number of expected mito-annotated transcripts
length(GCFdownsig$gene.id)*MApropF
# 10.89939


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMA, GCFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.01513

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCFdownF <- rbind(c(GCFdownsigMA, GCFdownsigNotMA), 
                     c(length(MAfswap7g$From)-GCFdownsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GCFdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GCFdownF) <- c("MA", "NotMA")
rownames(MA.GCFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCFdownF)
# P = 0.007772

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many mito-annotated transcripts are there?
GDFupsigMA <- sum(GDFupsig$gene.id %in% MAfswap7g$From)
# 18
# How many are not?
GDFupsigNotMA <- length(GDFupsig$gene.id)-GDFupsigMA
# 158

# The number of expected mito-annotated transcripts
length(GDFupsig$gene.id)*MApropF
# 5.255595


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMA, GDFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 1.661e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many mito-annotated transcripts are there?
GDFdownsigMA <- sum(GDFdownsig$gene.id %in% MAfswap7g$From)
# 4
# How many are not?
GDFdownsigNotMA <- length(GDFdownsig$gene.id)-GDFdownsigMA
# 279

# The number of expected mito-annotated transcripts
length(GDFdownsig$gene.id)*MApropF
# 8.450758


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMA, GDFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.1201




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many mito-annotated transcripts are there?
GEFupsigMA <- sum(GEFupsig$gene.id %in% MAfswap7g$From)
# 1
# How many are not?
GEFupsigNotMA <- length(GEFupsig$gene.id)-GEFupsigMA
# 90

# The number of expected mito-annotated transcripts
length(GEFupsig$gene.id)*MApropF
# 2.717381


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMA, GEFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.2902

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GEFupF <- rbind(c(GEFupsigMA, GEFupsigNotMA), 
                   c(length(MAfswap7g$From)-GEFupsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GEFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GEFupF) <- c("MA", "NotMA")
rownames(MA.GEFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GEFupF)
# P = 0.5285

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many mito-annotated transcripts are there?
GEFdownsigMA <- sum(GEFdownsig$gene.id %in% MAfswap7g$From)
# 17
# How many are not?
GEFdownsigNotMA <- length(GEFdownsig$gene.id)-GEFdownsigMA
# 246

# The number of expected mito-annotated transcripts
length(GEFdownsig$gene.id)*MApropF
# 7.853531


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMA, GEFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.0009209




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many mito-annotated transcripts are there?
GFFupsigMA <- sum(GFFupsig$gene.id %in% MAfswap7g$From)
# 2
# How many are not?
GFFupsigNotMA <- length(GFFupsig$gene.id)-GFFupsigMA
# 35

# The number of expected mito-annotated transcripts
length(GFFupsig$gene.id)*MApropF
# 1.104869


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMA, GFFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.3873

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GFFupF <- rbind(c(GFFupsigMA, GFFupsigNotMA), 
                   c(length(MAfswap7g$From)-GFFupsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GFFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GFFupF) <- c("MA", "NotMA")
rownames(MA.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GFFupF)
# P = 0.3033

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many mito-annotated transcripts are there?
GFFdownsigMA <- sum(GFFdownsig$gene.id %in% MAfswap7g$From)
# 3
# How many are not?
GFFdownsigNotMA <- length(GFFdownsig$gene.id)-GFFdownsigMA
# 260

# The number of expected mito-annotated transcripts
length(GFFdownsig$gene.id)*MApropF
# 7.853531


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMA, GFFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.07869




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many mito-annotated transcripts are there?
GGFupsigMA <- sum(GGFupsig$gene.id %in% MAfswap7g$From)
# 4
# How many are not?
GGFupsigNotMA <- length(GGFupsig$gene.id)-GGFupsigMA
# 91

# The number of expected mito-annotated transcripts
length(GGFupsig$gene.id)*MApropF
# 2.836827


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMA, GGFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.4832

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GGFupF <- rbind(c(GGFupsigMA, GGFupsigNotMA), 
                   c(length(MAfswap7g$From)-GGFupsigMA, (length(FLXswap7gF.data$gene.id)-length(MAfswap7g$From))-GGFupsigNotMA))
#Add names to the columns and rows
colnames(MA.GGFupF) <- c("MA", "NotMA")
rownames(MA.GGFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GGFupF)
# P = 0.3712

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many mito-annotated transcripts are there?
GGFdownsigMA <- sum(GGFdownsig$gene.id %in% MAfswap7g$From)
# 12
# How many are not?
GGFdownsigNotMA <- length(GGFdownsig$gene.id)-GGFdownsigMA
# 182

# The number of expected mito-annotated transcripts
length(GGFdownsig$gene.id)*MApropF
# 5.793099


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMA, GGFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.00884

