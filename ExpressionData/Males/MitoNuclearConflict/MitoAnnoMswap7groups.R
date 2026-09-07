################################################################################################
##################### MITO-ANNOTATED TRANSCRIPTS FLX SWAP MALES - 7 GROUPS #####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/MitoNuclearConflict")

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
MA.data <- read.table(file = "mitoannotated.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First I'll subset the mito-annotated data set to match the male data set
MAmwap7g <- subset(MA.data, MA.data$From %in% FLXswap7gM.data$gene.id)
# 463 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-annotated transcripts out of all my transcripts
MApropM <- length(MAmwap7g$From)/length(FLXswap7gM.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotMApropM <- 1-MApropM
# 0.97



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many mito-annotated transcripts are there?
GAMupsigMA <- sum(GAMupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GAMupsigNotMA <- length(GAMupsig$gene.id)-GAMupsigMA
# 5

# The number of expected mito-annotated transcripts
length(GAMupsig$gene.id)*MApropM
# 0.1493067


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMupsigMA, GAMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.6948

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GAMupF <- rbind(c(GAMupsigMA, GAMupsigNotMA), 
                   c(length(MAmwap7g$From)-GAMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GAMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GAMupF) <- c("MA", "NotMA")
rownames(MA.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many mito-annotated transcripts are there?
GAMdownsigMA <- sum(GAMdownsig$gene.id %in% MAmwap7g$From)
# 1
# How many are not?
GAMdownsigNotMA <- length(GAMdownsig$gene.id)-GAMdownsigMA
# 15

# The number of expected mito-annotated transcripts
length(GAMdownsig$gene.id)*MApropM
# 0.4461126


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigMA, GAMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.4431

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GAMdownF <- rbind(c(GAMdownsigMA, GAMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GAMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GAMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GAMdownF) <- c("MA", "NotMA")
rownames(MA.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GAMdownF)
# P = 0.3845

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many mito-annotated transcripts are there?
GBMupsigMA <- sum(GBMupsig$gene.id %in% MAmwap7g$From)
# 1
# How many are not?
GBMupsigNotMA <- length(GBMupsig$gene.id)-GBMupsigMA
# 50

# The number of expected mito-annotated transcripts
length(GBMupsig$gene.id)*MApropM
# 1.522928


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMupsigMA, GBMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.667

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GBMupF <- rbind(c(GBMupsigMA, GBMupsigNotMA), 
                   c(length(MAmwap7g$From)-GBMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GBMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GBMupF) <- c("MA", "NotMA")
rownames(MA.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GBMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many mito-annotated transcripts are there?
GBMdownsigMA <- sum(GBMdownsig$gene.id %in% MAmwap7g$From)
# 4
# How many are not?
GBMdownsigNotMA <- length(GBMdownsig$gene.id)-GBMdownsigMA
# 52

# The number of expected mito-annotated transcripts
length(GBMdownsig$gene.id)*MApropM
# 1.672235


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigMA, GBMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.06761

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GBMdownF <- rbind(c(GBMdownsigMA, GBMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GBMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GBMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GBMdownF) <- c("MA", "NotMA")
rownames(MA.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GBMdownF)
# P = 0.08575

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many mito-annotated transcripts are there?
GCMupsigMA <- sum(GCMupsig$gene.id %in% MAmwap7g$From)
# 1
# How many are not?
GCMupsigNotMA <- length(GCMupsig$gene.id)-GCMupsigMA
# 26

# The number of expected mito-annotated transcripts
length(GCMupsig$gene.id)*MApropM
# 0.806256


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMupsigMA, GCMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.8266

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCMupF <- rbind(c(GCMupsigMA, GCMupsigNotMA), 
                   c(length(MAmwap7g$From)-GCMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GCMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GCMupF) <- c("MA", "NotMA")
rownames(MA.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCMupF)
# P = 0.5592

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many mito-annotated transcripts are there?
GCMdownsigMA <- sum(GCMdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GCMdownsigNotMA <- length(GCMdownsig$gene.id)-GCMdownsigMA
# 38

# The number of expected mito-annotated transcripts
length(GCMdownsig$gene.id)*MApropM
# 1.134731


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigMA, GCMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.2795

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GCMdownF <- rbind(c(GCMdownsigMA, GCMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GCMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GCMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GCMdownF) <- c("MA", "NotMA")
rownames(MA.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GCMdownF)
# P = 0.629996

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many mito-annotated transcripts are there?
GDMupsigMA <- sum(GDMupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GDMupsigNotMA <- length(GDMupsig$gene.id)-GDMupsigMA
# 12

# The number of expected mito-annotated transcripts
length(GDMupsig$gene.id)*MApropM
# 0.358336


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMupsigMA, GDMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.5434

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GDMupF <- rbind(c(GDMupsigMA, GDMupsigNotMA), 
                   c(length(MAmwap7g$From)-GDMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GDMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GDMupF) <- c("MA", "NotMA")
rownames(MA.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many mito-annotated transcripts are there?
GDMdownsigMA <- sum(GDMdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GDMdownsigNotMA <- length(GDMdownsig$gene.id)-GDMdownsigMA
# 16

# The number of expected mito-annotated transcripts
length(GDMdownsig$gene.id)*MApropM
# 0.4777814


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigMA, GDMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.4828

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GDMdownF <- rbind(c(GDMdownsigMA, GDMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GDMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GDMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GDMdownF) <- c("MA", "NotMA")
rownames(MA.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GDMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many mito-annotated transcripts are there?
GEMupsigMA <- sum(GEMupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GEMupsigNotMA <- length(GEMupsig$gene.id)-GEMupsigMA
# 20

# The number of expected mito-annotated transcripts
length(GEMupsig$gene.id)*MApropM
# 0.5972267


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMupsigMA, GEMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.4327

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GEMupF <- rbind(c(GEMupsigMA, GEMupsigNotMA), 
                   c(length(MAmwap7g$From)-GEMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GEMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GEMupF) <- c("MA", "NotMA")
rownames(MA.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GEMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many mito-annotated transcripts are there?
GEMdownsigMA <- sum(GEMdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GEMdownsigNotMA <- length(GEMdownsig$gene.id)-GEMdownsigMA
# 26

# The number of expected mito-annotated transcripts
length(GEMdownsig$gene.id)*MApropM
# 0.7763947


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigMA, GEMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.371

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GEMdownF <- rbind(c(GEMdownsigMA, GEMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GEMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GEMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GEMdownF) <- c("MA", "NotMA")
rownames(MA.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GEMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many mito-annotated transcripts are there?
GFMupsigMA <- sum(GFMupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GFMupsigNotMA <- length(GFMupsig$gene.id)-GFMupsigMA
# 2

# The number of expected mito-annotated transcripts
length(GFMupsig$gene.id)*MApropM
# 0.05972267


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMupsigMA, GFMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.804

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GFMupF <- rbind(c(GFMupsigMA, GFMupsigNotMA), 
                   c(length(MAmwap7g$From)-GFMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GFMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GFMupF) <- c("MA", "NotMA")
rownames(MA.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many mito-annotated transcripts are there?
GFMdownsigMA <- sum(GFMdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GFMdownsigNotMA <- length(GFMdownsig$gene.id)-GFMdownsigMA
# 3

# The number of expected mito-annotated transcripts
length(GFMdownsig$gene.id)*MApropM
# 0.08958401


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigMA, GFMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.7612

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GFMdownF <- rbind(c(GFMdownsigMA, GFMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GFMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GFMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GFMdownF) <- c("MA", "NotMA")
rownames(MA.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GFMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many mito-annotated transcripts are there?
GGMupsigMA <- sum(GGMupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
GGMupsigNotMA <- length(GGMupsig$gene.id)-GGMupsigMA
# 10

# The number of expected mito-annotated transcripts
length(GGMupsig$gene.id)*MApropM
# 0.2986134


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMupsigMA, GGMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.579

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GGMupF <- rbind(c(GGMupsigMA, GGMupsigNotMA), 
                   c(length(MAmwap7g$From)-GGMupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GGMupsigNotMA))
#Add names to the columns and rows
colnames(MA.GGMupF) <- c("MA", "NotMA")
rownames(MA.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GGMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many mito-annotated transcripts are there?
GGMdownsigMA <- sum(GGMdownsig$gene.id %in% MAmwap7g$From)
# 2
# How many are not?
GGMdownsigNotMA <- length(GGMdownsig$gene.id)-GGMdownsigMA
# 20

# The number of expected mito-annotated transcripts
length(GGMdownsig$gene.id)*MApropM
# 0.6569494


# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigMA, GGMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.09251

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.GGMdownF <- rbind(c(GGMdownsigMA, GGMdownsigNotMA), 
                     c(length(MAmwap7g$From)-GGMdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-GGMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.GGMdownF) <- c("MA", "NotMA")
rownames(MA.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.GGMdownF)
# P = 0.1391

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group 2
G2Mupsig <- subset(FLXswap7gM.data, g2up == "Sig")
# 29

# How many mito-annotated transcripts are there?
G2MupsigMA <- sum(G2Mupsig$gene.id %in% MAmwap7g$From)
# 1
# How many are not?
G2MupsigNotMA <- length(G2Mupsig$gene.id)-G2MupsigMA
# 28

# The number of expected mito-annotated transcripts
length(G2Mupsig$gene.id)*MApropM
# 0.8085791

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MupsigMA, G2MupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.8291

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G2MupF <- rbind(c(G2MupsigMA, G2MupsigNotMA), 
                   c(length(MAmwap7g$From)-G2MupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-G2MupsigNotMA))
#Add names to the columns and rows
colnames(MA.G2MupF) <- c("MA", "NotMA")
rownames(MA.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G2MupF)
# P = 0.5599

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group 2
G2Mdownsig <- subset(FLXswap7gM.data, g2down == "Sig")
# 36

# How many mito-annotated transcripts are there?
G2MdownsigMA <- sum(G2Mdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
G2MdownsigNotMA <- length(G2Mdownsig$gene.id)-G2MdownsigMA
# 36

# The number of expected mito-annotated transcripts
length(G2Mdownsig$gene.id)*MApropM
# 1.003753

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MdownsigMA, G2MdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.3096

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G2MdownF <- rbind(c(G2MdownsigMA, G2MdownsigNotMA), 
                     c(length(MAmwap7g$From)-G2MdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-G2MdownsigNotMA))
#Add names to the columns and rows
colnames(MA.G2MdownF) <- c("MA", "NotMA")
rownames(MA.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G2MdownF)
# P = 0.6265

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in group 11
G11Mupsig <- subset(FLXswap7gM.data, g11up == "Sig")
# 21

# How many mito-annotated transcripts are there?
G11MupsigMA <- sum(G11Mupsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
G11MupsigNotMA <- length(G11Mupsig$gene.id)-G11MupsigMA
# 21

# The number of expected mito-annotated transcripts
length(G11Mupsig$gene.id)*MApropM
# 0.5855228

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MupsigMA, G11MupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.4377

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G11MupF <- rbind(c(G11MupsigMA, G11MupsigNotMA), 
                    c(length(MAmwap7g$From)-G11MupsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-G11MupsigNotMA))
#Add names to the columns and rows
colnames(MA.G11MupF) <- c("MA", "NotMA")
rownames(MA.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G11MupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in group 11
G11Mdownsig <- subset(FLXswap7gM.data, g11down == "Sig")
# 7

# How many mito-annotated transcripts are there?
G11MdownsigMA <- sum(G11Mdownsig$gene.id %in% MAmwap7g$From)
# 0
# How many are not?
G11MdownsigNotMA <- length(G11Mdownsig$gene.id)-G11MdownsigMA
# 7

# The number of expected mito-annotated transcripts
length(G11Mdownsig$gene.id)*MApropM
# 0.1951743

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MdownsigMA, G11MdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.6541

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.G11MdownF <- rbind(c(G11MdownsigMA, G11MdownsigNotMA), 
                      c(length(MAmwap7g$From)-G11MdownsigMA, (length(FLXswap7gM.data$gene.id)-length(MAmwap7g$From))-G11MdownsigNotMA))
#Add names to the columns and rows
colnames(MA.G11MdownF) <- c("MA", "NotMA")
rownames(MA.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.G11MdownF)
# P = 1

# We get the same results, so that is good

