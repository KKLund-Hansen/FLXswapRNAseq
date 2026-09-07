################################################################################################
##################### MITO-SENSITIVE TRANSCRIPTS FLX SWAP MALES - 7 GROUPS #####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the male data set
MSmswap7g <- subset(MS.data, MS.data$gene.id %in% FLXswap7gM.data$gene.id)
# 1210 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropM <- length(MSmswap7g$gene.id)/length(FLXswap7gM.data$gene.id)
# 0.08

# And for the chisq.test I also have to calculate
NotMSpropM <- 1-MSpropM
# 0.92




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many mito-sensitive transcripts are there?
GAMupsigMS <- sum(GAMupsig$gene.id %in% MSmswap7g$gene.id)
# 0
# How many are not?
GAMupsigNotMS <- length(GAMupsig$gene.id)-GAMupsigMS
# 5

# The number of expected mito-sensitive transcripts
length(GAMupsig$gene.id)*MSpropM
# 0.3901967


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMupsigMS, GAMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5153

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GAMupF <- rbind(c(GAMupsigMS, GAMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GAMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GAMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GAMupF) <- c("MS", "NotMS")
rownames(MS.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many mito-sensitive transcripts are there?
GAMdownsigMS <- sum(GAMdownsig$gene.id %in% MSmswap7g$gene.id)
# 2
# How many are not?
GAMdownsigNotMS <- length(GAMdownsig$gene.id)-GAMdownsigMS
# 14

# The number of expected mito-sensitive transcripts
length(GAMdownsig$gene.id)*MSpropM
# 1.248629


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigMS, GAMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4837

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GAMdownF <- rbind(c(GAMdownsigMS, GAMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GAMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GAMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GAMdownF) <- c("MS", "NotMS")
rownames(MS.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GAMdownF)
# P = 0.3584

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many mito-sensitive transcripts are there?
GBMupsigMS <- sum(GBMupsig$gene.id %in% MSmswap7g$gene.id)
# 3
# How many are not?
GBMupsigNotMS <- length(GBMupsig$gene.id)-GBMupsigMS
# 48

# The number of expected mito-sensitive transcripts
length(GBMupsig$gene.id)*MSpropM
# 3.980006


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMupsigMS, GBMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.6089

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GBMupF <- rbind(c(GBMupsigMS, GBMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GBMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GBMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GBMupF) <- c("MS", "NotMS")
rownames(MS.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GBMupF)
# P = 0.7962

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many mito-sensitive transcripts are there?
GBMdownsigMS <- sum(GBMdownsig$gene.id %in% MSmswap7g$gene.id)
# 3
# How many are not?
GBMdownsigNotMS <- length(GBMdownsig$gene.id)-GBMdownsigMS
# 53

# The number of expected mito-sensitive transcripts
length(GBMdownsig$gene.id)*MSpropM
# 4.370203


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigMS, GBMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4948

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GBMdownF <- rbind(c(GBMdownsigMS, GBMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GBMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GBMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GBMdownF) <- c("MS", "NotMS")
rownames(MS.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GBMdownF)
# P = 0.8005

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many mito-sensitive transcripts are there?
GCMupsigMS <- sum(GCMupsig$gene.id %in% MSmswap7g$gene.id)
# 2
# How many are not?
GCMupsigNotMS <- length(GCMupsig$gene.id)-GCMupsigMS
# 25

# The number of expected mito-sensitive transcripts
length(GCMupsig$gene.id)*MSpropM
# 2.107062


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMupsigMS, GCMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9388

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GCMupF <- rbind(c(GCMupsigMS, GCMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GCMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GCMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GCMupF) <- c("MS", "NotMS")
rownames(MS.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GCMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many mito-sensitive transcripts are there?
GCMdownsigMS <- sum(GCMdownsig$gene.id %in% MSmswap7g$gene.id)
# 3
# How many are not?
GCMdownsigNotMS <- length(GCMdownsig$gene.id)-GCMdownsigMS
# 35

# The number of expected mito-sensitive transcripts
length(GCMdownsig$gene.id)*MSpropM
# 2.965495


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigMS, GCMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9834

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GCMdownF <- rbind(c(GCMdownsigMS, GCMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GCMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GCMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GCMdownF) <- c("MS", "NotMS")
rownames(MS.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GCMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many mito-sensitive transcripts are there?
GDMupsigMS <- sum(GDMupsig$gene.id %in% MSmswap7g$gene.id)
# 1
# How many are not?
GDMupsigNotMS <- length(GDMupsig$gene.id)-GDMupsigMS
# 11

# The number of expected mito-sensitive transcripts
length(GDMupsig$gene.id)*MSpropM
# 0.9364721


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMupsigMS, GDMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.9455

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GDMupF <- rbind(c(GDMupsigMS, GDMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GDMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GDMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GDMupF) <- c("MS", "NotMS")
rownames(MS.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many mito-sensitive transcripts are there?
GDMdownsigMS <- sum(GDMdownsig$gene.id %in% MSmswap7g$gene.id)
# 2
# How many are not?
GDMdownsigNotMS <- length(GDMdownsig$gene.id)-GDMdownsigMS
# 14

# The number of expected mito-sensitive transcripts
length(GDMdownsig$gene.id)*MSpropM
# 1.248629


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigMS, GDMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4837

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GDMdownF <- rbind(c(GDMdownsigMS, GDMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GDMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GDMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GDMdownF) <- c("MS", "NotMS")
rownames(MS.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GDMdownF)
# P = 0.3584

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many mito-sensitive transcripts are there?
GEMupsigMS <- sum(GEMupsig$gene.id %in% MSmswap7g$gene.id)
# 0
# How many are not?
GEMupsigNotMS <- length(GEMupsig$gene.id)-GEMupsigMS
# 20

# The number of expected mito-sensitive transcripts
length(GEMupsig$gene.id)*MSpropM
# 1.560787


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMupsigMS, GEMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.1932

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GEMupF <- rbind(c(GEMupsigMS, GEMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GEMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GEMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GEMupF) <- c("MS", "NotMS")
rownames(MS.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GEMupF)
# P = 0.3983

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many mito-sensitive transcripts are there?
GEMdownsigMS <- sum(GEMdownsig$gene.id %in% MSmswap7g$gene.id)
# 1
# How many are not?
GEMdownsigNotMS <- length(GEMdownsig$gene.id)-GEMdownsigMS
# 25

# The number of expected mito-sensitive transcripts
length(GEMdownsig$gene.id)*MSpropM
# 2.029023


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigMS, GEMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.4518

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GEMdownF <- rbind(c(GEMdownsigMS, GEMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GEMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GEMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GEMdownF) <- c("MS", "NotMS")
rownames(MS.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GEMdownF)
# P = 0.7182

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many mito-sensitive transcripts are there?
GFMupsigMS <- sum(GFMupsig$gene.id %in% MSmswap7g$gene.id)
# 0
# How many are not?
GFMupsigNotMS <- length(GFMupsig$gene.id)-GFMupsigMS
# 2

# The number of expected mito-sensitive transcripts
length(GFMupsig$gene.id)*MSpropM
# 0.1560787


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMupsigMS, GFMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.6807

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFMupF <- rbind(c(GFMupsigMS, GFMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GFMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GFMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GFMupF) <- c("MS", "NotMS")
rownames(MS.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many mito-sensitive transcripts are there?
GFMdownsigMS <- sum(GFMdownsig$gene.id %in% MSmswap7g$gene.id)
# 2
# How many are not?
GFMdownsigNotMS <- length(GFMdownsig$gene.id)-GFMdownsigMS
# 1

# The number of expected mito-sensitive transcripts
length(GFMdownsig$gene.id)*MSpropM
# 0.234118


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigMS, GFMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.0001442

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFMdownF <- rbind(c(GFMdownsigMS, GFMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GFMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GFMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GFMdownF) <- c("MS", "NotMS")
rownames(MS.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFMdownF)
# P = 0.01731

# It goes from significant to non-significant




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many mito-sensitive transcripts are there?
GGMupsigMS <- sum(GGMupsig$gene.id %in% MSmswap7g$gene.id)
# 1
# How many are not?
GGMupsigNotMS <- length(GGMupsig$gene.id)-GGMupsigMS
# 9

# The number of expected mito-sensitive transcripts
length(GGMupsig$gene.id)*MSpropM
# 0.7803934


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMupsigMS, GGMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.7957

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GGMupF <- rbind(c(GGMupsigMS, GGMupsigNotMS), 
                   c(length(MSmswap7g$gene.id)-GGMupsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GGMupsigNotMS))
#Add names to the columns and rows
colnames(MS.GGMupF) <- c("MS", "NotMS")
rownames(MS.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GGMupF)
# P = 0.5564

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many mito-sensitive transcripts are there?
GGMdownsigMS <- sum(GGMdownsig$gene.id %in% MSmswap7g$gene.id)
# 1
# How many are not?
GGMdownsigNotMS <- length(GGMdownsig$gene.id)-GGMdownsigMS
# 21

# The number of expected mito-sensitive transcripts
length(GGMdownsig$gene.id)*MSpropM
# 1.716866


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigMS, GGMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5688

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GGMdownF <- rbind(c(GGMdownsigMS, GGMdownsigNotMS), 
                     c(length(MSmswap7g$gene.id)-GGMdownsigMS, (length(FLXswap7gM.data$gene.id)-length(MSmswap7g$gene.id))-GGMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.GGMdownF) <- c("MS", "NotMS")
rownames(MS.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GGMdownF)
# P = 1

# We get the same results, so that is good

