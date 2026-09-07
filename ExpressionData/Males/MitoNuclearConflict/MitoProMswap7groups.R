################################################################################################
###################### MITO-PROTEOME TRANSCRIPTS FLX SWAP MALES - 7 GROUPS #####################
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
MP.data <- read.table(file = "mitoproteome.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# There are some gene IDs that are multiple in the MP.data so first I'll remove the doubles
str(MP.data)
MP.data <- unique(MP.data)

# Next add the CG names as row names
rownames(MP.data) <- MP.data[,1]

# Next add a new column with the flybase names
MP.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MP.data), keytype = "FLYBASECG", column = "FLYBASE")

# Now I can finally subset the mito-proteome data set to match the male data set
MPmswap7g <- subset(MP.data, MP.data$gene.id %in% FLXswap7gM.data$gene.id)
# 762 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcriptsout of all my transcripts
MPpropM <- length(MPmswap7g$gene.id)/length(FLXswap7gM.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropM <- 1-MPpropM
# 0.95




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# How many mito-proteome transcripts are there?
GAMupsigMP <- sum(GAMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GAMupsigNotMP <- length(GAMupsig$gene.id)-GAMupsigMP
# 5

# The number of expected mito-proteome transcripts
length(GAMupsig$gene.id)*MPpropM
# 0.2457272


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMupsigMP, GAMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.6112

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GAMupF <- rbind(c(GAMupsigMP, GAMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GAMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GAMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GAMupF) <- c("MP", "NotMP")
rownames(MP.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# How many mito-proteome transcripts are there?
GAMdownsigMP <- sum(GAMdownsig$gene.id %in% MPmswap7g$gene.id)
# 1
# How many are not?
GAMdownsigNotMP <- length(GAMdownsig$gene.id)-GAMdownsigMP
# 15

# The number of expected mito-proteome transcripts
length(GAMdownsig$gene.id)*MPpropM
# 0.786327


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigMP, GAMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.8048

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GAMdownF <- rbind(c(GAMdownsigMP, GAMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GAMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GAMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GAMdownF) <- c("MP", "NotMP")
rownames(MP.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GAMdownF)
# P = 0.5537

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# How many mito-proteome transcripts are there?
GBMupsigMP <- sum(GBMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GBMupsigNotMP <- length(GBMupsig$gene.id)-GBMupsigMP
# 51

# The number of expected mito-proteome transcripts
length(GBMupsig$gene.id)*MPpropM
# 2.506417


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMupsigMP, GBMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.1045

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GBMupF <- rbind(c(GBMupsigMP, GBMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GBMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GBMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GBMupF) <- c("MP", "NotMP")
rownames(MP.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GBMupF)
# P = 0.1804

# We get the same results, so that is good




##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# How many mito-proteome transcripts are there?
GBMdownsigMP <- sum(GBMdownsig$gene.id %in% MPmswap7g$gene.id)
# 9
# How many are not?
GBMdownsigNotMP <- length(GBMdownsig$gene.id)-GBMdownsigMP
# 47

# The number of expected mito-proteome transcripts
length(GBMdownsig$gene.id)*MPpropM
# 2.752144


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigMP, GBMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.0001124

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GBMdownF <- rbind(c(GBMdownsigMP, GBMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GBMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GBMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GBMdownF) <- c("MP", "NotMP")
rownames(MP.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GBMdownF)
# P = 0.001514

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# How many mito-proteome transcripts are there?
GCMupsigMP <- sum(GCMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GCMupsigNotMP <- length(GCMupsig$gene.id)-GCMupsigMP
# 27

# The number of expected mito-proteome transcripts
length(GCMupsig$gene.id)*MPpropM
# 1.326927


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMupsigMP, GCMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.2375

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GCMupF <- rbind(c(GCMupsigMP, GCMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GCMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GCMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GCMupF) <- c("MP", "NotMP")
rownames(MP.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GCMupF)
# P = 0.6418

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# How many mito-proteome transcripts are there?
GCMdownsigMP <- sum(GCMdownsig$gene.id %in% MPmswap7g$gene.id)
# 1
# How many are not?
GCMdownsigNotMP <- length(GCMdownsig$gene.id)-GCMdownsigMP
# 37

# The number of expected mito-proteome transcripts
length(GCMdownsig$gene.id)*MPpropM
# 1.867527


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigMP, GCMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.515

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GCMdownF <- rbind(c(GCMdownsigMP, GCMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GCMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GCMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GCMdownF) <- c("MP", "NotMP")
rownames(MP.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GCMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# How many mito-proteome transcripts are there?
GDMupsigMP <- sum(GDMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GDMupsigNotMP <- length(GDMupsig$gene.id)-GDMupsigMP
# 12

# The number of expected mito-proteome transcripts
length(GDMupsig$gene.id)*MPpropM
# 0.5897452


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMupsigMP, GDMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.431

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GDMupF <- rbind(c(GDMupsigMP, GDMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GDMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GDMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GDMupF) <- c("MP", "NotMP")
rownames(MP.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# How many mito-proteome transcripts are there?
GDMdownsigMP <- sum(GDMdownsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GDMdownsigNotMP <- length(GDMdownsig$gene.id)-GDMdownsigMP
# 16

# The number of expected mito-proteome transcripts
length(GDMdownsig$gene.id)*MPpropM
# 0.786327


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigMP, GDMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.3632

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GDMdownF <- rbind(c(GDMdownsigMP, GDMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GDMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GDMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GDMdownF) <- c("MP", "NotMP")
rownames(MP.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GDMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# How many mito-proteome transcripts are there?
GEMupsigMP <- sum(GEMupsig$gene.id %in% MPmswap7g$gene.id)
# 1
# How many are not?
GEMupsigNotMP <- length(GEMupsig$gene.id)-GEMupsigMP
# 19

# The number of expected mito-proteome transcripts
length(GEMupsig$gene.id)*MPpropM
# 0.9829087


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMupsigMP, GEMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.9859

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEMupF <- rbind(c(GEMupsigMP, GEMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GEMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GEMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GEMupF) <- c("MP", "NotMP")
rownames(MP.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# How many mito-proteome transcripts are there?
GEMdownsigMP <- sum(GEMdownsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GEMdownsigNotMP <- length(GEMdownsig$gene.id)-GEMdownsigMP
# 26

# The number of expected mito-proteome transcripts
length(GEMdownsig$gene.id)*MPpropM
# 1.277781


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigMP, GEMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.2464

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEMdownF <- rbind(c(GEMdownsigMP, GEMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GEMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GEMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GEMdownF) <- c("MP", "NotMP")
rownames(MP.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEMdownF)
# P = 0.6373

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# How many mito-proteome transcripts are there?
GFMupsigMP <- sum(GFMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GFMupsigNotMP <- length(GFMupsig$gene.id)-GFMupsigMP
# 2

# The number of expected mito-proteome transcripts
length(GFMupsig$gene.id)*MPpropM
# 0.09829087


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMupsigMP, GFMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.7478

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFMupF <- rbind(c(GFMupsigMP, GFMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GFMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GFMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GFMupF) <- c("MP", "NotMP")
rownames(MP.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# How many mito-proteome transcripts are there?
GFMdownsigMP <- sum(GFMdownsig$gene.id %in% MPmswap7g$gene.id)
# 2
# How many are not?
GFMdownsigNotMP <- length(GFMdownsig$gene.id)-GFMdownsigMP
# 1

# The number of expected mito-proteome transcripts
length(GFMdownsig$gene.id)*MPpropM
# 0.1474363


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigMP, GFMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 7.505e-07

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFMdownF <- rbind(c(GFMdownsigMP, GFMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GFMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GFMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GFMdownF) <- c("MP", "NotMP")
rownames(MP.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFMdownF)
# P = 0.007

# It goes from significant to non-significant




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# How many mito-proteome transcripts are there?
GGMupsigMP <- sum(GGMupsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GGMupsigNotMP <- length(GGMupsig$gene.id)-GGMupsigMP
# 10

# The number of expected mito-proteome transcripts
length(GGMupsig$gene.id)*MPpropM
# 0.4914544


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMupsigMP, GGMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.4722

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGMupF <- rbind(c(GGMupsigMP, GGMupsigNotMP), 
                   c(length(MPmswap7g$gene.id)-GGMupsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GGMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GGMupF) <- c("MP", "NotMP")
rownames(MP.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# How many mito-proteome transcripts are there?
GGMdownsigMP <- sum(GGMdownsig$gene.id %in% MPmswap7g$gene.id)
# 0
# How many are not?
GGMdownsigNotMP <- length(GGMdownsig$gene.id)-GGMdownsigMP
# 22

# The number of expected mito-proteome transcripts
length(GGMdownsig$gene.id)*MPpropM
# 1.0812


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigMP, GGMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.2863

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGMdownF <- rbind(c(GGMdownsigMP, GGMdownsigNotMP), 
                     c(length(MPmswap7g$gene.id)-GGMdownsigMP, (length(FLXswap7gM.data$gene.id)-length(MPmswap7g$gene.id))-GGMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GGMdownF) <- c("MP", "NotMP")
rownames(MP.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGMdownF)
# P = 0.6245

# We get the same results, so that is good

