################################################################################################
########################### MITO-PROTEOME TRANSCRIPTS FLX SWAP MALES ###########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
MPmaleswap <- subset(MP.data, MP.data$gene.id %in% RNAswapmale.data$gene.id)
# 767 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcriptsout of all my transcripts
MPpropM <- length(MPmaleswap$gene.id)/length(RNAswapmale.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropM <- 1-MPpropM
# 0.95




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# How many mito-proteome transcripts are there?
GAMupsigMP <- sum(GAMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GAMupsigNotMP <- length(GAMupsig$gene.id)-GAMupsigMP
# 4

# The number of expected mito-proteome transcripts
length(GAMupsig$gene.id)*MPpropM
# 0.1827822

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMupsigMP, GAMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.6616

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GAMupF <- rbind(c(GAMupsigMP, GAMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GAMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GAMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GAMupF) <- c("MP", "NotMP")
rownames(MP.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# How many mito-proteome transcripts are there?
GAMdownsigMP <- sum(GAMdownsig$gene.id %in% MPmaleswap$gene.id)
# 1
# How many are not?
GAMdownsigNotMP <- length(GAMdownsig$gene.id)-GAMdownsigMP
# 15

# The number of expected mito-proteome transcripts
length(GAMdownsig$gene.id)*MPpropM
# 0.731129

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAMdownsigMP, GAMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.7475

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GAMdownF <- rbind(c(GAMdownsigMP, GAMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GAMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GAMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GAMdownF) <- c("MP", "NotMP")
rownames(MP.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GAMdownF)
# P = 0.527

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# How many mito-proteome transcripts are there?
GBMupsigMP <- sum(GBMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GBMupsigNotMP <- length(GBMupsig$gene.id)-GBMupsigMP
# 54

# The number of expected mito-proteome transcripts
length(GBMupsig$gene.id)*MPpropM
# 2.46756

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMupsigMP, GBMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.1078

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GBMupF <- rbind(c(GBMupsigMP, GBMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GBMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GBMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GBMupF) <- c("MP", "NotMP")
rownames(MP.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GBMupF)
# P = 0.1791

# We get the same results, so that is good




##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# How many mito-proteome transcripts are there?
GBMdownsigMP <- sum(GBMdownsig$gene.id %in% MPmaleswap$gene.id)
# 9
# How many are not?
GBMdownsigNotMP <- length(GBMdownsig$gene.id)-GBMdownsigMP
# 49

# The number of expected mito-proteome transcripts
length(GBMdownsig$gene.id)*MPpropM
# 2.650343

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBMdownsigMP, GBMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 6.535e-05

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GBMdownF <- rbind(c(GBMdownsigMP, GBMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GBMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GBMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GBMdownF) <- c("MP", "NotMP")
rownames(MP.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GBMdownF)
# P = 0.001181

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# How many mito-proteome transcripts are there?
GCMupsigMP <- sum(GCMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GCMupsigNotMP <- length(GCMupsig$gene.id)-GCMupsigMP
# 27

# The number of expected mito-proteome transcripts
length(GCMupsig$gene.id)*MPpropM
# 1.23378

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMupsigMP, GCMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.2555

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GCMupF <- rbind(c(GCMupsigMP, GCMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GCMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GCMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GCMupF) <- c("MP", "NotMP")
rownames(MP.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GCMupF)
# P = 0.6341

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# How many mito-proteome transcripts are there?
GCMdownsigMP <- sum(GCMdownsig$gene.id %in% MPmaleswap$gene.id)
# 1
# How many are not?
GCMdownsigNotMP <- length(GCMdownsig$gene.id)-GCMdownsigMP
# 36

# The number of expected mito-proteome transcripts
length(GCMdownsig$gene.id)*MPpropM
# 1.690736

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCMdownsigMP, GCMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.5866

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GCMdownF <- rbind(c(GCMdownsigMP, GCMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GCMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GCMdownsigNotMP))
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
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# How many mito-proteome transcripts are there?
GDMupsigMP <- sum(GDMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GDMupsigNotMP <- length(GDMupsig$gene.id)-GDMupsigMP
# 14

# The number of expected mito-proteome transcripts
length(GDMupsig$gene.id)*MPpropM
# 0.6397379

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMupsigMP, GDMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.4129

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GDMupF <- rbind(c(GDMupsigMP, GDMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GDMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GDMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GDMupF) <- c("MP", "NotMP")
rownames(MP.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GDMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# How many mito-proteome transcripts are there?
GDMdownsigMP <- sum(GDMdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GDMdownsigNotMP <- length(GDMdownsig$gene.id)-GDMdownsigMP
# 16

# The number of expected mito-proteome transcripts
length(GDMdownsig$gene.id)*MPpropM
# 0.731129

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDMdownsigMP, GDMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.3814

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GDMdownF <- rbind(c(GDMdownsigMP, GDMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GDMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GDMdownsigNotMP))
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
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# How many mito-proteome transcripts are there?
GEMupsigMP <- sum(GEMupsig$gene.id %in% MPmaleswap$gene.id)
# 1
# How many are not?
GEMupsigNotMP <- length(GEMupsig$gene.id)-GEMupsigMP
# 16

# The number of expected mito-proteome transcripts
length(GEMupsig$gene.id)*MPpropM
# 0.7768245

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMupsigMP, GEMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.7955

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEMupF <- rbind(c(GEMupsigMP, GEMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GEMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GEMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GEMupF) <- c("MP", "NotMP")
rownames(MP.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEMupF)
# P = 0.5487

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# How many mito-proteome transcripts are there?
GEMdownsigMP <- sum(GEMdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GEMdownsigNotMP <- length(GEMdownsig$gene.id)-GEMdownsigMP
# 26

# The number of expected mito-proteome transcripts
length(GEMdownsig$gene.id)*MPpropM
# 1.188085

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEMdownsigMP, GEMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.2645

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEMdownF <- rbind(c(GEMdownsigMP, GEMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GEMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GEMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GEMdownF) <- c("MP", "NotMP")
rownames(MP.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEMdownF)
# P = 0.6308

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# How many mito-proteome transcripts are there?
GFMupsigMP <- sum(GFMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GFMupsigNotMP <- length(GFMupsig$gene.id)-GFMupsigMP
# 2

# The number of expected mito-proteome transcripts
length(GFMupsig$gene.id)*MPpropM
# 0.09139112

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMupsigMP, GFMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.757

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFMupF <- rbind(c(GFMupsigMP, GFMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GFMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GFMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GFMupF) <- c("MP", "NotMP")
rownames(MP.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# How many mito-proteome transcripts are there?
GFMdownsigMP <- sum(GFMdownsig$gene.id %in% MPmaleswap$gene.id)
# 2
# How many are not?
GFMdownsigNotMP <- length(GFMdownsig$gene.id)-GFMdownsigMP
# 1

# The number of expected mito-proteome transcripts
length(GFMdownsig$gene.id)*MPpropM
# 0.1370867

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFMdownsigMP, GFMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 2.598e-07

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFMdownF <- rbind(c(GFMdownsigMP, GFMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GFMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GFMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GFMdownF) <- c("MP", "NotMP")
rownames(MP.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFMdownF)
# P = 0.006066

# It goes from significant to non-significant




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# How many mito-proteome transcripts are there?
GGMupsigMP <- sum(GGMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GGMupsigNotMP <- length(GGMupsig$gene.id)-GGMupsigMP
# 10

# The number of expected mito-proteome transcripts
length(GGMupsig$gene.id)*MPpropM
# 0.4569556

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMupsigMP, GGMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.4889

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGMupF <- rbind(c(GGMupsigMP, GGMupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-GGMupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GGMupsigNotMP))
#Add names to the columns and rows
colnames(MP.GGMupF) <- c("MP", "NotMP")
rownames(MP.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# How many mito-proteome transcripts are there?
GGMdownsigMP <- sum(GGMdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
GGMdownsigNotMP <- length(GGMdownsig$gene.id)-GGMdownsigMP
# 19

# The number of expected mito-proteome transcripts
length(GGMdownsig$gene.id)*MPpropM
# 0.8682157

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGMdownsigMP, GGMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.3402

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGMdownF <- rbind(c(GGMdownsigMP, GGMdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-GGMdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-GGMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.GGMdownF) <- c("MP", "NotMP")
rownames(MP.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGMdownF)
# P = 1

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# How many mito-proteome transcripts are there?
G2MupsigMP <- sum(G2Mupsig$gene.id %in% MPmaleswap$gene.id)
# 2
# How many are not?
G2MupsigNotMP <- length(G2Mupsig$gene.id)-G2MupsigMP
# 27

# The number of expected mito-proteome transcripts
length(G2Mupsig$gene.id)*MPpropM
# 1.325171

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MupsigMP, G2MupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.5484

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G2MupF <- rbind(c(G2MupsigMP, G2MupsigNotMP), 
                   c(length(MPmaleswap$gene.id)-G2MupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-G2MupsigNotMP))
#Add names to the columns and rows
colnames(MP.G2MupF) <- c("MP", "NotMP")
rownames(MP.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G2MupF)
# P = 0.3848

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# How many mito-proteome transcripts are there?
G2MdownsigMP <- sum(G2Mdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
G2MdownsigNotMP <- length(G2Mdownsig$gene.id)-G2MdownsigMP
# 36

# The number of expected mito-proteome transcripts
length(G2Mdownsig$gene.id)*MPpropM
# 1.64504

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G2MdownsigMP, G2MdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.1892

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G2MdownF <- rbind(c(G2MdownsigMP, G2MdownsigNotMP), 
                     c(length(MPmaleswap$gene.id)-G2MdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-G2MdownsigNotMP))
#Add names to the columns and rows
colnames(MP.G2MdownF) <- c("MP", "NotMP")
rownames(MP.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G2MdownF)
# P = 0.4113

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# How many mito-proteome transcripts are there?
G11MupsigMP <- sum(G11Mupsig$gene.id %in% MPmaleswap$gene.id)
# 3
# How many are not?
G11MupsigNotMP <- length(G11Mupsig$gene.id)-G11MupsigMP
# 18

# The number of expected mito-proteome transcripts
length(G11Mupsig$gene.id)*MPpropM
# 0.9596068

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MupsigMP, G11MupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.03299

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G11MupF <- rbind(c(G11MupsigMP, G11MupsigNotMP), 
                    c(length(MPmaleswap$gene.id)-G11MupsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-G11MupsigNotMP))
#Add names to the columns and rows
colnames(MP.G11MupF) <- c("MP", "NotMP")
rownames(MP.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G11MupF)
# P = 0.06854

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# How many mito-proteome transcripts are there?
G11MdownsigMP <- sum(G11Mdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
G11MdownsigNotMP <- length(G11Mdownsig$gene.id)-G11MdownsigMP
# 7

# The number of expected mito-proteome transcripts
length(G11Mdownsig$gene.id)*MPpropM
# 0.3198689

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G11MdownsigMP, G11MdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.5626

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G11MdownF <- rbind(c(G11MdownsigMP, G11MdownsigNotMP), 
                      c(length(MPmaleswap$gene.id)-G11MdownsigMP, (length(RNAswapmale.data$gene.id)-length(MPmaleswap$gene.id))-G11MdownsigNotMP))
#Add names to the columns and rows
colnames(MP.G11MdownF) <- c("MP", "NotMP")
rownames(MP.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G11MdownF)
# P = 1

# We get the same results, so that is good


