################################################################################################
########################## MITO-PROTEOME TRANSCRIPTS FLX SWAP FEMALES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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

# Now I can finally subset the mito-proteome data set to match the female data set
MPfemaleswap <- subset(MP.data, MP.data$gene.id %in% RNAswapfemale.data$gene.id)
# 767 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcriptsout of all my transcripts
MPpropF <- length(MPfemaleswap$gene.id)/length(RNAswapfemale.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropF <- 1-MPpropF
# 0.95




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# How many mito-proteome transcripts are there?
GAFupsigMP <- sum(GAFupsig$gene.id %in% MPfemaleswap$gene.id)
# 7
# How many are not?
GAFupsigNotMP <- length(GAFupsig$gene.id)-GAFupsigMP
# 122

# The number of expected mito-proteome transcripts
length(GAFupsig$gene.id)*MPpropF
# 5.894727

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMP, GAFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.6412



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# How many mito-proteome transcripts are there?
GAFdownsigMP <- sum(GAFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 91
# How many are not?
GAFdownsigNotMP <- length(GAFdownsig$gene.id)-GAFdownsigMP
# 1930

# The number of expected mito-proteome transcripts
length(GAFdownsig$gene.id)*MPpropF
# 92.35073

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMP, GAFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.8856




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# How many mito-proteome transcripts are there?
GBFupsigMP <- sum(GBFupsig$gene.id %in% MPfemaleswap$gene.id)
# 62
# How many are not?
GBFupsigNotMP <- length(GBFupsig$gene.id)-GBFupsigMP
# 508

# The number of expected mito-proteome transcripts
length(GBFupsig$gene.id)*MPpropF
# 26.04647

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMP, GBFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 5.535e-13



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# How many mito-proteome transcripts are there?
GBFdownsigMP <- sum(GBFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 32
# How many are not?
GBFdownsigNotMP <- length(GBFdownsig$gene.id)-GBFdownsigMP
# 350

# The number of expected mito-proteome transcripts
length(GBFdownsig$gene.id)*MPpropF
# 17.4557

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMP, GBFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.0003659




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# How many mito-proteome transcripts are there?
GCFupsigMP <- sum(GCFupsig$gene.id %in% MPfemaleswap$gene.id)
# 24
# How many are not?
GCFupsigNotMP <- length(GCFupsig$gene.id)-GCFupsigMP
# 153

# The number of expected mito-proteome transcripts
length(GCFupsig$gene.id)*MPpropF
# 8.088114

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMP, GCFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 1.02e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# How many mito-proteome transcripts are there?
GCFdownsigMP <- sum(GCFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 15
# How many are not?
GCFdownsigNotMP <- length(GCFdownsig$gene.id)-GCFdownsigMP
# 357

# The number of expected mito-proteome transcripts
length(GCFdownsig$gene.id)*MPpropF
# 16.99875

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMP, GCFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.6197




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# How many mito-proteome transcripts are there?
GDFupsigMP <- sum(GDFupsig$gene.id %in% MPfemaleswap$gene.id)
# 33
# How many are not?
GDFupsigNotMP <- length(GDFupsig$gene.id)-GDFupsigMP
# 158

# The number of expected mito-proteome transcripts
length(GDFupsig$gene.id)*MPpropF
# 8.727852

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMP, GDFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# How many mito-proteome transcripts are there?
GDFdownsigMP <- sum(GDFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 11
# How many are not?
GDFdownsigNotMP <- length(GDFdownsig$gene.id)-GDFdownsigMP
# 273

# The number of expected mito-proteome transcripts
length(GDFdownsig$gene.id)*MPpropF
# 12.97754

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMP, GDFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.5742




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# How many mito-proteome transcripts are there?
GEFupsigMP <- sum(GEFupsig$gene.id %in% MPfemaleswap$gene.id)
# 5
# How many are not?
GEFupsigNotMP <- length(GEFupsig$gene.id)-GEFupsigMP
# 93

# The number of expected mito-proteome transcripts
length(GEFupsig$gene.id)*MPpropF
# 4.478165

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMP, GEFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.8007

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEFupF <- rbind(c(GEFupsigMP, GEFupsigNotMP), 
                   c(length(MPfemaleswap$gene.id)-GEFupsigMP, (length(RNAswapfemale.data$gene.id)-length(MPfemaleswap$gene.id))-GEFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GEFupF) <- c("MP", "NotMP")
rownames(MP.GEFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEFupF)
# P = 0.8055

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# How many mito-proteome transcripts are there?
GEFdownsigMP <- sum(GEFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 18
# How many are not?
GEFdownsigNotMP <- length(GEFdownsig$gene.id)-GEFdownsigMP
# 246

# The number of expected mito-proteome transcripts
length(GEFdownsig$gene.id)*MPpropF
# 12.06363

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMP, GEFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.08019




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# How many mito-proteome transcripts are there?
GFFupsigMP <- sum(GFFupsig$gene.id %in% MPfemaleswap$gene.id)
# 4
# How many are not?
GFFupsigNotMP <- length(GFFupsig$gene.id)-GFFupsigMP
# 36

# The number of expected mito-proteome transcripts
length(GFFupsig$gene.id)*MPpropF
# 1.827822

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMP, GFFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.1

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFFupF <- rbind(c(GFFupsigMP, GFFupsigNotMP), 
                   c(length(MPfemaleswap$gene.id)-GFFupsigMP, (length(RNAswapfemale.data$gene.id)-length(MPfemaleswap$gene.id))-GFFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GFFupF) <- c("MP", "NotMP")
rownames(MP.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFFupF)
# P = 0.1084

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# How many mito-proteome transcripts are there?
GFFdownsigMP <- sum(GFFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 10
# How many are not?
GFFdownsigNotMP <- length(GFFdownsig$gene.id)-GFFdownsigMP
# 254

# The number of expected mito-proteome transcripts
length(GFFdownsig$gene.id)*MPpropF
# 12.06363

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMP, GFFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.5431




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# How many mito-proteome transcripts are there?
GGFupsigMP <- sum(GGFupsig$gene.id %in% MPfemaleswap$gene.id)
# 7
# How many are not?
GGFupsigNotMP <- length(GGFupsig$gene.id)-GGFupsigMP
# 85

# The number of expected mito-proteome transcripts
length(GGFupsig$gene.id)*MPpropF
# 4.203992

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMP, GGFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.1627

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGFupF <- rbind(c(GGFupsigMP, GGFupsigNotMP), 
                   c(length(MPfemaleswap$gene.id)-GGFupsigMP, (length(RNAswapfemale.data$gene.id)-length(MPfemaleswap$gene.id))-GGFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GGFupF) <- c("MP", "NotMP")
rownames(MP.GGFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGFupF)
# P = 0.2001

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# How many mito-proteome transcripts are there?
GGFdownsigMP <- sum(GGFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 22
# How many are not?
GGFdownsigNotMP <- length(GGFdownsig$gene.id)-GGFdownsigMP
# 177

# The number of expected mito-proteome transcripts
length(GGFdownsig$gene.id)*MPpropF
# 9.093417

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMP, GGFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 1.18e-05




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# How many mito-proteome transcripts are there?
G23FupsigMP <- sum(G23Fupsig$gene.id %in% MPfemaleswap$gene.id)
# 1
# How many are not?
G23FupsigNotMP <- length(G23Fupsig$gene.id)-G23FupsigMP
# 15

# The number of expected mito-proteome transcripts
length(G23Fupsig$gene.id)*MPpropF
# 0.731129

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FupsigMP, G23FupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.7475

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G23FupF <- rbind(c(G23FupsigMP, G23FupsigNotMP), 
                    c(length(MPfemaleswap$gene.id)-G23FupsigMP, (length(RNAswapfemale.data$gene.id)-length(MPfemaleswap$gene.id))-G23FupsigNotMP))
#Add names to the columns and rows
colnames(MP.G23FupF) <- c("MP", "NotMP")
rownames(MP.G23FupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G23FupF)
# P = 0.527

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# How many mito-proteome transcripts are there?
G23FdownsigMP <- sum(G23Fdownsig$gene.id %in% MPfemaleswap$gene.id)
# 0
# How many are not?
G23FdownsigNotMP <- length(G23Fdownsig$gene.id)-G23FdownsigMP
# 6

# The number of expected mito-proteome transcripts
length(G23Fdownsig$gene.id)*MPpropF
# 0.2741734

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(G23FdownsigMP, G23FdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.592

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.G23FdownF <- rbind(c(G23FdownsigMP, G23FdownsigNotMP),
                      c(length(MPfemaleswap$gene.id)-G23FdownsigMP, (length(RNAswapfemale.data$gene.id)-length(MPfemaleswap$gene.id))-G23FdownsigNotMP))
#Add names to the columns and rows
colnames(MP.G23FdownF) <- c("MP", "NotMP")
rownames(MP.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.G23FdownF)
# P = 1

# We get the same results, so that is good


