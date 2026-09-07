################################################################################################
########################### MITO-SENSITIVE TRANSCRIPTS FLX SWAP MALES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocMPnager", quietly = TRUE)) install.packages("BiocMPnager")
BiocMPnager::install("drosophila2.db")
BiocMPnager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the female data set
MSmaleswap <- subset(MS.data, MS.data$gene.id %in% FLXswapM.data$gene.id)
# 1210 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropM <- length(MSmaleswap$gene.id)/length(FLXswapM.data$gene.id)
# 0.08

# And for the chisq.test I also have to calculate
NotMSpropM <- 1-MSpropM
# 0.92



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many mito-sensitive transcripts are there?
FvCmupsigMS <- sum(FvCmupsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
FvCmupsigNotMS <- length(FvCmupsig$gene.id)-FvCmupsigMS
# 55

# The number of expected mito-sensitive transcripts
length(FvCmupsig$gene.id)*MSpropM
# 4.448243

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmupsigMS, FvCmupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.2267

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FvCmupF <- rbind(c(FvCmupsigMS, FvCmupsigNotMS), 
                    c(length(MSmaleswap$gene.id)-FvCmupsigMS, (length(FLXswapM.data$gene.id)-length(MSmaleswap$gene.id))-FvCmupsigNotMS))
#Add names to the columns and rows
colnames(MS.FvCmupF) <- c("MS", "NotMS")
rownames(MS.FvCmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FvCmupF)
# P = 0.3215

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many mito-sensitive transcripts are there?
FvCmdownsigMS <- sum(FvCmdownsig$gene.id %in% MSmaleswap$gene.id)
# 2
# How many are not?
FvCmdownsigNotMS <- length(FvCmdownsig$gene.id)-FvCmdownsigMS
# 55

# The number of expected mito-sensitive transcripts
length(FvCmdownsig$gene.id)*MSpropM
# 4.448243

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmdownsigMS, FvCmdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.2267

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FvCmdownF <- rbind(c(FvCmdownsigMS, FvCmdownsigNotMS), 
                      c(length(MSmaleswap$gene.id)-FvCmdownsigMS, (length(FLXswapM.data$gene.id)-length(MSmaleswap$gene.id))-FvCmdownsigNotMS))
#Add names to the columns and rows
colnames(MS.FvCmdownF) <- c("MS", "NotMS")
rownames(MS.FvCmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FvCmdownF)
# P = 0.3215

# We get the same results, so that is good



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many mito-sensitive transcripts are there?
AvSmupsigMS <- sum(AvSmupsig$gene.id %in% MSmaleswap$gene.id)
# 9
# How many are not?
AvSmupsigNotMS <- length(AvSmupsig$gene.id)-AvSmupsigMS
# 206

# The number of expected mito-sensitive transcripts
length(AvSmupsig$gene.id)*MSpropM
# 16.77846

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmupsigMS, AvSmupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.04796



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcriptsin the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many mito-sensitive transcripts are there?
AvSmdownsigMS <- sum(AvSmdownsig$gene.id %in% MSmaleswap$gene.id)
# 14
# How many are not?
AvSmdownsigNotMS <- length(AvSmdownsig$gene.id)-AvSmdownsigMS
# 248

# The number of expected mito-sensitive transcripts
length(AvSmdownsig$gene.id)*MSpropM
# 20.44631

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmdownsigMS, AvSmdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.1376




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many mito-sensitive transcripts are there?
FsvCsMupsigMS <- sum(FsvCsMupsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
FsvCsMupsigNotMS <- length(FsvCsMupsig$gene.id)-FsvCsMupsigMS
# 55

# The number of expected mito-sensitive transcripts
length(FsvCsMupsig$gene.id)*MSpropM
# 4.370203

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMupsigMS, FsvCsMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.09315

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FsvCsMupF <- rbind(c(FsvCsMupsigMS, FsvCsMupsigNotMS), 
                      c(length(MSmaleswap$gene.id)-FsvCsMupsigMS, (length(FLXswapM.data$gene.id)-length(MSmaleswap$gene.id))-FsvCsMupsigNotMS))
#Add names to the columns and rows
colnames(MS.FsvCsMupF) <- c("MS", "NotMS")
rownames(MS.FsvCsMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FsvCsMupF)
# P = 0.1281

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many mito-sensitive transcripts are there?
FsvCsMdownsigMS <- sum(FsvCsMdownsig$gene.id %in% MSmaleswap$gene.id)
# 1
# How many are not?
FsvCsMdownsigNotMS <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigMS
# 48

# The number of expected mito-sensitive transcripts
length(FsvCsMdownsig$gene.id)*MSpropM
# 3.823928

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMdownsigMS, FsvCsMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.1326

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FsvCsMdownF <- rbind(c(FsvCsMdownsigMS, FsvCsMdownsigNotMS), 
                        c(length(MSmaleswap$gene.id)-FsvCsMdownsigMS, (length(FLXswapM.data$gene.id)-length(MSmaleswap$gene.id))-FsvCsMdownsigNotMS))
#Add names to the columns and rows
colnames(MS.FsvCsMdownF) <- c("MS", "NotMS")
rownames(MS.FsvCsMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FsvCsMdownF)
# P = 0.1801

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many mito-sensitive transcripts are there?
FavCaMupsigMS <- sum(FavCaMupsig$gene.id %in% MSmaleswap$gene.id)
# 4
# How many are not?
FavCaMupsigNotMS <- length(FavCaMupsig$gene.id)-FavCaMupsigMS
# 61

# The number of expected mito-sensitive transcripts
length(FavCaMupsig$gene.id)*MSpropM
# 5.072557

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMupsigMS, FavCaMupsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.6199



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many mito-sensitive transcripts are there?
FavCaMdownsigMS <- sum(FavCaMdownsig$gene.id %in% MSmaleswap$gene.id)
# 11
# How many are not?
FavCaMdownsigNotMS <- length(FavCaMdownsig$gene.id)-FavCaMdownsigMS
# 110

# The number of expected mito-sensitive transcripts
length(FavCaMdownsig$gene.id)*MSpropM
# 9.44276

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMdownsigMS, FavCaMdownsigNotMS), p = c(MSpropM, NotMSpropM))
# P = 0.5977

