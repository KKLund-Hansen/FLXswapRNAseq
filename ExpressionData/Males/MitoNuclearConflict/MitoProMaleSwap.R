################################################################################################
########################### MITO-PROTEOME TRANSCRIPTS FLX SWAP MALES ###########################
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
MPmaleswap <- subset(MP.data, MP.data$gene.id %in% FLXswapM.data$gene.id)
# 762 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcripts out of all my transcripts
MPpropM <- length(MPmaleswap$gene.id)/length(FLXswapM.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropM <- 1-MPpropM
# 0.95



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many mito-proteome transcripts are there?
FvCmupsigMP <- sum(FvCmupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
FvCmupsigNotMP <- length(FvCmupsig$gene.id)-FvCmupsigMP
# 57

# The number of expected mito-proteome transcripts
length(FvCmupsig$gene.id)*MPpropM
# 2.80129

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmupsigMP, FvCmupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.08609

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FvCmupF <- rbind(c(FvCmupsigMP, FvCmupsigNotMP), 
                    c(length(MPmaleswap$gene.id)-FvCmupsigMP, (length(FLXswapM.data$gene.id)-length(MPmaleswap$gene.id))-FvCmupsigNotMP))
#Add names to the columns and rows
colnames(MP.FvCmupF) <- c("MP", "NotMP")
rownames(MP.FvCmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FvCmupF)
# P = 0.1163

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many mito-proteome transcripts are there?
FvCmdownsigMP <- sum(FvCmdownsig$gene.id %in% MPmaleswap$gene.id)
# 7
# How many are not?
FvCmdownsigNotMP <- length(FvCmdownsig$gene.id)-FvCmdownsigMP
# 50

# The number of expected mito-proteome transcripts
length(FvCmdownsig$gene.id)*MPpropM
# 2.80129

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmdownsigMP, FvCmdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.01009

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FvCmdownF <- rbind(c(FvCmdownsigMP, FvCmdownsigNotMP), 
                      c(length(MPmaleswap$gene.id)-FvCmdownsigMP, (length(FLXswapM.data$gene.id)-length(MPmaleswap$gene.id))-FvCmdownsigNotMP))
#Add names to the columns and rows
colnames(MP.FvCmdownF) <- c("MP", "NotMP")
rownames(MP.FvCmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FvCmdownF)
# P = 0.02104

# We get the same results, so that is good



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many mito-proteome transcripts are there?
AvSmupsigMP <- sum(AvSmupsig$gene.id %in% MPmaleswap$gene.id)
# 12
# How many are not?
AvSmupsigNotMP <- length(AvSmupsig$gene.id)-AvSmupsigMP
# 203

# The number of expected mito-proteome transcripts
length(AvSmupsig$gene.id)*MPpropM
# 10.56627

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmupsigMP, AvSmupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.651



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcriptsin the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many mito-proteome transcripts are there?
AvSmdownsigMP <- sum(AvSmdownsig$gene.id %in% MPmaleswap$gene.id)
# 11
# How many are not?
AvSmdownsigNotMP <- length(AvSmdownsig$gene.id)-AvSmdownsigMP
# 251

# The number of expected mito-proteome transcripts
length(AvSmdownsig$gene.id)*MPpropM
# 12.8761

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmdownsigMP, AvSmdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.5918




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many mito-proteome transcripts are there?
FsvCsMupsigMP <- sum(FsvCsMupsig$gene.id %in% MPmaleswap$gene.id)
# 4
# How many are not?
FsvCsMupsigNotMP <- length(FsvCsMupsig$gene.id)-FsvCsMupsigMP
# 52

# The number of expected mito-proteome transcripts
length(FsvCsMupsig$gene.id)*MPpropM
# 2.752144

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMupsigMP, FsvCsMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.4405

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FsvCsMupF <- rbind(c(FsvCsMupsigMP, FsvCsMupsigNotMP), 
                      c(length(MPmaleswap$gene.id)-FsvCsMupsigMP, (length(FLXswapM.data$gene.id)-length(MPmaleswap$gene.id))-FsvCsMupsigNotMP))
#Add names to the columns and rows
colnames(MP.FsvCsMupF) <- c("MP", "NotMP")
rownames(MP.FsvCsMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FsvCsMupF)
# P = 0.355

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many mito-proteome transcripts are there?
FsvCsMdownsigMP <- sum(FsvCsMdownsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
FsvCsMdownsigNotMP <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigMP
# 49

# The number of expected mito-proteome transcripts
length(FsvCsMdownsig$gene.id)*MPpropM
# 2.408126

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMdownsigMP, FsvCsMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.1115

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FsvCsMdownF <- rbind(c(FsvCsMdownsigMP, FsvCsMdownsigNotMP), 
                        c(length(MPmaleswap$gene.id)-FsvCsMdownsigMP, (length(FLXswapM.data$gene.id)-length(MPmaleswap$gene.id))-FsvCsMdownsigNotMP))
#Add names to the columns and rows
colnames(MP.FsvCsMdownF) <- c("MP", "NotMP")
rownames(MP.FsvCsMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FsvCsMdownF)
# P = 0.1758

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many mito-proteome transcripts are there?
FavCaMupsigMP <- sum(FavCaMupsig$gene.id %in% MPmaleswap$gene.id)
# 0
# How many are not?
FavCaMupsigNotMP <- length(FavCaMupsig$gene.id)-FavCaMupsigMP
# 65

# The number of expected mito-proteome transcripts
length(FavCaMupsig$gene.id)*MPpropM
# 3.194453

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMupsigMP, FavCaMupsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.06682

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FavCaMupF <- rbind(c(FavCaMupsigMP, FavCaMupsigNotMP), 
                      c(length(MPmaleswap$gene.id)-FavCaMupsigMP, (length(FLXswapM.data$gene.id)-length(MPmaleswap$gene.id))-FavCaMupsigNotMP))
#Add names to the columns and rows
colnames(MP.FavCaMupF) <- c("MP", "NotMP")
rownames(MP.FavCaMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FavCaMupF)
# P = 0.07724

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many mito-proteome transcripts are there?
FavCaMdownsigMP <- sum(FavCaMdownsig$gene.id %in% MPmaleswap$gene.id)
# 8
# How many are not?
FavCaMdownsigNotMP <- length(FavCaMdownsig$gene.id)-FavCaMdownsigMP
# 113

# The number of expected mito-proteome transcripts
length(FavCaMdownsig$gene.id)*MPpropM
# 5.946598

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMdownsigMP, FavCaMdownsigNotMP), p = c(MPpropM, NotMPpropM))
# P = 0.3878

