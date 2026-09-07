################################################################################################
########################### MITO-ANNOTATED TRANSCRIPTS FLX SWAP MALES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/MitoNuclearConflict")

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
MA.data <- read.table(file = "mitoannotated.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First I'll subset the mito-annotated data set to match the male data set
MAmaleswap <- subset(MA.data, MA.data$From %in% FLXswapM.data$gene.id)
# 463 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-annotated transcripts out of all my transcripts
MApropM <- length(MAmaleswap$From)/length(FLXswapM.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotMApropM <- 1-MApropM
# 0.97



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many mito-annotated transcripts are there?
FvCmupsigMA <- sum(FvCmupsig$gene.id %in% MAmaleswap$From)
# 0
# How many are not?
FvCmupsigNotMA <- length(FvCmupsig$gene.id)-FvCmupsigMA
# 57

# The number of expected mito-annotated transcripts
length(FvCmupsig$gene.id)*MApropM
# 1.702096

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmupsigMA, FvCmupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.1853

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FvCmupF <- rbind(c(FvCmupsigMA, FvCmupsigNotMA), 
                    c(length(MAmaleswap$From)-FvCmupsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FvCmupsigNotMA))
#Add names to the columns and rows
colnames(MA.FvCmupF) <- c("MA", "NotMA")
rownames(MA.FvCmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FvCmupF)
# P = 0.419

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in FLX vs CFM 
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many mito-annotated transcripts are there?
FvCmdownsigMA <- sum(FvCmdownsig$gene.id %in% MAmaleswap$From)
# 3
# How many are not?
FvCmdownsigNotMA <- length(FvCmdownsig$gene.id)-FvCmdownsigMA
# 54

# The number of expected mito-annotated transcripts
length(FvCmdownsig$gene.id)*MApropM
# 1.702096

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCmdownsigMA, FvCmdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.3125

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FvCmdownF <- rbind(c(FvCmdownsigMA, FvCmdownsigNotMA), 
                      c(length(MAmaleswap$From)-FvCmdownsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FvCmdownsigNotMA))
#Add names to the columns and rows
colnames(MA.FvCmdownF) <- c("MA", "NotMA")
rownames(MA.FvCmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FvCmdownF)
# P = 0.2419

# We get the same results, so that is good



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many mito-annotated transcripts are there?
AvSmupsigMA <- sum(AvSmupsig$gene.id %in% MAmaleswap$From)
# 9
# How many are not?
AvSmupsigNotMA <- length(AvSmupsig$gene.id)-AvSmupsigMA
# 206

# The number of expected mito-annotated transcripts
length(AvSmupsig$gene.id)*MApropM
# 6.420187

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmupsigMA, AvSmupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.3013



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many mito-annotated transcripts are there?
AvSmdownsigMA <- sum(AvSmdownsig$gene.id %in% MAmaleswap$From)
# 3
# How many are not?
AvSmdownsigNotMA <- length(AvSmdownsig$gene.id)-AvSmdownsigMA
# 259

# The number of expected mito-annotated transcripts
length(AvSmdownsig$gene.id)*MApropM
# 7.82367

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSmdownsigMA, AvSmdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.07997




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many mito-annotated transcripts are there?
FsvCsMupsigMA <- sum(FsvCsMupsig$gene.id %in% MAmaleswap$From)
# 1
# How many are not?
FsvCsMupsigNotMA <- length(FsvCsMupsig$gene.id)-FsvCsMupsigMA
# 55

# The number of expected mito-annotated transcripts
length(FsvCsMupsig$gene.id)*MApropM
# 1.672235

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMupsigMA, FsvCsMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.5976

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FsvCsMupF <- rbind(c(FsvCsMupsigMA, FsvCsMupsigNotMA), 
                      c(length(MAmaleswap$From)-FsvCsMupsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FsvCsMupsigNotMA))
#Add names to the columns and rows
colnames(MA.FsvCsMupF) <- c("MA", "NotMA")
rownames(MA.FsvCsMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FsvCsMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many mito-annotated transcripts are there?
FsvCsMdownsigMA <- sum(FsvCsMdownsig$gene.id %in% MAmaleswap$From)
# 1
# How many are not?
FsvCsMdownsigNotMA <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigMA
# 48

# The number of expected mito-annotated transcripts
length(FsvCsMdownsig$gene.id)*MApropM
# 1.463205

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsMdownsigMA, FsvCsMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.6974

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FsvCsMdownF <- rbind(c(FsvCsMdownsigMA, FsvCsMdownsigNotMA), 
                        c(length(MAmaleswap$From)-FsvCsMdownsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FsvCsMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.FsvCsMdownF) <- c("MA", "NotMA")
rownames(MA.FsvCsMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FsvCsMdownF)
# P = 1

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many mito-annotated transcripts are there?
FavCaMupsigMA <- sum(FavCaMupsig$gene.id %in% MAmaleswap$From)
# 1
# How many are not?
FavCaMupsigNotMA <- length(FavCaMupsig$gene.id)-FavCaMupsigMA
# 64

# The number of expected mito-annotated transcripts
length(FavCaMupsig$gene.id)*MApropM
# 1.940987

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMupsigMA, FavCaMupsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.5405

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FavCaMupF <- rbind(c(FavCaMupsigMA, FavCaMupsigNotMA), 
                      c(length(MAmaleswap$From)-FavCaMupsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FavCaMupsigNotMA))
#Add names to the columns and rows
colnames(MA.FavCaMupF) <- c("MA", "NotMA")
rownames(MA.FavCaMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FavCaMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many mito-annotated transcripts are there?
FavCaMdownsigMA <- sum(FavCaMdownsig$gene.id %in% MAmaleswap$From)
# 1
# How many are not?
FavCaMdownsigNotMA <- length(FavCaMdownsig$gene.id)-FavCaMdownsigMA
# 120

# The number of expected mito-annotated transcripts
length(FavCaMdownsig$gene.id)*MApropM
# 3.613222

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaMdownsigMA, FavCaMdownsigNotMA), p = c(MApropM, NotMApropM))
# P = 0.1628

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FavCaMdownF <- rbind(c(FavCaMdownsigMA, FavCaMdownsigNotMA), 
                        c(length(MAmaleswap$From)-FavCaMdownsigMA, (length(FLXswapM.data$gene.id)-length(MAmaleswap$From))-FavCaMdownsigNotMA))
#Add names to the columns and rows
colnames(MA.FavCaMdownF) <- c("MA", "NotMA")
rownames(MA.FavCaMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FavCaMdownF)
# P = 0.2736

# We get the same results, so that is good

