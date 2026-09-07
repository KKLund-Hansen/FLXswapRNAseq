################################################################################################
########### SEXUAL ANTAGONISTIC TRANSCRIPTS USING INNOCENTI AND MORROW FLX SWAP MALES ##########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexualAntagonisticTranscripts")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
IMsa.data <- read.table(file = "InnocentiAndMorrowAntagonisticTranscripts.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(IMsa.data) <- IMsa.data[,1]

# Next add a new column with the flybase names
IMsa.data$gene.id <- mapIds(drosophila2.db, keys = row.names(IMsa.data), keytype = "PROBEID", column = "FLYBASE")

# Next subset the Innocenti and Morrow transcript set to match the female data set
IMsaM.data <- subset(IMsa.data, IMsa.data$gene.id %in% FLXswapM.data$gene.id)
# 1337

# Next I'll calculate the proportion of the number of SA transcripts out of all transcripts
SApropIMM <- length(IMsaM.data$gene.id)/length(FLXswapM.data$gene.id)
# 0.09
# And for the chisq.test I also have to calculate
NotSApropIMM <- 1-SApropIMM
# 0.91



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# How many transcripts are sexual antagonistic?
FvCmupsigSAIM <- sum(FvCmupsig$gene.id %in% IMsaM.data$gene.id)
# 5
# How many are not?
FvCmupsigNotSAIM <- length(FvCmupsig$gene.id)-FvCmupsigSAIM
# 52

# The number of expected transcripts which are SA
length(FvCmupsig$gene.id)*SApropIMM
# 4.915124

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCmupsigSAIM, FvCmupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.9681

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FvCmupF <- rbind(c(FvCmupsigSAIM, FvCmupsigNotSAIM),
                      c(length(IMsaM.data$gene.id)-FvCmupsigSAIM, (length(FLXswapM.data$gene.id)-length(IMsaM.data$gene.id))-FvCmupsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FvCmupF) <- c("SA", "NotSA")
rownames(SAIM.FvCmupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FvCmupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# How many transcripts are sexual antagonistic?
FvCmdownsigSAIM <- sum(FvCmdownsig$gene.id %in% IMsaM.data$gene.id)
# 10
# How many are not?
FvCmdownsigNotSAIM <- length(FvCmdownsig$gene.id)-FvCmdownsigSAIM
# 47

# The number of expected transcripts which are SA
length(FvCmdownsig$gene.id)*SApropIMM
# 4.915124

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FvCmdownsigSAIM, FvCmdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.01642

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FvCmdownF <- rbind(c(FvCmdownsigSAIM, FvCmdownsigNotSAIM),
                        c(length(IMsaM.data$gene.id)-FvCmdownsigSAIM, (length(FLXswapM.data$gene.id)-length(IMsaM.data$gene.id))-FvCmdownsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FvCmdownF) <- c("SA", "NotSA")
rownames(SAIM.FvCmdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FvCmdownF)
# P = 0.02847

# It goes from significant to none significant




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# How many transcripts are sexual antagonistic?
AvSmupsigSAIM <- sum(AvSmupsig$gene.id %in% IMsaM.data$gene.id)
# 29
# How many are not?
AvSmupsigNotSAIM <- length(AvSmupsig$gene.id)-AvSmupsigSAIM
# 186

# The number of expected transcripts which are SA
length(AvSmupsig$gene.id)*SApropIMM
# 18.5395

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSmupsigSAIM, AvSmupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.01104



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# How many transcripts are sexual antagonistic?
AvSmdownsigSAIM <- sum(AvSmdownsig$gene.id %in% IMsaM.data$gene.id)
# 30
# How many are not?
AvSmdownsigNotSAIM <- length(AvSmdownsig$gene.id)-AvSmdownsigSAIM
# 232

# The number of expected transcripts which are SA
length(AvSmdownsig$gene.id)*SApropIMM
# 22.59233

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(AvSmdownsigSAIM, AvSmdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.103




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# How many transcripts are sexual antagonistic?
FsvCsMupsigSAIM <- sum(FsvCsMupsig$gene.id %in% IMsaM.data$gene.id)
# 13
# How many are not?
FsvCsMupsigNotSAIM <- length(FsvCsMupsig$gene.id)-FsvCsMupsigSAIM
# 43

# The number of expected transcripts which are SA
length(FsvCsMupsig$gene.id)*SApropIMM
# 4.828894

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsMupsigSAIM, FsvCsMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.0001003

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FsvCsMupF <- rbind(c(FsvCsMupsigSAIM, FsvCsMupsigNotSAIM),
                        c(length(IMsaM.data$gene.id)-FsvCsMupsigSAIM, (length(FLXswapM.data$gene.id)-length(IMsaM.data$gene.id))-FsvCsMupsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FsvCsMupF) <- c("SA", "NotSA")
rownames(SAIM.FsvCsMupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FsvCsMupF)
# P = 0.0007726

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# How many transcripts are sexual antagonistic?
FsvCsMdownsigSAIM <- sum(FsvCsMdownsig$gene.id %in% IMsaM.data$gene.id)
# 5
# How many are not?
FsvCsMdownsigNotSAIM <- length(FsvCsMdownsig$gene.id)-FsvCsMdownsigSAIM
# 44

# The number of expected transcripts which are SA
length(FsvCsMdownsig$gene.id)*SApropIMM
# 4.225282

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FsvCsMdownsigSAIM, FsvCsMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.6934

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
SAIM.FsvCsMdownF <- rbind(c(FsvCsMdownsigSAIM, FsvCsMdownsigNotSAIM),
                          c(length(IMsaM.data$gene.id)-FsvCsMdownsigSAIM, (length(FLXswapM.data$gene.id)-length(IMsaM.data$gene.id))-FsvCsMdownsigNotSAIM))
# Add names to the columns and rows
colnames(SAIM.FsvCsMdownF) <- c("SA", "NotSA")
rownames(SAIM.FsvCsMdownF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SAIM.FsvCsMdownF)
# P = 0.6105

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# How many transcripts are sexual antagonistic?
FavCaMupsigSAIM <- sum(FavCaMupsig$gene.id %in% IMsaM.data$gene.id)
# 9
# How many are not?
FavCaMupsigNotSAIM <- length(FavCaMupsig$gene.id)-FavCaMupsigSAIM
# 56

# The number of expected transcripts which are SA
length(FavCaMupsig$gene.id)*SApropIMM
# 5.604966

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaMupsigSAIM, FavCaMupsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.1336



##### DOWN-REGULATED TRANSCRIPTS #####

## Proportion of sexual antagonistic transcripts

# First subset the data to be able to identify sexual antagonistic transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# How many transcripts are sexual antagonistic?
FavCaMdownsigSAIM <- sum(FavCaMdownsig$gene.id %in% IMsaM.data$gene.id)
# 20
# How many are not?
FavCaMdownsigNotSAIM <- length(FavCaMdownsig$gene.id)-FavCaMdownsigSAIM
# 101

# The number of expected transcripts which are SA
length(FavCaMdownsig$gene.id)*SApropIMM
# 10.43386

# Do the test, with the number of SA transcripts and not SA transcripts and the proportion of overall SA transcripts 
chisq.test(c(FavCaMdownsigSAIM, FavCaMdownsigNotSAIM), p = c(SApropIMM, NotSApropIMM))
# P = 0.001948

