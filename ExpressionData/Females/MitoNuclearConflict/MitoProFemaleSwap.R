################################################################################################
########################## MITO-PROTEOME TRANSCRIPTS FLX SWAP FEMALES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
MPfemaleswap <- subset(MP.data, MP.data$gene.id %in% FLXswapF.data$gene.id)
# 762 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcriptsout of all my transcripts
MPpropF <- length(MPfemaleswap$gene.id)/length(FLXswapF.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropF <- 1-MPpropF
# 0.95


######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many mito-proteome transcripts are there?
FvCfupsigMP <- sum(FvCfupsig$gene.id %in% MPfemaleswap$gene.id)
# 1
# How many are not?
FvCfupsigNotMP <- length(FvCfupsig$gene.id)-FvCfupsigMP
# 34

# The number of expected mito-proteome transcripts
length(FvCfupsig$gene.id)*MPpropF
# 1.72009

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfupsigMP, FvCfupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.5734

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FvCfupF <- rbind(c(FvCfupsigMP, FvCfupsigNotMP), 
                    c(length(MPfemaleswap$gene.id)-FvCfupsigMP, (length(FLXswapF.data$gene.id)-length(MPfemaleswap$gene.id))-FvCfupsigNotMP))
#Add names to the columns and rows
colnames(MP.FvCfupF) <- c("MP", "NotMP")
rownames(MP.FvCfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FvCfupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts 
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many mito-proteome transcripts are there?
FvCfdownsigMP <- sum(FvCfdownsig$gene.id %in% MPfemaleswap$gene.id)
# 1
# How many are not?
FvCfdownsigNotMP <- length(FvCfdownsig$gene.id)-FvCfdownsigMP
# 39

# The number of expected mito-proteome transcripts
length(FvCfdownsig$gene.id)*MPpropF
# 1.965817

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfdownsigMP, FvCfdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.4799

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FvCfdownF <- rbind(c(FvCfdownsigMP, FvCfdownsigNotMP), 
                      c(length(MPfemaleswap$gene.id)-FvCfdownsigMP, (length(FLXswapF.data$gene.id)-length(MPfemaleswap$gene.id))-FvCfdownsigNotMP))
#Add names to the columns and rows
colnames(MP.FvCfdownF) <- c("MP", "NotMP")
rownames(MP.FvCfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FvCfdownF)
# P = 0.7221

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many mito-proteome transcripts are there?
AvSfupsigMP <- sum(AvSfupsig$gene.id %in% MPfemaleswap$gene.id)
# 8
# How many are not?
AvSfupsigNotMP <- length(AvSfupsig$gene.id)-AvSfupsigMP
# 236

# The number of expected mito-proteome transcripts
length(AvSfupsig$gene.id)*MPpropF
# 11.99149

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfupsigMP, AvSfupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.2372



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many mito-proteome transcripts are there?
AvSfdownsigMP <- sum(AvSfdownsig$gene.id %in% MPfemaleswap$gene.id)
# 7
# How many are not?
AvSfdownsigNotMP <- length(AvSfdownsig$gene.id)-AvSfdownsigMP
# 154

# The number of expected mito-proteome transcripts
length(AvSfdownsig$gene.id)*MPpropF
# 7.912415

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfdownsigMP, AvSfdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.7394




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many mito-proteome transcripts are there?
FsvCsFupsigMP <- sum(FsvCsFupsig$gene.id %in% MPfemaleswap$gene.id)
# 1
# How many are not?
FsvCsFupsigNotMP <- length(FsvCsFupsig$gene.id)-FsvCsFupsigMP
# 136

# The number of expected mito-proteome transcripts
length(FsvCsFupsig$gene.id)*MPpropF
# 6.732925

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFupsigMP, FsvCsFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.02346



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many mito-proteome transcripts are there?
FsvCsFdownsigMP <- sum(FsvCsFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 4
# How many are not?
FsvCsFdownsigNotMP <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigMP
# 51

# The number of expected mito-proteome transcripts
length(FsvCsFdownsig$gene.id)*MPpropF
# 2.702999

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFdownsigMP, FsvCsFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.4185

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.FsvCsFdownF <- rbind(c(FsvCsFdownsigMP, FsvCsFdownsigNotMP), 
                        c(length(MPfemaleswap$gene.id)-FsvCsFdownsigMP, (length(FLXswapF.data$gene.id)-length(MPfemaleswap$gene.id))-FsvCsFdownsigNotMP))
#Add names to the columns and rows
colnames(MP.FsvCsFdownF) <- c("MP", "NotMP")
rownames(MP.FsvCsFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.FsvCsFdownF)
# P = 0.347

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many mito-proteome transcripts are there?
FavCaFupsigMP <- sum(FavCaFupsig$gene.id %in% MPfemaleswap$gene.id)
# 69
# How many are not?
FavCaFupsigNotMP <- length(FavCaFupsig$gene.id)-FavCaFupsigMP
# 456

# The number of expected mito-proteome transcripts
length(FavCaFupsig$gene.id)*MPpropF
# 25.80135

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFupsigMP, FavCaFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many mito-proteome transcripts are there?
FavCaFdownsigMP <- sum(FavCaFdownsig$gene.id %in% MPfemaleswap$gene.id)
# 126
# How many are not?
FavCaFdownsigNotMP <- length(FavCaFdownsig$gene.id)-FavCaFdownsigMP
# 2868

# The number of expected mito-proteome transcripts
length(FavCaFdownsig$gene.id)*MPpropF
# 147.1414

# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFdownsigMP, FavCaFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.07388


