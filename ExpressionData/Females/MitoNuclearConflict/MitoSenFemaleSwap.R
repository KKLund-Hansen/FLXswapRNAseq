################################################################################################
########################## MITO-SENSITIVE TRANSCRIPTS FLX SWAP FEMALES ######################### 
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
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the female data set
MSfemaleswap <- subset(MS.data, MS.data$gene.id %in% FLXswapF.data$gene.id)
# 1210 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropF <- length(MSfemaleswap$gene.id)/length(FLXswapF.data$gene.id)
# 0.08

# And for the chisq.test I also have to calculate
NotMSpropF <- 1-MSpropF
# 0.92



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many mito-sensitive transcripts are there?
FvCfupsigMS <- sum(FvCfupsig$gene.id %in% MSfemaleswap$gene.id)
# 2
# How many are not?
FvCfupsigNotMS <- length(FvCfupsig$gene.id)-FvCfupsigMS
# 33

# The number of expected mito-sensitive transcripts
length(FvCfupsig$gene.id)*MSpropF
# 2.731377

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfupsigMS, FvCfupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.6449

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FvCfupF <- rbind(c(FvCfupsigMS, FvCfupsigNotMS), 
                    c(length(MSfemaleswap$gene.id)-FvCfupsigMS, (length(FLXswapF.data$gene.id)-length(MSfemaleswap$gene.id))-FvCfupsigNotMS))
#Add names to the columns and rows
colnames(MS.FvCfupF) <- c("MS", "NotMS")
rownames(MS.FvCfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FvCfupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the FLXvCFM group 
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many mito-sensitive transcripts are there?
FvCfdownsigMS <- sum(FvCfdownsig$gene.id %in% MSfemaleswap$gene.id)
# 1
# How many are not?
FvCfdownsigNotMS <- length(FvCfdownsig$gene.id)-FvCfdownsigMS
# 39

# The number of expected mito-sensitive transcripts
length(FvCfdownsig$gene.id)*MSpropF
# 3.121574

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfdownsigMS, FvCfdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2111

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FvCfdownF <- rbind(c(FvCfdownsigMS, FvCfdownsigNotMS), 
                      c(length(MSfemaleswap$gene.id)-FvCfdownsigMS, (length(FLXswapF.data$gene.id)-length(MSfemaleswap$gene.id))-FvCfdownsigNotMS))
#Add names to the columns and rows
colnames(MS.FvCfdownF) <- c("MS", "NotMS")
rownames(MS.FvCfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FvCfdownF)
# P = 0.3686

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many mito-sensitive transcripts are there?
AvSfupsigMS <- sum(AvSfupsig$gene.id %in% MSfemaleswap$gene.id)
# 16
# How many are not?
AvSfupsigNotMS <- length(AvSfupsig$gene.id)-AvSfupsigMS
# 228

# The number of expected mito-sensitive transcripts
length(AvSfupsig$gene.id)*MSpropF
# 19.0416

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfupsigMS, AvSfupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.4679



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many mito-sensitive transcripts are there?
AvSfdownsigMS <- sum(AvSfdownsig$gene.id %in% MSfemaleswap$gene.id)
# 3
# How many are not?
AvSfdownsigNotMS <- length(AvSfdownsig$gene.id)-AvSfdownsigMS
# 158

# The number of expected mito-sensitive transcripts
length(AvSfdownsig$gene.id)*MSpropF
# 12.56433

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfdownsigMS, AvSfdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.004952




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many mito-sensitive transcripts are there?
FsvCsFupsigMS <- sum(FsvCsFupsig$gene.id %in% MSfemaleswap$gene.id)
# 5
# How many are not?
FsvCsFupsigNotMS <- length(FsvCsFupsig$gene.id)-FsvCsFupsigMS
# 132

# The number of expected mito-sensitive transcripts
length(FsvCsFupsig$gene.id)*MSpropF
# 10.69139

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFupsigMS, FsvCsFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.06987



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many mito-sensitive transcripts are there?
FsvCsFdownsigMS <- sum(FsvCsFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 2
# How many are not?
FsvCsFdownsigNotMS <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigMS
# 53

# The number of expected mito-sensitive transcripts
length(FsvCsFdownsig$gene.id)*MSpropF
# 4.292164

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFdownsigMS, FsvCsFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2492

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.FsvCsFdownF <- rbind(c(FsvCsFdownsigMS, FsvCsFdownsigNotMS), 
                        c(length(MSfemaleswap$gene.id)-FsvCsFdownsigMS, (length(FLXswapF.data$gene.id)-length(MSfemaleswap$gene.id))-FsvCsFdownsigNotMS))
#Add names to the columns and rows
colnames(MS.FsvCsFdownF) <- c("MS", "NotMS")
rownames(MS.FsvCsFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.FsvCsFdownF)
# P = 0.3207

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many mito-sensitive transcripts are there?
FavCaFupsigMS <- sum(FavCaFupsig$gene.id %in% MSfemaleswap$gene.id)
# 41
# How many are not?
FavCaFupsigNotMS <- length(FavCaFupsig$gene.id)-FavCaFupsigMS
# 484

# The number of expected mito-sensitive transcripts
length(FavCaFupsig$gene.id)*MSpropF
# 40.97065

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFupsigMS, FavCaFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.9962



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many mito-sensitive transcripts are there?
FavCaFdownsigMS <- sum(FavCaFdownsig$gene.id %in% MSfemaleswap$gene.id)
# 268
# How many are not?
FavCaFdownsigNotMS <- length(FavCaFdownsig$gene.id)-FavCaFdownsigMS
# 2726

# The number of expected mito-sensitive transcripts
length(FavCaFdownsig$gene.id)*MSpropF
# 233.6498

# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFdownsigMS, FavCaFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.01926


