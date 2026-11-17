################################################################################################
############################ SHARED INTERSECTIONS BETWEEN THE SEXES ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Both/Dataprep")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

# Set up environment
library(DESeq2)

# Read in files with data
FLXswapsampleinfo <- read.table(file = "FLXswapSampleInfo.csv", h = T, sep = ",", stringsAsFactors = T)
FLXswapreadcountUnfilt <- read.table(file = "FLXswapcount.txt", h = T, sep = "", skip = 1, row.names = "Geneid", stringsAsFactors = T)


########################################## DIFFIRENTIAL EXPRESSION MODEL #########################################

##### PREPARE THE DATA #####

# Before we run the model, we'll tidy up the count data a bit. First, we remove the first column as they contain information that is not necessary to create a DESeqDataSet
FLXswapreadcount <- FLXswapreadcountUnfilt[,6:101]

# And then we remove the ".bam" and X from the column names so they are similar to the sample file row names
colnames(FLXswapreadcount) <- sub(".bam", "", colnames(FLXswapreadcount))
colnames(FLXswapreadcount) <- sub("X", "", colnames(FLXswapreadcount))

# We also need to make sure that the nested factor is a factor
FLXswapsampleinfo$nested.RP <- as.factor(FLXswapsampleinfo$nested.RP)


######################## FULL MODEL WITH BOTH SEXES ########################

# Create the DESeqDataSet object
FLXswapsex <- DESeqDataSetFromMatrix(countData = FLXswapreadcount, colData = FLXswapsampleinfo, design = ~  type + sex + type:sex + type:nested.RP)
# 17609 transcripts

# Next we do a prefilter and remove all transcripts with a count less than 10 across the three biological replicates
smallestGroupSize <- 3
keep <- rowSums(counts(FLXswapsex) >= 10) >= smallestGroupSize
FLXswapsex <- FLXswapsex[keep,]

# Check how many transcripts are left
nrow(FLXswapsex)
# 15505

# Now are we ready to do a differential expression analysis
FLXswap <- DESeq(FLXswapsex)


######################## SEPERATE SEX MODEL ########################

# As we want to analyse the data separately per sex, we start by sub-setting the data

# First, separate the data into females
FLXswapF <- FLXswap[,FLXswap$sex %in% "f"]
# And males
FLXswapM <- FLXswap[,FLXswap$sex %in% "m"]

# Next drop sex as a level, as it is not need anymore
FLXswapF$sex <- droplevels(FLXswapF$sex)
FLXswapM$sex <- droplevels(FLXswapM$sex)

# Last thing is to update the model to not contain sex
design(FLXswapF) <- formula(~ type + type:nested.RP)
design(FLXswapM) <- formula(~ type + type:nested.RP)

# And the run the DESeq again
FLXswapF <- DESeq(FLXswapF)
FLXswapM <- DESeq(FLXswapM)


########################################## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########################################


######################## FEMALE DATA ########################

# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaFsig <- subset(FLXsexFLXaF, padj < 0.05)

# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexFsig <- subset(FLXsexCFMsexF, padj < 0.05)

# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaFsig <- subset(FLXsexCFMaF, padj < 0.05)

# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXaCFMsexFsig <- subset(FLXaCFMsexF, padj < 0.05)

# FLX-autosomes and CFM-autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaFsig <- subset(FLXaCFMaF, padj < 0.05)

# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaF <- results(FLXswapF, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
CFMsexCFMaFsig <- subset(CFMsexCFMaF, padj < 0.05)


########## COLLECT THE DATA ##########

# Make a new vector with the shared significant transcripts in females
# First collect the transcripts in the female only groups

# Group 35
G35.1 <- intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig)), rownames(CFMsexCFMaFsig))
Group35 <- setdiff(setdiff(setdiff(G35.1, rownames(FLXsexFLXaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))

# Group 39
G39.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))
Group39 <- setdiff(setdiff(setdiff(G39.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(CFMsexCFMaFsig))

# Group 40
G40.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(CFMsexCFMaFsig))
Group40 <- setdiff(setdiff(setdiff(G40.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXaCFMaFsig))

# Group 42
G42.1 <- intersect(intersect(rownames(FLXaCFMsexFsig), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
Group42 <- setdiff(setdiff(setdiff(G42.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig))

# Group 44
G44.1 <- intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMaFsig))
Group44 <- setdiff(setdiff(G44.1, rownames(FLXaCFMsexFsig)), rownames(CFMsexCFMaFsig))

# Group 49
G49.1 <- intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))
Group49 <- setdiff(setdiff(G49.1, rownames(FLXsexCFMsexFsig)), rownames(CFMsexCFMaFsig))

# Group 56
G56.1 <- intersect(intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
Group56 <- setdiff(setdiff(G56.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMaFsig))

# Group 57
G57.1 <- intersect(intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
Group57 <- setdiff(setdiff(G57.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig))

# Group 59
G59.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), 
                   rownames(CFMsexCFMaFsig))
Group59 <- setdiff(G59.1, rownames(FLXaCFMaFsig))

# Group 62
G62.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), 
                   rownames(CFMsexCFMaFsig))
Group62 <- setdiff(G62.1, rownames(FLXsexCFMsexFsig))

# Group 63
G63.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), 
                   rownames(CFMsexCFMaFsig))
Group63 <- setdiff(G63.1, rownames(FLXsexFLXaFsig))

# Now collect the female only groups
femaleonly <- c(Group35, Group39, Group40, Group42, Group44, Group49, Group56, Group57, Group59, Group62, Group63)

# Collect all significant female transcripts
AllsigFemale <- unique(cbind(c(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig),
                               rownames(FLXaCFMsexFsig), rownames(FLXaCFMaFsig), rownames(CFMsexCFMaFsig))))

# For use in other analyses I save all significant female transcripts as a csv file
write.csv(AllsigFemale, "AllsigFemaleres.csv", row.names = F)



######################## MALE DATA ########################

# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaMsig <- subset(FLXsexFLXaM, padj < 0.05)

# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexMsig <- subset(FLXsexCFMsexM, padj < 0.05)

# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaMsig <- subset(FLXsexCFMaM, padj < 0.05)

# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXaCFMsexMsig <- subset(FLXaCFMsexM, padj < 0.05)

# FLX-autosomes and CFM-autosomes
FLXaCFMaM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaMsig <- subset(FLXaCFMaM, padj < 0.05)

# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaM <- results(FLXswapM, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
CFMsexCFMaMsig <- subset(CFMsexCFMaM, padj < 0.05)


########## COLLECT THE DATA ##########

# Make a new vector with the shared significant transcripts in males
# First collect the transcripts in the male only groups

# Group 33
G33 <- intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig)), rownames(FLXaCFMsexMsig))
Group33 <- setdiff(setdiff(setdiff(G33, rownames(FLXsexFLXaMsig)), rownames(FLXaCFMaMsig)), rownames(CFMsexCFMaMsig))

# Group 36
G36 <- intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXaCFMsexMsig)), rownames(FLXaCFMaMsig))
Group36 <- setdiff(setdiff(setdiff(G36, rownames(FLXsexFLXaMsig)), rownames(FLXsexCFMaMsig)), rownames(CFMsexCFMaMsig))

# Group 55
G55 <- intersect(intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig)), rownames(FLXaCFMaMsig)),
                 rownames(CFMsexCFMaMsig))
Group55 <- setdiff(setdiff(G55, rownames(FLXsexFLXaMsig)), rownames(FLXaCFMsexMsig))

# Now collect the male only groups
maleonly <- c(Group33, Group36, Group55)

# Collect all significant male transcripts
AllsigMale <- unique(cbind(c(rownames(FLXsexFLXaMsig), rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig),
                             rownames(FLXaCFMsexMsig), rownames(FLXaCFMaMsig), rownames(CFMsexCFMaMsig))))

# For use in other analyses I save all significant female transcripts as a csv file
write.csv(AllsigMale, "AllsigMaleres.csv", row.names = F)



######################## OVERLAP ########################

# Now it is time to test for how many significant transcripts overlap between the two sexes

# Collect all transcripts in the shared intersections in females
SharedFemale <- setdiff(AllsigFemale, femaleonly)

# Collect the transcripts in the shared intersections in males
SharedMale <- setdiff(AllsigMale, maleonly)

# How many overlap?
Shared <- intersect(SharedMale, SharedFemale)
# 264


# Percentage of overlap in males
(length(Shared)/length(SharedMale))*100
# 44.59% of male transcripts

# Percentage of overlap in females
(length(Shared)/length(SharedFemale))*100
# 4.80262

