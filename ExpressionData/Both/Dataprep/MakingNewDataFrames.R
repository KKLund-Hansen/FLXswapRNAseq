################################################################################################
################################# DIFFERENCES BETWEEN THE SEXES ################################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Both/Dataprep")

# First install the software from Bioconductor
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

# First collect the gene IDs from the full model
genesswap <- rownames(FLXswapF)

# Then to create a multiple new data sets
FLXswapSE <- data.frame(gene.id = genesswap)
FLXswapL2F <- data.frame(gene.id = genesswap)
FLXswapL2Fabs <- data.frame(gene.id = genesswap)


########## FEMALE DATA ##########

# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Add the data to the new data frames
FLXswapSE$FLXsFLXaf <- FLXsexFLXaF$lfcSE
FLXswapL2F$FLXsFLXaf <- FLXsexFLXaF$log2FoldChange
FLXswapL2Fabs$FLXsFLXaf <- abs(FLXsexFLXaF$log2FoldChange)

# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Add the data to the new data frames
FLXswapSE$FLXsCFMsf <- FLXsexCFMsexF$lfcSE
FLXswapL2F$FLXsCFMsf <- FLXsexCFMsexF$log2FoldChange
FLXswapL2Fabs$FLXsCFMsf <- abs(FLXsexCFMsexF$log2FoldChange)

# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$FLXsCFMaf <- FLXsexCFMaF$lfcSE
FLXswapL2F$FLXsCFMaf <- FLXsexCFMaF$log2FoldChange
FLXswapL2Fabs$FLXsCFMaf <- abs(FLXsexCFMaF$log2FoldChange)

# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Add the data to the new data frames
FLXswapSE$FLXaCFMsf <- FLXaCFMsexF$lfcSE
FLXswapL2F$FLXaCFMsf <- FLXaCFMsexF$log2FoldChange
FLXswapL2Fabs$FLXaCFMsf <- abs(FLXaCFMsexF$log2FoldChange)

# FLX-autosomes and CFM-autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$FLXaCFMaf <- FLXaCFMaF$lfcSE
FLXswapL2F$FLXaCFMaf <- FLXaCFMaF$log2FoldChange
FLXswapL2Fabs$FLXaCFMaf <- abs(FLXaCFMaF$log2FoldChange)

# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaF <- results(FLXswapF, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$CFMsCFMaf <- CFMsexCFMaF$lfcSE
FLXswapL2F$CFMsCFMaf <- CFMsexCFMaF$log2FoldChange
FLXswapL2Fabs$CFMsCFMaf <- abs(CFMsexCFMaF$log2FoldChange)


########## MALE DATA ##########

# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Add the data to the new data frames
FLXswapSE$FLXsFLXam <- FLXsexFLXaM$lfcSE
FLXswapL2F$FLXsFLXam <- FLXsexFLXaM$log2FoldChange
FLXswapL2Fabs$FLXsFLXam <- abs(FLXsexFLXaM$log2FoldChange)

# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Add the data to the new data frames
FLXswapSE$FLXsCFMsm <- FLXsexCFMsexM$lfcSE
FLXswapL2F$FLXsCFMsm <- FLXsexCFMsexM$log2FoldChange
FLXswapL2Fabs$FLXsCFMsm <- abs(FLXsexCFMsexM$log2FoldChange)

# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$FLXsCFMam <- FLXsexCFMaM$lfcSE
FLXswapL2F$FLXsCFMam <- FLXsexCFMaM$log2FoldChange
FLXswapL2Fabs$FLXsCFMam <- abs(FLXsexCFMaM$log2FoldChange)

# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Add the data to the new data frames
FLXswapSE$FLXaCFMsm <- FLXaCFMsexM$lfcSE
FLXswapL2F$FLXaCFMsm <- FLXaCFMsexM$log2FoldChange
FLXswapL2Fabs$FLXaCFMsm <- abs(FLXaCFMsexM$log2FoldChange)

# FLX-autosomes and CFM-autosomes
FLXaCFMaM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$FLXaCFMam <- FLXaCFMaM$lfcSE
FLXswapL2F$FLXaCFMam <- FLXaCFMaM$log2FoldChange
FLXswapL2Fabs$FLXaCFMam <- abs(FLXaCFMaM$log2FoldChange)

# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaM <- results(FLXswapM, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Add the data to the new data frames
FLXswapSE$CFMsCFMam <- CFMsexCFMaM$lfcSE
FLXswapL2F$CFMsCFMam <- CFMsexCFMaM$log2FoldChange
FLXswapL2Fabs$CFMsCFMam <- abs(CFMsexCFMaM$log2FoldChange)

# Now save the full new data frames as a csv file
write.csv(FLXswapSE, "FLXswapSEres.csv", row.names = F)
write.csv(FLXswapL2F, "FLXswapL2Fres.csv", row.names = F)
write.csv(FLXswapL2Fabs, "FLXswapL2Fabsres.csv", row.names = F)

