################################################################################################
#################### DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP MALES ####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/DifferentialExpression")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

# Set up environment
library(DESeq2)

# Read in files with data
FLXswapsampleinfo <- read.table(file = "FLXswapSampleInfo.csv", h = T, sep = ",", stringsAsFactors = T)
FLXswapreadcountUnfilt <- read.table(file = "FLXswapcount.txt", h = T, sep = "", skip = 1, row.names = "Geneid", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# Before we run the model, we'll tidy up the count data a bit. First, we remove the first column as they contain information that is not necessary to create a DESeqDataSet
FLXswapreadcount <- FLXswapreadcountUnfilt[,6:101]

# And then we remove the ".bam" and X from the column names so they are similar to the sample file row names
colnames(FLXswapreadcount) <- sub(".bam", "", colnames(FLXswapreadcount))
colnames(FLXswapreadcount) <- sub("X", "", colnames(FLXswapreadcount))

# We also need to make sure that the nested factor is a factor
FLXswapsampleinfo$nested.RP <- as.factor(FLXswapsampleinfo$nested.RP)



#################################################
########### FULL MODEL WITH BOTH SEXES ##########
#################################################

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


######################## MALES ########################

# As we want to analyse the data separately per sex, we start by sub-setting the data

# First, separate the data into males
FLXswapM <- FLXswap[,FLXswap$sex %in% "m"]

# Next drop sex as a level, as it is not need anymore
FLXswapM$sex <- droplevels(FLXswapM$sex)

# Last thing is to update the model to not contain sex
design(FLXswapM) <- formula(~ type + type:nested.RP)

# And the run the DESeq again
FLXswapM <- DESeq(FLXswapM)


########## MALE VS MALE ##########

# Create a function to be able to combine different types in one group. From https://www.atakanekiz.com/technical/a-guide-to-designs-and-contrasts-in-DESeq2/
contraster <- function(dds, group1, group2, weighted = F
){
  mod_mat <- model.matrix(design(dds), colData(dds))
  grp1_rows <- list()
  grp2_rows <- list()
  for(i in 1:length(group1)){
    grp1_rows[[i]] <- colData(dds)[[group1[[i]][1]]] %in% group1[[i]][2:length(group1[[i]])]
  }
  for(i in 1:length(group2)){
    grp2_rows[[i]] <- colData(dds)[[group2[[i]][1]]] %in% group2[[i]][2:length(group2[[i]])]
  }
  grp1_rows <- Reduce(function(x, y) x & y, grp1_rows)
  grp2_rows <- Reduce(function(x, y) x & y, grp2_rows)
  mod_mat1 <- mod_mat[grp1_rows, ,drop=F]
  mod_mat2 <- mod_mat[grp2_rows, ,drop=F]
  if(!weighted){
    mod_mat1 <- mod_mat1[!duplicated(mod_mat1),,drop=F]
    mod_mat2 <- mod_mat2[!duplicated(mod_mat2),,drop=F]
  }
  return(colMeans(mod_mat1)-colMeans(mod_mat2))
}


# FLX and CFM
FLXCFMm <- results(FLXswapM, contrast = contraster(FLXswapM, group1 = list(c("type", "aFLXsCwta","bCwtsFLXa")), 
                                                   group2 = list(c("type", "cCFMsCwta", "dCwtsCFMa"))))
# Then select the transcripts that are significantly different
FLXCFMmsig <- subset(FLXCFMm, padj < 0.05)
summary(FLXCFMmsig)
# There are 114 transcripts that are different, 57 up-regulated and 57 down-regulated
# Up-regulated transcripts
FLXCFMmsigup <- subset(FLXCFMmsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXCFMmsigdown <- subset(FLXCFMmsig, log2FoldChange < 0)


# Autosomes and sex chromosomes
AutoSexM <- results(FLXswapM, contrast = contraster(FLXswapM, group1 = list(c("type", "bCwtsFLXa", "dCwtsCFMa")),
                                                    group2 = list(c("type", "aFLXsCwta","cCFMsCwta"))))
# Then select the transcripts that are significantly different
AutoSexMsig <- subset(AutoSexM, padj < 0.05)
summary(AutoSexMsig)
# There are 477 transcripts that are different, 215 up-regulated and 262 down-regulated
# Up-regulated transcripts
AutoSexMsigup <- subset(AutoSexMsig, log2FoldChange > 0)
# Down-regulated transcripts
AutoSexMsigdown <- subset(AutoSexMsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexMsig <- subset(FLXsexCFMsexM, padj < 0.05)
summary(FLXsexCFMsexMsig)
# There are 105 transcripts that are different, 56 up-regulated and 49 down-regulated
# Up-regulated transcripts
FLXsexCFMsexMsigup <- subset(FLXsexCFMsexMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMsexMsigdown <- subset(FLXsexCFMsexMsig, log2FoldChange < 0)


# FLX autosomes and CFM autosomes
FLXaCFMaM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaMsig <- subset(FLXaCFMaM, padj < 0.05)
summary(FLXaCFMaMsig)
# There are 186 transcripts that are different, 65 up-regulated and 121 down-regulated
# Up-regulated transcripts
FLXaCFMaMsigup <- subset(FLXaCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaMsigdown <- subset(FLXaCFMaMsig, log2FoldChange < 0)



########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########


### NEW GROUPS ###

# First I make a data frame to include all transcripts that are specifically only found in females

# First collect the transcripts from the last DGElist
genesMswap <- rownames(FLXswapM)

# Now create the new data frame
FLXswapResM <- data.frame(gene.id = genesMswap, FLX.CFM = numeric(15505), FLX.CFMup = numeric(15505), FLX.CFMdown = numeric(15505), 
                          auto.sex = numeric(15505), auto.sexup = numeric(15505), auto.sexdown = numeric(15505),
                          FLXs.CFMs = numeric(15505), FLXs.CFMsup = numeric(15505), FLXs.CFMsdown = numeric(15505),
                          FLXa.CFMa = numeric(15505), FLXa.CFMaup = numeric(15505), FLXa.CFMadown = numeric(15505))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts from FLX vs CFM
FLXswapResM$FLX.CFM <- "NotSig"
# To be able to add sig make a list of transcripts
genesFLXCFMm <- rownames(FLXCFMmsig)
# Then add Sig for transcripts that are significant
FLXswapResM$FLX.CFM[which(FLXswapResM$gene.id %in% genesFLXCFMm)] <- "Sig"
# Make this a factor
FLXswapResM$FLX.CFM <- factor(FLXswapResM$FLX.CFM)
# Next add the up- and down-regulated transcripts separately
FLXswapResM$FLX.CFMup <- "NotSig"
FLXswapResM$FLX.CFMdown <- "NotSig"
genesFLXCFMmup <- rownames(FLXCFMmsigup)
genesFLXCFMmdown <- rownames(FLXCFMmsigdown)
FLXswapResM$FLX.CFMup[which(FLXswapResM$gene.id %in% genesFLXCFMmup)] <- "Sig"
FLXswapResM$FLX.CFMdown[which(FLXswapResM$gene.id %in% genesFLXCFMmdown)] <- "Sig"
FLXswapResM$FLX.CFMup <- factor(FLXswapResM$FLX.CFMup)
FLXswapResM$FLX.CFMdown <- factor(FLXswapResM$FLX.CFMdown)

# Then I add the transcripts from sex chromosomes vs autosomes
FLXswapResM$auto.sex <- "NotSig"
genesAutoSexM <- rownames(AutoSexMsig)
FLXswapResM$auto.sex[which(FLXswapResM$gene.id %in% genesAutoSexM)] <- "Sig"
FLXswapResM$auto.sex <- factor(FLXswapResM$auto.sex)
FLXswapResM$auto.sexup <- "NotSig"
FLXswapResM$auto.sexdown <- "NotSig"
genesAutoSexMup <- rownames(AutoSexMsigup)
genesAutoSexMdown <- rownames(AutoSexMsigdown)
FLXswapResM$auto.sexup[which(FLXswapResM$gene.id %in% genesAutoSexMup)] <- "Sig"
FLXswapResM$auto.sexdown[which(FLXswapResM$gene.id %in% genesAutoSexMdown)] <- "Sig"
FLXswapResM$auto.sexup <- factor(FLXswapResM$auto.sexup)
FLXswapResM$auto.sexdown <- factor(FLXswapResM$auto.sexdown)

# Then I add the transcripts from FLX sex chromosomes vs CFM sex chromosomes
FLXswapResM$FLXs.CFMs <- "NotSig"
genesFLXsexCFMsexM <- rownames(FLXsexCFMsexMsig)
FLXswapResM$FLXs.CFMs[which(FLXswapResM$gene.id %in% genesFLXsexCFMsexM)] <- "Sig"
FLXswapResM$FLXs.CFMs <- factor(FLXswapResM$FLXs.CFMs)
FLXswapResM$FLXs.CFMsup <- "NotSig"
FLXswapResM$FLXs.CFMsdown <- "NotSig"
genesFLXsexCFMsexMup <- rownames(FLXsexCFMsexMsigup)
genesFLXsexCFMsexMdown <- rownames(FLXsexCFMsexMsigdown)
FLXswapResM$FLXs.CFMsup[which(FLXswapResM$gene.id %in% genesFLXsexCFMsexMup)] <- "Sig"
FLXswapResM$FLXs.CFMsdown[which(FLXswapResM$gene.id %in% genesFLXsexCFMsexMdown)] <- "Sig"
FLXswapResM$FLXs.CFMsup <- factor(FLXswapResM$FLXs.CFMsup)
FLXswapResM$FLXs.CFMsdown <- factor(FLXswapResM$FLXs.CFMsdown)

# Then I add the transcripts from FLX autosomes vs CFM autosomes
FLXswapResM$FLXa.CFMa <- "NotSig"
genesFLXaCFMaM <- rownames(FLXaCFMaMsig)
FLXswapResM$FLXa.CFMa[which(FLXswapResM$gene.id %in% genesFLXaCFMaM)] <- "Sig"
FLXswapResM$FLXa.CFMa <- factor(FLXswapResM$FLXa.CFMa)
FLXswapResM$FLXa.CFMaup <- "NotSig"
FLXswapResM$FLXa.CFMadown <- "NotSig"
genesFLXaCFMaMup <- rownames(FLXaCFMaMsigup)
genesFLXaCFMaMdown <- rownames(FLXaCFMaMsigdown)
FLXswapResM$FLXa.CFMaup[which(FLXswapResM$gene.id %in% genesFLXaCFMaMup)] <- "Sig"
FLXswapResM$FLXa.CFMadown[which(FLXswapResM$gene.id %in% genesFLXaCFMaMdown)] <- "Sig"
FLXswapResM$FLXa.CFMaup <- factor(FLXswapResM$FLXa.CFMaup)
FLXswapResM$FLXa.CFMadown <- factor(FLXswapResM$FLXa.CFMadown)

# Add chromosome location to the data set as well
# First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
# Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
# Now make a new vector with the chromosomes matching the FLXswapResM data frame
chrdataMswap <- cbind(genesMswap, chrIDswap[match(FLXswapResM$gene.id, rownames(FLXswapreadcount)), c(1,2)])

# Last add the chromosome location to the FLXswapResM data frame
FLXswapResM$chr <- character(15505)
for(i in 1:15505) { FLXswapResM$chr[FLXswapResM$gene.id == chrdataMswap[i, 1]] <- chrdataMswap[i, 2]}

# Clean up the chromosome location column
FLXswapResM$chr <- gsub("\\;.*", "", as.character(FLXswapResM$chr))


# Now save the full new data frame as a csv file
write.csv(FLXswapResM, "FLXswapMres.csv", row.names = F)

