################################################################################################
################### DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP FEMALES ###################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/DifferentialExpression")

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


######################## FEMALES ########################

# As we want to analyse the data separately per sex, we start by sub-setting the data

# First, separate the data into females
FLXswapF <- FLXswap[,FLXswap$sex %in% "f"]

# Next drop sex as a level, as it is not need anymore
FLXswapF$sex <- droplevels(FLXswapF$sex)

# Last thing is to update the model to not contain sex
design(FLXswapF) <- formula(~ type + type:nested.RP)

# And the run the DESeq again
FLXswapF <- DESeq(FLXswapF)


########## FEMALE VS FEMALE ##########

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
FLXCFMf <- results(FLXswapF, contrast = contraster(FLXswapF, group1 = list(c("type", "aFLXsCwta","bCwtsFLXa")), 
                                                   group2 = list(c("type", "cCFMsCwta", "dCwtsCFMa"))))
# Then select the transcripts that are significantly different
FLXCFMfsig <- subset(FLXCFMf, padj < 0.05)
summary(FLXCFMfsig)
# There are 75 transcripts that are different, 35 up-regulated and 40 down-regulated
# Up-regulated transcripts
FLXCFMfsigup <- subset(FLXCFMfsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXCFMfsigdown <- subset(FLXCFMfsig, log2FoldChange < 0)


# Sex chromosomes and autosomes
AutoSexF <- results(FLXswapF, contrast = contraster(FLXswapF, group1 = list(c("type", "bCwtsFLXa", "dCwtsCFMa")),
                                                    group2 = list(c("type", "aFLXsCwta","cCFMsCwta"), weighted = T)))
# Then select the transcripts that are significantly different
AutoSexFsig <- subset(AutoSexF, padj < 0.05)
summary(AutoSexFsig)
# There are 405 transcripts that are different, 161 up-regulated and 244 down-regulated
# Up-regulated transcripts
AutoSexFsigup <- subset(AutoSexFsig, log2FoldChange > 0)
# Down-regulated transcripts
AutoSexFsigdown <- subset(AutoSexFsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexFsig <- subset(FLXsexCFMsexF, padj < 0.05)
summary(FLXsexCFMsexFsig)
# There are 192 transcripts that are different, 137 up-regulated and 55 down-regulated
# Up-regulated transcripts
FLXsexCFMsexFsigup <- subset(FLXsexCFMsexFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMsexFsigdown <- subset(FLXsexCFMsexFsig, log2FoldChange < 0)


# FLX autosomes and CFM autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaFsig <- subset(FLXaCFMaF, padj < 0.05)
summary(FLXaCFMaFsig)
# There are 3519 transcripts that are different, 525 up-regulated and 2994 down-regulated
# Up-regulated transcripts
FLXaCFMaFsigup <- subset(FLXaCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaFsigdown <- subset(FLXaCFMaFsig, log2FoldChange < 0)



########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########

# First collect the transcripts from the last DGElist
genesFswap <- rownames(FLXswapF)

# Now create the new data frame
FLXswapResF <- data.frame(gene.id = genesFswap, FLX.CFM = numeric(15505), FLX.CFMup = numeric(15505), FLX.CFMdown = numeric(15505), 
                          auto.sex = numeric(15505), auto.sexup = numeric(15505), auto.sexdown = numeric(15505),
                          FLXs.CFMs = numeric(15505), FLXs.CFMsup = numeric(15505), FLXs.CFMsdown = numeric(15505),
                          FLXa.CFMa = numeric(15505), FLXa.CFMaup = numeric(15505), FLXa.CFMadown = numeric(15505))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts from FLX vs CFM
FLXswapResF$FLX.CFM <- "NotSig"
# To be able to add sig make a list of transcripts
genesFLXCFMf <- rownames(FLXCFMfsig)
# Then add Sig for transcripts that are significant
FLXswapResF$FLX.CFM[which(FLXswapResF$gene.id %in% genesFLXCFMf)] <- "Sig"
# Make this a factor
FLXswapResF$FLX.CFM <- factor(FLXswapResF$FLX.CFM)
# Next add the up- and down-regulated transcripts separately
FLXswapResF$FLX.CFMup <- "NotSig"
FLXswapResF$FLX.CFMdown <- "NotSig"
genesFLXCFMfup <- rownames(FLXCFMfsigup)
genesFLXCFMfdown <- rownames(FLXCFMfsigdown)
FLXswapResF$FLX.CFMup[which(FLXswapResF$gene.id %in% genesFLXCFMfup)] <- "Sig"
FLXswapResF$FLX.CFMdown[which(FLXswapResF$gene.id %in% genesFLXCFMfdown)] <- "Sig"
FLXswapResF$FLX.CFMup <- factor(FLXswapResF$FLX.CFMup)
FLXswapResF$FLX.CFMdown <- factor(FLXswapResF$FLX.CFMdown)

# Then I add the transcripts from autosomes vs sex chromosomes 
FLXswapResF$auto.sex <- "NotSig"
genesAutoSexF <- rownames(AutoSexFsig)
FLXswapResF$auto.sex[which(FLXswapResF$gene.id %in% genesAutoSexF)] <- "Sig"
FLXswapResF$auto.sex <- factor(FLXswapResF$auto.sex)
FLXswapResF$auto.sexup <- "NotSig"
FLXswapResF$auto.sexdown <- "NotSig"
genesAutoSexFup <- rownames(AutoSexFsigup)
genesAutoSexFdown <- rownames(AutoSexFsigdown)
FLXswapResF$auto.sexup[which(FLXswapResF$gene.id %in% genesAutoSexFup)] <- "Sig"
FLXswapResF$auto.sexdown[which(FLXswapResF$gene.id %in% genesAutoSexFdown)] <- "Sig"
FLXswapResF$auto.sexup <- factor(FLXswapResF$auto.sexup)
FLXswapResF$auto.sexdown <- factor(FLXswapResF$auto.sexdown)

# Then I add the transcripts from FLX sex chromosomes vs CFM sex chromosomes
FLXswapResF$FLXs.CFMs <- "NotSig"
genesFLXsexCFMsexF <- rownames(FLXsexCFMsexFsig)
FLXswapResF$FLXs.CFMs[which(FLXswapResF$gene.id %in% genesFLXsexCFMsexF)] <- "Sig"
FLXswapResF$FLXs.CFMs <- factor(FLXswapResF$FLXs.CFMs)
FLXswapResF$FLXs.CFMsup <- "NotSig"
FLXswapResF$FLXs.CFMsdown <- "NotSig"
genesFLXsexCFMsexFup <- rownames(FLXsexCFMsexFsigup)
genesFLXsexCFMsexFdown <- rownames(FLXsexCFMsexFsigdown)
FLXswapResF$FLXs.CFMsup[which(FLXswapResF$gene.id %in% genesFLXsexCFMsexFup)] <- "Sig"
FLXswapResF$FLXs.CFMsdown[which(FLXswapResF$gene.id %in% genesFLXsexCFMsexFdown)] <- "Sig"
FLXswapResF$FLXs.CFMsup <- factor(FLXswapResF$FLXs.CFMsup)
FLXswapResF$FLXs.CFMsdown <- factor(FLXswapResF$FLXs.CFMsdown)

# Then I add the transcripts from FLX autosomes vs CFM autosomes
FLXswapResF$FLXa.CFMa <- "NotSig"
genesFLXaCFMaF <- rownames(FLXaCFMaFsig)
FLXswapResF$FLXa.CFMa[which(FLXswapResF$gene.id %in% genesFLXaCFMaF)] <- "Sig"
FLXswapResF$FLXa.CFMa <- factor(FLXswapResF$FLXa.CFMa)
FLXswapResF$FLXa.CFMaup <- "NotSig"
FLXswapResF$FLXa.CFMadown <- "NotSig"
genesFLXaCFMaFup <- rownames(FLXaCFMaFsigup)
genesFLXaCFMaFdown <- rownames(FLXaCFMaFsigdown)
FLXswapResF$FLXa.CFMaup[which(FLXswapResF$gene.id %in% genesFLXaCFMaFup)] <- "Sig"
FLXswapResF$FLXa.CFMadown[which(FLXswapResF$gene.id %in% genesFLXaCFMaFdown)] <- "Sig"
FLXswapResF$FLXa.CFMaup <- factor(FLXswapResF$FLXa.CFMaup)
FLXswapResF$FLXa.CFMadown <- factor(FLXswapResF$FLXa.CFMadown)

# Add chromosome location to the data set as well
# First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
# Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
# Now make a new vector with the chromosomes matching the FLXswapResF data frame
chrdataFswap <- cbind(genesFswap, chrIDswap[match(FLXswapResF$gene.id, rownames(FLXswapreadcount)), c(1,2)])

# Last add the chromosome location to the FLXswapResF data frame
FLXswapResF$chr <- character(15505)
for(i in 1:15505) { FLXswapResF$chr[FLXswapResF$gene.id == chrdataFswap[i, 1]] <- chrdataFswap[i, 2]}

# Clean up the chromosome location column
FLXswapResF$chr <- gsub("\\;.*", "", as.character(FLXswapResF$chr))


# Now save the full new data frame as a csv file
write.csv(FLXswapResF, "FLXswapFres.csv", row.names = F)

