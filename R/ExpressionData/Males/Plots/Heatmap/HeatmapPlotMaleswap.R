################################################################################################
############################### HEATMAP FLX CHROMOSOME SWAP MALES ##############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Flies/FLX/5.FLXchromosomeSwap/R/Males/Plots/Heatmap")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

# Set up environment
library(DESeq2)
library(RColorBrewer)

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

# Next we do a very simple prefiltering and remove all transcripts with a count less than 1
keep <- rowSums(counts(FLXswapsex)) >= 1
FLXswapsex <- FLXswapsex[keep,]

# Check how many transcripts are left
nrow(FLXswapsex)
# 16785

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


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaMsig <- subset(FLXsexFLXaM, padj < 0.05)
summary(FLXsexFLXaMsig)
# There are 191 transcripts that are different, 106 up-regulated and 85 down-regulated
# Up-regulated transcripts
FLXsexFLXaMsigup <- subset(FLXsexFLXaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexFLXaMsigdown <- subset(FLXsexFLXaMsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexMsig <- subset(FLXsexCFMsexM, padj < 0.05)
summary(FLXsexCFMsexMsig)
# There are 100 transcripts that are different, 53 up-regulated and 47 down-regulated
# Up-regulated transcripts
FLXsexCFMsexMsigup <- subset(FLXsexCFMsexMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMsexMsigdown <- subset(FLXsexCFMsexMsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaMsig <- subset(FLXsexCFMaM, padj < 0.05)
summary(FLXsexCFMaMsig)
# There are 161 transcripts that are different, 68 up-regulated and 93 down-regulated
# Up-regulated transcripts
FLXsexCFMaMsigup <- subset(FLXsexCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMaMsigdown <- subset(FLXsexCFMaMsig, log2FoldChange < 0)


# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXaCFMsexMsig <- subset(FLXaCFMsexM, padj < 0.05)
summary(FLXaCFMsexMsig)
# There are 80 transcripts that are different, 35 up-regulated and 45 down-regulated
# Up-regulated transcripts
FLXaCFMsexMsigup <- subset(FLXaCFMsexMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMsexMsigdown <- subset(FLXaCFMsexMsig, log2FoldChange < 0)


# FLX-autosomes and CFM-autosomes
FLXaCFMaM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaMsig <- subset(FLXaCFMaM, padj < 0.05)
summary(FLXaCFMaMsig)
# There are 185 transcripts that are different, 65 up-regulated and 120 down-regulated
# Up-regulated transcripts
FLXaCFMaMsigup <- subset(FLXaCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaMsigdown <- subset(FLXaCFMaMsig, log2FoldChange < 0)


# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaM <- results(FLXswapM, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
CFMsexCFMaMsig <- subset(CFMsexCFMaM, padj < 0.05)
summary(CFMsexCFMaMsig)
# There are 266 transcripts that are different, 115 up-regulated and 151 down-regulated
# Up-regulated transcripts
CFMsexCFMaMsigup <- subset(CFMsexCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
CFMsexCFMaMsigdown <- subset(CFMsexCFMaMsig, log2FoldChange < 0)




#########################################  PLOT DATA  #########################################


# First make a new data set with all the significant transcripts for the males

# First collect all the significant transcripts in one vector
genesMswapOverall <- unique(cbind(c(rownames(FLXsexFLXaMsig), rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig),
                                    rownames(FLXaCFMsexMsig), rownames(FLXaCFMaMsig), rownames(CFMsexCFMaMsig))))

# Next make two new data sets with the counts for the significant transcripts
heatmap.dataMswap <- cbind(type = FLXswapsampleinfo$type, sex = FLXswapsampleinfo$sex, as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesMswapOverall, ])))
heatmap.dataMswaprp <- cbind(typeRP = paste(FLXswapsampleinfo$type, FLXswapsampleinfo$rep_pop), sex = FLXswapsampleinfo$sex,
                             as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesMswapOverall, ])))

# Last subset the data to only include males
heatmapMswap <- subset(heatmap.dataMswap, sex == "m")
heatmapMswaprp <- subset(heatmap.dataMswaprp, sex == "m")

# Save the new data sets
write.csv(heatmapMswap, "heatmapMswap.csv", row.names = F)
write.csv(heatmapMswaprp, "heatmapMswaprp.csv", row.names = F)


# Now it is time to make the heatmap, I'll make two one with types and one with the rep pop as well


# Read in files with data
heatmapMswap <- read.table(file = "heatmapMswap.csv", h = T, sep = ",", stringsAsFactors = T)
heatmapMswaprp <- read.table(file = "heatmapMswaprp.csv", h = T, sep = ",", stringsAsFactors = T)


### Type only
# Make a matrix to add the mean data into
M.heatmapMswap <- matrix(0, 588, 4)
# Add the names of the types and rep pop
colnames(M.heatmapMswap) <- c("FLXsex", "FLXauto", "CFMsex", "CFMauto")
# and the transcripts names
rownames(M.heatmapMswap) <- colnames(heatmapMswap[,3:590])
# Calculate the means for each transcripts
for (i in 1:588) M.heatmapMswap[i,] <- as.numeric(tapply(heatmapMswap[,i+2], heatmapMswap$type, mean))


### Type and rep pop
# Make a matrix to add the mean data into
M.heatmapMswaprp <- matrix(0, 588, 16)
# Add the names of the types and rep pop
colnames(M.heatmapMswaprp) <- c("FLXsex-RP1", "FLXsex-RP2", "FLXsex-RP3", "FLXsex-RP4",
                                "FLXauto-RP1", "FLXauto-RP2", "FLXauto-RP3", "FLXauto-RP4",
                                "CFMsex-RP1", "CFMsex-RP2", "CFMsex-RP3", "CFMsex-RP4",
                                "CFMauto-RP1", "CFMauto-RP2", "CFMauto-RP3", "CFMauto-RP4")
# and the transcripts names
rownames(M.heatmapMswaprp) <- colnames(heatmapMswaprp[,3:590])
# Calculate the means for each transcripts
for (i in 1:588) M.heatmapMswaprp[i,] <- as.numeric(tapply(heatmapMswaprp[,i+2], heatmapMswaprp$typeRP, mean))


#### PLOTS ####

### Type only

pdf("HMmswap.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapMswap, Rowv = NA, col =  brewer.pal(9, "Greys"), labRow = "", margins = c(6, 2), scale = "row", cexCol = 1.4)

dev.off()



### Type and rep pop

pdf("HMmswaprp.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapMswaprp, Rowv = NA, col = colM, labRow = "", margins = c(6, 2), scale = "row")

dev.off()


