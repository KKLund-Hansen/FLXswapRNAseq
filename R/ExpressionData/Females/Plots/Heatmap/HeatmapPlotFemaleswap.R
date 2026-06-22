################################################################################################
############################## HEATMAP FLX CHROMOSOME SWAP FEMALES #############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Flies/FLX/5.FLXchromosomeSwap/R/Females/Plots/Heatmap")

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


######################## FEMALES ########################

# As we want to analyse the data separately per sex, we start by sub-setting the data

# First, separate the data into males
FLXswapF <- FLXswap[,FLXswap$sex %in% "f"]

# Next drop sex as a level, as it is not need anymore
FLXswapF$sex <- droplevels(FLXswapF$sex)

# Last thing is to update the model to not contain sex
design(FLXswapF) <- formula(~ type + type:nested.RP)

# And the run the DESeq again
FLXswapF <- DESeq(FLXswapF)


########## FEMALE VS FEMALE ##########


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaFsig <- subset(FLXsexFLXaF, padj < 0.05)
summary(FLXsexFLXaFsig)
# There are 287 transcripts that are different, 199 up-regulated and 88, down-regulated
# Up-regulated transcripts
FLXsexFLXaFsigup <- subset(FLXsexFLXaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexFLXaFsigdown <- subset(FLXsexFLXaFsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexFsig <- subset(FLXsexCFMsexF, padj < 0.05)
summary(FLXsexCFMsexFsig)
# There are 204 transcripts that are different, 143 up-regulated and 61 down-regulated
# Up-regulated transcripts
FLXsexCFMsexFsigup <- subset(FLXsexCFMsexFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMsexFsigdown <- subset(FLXsexCFMsexFsig, log2FoldChange < 0)


# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaFsig <- subset(FLXsexCFMaF, padj < 0.05)
summary(FLXsexCFMaFsig)
# There are 3388 transcripts that are different, 531 up-regulated and 2857 down-regulated
# Up-regulated transcripts
FLXsexCFMaFsigup <- subset(FLXsexCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMaFsigdown <- subset(FLXsexCFMaFsig, log2FoldChange < 0)


# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
#Then select the transcripts that are significantly different
FLXaCFMsexFsig <- subset(FLXaCFMsexF, padj < 0.05)
summary(FLXaCFMsexFsig)
# There are 434 transcripts that are different, 179 up-regulated and 255 down-regulated
# Up-regulated transcripts
FLXaCFMsexFsigup <- subset(FLXaCFMsexFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMsexFsigdown <- subset(FLXaCFMsexFsig, log2FoldChange < 0)


# FLX-autosomes and CFM-autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaFsig <- subset(FLXaCFMaF, padj < 0.05)
summary(FLXaCFMaFsig)
# There are 3631 transcripts that are different, 569 up-regulated and 3062 down-regulated
# Up-regulated transcripts
FLXaCFMaFsigup <- subset(FLXaCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaFsigdown <- subset(FLXaCFMaFsig, log2FoldChange < 0)


# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaF <- results(FLXswapF, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
#Then select the transcripts that are significantly different
CFMsexCFMaFsig <- subset(CFMsexCFMaF, padj < 0.05)
summary(CFMsexCFMaFsig)
# There are 4234 transcripts that are different, 1162 up-regulated and 3072 down-regulated
# Up-regulated transcripts
CFMsexCFMaFsigup <- subset(CFMsexCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
CFMsexCFMaFsigdown <- subset(CFMsexCFMaFsig, log2FoldChange < 0)




#########################################  PLOT DATA  #########################################


# First make a new data set with all the significant transcripts for the females

# First collect all the significant transcripts in one vector
genesFswapOverall <- unique(cbind(c(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig),
                                    rownames(FLXaCFMsexFsig), rownames(FLXaCFMaFsig), rownames(CFMsexCFMaFsig))))
# 5725

# Next make two new data sets with the counts for the significant transcripts
heatmap.dataFswap <- cbind(type = FLXswapsampleinfo$type, sex = FLXswapsampleinfo$sex, as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesFswapOverall, ])))
heatmap.dataFswaprp <- cbind(typeRP = paste(FLXswapsampleinfo$type, FLXswapsampleinfo$rep_pop), sex = FLXswapsampleinfo$sex,
                             as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesFswapOverall, ])))

# Last subset the data to only include females
heatmapFswap <- subset(heatmap.dataFswap, sex == "f")
heatmapFswaprp <- subset(heatmap.dataFswaprp, sex == "f")

# Save the new data sets
write.csv(heatmapFswap, "heatmapFswap.csv", row.names = F)
write.csv(heatmapFswaprp, "heatmapFswaprp.csv", row.names = F)


# Now it is time to make the heatmap, I'll make two one with types and one with the rep pop as well


# Read in files with data
heatmapFswap <- read.table(file = "heatmapFswap.csv", h = T, sep = ",", stringsAsFactors = T)
heatmapFswaprp <- read.table(file = "heatmapFswaprp.csv", h = T, sep = ",", stringsAsFactors = T)


### Type only
# Make a matrix to add the mean data into
M.heatmapFswap <- matrix(0, 5725, 4)
# Add the names of the types and rep pop
colnames(M.heatmapFswap) <- c("FLXsex", "FLXauto", "CFMsex", "CFMauto")
# and the transcripts names
rownames(M.heatmapFswap) <- colnames(heatmapFswap[,3:5727])
# Calculate the means for each transcripts
for (i in 1:5725) M.heatmapFswap[i,] <- as.numeric(tapply(heatmapFswap[,i+2], heatmapFswap$type, mean))


### Type and rep pop
# Make a matrix to add the mean data into
M.heatmapFswaprp <- matrix(0, 5725, 16)
# Add the names of the types and rep pop
colnames(M.heatmapFswaprp) <- c("FLXsex-RP1", "FLXsex-RP2", "FLXsex-RP3", "FLXsex-RP4",
                                "FLXauto-RP1", "FLXauto-RP2", "FLXauto-RP3", "FLXauto-RP4",
                                "CFMsex-RP1", "CFMsex-RP2", "CFMsex-RP3", "CFMsex-RP4",
                                "CFMauto-RP1", "CFMauto-RP2", "CFMauto-RP3", "CFMauto-RP4")
# and the transcripts names
rownames(M.heatmapFswaprp) <- colnames(heatmapFswaprp[,3:5727])
# Calculate the means for each transcripts
for (i in 1:5725) M.heatmapFswaprp[i,] <- as.numeric(tapply(heatmapFswaprp[,i+2], heatmapFswaprp$typeRP, mean))


#### PLOTS ####


### Type only

pdf("HMfswap.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapFswap, Rowv = NA, col = brewer.pal(9, "Greys"), labRow = "", margins = c(6, 2), scale = "row", cexCol = 1.4)

dev.off()



### Type and rep pop

pdf("HMfswaprp.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapFswaprp, Rowv = NA, col = colF, labRow = "", margins = c(6, 2), scale = "row")

dev.off()


