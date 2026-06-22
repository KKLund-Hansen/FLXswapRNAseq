################################################################################################
############################ HEATMAP FLX CHROMOSOME SWAP BOTH SEXES ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Flies/FLX/5.FLXchromosomeSwap/R/Both/Heatmap")

#First install the software from Bioconductor
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


########## GENOTYPES ##########

# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXa <- results(FLXswap, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXasig <- subset(FLXsexFLXa, padj < 0.05)
summary(FLXsexFLXasig)
# There are 571 transcripts that are different, 194 up-regulated and 377 down-regulated

# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsex <- results(FLXswap, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexsig <- subset(FLXsexCFMsex, padj < 0.05)
summary(FLXsexCFMsexsig)
# There are 691 transcripts that are different, 482 up-regulated and 209 down-regulated

# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMa <- results(FLXswap, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMasig <- subset(FLXsexCFMa, padj < 0.05)
summary(FLXsexCFMasig)
# There are 5239 transcripts that are different, 300 up-regulated and 4939 down-regulated

# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsex <- results(FLXswap, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXaCFMsexsig <- subset(FLXaCFMsex, padj < 0.05)
summary(FLXaCFMsexsig)
# There are 633 transcripts that are different, 447 up-regulated and 186 down-regulated

# FLX-autosomes and CFM-autosomes
FLXaCFMa <- results(FLXswap, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMasig <- subset(FLXaCFMa, padj < 0.05)
summary(FLXaCFMasig)
# There are 4637 transcripts that are different, 166 up-regulated and 4471 down-regulated

# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMa <- results(FLXswap, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
CFMsexCFMasig <- subset(CFMsexCFMa, padj < 0.05)
summary(CFMsexCFMasig)
# There are 5715 transcripts that are different, 728 up-regulated and 4987 down-regulated




#########################################  PLOT DATA  #########################################


# First make a new data set with all the significant transcripts

# First collect all the significant transcripts in one vector
genesSwapOverall <- unique(cbind(c(rownames(FLXsexFLXasig), rownames(FLXsexCFMsexsig), rownames(FLXsexCFMasig),
                                   rownames(FLXaCFMsexsig), rownames(FLXaCFMasig), rownames(CFMsexCFMasig))))

# Next make two new data sets with the counts for the significant transcripts
heatmap.dataswap <- cbind(typesex = paste(FLXswapsampleinfo$type, FLXswapsampleinfo$sex), as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesSwapOverall, ])))
heatmap.dataswaprp <- cbind(typesexRP = paste(FLXswapsampleinfo$type, FLXswapsampleinfo$sex, FLXswapsampleinfo$rep_pop),
                            as.data.frame(t(FLXswapreadcount[rownames(FLXswapreadcount) %in% genesSwapOverall, ])))

# Save the new data sets
write.csv(heatmap.dataswap, "heatmapswap.csv", row.names = F)
write.csv(heatmap.dataswaprp, "heatmapswaprp.csv", row.names = F)



# Now it is time to make the heat map, I'll make two one with types and one with the rep pop as well

# Read in files with data
heatmapswap <- read.table(file = "heatmapswap.csv", h = T, sep = ",", stringsAsFactors = T)
heatmapswaprp <- read.table(file = "heatmapswaprp.csv", h = T, sep = ",", stringsAsFactors = T)


### Type only
# Make a matrix to add the mean data into
M.heatmapswap <- matrix(0, 6961, 8)
# Add the names of the types and rep pop
colnames(M.heatmapswap) <- c("FLXsex F", "FLXsex M", "FLXauto F",  "FLXauto M",
                             "CFMsex F", "CFMsex M", "CFMauto F", "CFMauto M")
# and the transcripts names
rownames(M.heatmapswap) <- colnames(heatmapswap[,2:6962])
# Calculate the means for each transcripts
for (i in 1:6961) M.heatmapswap[i,] <- as.numeric(tapply(heatmapswap[,i+1], heatmapswap$typesex, mean))


### Type and rep pop
# Make a matrix to add the mean data into
M.heatmapswaprp <- matrix(0, 6961, 32)
# Add the names of the types and rep pop
colnames(M.heatmapswaprp) <- c("FLXsex-F-RP1", "FLXsex-F-RP2", "FLXsex-F-RP3", "FLXsex-F-RP4",
                               "FLXsex-M-RP1", "FLXsex-M-RP2", "FLXsex-M-RP3", "FLXsex-M-RP4",
                               "FLXauto-F-RP1", "FLXauto-F-RP2", "FLXauto-F-RP3", "FLXauto-F-RP4",
                               "FLXauto-M-RP1", "FLXauto-M-RP2", "FLXauto-M-RP3", "FLXauto-M-RP4",
                               "CFMsex-F-RP1", "CFMsex-F-RP2", "CFMsex-F-RP3", "CFMsex-F-RP4",
                               "CFMsex-M-RP1", "CFMsex-M-RP2", "CFMsex-M-RP3", "CFMsex-M-RP4",
                               "CFMauto-F-RP1", "CFMauto-F-RP2", "CFMauto-F-RP3", "CFMauto-F-RP4",
                               "CFMauto-M-RP1", "CFMauto-M-RP2", "CFMauto-M-RP3", "CFMauto-M-RP4")
#and the transcripts names
rownames(M.heatmapswaprp) <- colnames(heatmapswaprp[,2:6962])
#Calculate the means for each transcripts
for (i in 1:6961) M.heatmapswaprp[i,] <- as.numeric(tapply(heatmapswaprp[,i+1], heatmapswaprp$typesexRP, mean))


#### PLOTS ####


### Type only

pdf("HMswap.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapswap, Rowv = NA, col = brewer.pal(9, "Greys"), labRow = "", margins = c(6, 2), scale = "row", cexCol = 1.4)

dev.off()


### Type and rep pop

pdf("HMswaprp.pdf", family = "serif", width = 5, height = 8)

par(mar = c(5, 2, 2, 2) + 0.1)
heatmap(M.heatmapswaprp, Rowv = NA, col = brewer.pal(9, "Greys"), labRow = "", margins = c(6, 2), scale = "row")

dev.off()


