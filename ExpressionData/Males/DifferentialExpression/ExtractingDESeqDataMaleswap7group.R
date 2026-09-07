################################################################################################
########### DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP MALES - THE SEVEN GROUPS ##########
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


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
#Then select the transcripts that are significantly different
FLXsexFLXaMsig <- subset(FLXsexFLXaM, padj < 0.05)
summary(FLXsexFLXaMsig)
# There are 194 transcripts that are different, 107 up-regulated and 87 down-regulated
# Up-regulated transcripts
FLXsexFLXaMsigup <- subset(FLXsexFLXaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexFLXaMsigdown <- subset(FLXsexFLXaMsig, log2FoldChange < 0)


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


# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaMsig <- subset(FLXsexCFMaM, padj < 0.05)
summary(FLXsexCFMaMsig)
# There are 168 transcripts that are different, 71 up-regulated and 97 down-regulated
# Up-regulated transcripts
FLXsexCFMaMsigup <- subset(FLXsexCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMaMsigdown <- subset(FLXsexCFMaMsig, log2FoldChange < 0)


# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXaCFMsexMsig <- subset(FLXaCFMsexM, padj < 0.05)
summary(FLXaCFMsexMsig)
# There are 83 transcripts that are different, 35 up-regulated and 48 down-regulated
# Up-regulated transcripts
FLXaCFMsexMsigup <- subset(FLXaCFMsexMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMsexMsigdown <- subset(FLXaCFMsexMsig, log2FoldChange < 0)


# FLX-autosomes and CFM-autosomes
FLXaCFMaM <- results(FLXswapM, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaMsig <- subset(FLXaCFMaM, padj < 0.05)
summary(FLXaCFMaMsig)
# There are 186 transcripts that are different, 65 up-regulated and 121 down-regulated
# Up-regulated transcripts
FLXaCFMaMsigup <- subset(FLXaCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaMsigdown <- subset(FLXaCFMaMsig, log2FoldChange < 0)


# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaM <- results(FLXswapM, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
CFMsexCFMaMsig <- subset(CFMsexCFMaM, padj < 0.05)
summary(CFMsexCFMaMsig)
# There are 265 transcripts that are different, 113 up-regulated and 152 down-regulated
# Up-regulated transcripts
CFMsexCFMaMsigup <- subset(CFMsexCFMaMsig, log2FoldChange > 0)
# Down-regulated transcripts
CFMsexCFMaMsigdown <- subset(CFMsexCFMaMsig, log2FoldChange < 0)




########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########


### FIRST SEVEN FEMALE GROUPS IN THE MALES ###

# First collect the transcripts from the last DGElist
genesMswap <- rownames(FLXswapM)

# Time to create a new data set
FLXswap7groupsM <- data.frame(gene.id = genesMswap, gA = numeric(15505), gAup = numeric(15505), gAdown = numeric(15505),
                              gB = numeric(15505), gBup = numeric(15505), gBdown = numeric(15505),
                              gC = numeric(15505), gCup = numeric(15505), gCdown = numeric(15505),
                              gD = numeric(15505), gDup = numeric(15505), gDdown = numeric(15505),
                              gE = numeric(15505), gEup = numeric(15505), gEdown = numeric(15505),
                              gF = numeric(15505), gFup = numeric(15505), gFdown = numeric(15505),
                              gG = numeric(15505), gGup = numeric(15505), gGdown = numeric(15505))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts for group A
FLXswap7groupsM$gA <- "NotSig"
# To be able to add sig make a list of transcripts
GAm.1 <- intersect(intersect(rownames(FLXsexCFMaMsig), rownames(FLXaCFMaMsig)), rownames(CFMsexCFMaMsig))
genesMGA <- setdiff(setdiff(setdiff(GAm.1, rownames(FLXsexFLXaMsig)), rownames(FLXsexCFMsexMsig)), rownames(FLXaCFMsexMsig))
# Then add Sig for transcripts that are significant
FLXswap7groupsM$gA[which(FLXswap7groupsM$gene.id %in% genesMGA)] <- "Sig"
# Make this a factor
FLXswap7groupsM$gA <- factor(FLXswap7groupsM$gA)
# Next add the up- and down-regulated transcripts separately
FLXswap7groupsM$gAup <- "NotSig"
FLXswap7groupsM$gAdown <- "NotSig"
genesMGAup <- intersect(genesMGA, rownames(FLXsexCFMaMsigup))
genesMGAdown <- intersect(genesMGA, rownames(FLXsexCFMaMsigdown))
FLXswap7groupsM$gAup[which(FLXswap7groupsM$gene.id %in% genesMGAup)] <- "Sig"
FLXswap7groupsM$gAdown[which(FLXswap7groupsM$gene.id %in% genesMGAdown)] <- "Sig"
FLXswap7groupsM$gAup <- factor(FLXswap7groupsM$gAup)
FLXswap7groupsM$gAdown <- factor(FLXswap7groupsM$gAdown)

# Then I add the transcripts from group B
FLXswap7groupsM$gB <- "NotSig"
genesMGB <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(CFMsexCFMaMsig), rownames(FLXsexFLXaMsig)),
                                            rownames(FLXsexCFMsexMsig)),
                                    rownames(FLXsexCFMaMsig)),
                            rownames(FLXaCFMsexMsig)),
                    rownames(FLXaCFMaMsig))
FLXswap7groupsM$gB[which(FLXswap7groupsM$gene.id %in% genesMGB)] <- "Sig"
FLXswap7groupsM$gB <- factor(FLXswap7groupsM$gB)
FLXswap7groupsM$gBup <- "NotSig"
FLXswap7groupsM$gBdown <- "NotSig"
genesMGBup <- intersect(genesMGB, rownames(CFMsexCFMaMsigup))
genesMGBdown <- intersect(genesMGB, rownames(CFMsexCFMaMsigdown))
FLXswap7groupsM$gBup[which(FLXswap7groupsM$gene.id %in% genesMGBup)] <- "Sig"
FLXswap7groupsM$gBdown[which(FLXswap7groupsM$gene.id %in% genesMGBdown)] <- "Sig"
FLXswap7groupsM$gBup <- factor(FLXswap7groupsM$gBup)
FLXswap7groupsM$gBdown <- factor(FLXswap7groupsM$gBdown)

# Group C
FLXswap7groupsM$gC <- "NotSig"
genesMGC <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXaCFMaMsig), rownames(FLXsexFLXaMsig)),
                                            rownames(FLXsexCFMsexMsig)),
                                    rownames(FLXsexCFMaMsig)),
                            rownames(FLXaCFMsexMsig)),
                    rownames(CFMsexCFMaMsig))
FLXswap7groupsM$gC[which(FLXswap7groupsM$gene.id %in% genesMGC)] <- "Sig"
FLXswap7groupsM$gC <- factor(FLXswap7groupsM$gC)
FLXswap7groupsM$gCup <- "NotSig"
FLXswap7groupsM$gCdown <- "NotSig"
genesMGCup <- intersect(genesMGC, rownames(FLXaCFMaMsigup))
genesMGCdown <- intersect(genesMGC, rownames(FLXaCFMaMsigdown))
FLXswap7groupsM$gCup[which(FLXswap7groupsM$gene.id %in% genesMGCup)] <- "Sig"
FLXswap7groupsM$gCdown[which(FLXswap7groupsM$gene.id %in% genesMGCdown)] <- "Sig"
FLXswap7groupsM$gCup <- factor(FLXswap7groupsM$gCup)
FLXswap7groupsM$gCdown <- factor(FLXswap7groupsM$gCdown)

# Group D
FLXswap7groupsM$gD <- "NotSig"
GDm.1 <- intersect(rownames(FLXaCFMaMsig), rownames(CFMsexCFMaMsig))
genesMGD <- setdiff(setdiff(setdiff(setdiff(GDm.1, rownames(FLXsexFLXaMsig)),
                                    rownames(FLXsexCFMsexMsig)),
                            rownames(FLXsexCFMaMsig)),
                    rownames(FLXaCFMsexMsig))
FLXswap7groupsM$gD[which(FLXswap7groupsM$gene.id %in% genesMGD)] <- "Sig"
FLXswap7groupsM$gD <- factor(FLXswap7groupsM$gD)
FLXswap7groupsM$gDup <- "NotSig"
FLXswap7groupsM$gDdown <- "NotSig"
genesMGDup <- intersect(genesMGD, rownames(FLXaCFMaMsigup))
genesMGDdown <- intersect(genesMGD, rownames(FLXaCFMaMsigdown))
FLXswap7groupsM$gDup[which(FLXswap7groupsM$gene.id %in% genesMGDup)] <- "Sig"
FLXswap7groupsM$gDdown[which(FLXswap7groupsM$gene.id %in% genesMGDdown)] <- "Sig"
FLXswap7groupsM$gDup <- factor(FLXswap7groupsM$gDup)
FLXswap7groupsM$gDdown <- factor(FLXswap7groupsM$gDdown)

# Group E
FLXswap7groupsM$gE <- "NotSig"
genesMGE <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXsexCFMaMsig), rownames(FLXsexFLXaMsig)),
                                            rownames(FLXsexCFMsexMsig)),
                                    rownames(FLXaCFMsexMsig)),
                            rownames(FLXaCFMaMsig)),
                    rownames(CFMsexCFMaMsig))
FLXswap7groupsM$gE[which(FLXswap7groupsM$gene.id %in% genesMGE)] <- "Sig"
FLXswap7groupsM$gE <- factor(FLXswap7groupsM$gE)
FLXswap7groupsM$gEup <- "NotSig"
FLXswap7groupsM$gEdown <- "NotSig"
genesMGEup <- intersect(genesMGE, rownames(FLXsexCFMaMsigup))
genesMGEdown <- intersect(genesMGE, rownames(FLXsexCFMaMsigdown))
FLXswap7groupsM$gEup[which(FLXswap7groupsM$gene.id %in% genesMGEup)] <- "Sig"
FLXswap7groupsM$gEdown[which(FLXswap7groupsM$gene.id %in% genesMGEdown)] <- "Sig"
FLXswap7groupsM$gEup <- factor(FLXswap7groupsM$gEup)
FLXswap7groupsM$gEdown <- factor(FLXswap7groupsM$gEdown)

# Group F
FLXswap7groupsM$gF <- "NotSig"
GFm.1 <- intersect(rownames(FLXsexCFMaMsig), rownames(FLXaCFMaMsig))
genesMGF <- setdiff(setdiff(setdiff(setdiff(GFm.1, rownames(FLXsexFLXaMsig)),
                                    rownames(FLXsexCFMsexMsig)),
                            rownames(FLXaCFMsexMsig)),
                    rownames(CFMsexCFMaMsig))
FLXswap7groupsM$gF[which(FLXswap7groupsM$gene.id %in% genesMGF)] <- "Sig"
FLXswap7groupsM$gF <- factor(FLXswap7groupsM$gF)
FLXswap7groupsM$gFup <- "NotSig"
FLXswap7groupsM$gFdown <- "NotSig"
genesMGFup <- intersect(genesMGF, rownames(FLXsexCFMaMsigup))
genesMGFdown <- intersect(genesMGF, rownames(FLXsexCFMaMsigdown))
FLXswap7groupsM$gFup[which(FLXswap7groupsM$gene.id %in% genesMGFup)] <- "Sig"
FLXswap7groupsM$gFdown[which(FLXswap7groupsM$gene.id %in% genesMGFdown)] <- "Sig"
FLXswap7groupsM$gFup <- factor(FLXswap7groupsM$gFup)
FLXswap7groupsM$gFdown <- factor(FLXswap7groupsM$gFdown)

# Group G
FLXswap7groupsM$gG <- "NotSig"
GGm.1 <- intersect(rownames(FLXsexCFMaMsig), rownames(CFMsexCFMaMsig))
genesMGG <- setdiff(setdiff(setdiff(setdiff(GGm.1, rownames(FLXsexFLXaMsig)),
                                    rownames(FLXsexCFMsexMsig)),
                            rownames(FLXaCFMsexMsig)),
                    rownames(FLXaCFMaMsig))
FLXswap7groupsM$gG[which(FLXswap7groupsM$gene.id %in% genesMGG)] <- "Sig"
FLXswap7groupsM$gG <- factor(FLXswap7groupsM$gG)
FLXswap7groupsM$gGup <- "NotSig"
FLXswap7groupsM$gGdown <- "NotSig"
genesMGGup <- intersect(genesMGG, rownames(FLXsexCFMaMsigup))
genesMGGdown <- intersect(genesMGG, rownames(FLXsexCFMaMsigdown))
FLXswap7groupsM$gGup[which(FLXswap7groupsM$gene.id %in% genesMGGup)] <- "Sig"
FLXswap7groupsM$gGdown[which(FLXswap7groupsM$gene.id %in% genesMGGdown)] <- "Sig"
FLXswap7groupsM$gGup <- factor(FLXswap7groupsM$gGup)
FLXswap7groupsM$gGdown <- factor(FLXswap7groupsM$gGdown)


# Add chromosome location to the data set as well
# First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
# Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
# Now make a new vector with the chromosomes matching the FLXswap7groupsM data frame
chrdataMswap <- cbind(genesMswap, chrIDswap[match(FLXswap7groupsM$gene.id, rownames(FLXswapreadcount)), c(1,2)])

# Last add the chromosome location to the FLXswap7groupsM data frame
FLXswap7groupsM$chr <- character(15505)
for(i in 1:15505) { FLXswap7groupsM$chr[FLXswap7groupsM$gene.id == chrdataMswap[i, 1]] <- chrdataMswap[i, 2]}

# Clean up the chromosome location column
FLXswap7groupsM$chr <- gsub("\\;.*", "", as.character(FLXswap7groupsM$chr))


#Now save the full new data frame as a csv file
write.csv(FLXswap7groupsM, "FLXswapMres7G.csv", row.names = F)

