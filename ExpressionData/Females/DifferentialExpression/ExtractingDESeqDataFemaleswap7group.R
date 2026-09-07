################################################################################################
########## DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP FEMALES - THE SEVEN GROUPS #########
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


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaFsig <- subset(FLXsexFLXaF, padj < 0.05)
summary(FLXsexFLXaFsig)
# There are 277 transcripts that are different, 193 up-regulated and 84 down-regulated
# Up-regulated transcripts
FLXsexFLXaFsigup <- subset(FLXsexFLXaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexFLXaFsigdown <- subset(FLXsexFLXaFsig, log2FoldChange < 0)


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


# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaFsig <- subset(FLXsexCFMaF, padj < 0.05)
summary(FLXsexCFMaFsig)
# There are 3293 transcripts that are different, 508 up-regulated and 2785 down-regulated
# Up-regulated transcripts
FLXsexCFMaFsigup <- subset(FLXsexCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXsexCFMaFsigdown <- subset(FLXsexCFMaFsig, log2FoldChange < 0)


# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
#Then select the transcripts that are significantly different
FLXaCFMsexFsig <- subset(FLXaCFMsexF, padj < 0.05)
summary(FLXaCFMsexFsig)
# There are 411 transcripts that are different, 166 up-regulated and 245 down-regulated
# Up-regulated transcripts
FLXaCFMsexFsigup <- subset(FLXaCFMsexFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMsexFsigdown <- subset(FLXaCFMsexFsig, log2FoldChange < 0)


# FLX-autosomes and CFM-autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaFsig <- subset(FLXaCFMaF, padj < 0.05)
summary(FLXaCFMaFsig)
# There are 3519 transcripts that are different, 525 up-regulated and 2994 down-regulated
# Up-regulated transcripts
FLXaCFMaFsigup <- subset(FLXaCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
FLXaCFMaFsigdown <- subset(FLXaCFMaFsig, log2FoldChange < 0)


# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaF <- results(FLXswapF, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
#Then select the transcripts that are significantly different
CFMsexCFMaFsig <- subset(CFMsexCFMaF, padj < 0.05)
summary(CFMsexCFMaFsig)
# There are 4137  transcripts that are different, 1130 up-regulated and 3007 down-regulated
# Up-regulated transcripts
CFMsexCFMaFsigup <- subset(CFMsexCFMaFsig, log2FoldChange > 0)
# Down-regulated transcripts
CFMsexCFMaFsigdown <- subset(CFMsexCFMaFsig, log2FoldChange < 0)



########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########


### FIRST SEVEN FEMALE GROUPS ###

# First collect the transcripts from the last DGElist
genesFswap <- rownames(FLXswapF)

# Time to create a new data set
FLXswap7groupsF <- data.frame(gene.id = genesFswap, gA = numeric(15505), gAup = numeric(15505), gAdown = numeric(15505),
                              gB = numeric(15505), gBup = numeric(15505), gBdown = numeric(15505),
                              gC = numeric(15505), gCup = numeric(15505), gCdown = numeric(15505),
                              gD = numeric(15505), gDup = numeric(15505), gDdown = numeric(15505),
                              gE = numeric(15505), gEup = numeric(15505), gEdown = numeric(15505),
                              gF = numeric(15505), gFup = numeric(15505), gFdown = numeric(15505),
                              gG = numeric(15505), gGup = numeric(15505), gGdown = numeric(15505))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts for group A
FLXswap7groupsF$gA <- "NotSig"
# To be able to add sig make a list of transcripts
GAf.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFGA <- setdiff(setdiff(setdiff(GAf.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXaCFMsexFsig))
# Then add Sig for transcripts that are significant
FLXswap7groupsF$gA[which(FLXswap7groupsF$gene.id %in% genesFGA)] <- "Sig"
# Make this a factor
FLXswap7groupsF$gA <- factor(FLXswap7groupsF$gA)
# Next add the up- and down-regulated transcripts separately
FLXswap7groupsF$gAup <- "NotSig"
FLXswap7groupsF$gAdown <- "NotSig"
genesFGAup <- intersect(genesFGA, rownames(FLXsexCFMaFsigup))
genesFGAdown <- intersect(genesFGA, rownames(FLXsexCFMaFsigdown))
FLXswap7groupsF$gAup[which(FLXswap7groupsF$gene.id %in% genesFGAup)] <- "Sig"
FLXswap7groupsF$gAdown[which(FLXswap7groupsF$gene.id %in% genesFGAdown)] <- "Sig"
FLXswap7groupsF$gAup <- factor(FLXswap7groupsF$gAup)
FLXswap7groupsF$gAdown <- factor(FLXswap7groupsF$gAdown)

# Then I add the transcripts from group B
FLXswap7groupsF$gB <- "NotSig"
genesFGB <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(CFMsexCFMaFsig), rownames(FLXsexFLXaFsig)),
                                            rownames(FLXsexCFMsexFsig)),
                                    rownames(FLXsexCFMaFsig)),
                            rownames(FLXaCFMsexFsig)),
                    rownames(FLXaCFMaFsig))
FLXswap7groupsF$gB[which(FLXswap7groupsF$gene.id %in% genesFGB)] <- "Sig"
FLXswap7groupsF$gB <- factor(FLXswap7groupsF$gB)
FLXswap7groupsF$gBup <- "NotSig"
FLXswap7groupsF$gBdown <- "NotSig"
genesFGBup <- intersect(genesFGB, rownames(CFMsexCFMaFsigup))
genesFGBdown <- intersect(genesFGB, rownames(CFMsexCFMaFsigdown))
FLXswap7groupsF$gBup[which(FLXswap7groupsF$gene.id %in% genesFGBup)] <- "Sig"
FLXswap7groupsF$gBdown[which(FLXswap7groupsF$gene.id %in% genesFGBdown)] <- "Sig"
FLXswap7groupsF$gBup <- factor(FLXswap7groupsF$gBup)
FLXswap7groupsF$gBdown <- factor(FLXswap7groupsF$gBdown)

# Group C
FLXswap7groupsF$gC <- "NotSig"
genesFGC <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXaCFMaFsig), rownames(FLXsexFLXaFsig)),
                                            rownames(FLXsexCFMsexFsig)),
                                    rownames(FLXsexCFMaFsig)),
                            rownames(FLXaCFMsexFsig)),
                    rownames(CFMsexCFMaFsig))
FLXswap7groupsF$gC[which(FLXswap7groupsF$gene.id %in% genesFGC)] <- "Sig"
FLXswap7groupsF$gC <- factor(FLXswap7groupsF$gC)
FLXswap7groupsF$gCup <- "NotSig"
FLXswap7groupsF$gCdown <- "NotSig"
genesFGCup <- intersect(genesFGC, rownames(FLXaCFMaFsigup))
genesFGCdown <- intersect(genesFGC, rownames(FLXaCFMaFsigdown))
FLXswap7groupsF$gCup[which(FLXswap7groupsF$gene.id %in% genesFGCup)] <- "Sig"
FLXswap7groupsF$gCdown[which(FLXswap7groupsF$gene.id %in% genesFGCdown)] <- "Sig"
FLXswap7groupsF$gCup <- factor(FLXswap7groupsF$gCup)
FLXswap7groupsF$gCdown <- factor(FLXswap7groupsF$gCdown)

# Group D
FLXswap7groupsF$gD <- "NotSig"
GDf.1 <- intersect(rownames(FLXaCFMaFsig), rownames(CFMsexCFMaFsig))
genesFGD <- setdiff(setdiff(setdiff(setdiff(GDf.1, rownames(FLXsexFLXaFsig)),
                                    rownames(FLXsexCFMsexFsig)),
                            rownames(FLXsexCFMaFsig)),
                    rownames(FLXaCFMsexFsig))
FLXswap7groupsF$gD[which(FLXswap7groupsF$gene.id %in% genesFGD)] <- "Sig"
FLXswap7groupsF$gD <- factor(FLXswap7groupsF$gD)
FLXswap7groupsF$gDup <- "NotSig"
FLXswap7groupsF$gDdown <- "NotSig"
genesFGDup <- intersect(genesFGD, rownames(FLXaCFMaFsigup))
genesFGDdown <- intersect(genesFGD, rownames(FLXaCFMaFsigdown))
FLXswap7groupsF$gDup[which(FLXswap7groupsF$gene.id %in% genesFGDup)] <- "Sig"
FLXswap7groupsF$gDdown[which(FLXswap7groupsF$gene.id %in% genesFGDdown)] <- "Sig"
FLXswap7groupsF$gDup <- factor(FLXswap7groupsF$gDup)
FLXswap7groupsF$gDdown <- factor(FLXswap7groupsF$gDdown)

# Group E
FLXswap7groupsF$gE <- "NotSig"
genesFGE <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXsexCFMaFsig), rownames(FLXsexFLXaFsig)),
                                            rownames(FLXsexCFMsexFsig)),
                                    rownames(FLXaCFMsexFsig)),
                            rownames(FLXaCFMaFsig)),
                    rownames(CFMsexCFMaFsig))
FLXswap7groupsF$gE[which(FLXswap7groupsF$gene.id %in% genesFGE)] <- "Sig"
FLXswap7groupsF$gE <- factor(FLXswap7groupsF$gE)
FLXswap7groupsF$gEup <- "NotSig"
FLXswap7groupsF$gEdown <- "NotSig"
genesFGEup <- intersect(genesFGE, rownames(FLXsexCFMaFsigup))
genesFGEdown <- intersect(genesFGE, rownames(FLXsexCFMaFsigdown))
FLXswap7groupsF$gEup[which(FLXswap7groupsF$gene.id %in% genesFGEup)] <- "Sig"
FLXswap7groupsF$gEdown[which(FLXswap7groupsF$gene.id %in% genesFGEdown)] <- "Sig"
FLXswap7groupsF$gEup <- factor(FLXswap7groupsF$gEup)
FLXswap7groupsF$gEdown <- factor(FLXswap7groupsF$gEdown)

# Group F
FLXswap7groupsF$gF <- "NotSig"
GFf.1 <- intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMaFsig))
genesFGF <- setdiff(setdiff(setdiff(setdiff(GFf.1, rownames(FLXsexFLXaFsig)),
                                    rownames(FLXsexCFMsexFsig)),
                            rownames(FLXaCFMsexFsig)),
                    rownames(CFMsexCFMaFsig))
FLXswap7groupsF$gF[which(FLXswap7groupsF$gene.id %in% genesFGF)] <- "Sig"
FLXswap7groupsF$gF <- factor(FLXswap7groupsF$gF)
FLXswap7groupsF$gFup <- "NotSig"
FLXswap7groupsF$gFdown <- "NotSig"
genesFGFup <- intersect(genesFGF, rownames(FLXsexCFMaFsigup))
genesFGFdown <- intersect(genesFGF, rownames(FLXsexCFMaFsigdown))
FLXswap7groupsF$gFup[which(FLXswap7groupsF$gene.id %in% genesFGFup)] <- "Sig"
FLXswap7groupsF$gFdown[which(FLXswap7groupsF$gene.id %in% genesFGFdown)] <- "Sig"
FLXswap7groupsF$gFup <- factor(FLXswap7groupsF$gFup)
FLXswap7groupsF$gFdown <- factor(FLXswap7groupsF$gFdown)

# Group G
FLXswap7groupsF$gG <- "NotSig"
GGf.1 <- intersect(rownames(FLXsexCFMaFsig), rownames(CFMsexCFMaFsig))
genesFGG <- setdiff(setdiff(setdiff(setdiff(GGf.1, rownames(FLXsexFLXaFsig)),
                                    rownames(FLXsexCFMsexFsig)),
                            rownames(FLXaCFMsexFsig)),
                    rownames(FLXaCFMaFsig))
FLXswap7groupsF$gG[which(FLXswap7groupsF$gene.id %in% genesFGG)] <- "Sig"
FLXswap7groupsF$gG <- factor(FLXswap7groupsF$gG)
FLXswap7groupsF$gGup <- "NotSig"
FLXswap7groupsF$gGdown <- "NotSig"
genesFGGup <- intersect(genesFGG, rownames(FLXsexCFMaFsigup))
genesFGGdown <- intersect(genesFGG, rownames(FLXsexCFMaFsigdown))
FLXswap7groupsF$gGup[which(FLXswap7groupsF$gene.id %in% genesFGGup)] <- "Sig"
FLXswap7groupsF$gGdown[which(FLXswap7groupsF$gene.id %in% genesFGGdown)] <- "Sig"
FLXswap7groupsF$gGup <- factor(FLXswap7groupsF$gGup)
FLXswap7groupsF$gGdown <- factor(FLXswap7groupsF$gGdown)


# Add chromosome location to the data set as well
# First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
# Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
# Now make a new vector with the chromosomes matching the FLXswap7groupsF data frame
chrdataFswap <- cbind(genesFswap, chrIDswap[match(FLXswap7groupsF$gene.id, rownames(FLXswapreadcount)), c(1,2)])

# Last add the chromosome location to the FLXswap7groupsF data frame
FLXswap7groupsF$chr <- character(15505)
for(i in 1:15505) { FLXswap7groupsF$chr[FLXswap7groupsF$gene.id == chrdataFswap[i, 1]] <- chrdataFswap[i, 2]}

# Clean up the chromosome location column
FLXswap7groupsF$chr <- gsub("\\;.*", "", as.character(FLXswap7groupsF$chr))


#Now save the full new data frame as a csv file
write.csv(FLXswap7groupsF, "FLXswapFres7G.csv", row.names = F)
