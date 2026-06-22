################################################################################################
#################### DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP MALES ####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/DifferentialExpression")

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

# Next we do a very simple prefiltering and remove all transcripts with a count less than 1
keep <- rowSums(counts(FLXswapsex)) >= 1
FLXswapsex <- FLXswapsex[keep,]

# Check how many transcripts are left
nrow(FLXswapsex)
# 16785

# Now are we ready to do a differential expression analysis
FLXswap <- DESeq(FLXswapsex)



######################## MALES ########################

#As we want to analyse the data separately per sex, we start by sub-setting the data

#First, separate the data into males
FLXswapM <- FLXswap[,FLXswap$sex %in% "m"]

#Next drop sex as a level, as it is not need anymore
FLXswapM$sex <- droplevels(FLXswapM$sex)

#Last thing is to update the model to not contain sex
design(FLXswapM) <- formula(~ type + type:nested.RP)

#And the run the DESeq again
FLXswapM <- DESeq(FLXswapM)



########## MALE VS MALE ##########


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaM <- results(FLXswapM, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
#Then select the transcripts that are significantly different
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




########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########


### MALE ONLY GROUPS ###

# First I make a data frame to include all transcripts that are specifically only found in males

# First collect the transcripts from the last DGElist
genesMswap <- rownames(FLXswapM)

# Now create the new data frame
FLXswapMG <- data.frame(gene.id = genesMswap, g26 = numeric(16785), g26up = numeric(16785), g26down = numeric(16785),
                        g33 = numeric(16785), g33up = numeric(16785), g33down = numeric(16785), 
                        g36 = numeric(16785), g36up = numeric(16785), g36down = numeric(16785),
                        g55 = numeric(16785), g55up = numeric(16785), g55down = numeric(16785), Overall = numeric(16785))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts from group 26 
FLXswapMG$g26 <- "NotSig"
# To be able to add sig make a list of transcripts
G26.1 <- intersect(intersect(rownames(FLXsexFLXaMsig), rownames(FLXsexCFMsexMsig)), rownames(CFMsexCFMaMsig))
genesMG26 <- setdiff(setdiff(setdiff(G26.1, rownames(FLXsexCFMaMsig)), rownames(FLXaCFMsexMsig)), rownames(FLXaCFMaMsig))
# Then add Sig for transcripts that are significant
FLXswapMG$g26[which(FLXswapMG$gene.id %in% genesMG26)] <- "Sig"
# Make this a factor
FLXswapMG$g26 <- factor(FLXswapMG$g26)
# Next add the up- and down-regulated transcripts separately
FLXswapMG$g26up <- "NotSig"
FLXswapMG$g26down <- "NotSig"
genesMG26up <- intersect(genesMG26, rownames(FLXsexFLXaMsigup))
genesMG26down <- intersect(genesMG26, rownames(FLXsexFLXaMsigdown))
FLXswapMG$g26up[which(FLXswapMG$gene.id %in% genesMG26up)] <- "Sig"
FLXswapMG$g26down[which(FLXswapMG$gene.id %in% genesMG26down)] <- "Sig"
FLXswapMG$g26up <- factor(FLXswapMG$g26up)
FLXswapMG$g26down <- factor(FLXswapMG$g26down)

# Then I add the transcripts from group 33
FLXswapMG$g33 <- "NotSig"
G33.1 <- intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig)), rownames(FLXaCFMsexMsig))
genesMG33 <- setdiff(setdiff(setdiff(G33.1, rownames(FLXsexFLXaMsig)), rownames(FLXaCFMaMsig)), rownames(CFMsexCFMaMsig))
FLXswapMG$g33[which(FLXswapMG$gene.id %in% genesMG33)] <- "Sig"
FLXswapMG$g33 <- factor(FLXswapMG$g33)
FLXswapMG$g33up <- "NotSig"
FLXswapMG$g33down <- "NotSig"
genesMG33up <- intersect(genesMG33, rownames(FLXsexCFMsexMsigup))
genesMG33down <- intersect(genesMG33, rownames(FLXsexCFMsexMsigdown))
FLXswapMG$g33up[which(FLXswapMG$gene.id %in% genesMG33up)] <- "Sig"
FLXswapMG$g33down[which(FLXswapMG$gene.id %in% genesMG33down)] <- "Sig"
FLXswapMG$g33up <- factor(FLXswapMG$g33up)
FLXswapMG$g33down <- factor(FLXswapMG$g33down)

# Group 36
FLXswapMG$g36 <- "NotSig"
G36.1 <- intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXaCFMsexMsig)), rownames(FLXaCFMaMsig))
genesMG36 <- setdiff(setdiff(setdiff(G36.1, rownames(FLXsexFLXaMsig)), rownames(FLXsexCFMaMsig)), rownames(CFMsexCFMaMsig))
FLXswapMG$g36[which(FLXswapMG$gene.id %in% genesMG36)] <- "Sig"
FLXswapMG$g36 <- factor(FLXswapMG$g36)
FLXswapMG$g36up <- "NotSig"
FLXswapMG$g36down <- "NotSig"
genesMG36up <- intersect(genesMG36, rownames(FLXsexCFMsexMsigup))
genesMG36down <- intersect(genesMG36, rownames(FLXsexCFMsexMsigdown))
FLXswapMG$g36up[which(FLXswapMG$gene.id %in% genesMG36up)] <- "Sig"
FLXswapMG$g36down[which(FLXswapMG$gene.id %in% genesMG36down)] <- "Sig"
FLXswapMG$g36up <- factor(FLXswapMG$g36up)
FLXswapMG$g36down <- factor(FLXswapMG$g36down)

# Group 55
FLXswapMG$g55 <- "NotSig"
G55.1 <- intersect(intersect(intersect(rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig)),
                             rownames(FLXaCFMaMsig)), 
                   rownames(CFMsexCFMaMsig))
genesMG55 <- setdiff(setdiff(G55.1, rownames(FLXsexFLXaMsig)), rownames(FLXaCFMsexMsig))
FLXswapMG$g55[which(FLXswapMG$gene.id %in% genesMG55)] <- "Sig"
FLXswapMG$g55 <- factor(FLXswapMG$g55)
FLXswapMG$g55up <- "NotSig"
FLXswapMG$g55down <- "NotSig"
genesMG55up <- intersect(genesMG55, rownames(FLXsexCFMsexMsigup))
genesMG55down <- intersect(genesMG55, rownames(FLXsexCFMsexMsigdown))
FLXswapMG$g55up[which(FLXswapMG$gene.id %in% genesMG55up)] <- "Sig"
FLXswapMG$g55down[which(FLXswapMG$gene.id %in% genesMG55down)] <- "Sig"
FLXswapMG$g55up <- factor(FLXswapMG$g55up)
FLXswapMG$g55down <- factor(FLXswapMG$g55down)

# Last I'll add the overall significant transcripts for all 
FLXswapMG$Overall <- "NotSig"
# To be able to add sig make a list of transcripts significant in all comparisons
genesMGoverall <- unique(cbind(c(rownames(FLXsexFLXaMsig), rownames(FLXsexCFMsexMsig), rownames(FLXsexCFMaMsig),
                                 rownames(FLXaCFMsexMsig), rownames(FLXaCFMaMsig), rownames(CFMsexCFMaMsig))))
# Then add Sig for transcripts that are significant
FLXswapMG$Overall[which(FLXswapMG$gene.id %in% genesMGoverall)] <- "Sig"
# Make this a factor
FLXswapMG$Overall <- factor(FLXswapMG$Overall)

# Add chromosome location to the data set as well
# First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
# Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
# Now make a new vector with the chromosomes matching the FLXswapMG data frame
chrdataMswap <- cbind(genesMswap, chrIDswap[match(FLXswapMG$gene.id, rownames(FLXswapreadcount)), c(1,2)])

# Last add the chromosome location to the FLXswapMG data frame
FLXswapMG$chr <- character(16785)
for(i in 1:16785) { FLXswapMG$chr[FLXswapMG$gene.id == chrdataMswap[i, 1]] <- chrdataMswap[i, 2]}

# Clean up the chromosome location column
FLXswapMG$chr <- gsub("\\;.*", "", as.character(FLXswapMG$chr))


# Now save the full new data frame as a csv file
write.csv(FLXswapMG, "FLXswapMaleGroups.csv", row.names = F)




### FIRST SEVEN FEMALE GROUPS ###

# Time to create another new data set
FLXswapResM <- data.frame(gene.id = genesMswap, gA41 = numeric(16785), gA41up = numeric(16785), gA41down = numeric(16785),
                          gB7 = numeric(16785), gB7up = numeric(16785), gB7down = numeric(16785), 
                          gC6 = numeric(16785), gC6up = numeric(16785), gC6down = numeric(16785),
                          gD22 = numeric(16785), gD22up = numeric(16785), gD22down = numeric(16785),
                          gE4 = numeric(16785), gE4up = numeric(16785), gE4down = numeric(16785),
                          gF18 = numeric(16785), gF18up = numeric(16785), gF18down = numeric(16785),
                          gG19 = numeric(16785), gG19up = numeric(16785), gG19down = numeric(16785),
                          g2 = numeric(16785), g2up = numeric(16785), g2down = numeric(16785),
                          g11 = numeric(16785), g11up = numeric(16785), g11down = numeric(16785))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts for group A/41 
FLXswapResM$gA41 <- "NotSig"
# To be able to add sig make a list of transcripts
GA41m.1 <- intersect(intersect(rownames(FLXsexCFMaMsig), rownames(FLXaCFMaMsig)), rownames(CFMsexCFMaMsig))
genesMGA41 <- setdiff(setdiff(setdiff(GA41m.1, rownames(FLXsexFLXaMsig)), rownames(FLXsexCFMsexMsig)), rownames(FLXaCFMsexMsig))
# Then add Sig for transcripts that are significant
FLXswapResM$gA41[which(FLXswapResM$gene.id %in% genesMGA41)] <- "Sig"
# Make this a factor
FLXswapResM$gA41 <- factor(FLXswapResM$gA41)
# Next add the up- and down-regulated transcripts separately
FLXswapResM$gA41up <- "NotSig"
FLXswapResM$gA41down <- "NotSig"
genesMGA41up <- intersect(genesMGA41, rownames(FLXsexCFMaMsigup))
genesMGA41down <- intersect(genesMGA41, rownames(FLXsexCFMaMsigdown))
FLXswapResM$gA41up[which(FLXswapResM$gene.id %in% genesMGA41up)] <- "Sig"
FLXswapResM$gA41down[which(FLXswapResM$gene.id %in% genesMGA41down)] <- "Sig"
FLXswapResM$gA41up <- factor(FLXswapResM$gA41up)
FLXswapResM$gA41down <- factor(FLXswapResM$gA41down)

# Then I add the transcripts from group B/7
FLXswapResM$gB7 <- "NotSig"
genesMGB7 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(CFMsexCFMaMsig), rownames(FLXsexFLXaMsig)), 
                                             rownames(FLXsexCFMsexMsig)), 
                                     rownames(FLXsexCFMaMsig)), 
                             rownames(FLXaCFMsexMsig)), 
                     rownames(FLXaCFMaMsig))
FLXswapResM$gB7[which(FLXswapResM$gene.id %in% genesMGB7)] <- "Sig"
FLXswapResM$gB7 <- factor(FLXswapResM$gB7)
FLXswapResM$gB7up <- "NotSig"
FLXswapResM$gB7down <- "NotSig"
genesMGB7up <- intersect(genesMGB7, rownames(CFMsexCFMaMsigup))
genesMGB7down <- intersect(genesMGB7, rownames(CFMsexCFMaMsigdown))
FLXswapResM$gB7up[which(FLXswapResM$gene.id %in% genesMGB7up)] <- "Sig"
FLXswapResM$gB7down[which(FLXswapResM$gene.id %in% genesMGB7down)] <- "Sig"
FLXswapResM$gB7up <- factor(FLXswapResM$gB7up)
FLXswapResM$gB7down <- factor(FLXswapResM$gB7down)

# Group C/6
FLXswapResM$gC6 <- "NotSig"
genesMGC6 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXaCFMaMsig), rownames(FLXsexFLXaMsig)), 
                                             rownames(FLXsexCFMsexMsig)), 
                                     rownames(FLXsexCFMaMsig)), 
                             rownames(FLXaCFMsexMsig)), 
                     rownames(CFMsexCFMaMsig))
FLXswapResM$gC6[which(FLXswapResM$gene.id %in% genesMGC6)] <- "Sig"
FLXswapResM$gC6 <- factor(FLXswapResM$gC6)
FLXswapResM$gC6up <- "NotSig"
FLXswapResM$gC6down <- "NotSig"
genesMGC6up <- intersect(genesMGC6, rownames(FLXaCFMaMsigup))
genesMGC6down <- intersect(genesMGC6, rownames(FLXaCFMaMsigdown))
FLXswapResM$gC6up[which(FLXswapResM$gene.id %in% genesMGC6up)] <- "Sig"
FLXswapResM$gC6down[which(FLXswapResM$gene.id %in% genesMGC6down)] <- "Sig"
FLXswapResM$gC6up <- factor(FLXswapResM$gC6up)
FLXswapResM$gC6down <- factor(FLXswapResM$gC6down)

# Group D/22
FLXswapResM$gD22 <- "NotSig"
GD22m.1 <- intersect(rownames(FLXaCFMaMsig), rownames(CFMsexCFMaMsig))
genesMGD22 <- setdiff(setdiff(setdiff(setdiff(GD22m.1, rownames(FLXsexFLXaMsig)),
                                      rownames(FLXsexCFMsexMsig)),
                              rownames(FLXsexCFMaMsig)),
                      rownames(FLXaCFMsexMsig))
FLXswapResM$gD22[which(FLXswapResM$gene.id %in% genesMGD22)] <- "Sig"
FLXswapResM$gD22 <- factor(FLXswapResM$gD22)
FLXswapResM$gD22up <- "NotSig"
FLXswapResM$gD22down <- "NotSig"
genesMGD22up <- intersect(genesMGD22, rownames(FLXaCFMaMsigup))
genesMGD22down <- intersect(genesMGD22, rownames(FLXaCFMaMsigdown))
FLXswapResM$gD22up[which(FLXswapResM$gene.id %in% genesMGD22up)] <- "Sig"
FLXswapResM$gD22down[which(FLXswapResM$gene.id %in% genesMGD22down)] <- "Sig"
FLXswapResM$gD22up <- factor(FLXswapResM$gD22up)
FLXswapResM$gD22down <- factor(FLXswapResM$gD22down)

# Group E/4
FLXswapResM$gE4 <- "NotSig"
genesMGE4 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXsexCFMaMsig), rownames(FLXsexFLXaMsig)),
                                             rownames(FLXsexCFMsexMsig)),
                                     rownames(FLXaCFMsexMsig)),
                             rownames(FLXaCFMaMsig)),
                     rownames(CFMsexCFMaMsig))
FLXswapResM$gE4[which(FLXswapResM$gene.id %in% genesMGE4)] <- "Sig"
FLXswapResM$gE4 <- factor(FLXswapResM$gE4)
FLXswapResM$gE4up <- "NotSig"
FLXswapResM$gE4down <- "NotSig"
genesMGE4up <- intersect(genesMGE4, rownames(FLXsexCFMaMsigup))
genesMGE4down <- intersect(genesMGE4, rownames(FLXsexCFMaMsigdown))
FLXswapResM$gE4up[which(FLXswapResM$gene.id %in% genesMGE4up)] <- "Sig"
FLXswapResM$gE4down[which(FLXswapResM$gene.id %in% genesMGE4down)] <- "Sig"
FLXswapResM$gE4up <- factor(FLXswapResM$gE4up)
FLXswapResM$gE4down <- factor(FLXswapResM$gE4down)

# Group F/18
FLXswapResM$gF18 <- "NotSig"
GF18m.1 <- intersect(rownames(FLXsexCFMaMsig), rownames(FLXaCFMaMsig))
genesMGF18 <- setdiff(setdiff(setdiff(setdiff(GF18m.1, rownames(FLXsexFLXaMsig)),
                                      rownames(FLXsexCFMsexMsig)),
                              rownames(FLXaCFMsexMsig)),
                      rownames(CFMsexCFMaMsig))
FLXswapResM$gF18[which(FLXswapResM$gene.id %in% genesMGF18)] <- "Sig"
FLXswapResM$gF18 <- factor(FLXswapResM$gF18)
FLXswapResM$gF18up <- "NotSig"
FLXswapResM$gF18down <- "NotSig"
genesMGF18up <- intersect(genesMGF18, rownames(FLXsexCFMaMsigup))
genesMGF18down <- intersect(genesMGF18, rownames(FLXsexCFMaMsigdown))
FLXswapResM$gF18up[which(FLXswapResM$gene.id %in% genesMGF18up)] <- "Sig"
FLXswapResM$gF18down[which(FLXswapResM$gene.id %in% genesMGF18down)] <- "Sig"
FLXswapResM$gF18up <- factor(FLXswapResM$gF18up)
FLXswapResM$gF18down <- factor(FLXswapResM$gF18down)

# Group G/19
FLXswapResM$gG19 <- "NotSig"
GG19m.1 <- intersect(rownames(FLXsexCFMaMsig), rownames(CFMsexCFMaMsig))
genesMGG19 <- setdiff(setdiff(setdiff(setdiff(GG19m.1, rownames(FLXsexFLXaMsig)),
                                      rownames(FLXsexCFMsexMsig)),
                              rownames(FLXaCFMsexMsig)),
                      rownames(FLXaCFMaMsig))
FLXswapResM$gG19[which(FLXswapResM$gene.id %in% genesMGG19)] <- "Sig"
FLXswapResM$gG19 <- factor(FLXswapResM$gG19)
FLXswapResM$gG19up <- "NotSig"
FLXswapResM$gG19down <- "NotSig"
genesMGG19up <- intersect(genesMGG19, rownames(FLXsexCFMaMsigup))
genesMGG19down <- intersect(genesMGG19, rownames(FLXsexCFMaMsigdown))
FLXswapResM$gG19up[which(FLXswapResM$gene.id %in% genesMGG19up)] <- "Sig"
FLXswapResM$gG19down[which(FLXswapResM$gene.id %in% genesMGG19down)] <- "Sig"
FLXswapResM$gG19up <- factor(FLXswapResM$gG19up)
FLXswapResM$gG19down <- factor(FLXswapResM$gG19down)

# Group 2
FLXswapResM$g2 <- "NotSig"
genesMG2 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXsexFLXaMsig), rownames(FLXsexCFMsexMsig)),
                                             rownames(FLXsexCFMaMsig)),
                                     rownames(FLXaCFMsexMsig)),
                             rownames(FLXaCFMaMsig)),
                     rownames(CFMsexCFMaMsig))
FLXswapResM$g2[which(FLXswapResM$gene.id %in% genesMG2)] <- "Sig"
FLXswapResM$g2 <- factor(FLXswapResM$g2)
FLXswapResM$g2up <- "NotSig"
FLXswapResM$g2down <- "NotSig"
genesMG2up <- intersect(genesMG2, rownames(FLXsexFLXaMsigup))
genesMG2down <- intersect(genesMG2, rownames(FLXsexFLXaMsigdown))
FLXswapResM$g2up[which(FLXswapResM$gene.id %in% genesMG2up)] <- "Sig"
FLXswapResM$g2down[which(FLXswapResM$gene.id %in% genesMG2down)] <- "Sig"
FLXswapResM$g2up <- factor(FLXswapResM$g2up)
FLXswapResM$g2down <- factor(FLXswapResM$g2down)

# Group 11
FLXswapResM$g11 <- "NotSig"
G11m.1 <- intersect(rownames(FLXsexFLXaMsig), rownames(FLXaCFMaMsig))
genesMG11 <- setdiff(setdiff(setdiff(setdiff(G11m.1, rownames(FLXsexCFMsexMsig)),
                                     rownames(FLXsexCFMaMsig)),
                             rownames(FLXaCFMsexMsig)),
                     rownames(CFMsexCFMaMsig))
FLXswapResM$g11[which(FLXswapResM$gene.id %in% genesMG11)] <- "Sig"
FLXswapResM$g11 <- factor(FLXswapResM$g11)
FLXswapResM$g11up <- "NotSig"
FLXswapResM$g11down <- "NotSig"
genesMG11up <- intersect(genesMG11, rownames(FLXsexFLXaMsigup))
genesMG11down <- intersect(genesMG11, rownames(FLXsexFLXaMsigdown))
FLXswapResM$g11up[which(FLXswapResM$gene.id %in% genesMG11up)] <- "Sig"
FLXswapResM$g11down[which(FLXswapResM$gene.id %in% genesMG11down)] <- "Sig"
FLXswapResM$g11up <- factor(FLXswapResM$g11up)
FLXswapResM$g11down <- factor(FLXswapResM$g11down)


# Last add the chromosome location to the FLXswapResM data frame
FLXswapResM$chr <- character(16785)
for(i in 1:16785) { FLXswapResM$chr[FLXswapResM$gene.id == chrdataMswap[i, 1]] <- chrdataMswap[i, 2]}

# Clean up the chromosome location column
FLXswapResM$chr <- gsub("\\;.*", "", as.character(FLXswapResM$chr))


# Now save the full new data frame as a csv file
write.csv(FLXswapResM, "FLXswapMres.csv", row.names = F)

