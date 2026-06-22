################################################################################################
################### DIFFERENTIAL GENE EXPRESSION FLX CHROMOSOME SWAP FEMALES ###################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/DifferentialExpression")

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


######################## FEMALES ########################

#As we want to analyse the data separately per sex, we start by sub-setting the data

#First, separate the data into males
FLXswapF <- FLXswap[,FLXswap$sex %in% "f"]

#Next drop sex as a level, as it is not need anymore
FLXswapF$sex <- droplevels(FLXswapF$sex)

#Last thing is to update the model to not contain sex
design(FLXswapF) <- formula(~ type + type:nested.RP)

#And the run the DESeq again
FLXswapF <- DESeq(FLXswapF)


########## FEMALE VS FEMALE ##########


# FLX-sex chromosomes and FLX-autosomes
FLXsexFLXaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "bCwtsFLXa"))
# Then select the transcripts that are significantly different
FLXsexFLXaFsig <- subset(FLXsexFLXaF, padj < 0.05)
summary(FLXsexFLXaFsig)
# There are 287 transcripts that are different, 199 up-regulated and 88 down-regulated
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



########## MAKE THE DATA USEFULL FOR FURTHER ANALYSIS ##########


### FEMALE ONLY GROUPS ###

# First I make a data frame to include all transcripts that are specifically only found in females

# First collect the transcripts from the last DGElist
genesFswap <- rownames(FLXswapF)

# Now create the new data frame
FLXswapFG <- data.frame(gene.id = genesFswap, g35 = numeric(16785), g35up = numeric(16785), g35down = numeric(16785),
                        g39 = numeric(16785), g39up = numeric(16785), g39down = numeric(16785), g40 = numeric(16785), g40up = numeric(16785), g40down = numeric(16785),
                        g42 = numeric(16785), g42up = numeric(16785), g42down = numeric(16785), g44 = numeric(16785), g44up = numeric(16785), g44down = numeric(16785),
                        g49 = numeric(16785), g49up = numeric(16785), g49down = numeric(16785), g56 = numeric(16785), g56up = numeric(16785), g56down = numeric(16785),
                        g57 = numeric(16785), g57up = numeric(16785), g57down = numeric(16785), g59 = numeric(16785), g59up = numeric(16785), g59down = numeric(16785),
                        g62 = numeric(16785), g62up = numeric(16785), g62down = numeric(16785), g63 = numeric(16785), g63up = numeric(16785), g63down = numeric(16785),
                        Overall = numeric(16785))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts from group 35 
FLXswapFG$g35<- "NotSig"
# To be able to add sig make a list of transcripts
G35.1 <- intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFG35 <- setdiff(setdiff(setdiff(G35.1, rownames(FLXsexFLXaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))
# Then add Sig for transcripts that are significant
FLXswapFG$g35[which(FLXswapFG$gene.id %in% genesFG35)] <- "Sig"
# Make this a factor
FLXswapFG$g35 <- factor(FLXswapFG$g35)
# Next add the up- and down-regulated transcripts separately
FLXswapFG$g35up <- "NotSig"
FLXswapFG$g35down <- "NotSig"
genesFG35up <- intersect(genesFG35, rownames(FLXsexCFMsexFsigup))
genesFG35down <- intersect(genesFG35, rownames(FLXsexCFMsexFsigdown))
FLXswapFG$g35up[which(FLXswapFG$gene.id %in% genesFG35up)] <- "Sig"
FLXswapFG$g35down[which(FLXswapFG$gene.id %in% genesFG35down)] <- "Sig"
FLXswapFG$g35up <- factor(FLXswapFG$g35up)
FLXswapFG$g35down <- factor(FLXswapFG$g35down)

# Then I add the transcripts from group 39
FLXswapFG$g39 <- "NotSig"
G39.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))
genesFG39 <- setdiff(setdiff(setdiff(G39.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(CFMsexCFMaFsig))
FLXswapFG$g39[which(FLXswapFG$gene.id %in% genesFG39)] <- "Sig"
FLXswapFG$g39 <- factor(FLXswapFG$g39)
FLXswapFG$g39up <- "NotSig"
FLXswapFG$g39down <- "NotSig"
genesFG39up <- intersect(genesFG39, rownames(FLXsexCFMaFsigup))
genesFG39down <- intersect(genesFG39, rownames(FLXsexCFMaFsigdown))
FLXswapFG$g39up[which(FLXswapFG$gene.id %in% genesFG39up)] <- "Sig"
FLXswapFG$g39down[which(FLXswapFG$gene.id %in% genesFG39down)] <- "Sig"
FLXswapFG$g39up <- factor(FLXswapFG$g39up)
FLXswapFG$g39down <- factor(FLXswapFG$g39down)

# Group 40
FLXswapFG$g40 <- "NotSig"
G40.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(CFMsexCFMaFsig))
genesFG40 <- setdiff(setdiff(setdiff(G40.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXaCFMaFsig))
FLXswapFG$g40[which(FLXswapFG$gene.id %in% genesFG40)] <- "Sig"
FLXswapFG$g40 <- factor(FLXswapFG$g40)
FLXswapFG$g40up <- "NotSig"
FLXswapFG$g40down <- "NotSig"
genesFG40up <- intersect(genesFG40, rownames(FLXsexCFMaFsigup))
genesFG40down <- intersect(genesFG40, rownames(FLXsexCFMaFsigdown))
FLXswapFG$g40up[which(FLXswapFG$gene.id %in% genesFG40up)] <- "Sig"
FLXswapFG$g40down[which(FLXswapFG$gene.id %in% genesFG40down)] <- "Sig"
FLXswapFG$g40up <- factor(FLXswapFG$g40up)
FLXswapFG$g40down <- factor(FLXswapFG$g40down)

# Group 42
FLXswapFG$g42 <- "NotSig"
G42.1 <- intersect(intersect(rownames(FLXaCFMsexFsig), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFG42 <- setdiff(setdiff(setdiff(G42.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig))
FLXswapFG$g42[which(FLXswapFG$gene.id %in% genesFG42)] <- "Sig"
FLXswapFG$g42 <- factor(FLXswapFG$g42)
FLXswapFG$g42up <- "NotSig"
FLXswapFG$g42down <- "NotSig"
genesFG42up <- intersect(genesFG42, rownames(FLXaCFMsexFsigup))
genesFG42down <- intersect(genesFG42, rownames(FLXaCFMsexFsigdown))
FLXswapFG$g42up[which(FLXswapFG$gene.id %in% genesFG42up)] <- "Sig"
FLXswapFG$g42down[which(FLXswapFG$gene.id %in% genesFG42down)] <- "Sig"
FLXswapFG$g42up <- factor(FLXswapFG$g42up)
FLXswapFG$g42down <- factor(FLXswapFG$g42down)

# Group 44
FLXswapFG$g44 <- "NotSig"
G44.1 <- intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMaFsig))
genesFG44 <- setdiff(setdiff(G44.1, rownames(FLXaCFMsexFsig)), rownames(CFMsexCFMaFsig))
FLXswapFG$g44[which(FLXswapFG$gene.id %in% genesFG44)] <- "Sig"
FLXswapFG$g44 <- factor(FLXswapFG$g44)
FLXswapFG$g44up <- "NotSig"
FLXswapFG$g44down <- "NotSig"
genesFG44up <- intersect(genesFG44, rownames(FLXsexFLXaFsigup))
genesFG44down <- intersect(genesFG44, rownames(FLXsexFLXaFsigdown))
FLXswapFG$g44up[which(FLXswapFG$gene.id %in% genesFG44up)] <- "Sig"
FLXswapFG$g44down[which(FLXswapFG$gene.id %in% genesFG44down)] <- "Sig"
FLXswapFG$g44up <- factor(FLXswapFG$g44up)
FLXswapFG$g44down <- factor(FLXswapFG$g44down)

# Group 49
FLXswapFG$g49 <- "NotSig"
G49.1 <- intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig))
genesFG49 <- setdiff(setdiff(G49.1, rownames(FLXsexCFMsexFsig)), rownames(CFMsexCFMaFsig))
FLXswapFG$g49[which(FLXswapFG$gene.id %in% genesFG49)] <- "Sig"
FLXswapFG$g49 <- factor(FLXswapFG$g49)
FLXswapFG$g49up <- "NotSig"
FLXswapFG$g49down <- "NotSig"
genesFG49up <- intersect(genesFG49, rownames(FLXsexFLXaFsigup))
genesFG49down <- intersect(genesFG49, rownames(FLXsexFLXaFsigdown))
FLXswapFG$g49up[which(FLXswapFG$gene.id %in% genesFG49up)] <- "Sig"
FLXswapFG$g49down[which(FLXswapFG$gene.id %in% genesFG49down)] <- "Sig"
FLXswapFG$g49up <- factor(FLXswapFG$g49up)
FLXswapFG$g49down <- factor(FLXswapFG$g49down)

# Group 56
FLXswapFG$g56 <- "NotSig"
G56.1 <- intersect(intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFG56 <- setdiff(setdiff(G56.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMaFsig))
FLXswapFG$g56[which(FLXswapFG$gene.id %in% genesFG56)] <- "Sig"
FLXswapFG$g56 <- factor(FLXswapFG$g56)
FLXswapFG$g56up <- "NotSig"
FLXswapFG$g56down <- "NotSig"
genesFG56up <- intersect(genesFG56, rownames(FLXsexCFMsexFsigup))
genesFG56down <- intersect(genesFG56, rownames(FLXsexCFMsexFsigdown))
FLXswapFG$g56up[which(FLXswapFG$gene.id %in% genesFG56up)] <- "Sig"
FLXswapFG$g56down[which(FLXswapFG$gene.id %in% genesFG56down)] <- "Sig"
FLXswapFG$g56up <- factor(FLXswapFG$g56up)
FLXswapFG$g56down <- factor(FLXswapFG$g56down)

# Group 57
FLXswapFG$g57 <- "NotSig"
G57.1 <- intersect(intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFG57 <- setdiff(setdiff(G57.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig))
FLXswapFG$g57[which(FLXswapFG$gene.id %in% genesFG57)] <- "Sig"
FLXswapFG$g57 <- factor(FLXswapFG$g57)
FLXswapFG$g57up <- "NotSig"
FLXswapFG$g57down <- "NotSig"
genesFG57up <- intersect(genesFG57, rownames(FLXsexCFMaFsigup))
genesFG57down <- intersect(genesFG57, rownames(FLXsexCFMaFsigdown))
FLXswapFG$g57up[which(FLXswapFG$gene.id %in% genesFG57up)] <- "Sig"
FLXswapFG$g57down[which(FLXswapFG$gene.id %in% genesFG57down)] <- "Sig"
FLXswapFG$g57up <- factor(FLXswapFG$g57up)
FLXswapFG$g57down <- factor(FLXswapFG$g57down)

# Group 59
FLXswapFG$g59 <- "NotSig"
G59.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), 
                   rownames(CFMsexCFMaFsig))
genesFG59 <- setdiff(G59.1, rownames(FLXaCFMaFsig))
FLXswapFG$g59[which(FLXswapFG$gene.id %in% genesFG59)] <- "Sig"
FLXswapFG$g59 <- factor(FLXswapFG$g59)
FLXswapFG$g59up <- "NotSig"
FLXswapFG$g59down <- "NotSig"
genesFG59up <- intersect(genesFG59, rownames(FLXsexFLXaFsigup))
genesFG59down <- intersect(genesFG59, rownames(FLXsexFLXaFsigdown))
FLXswapFG$g59up[which(FLXswapFG$gene.id %in% genesFG59up)] <- "Sig"
FLXswapFG$g59down[which(FLXswapFG$gene.id %in% genesFG59down)] <- "Sig"
FLXswapFG$g59up <- factor(FLXswapFG$g59up)
FLXswapFG$g59down <- factor(FLXswapFG$g59down)

# Group 62
FLXswapFG$g62 <- "NotSig"
G62.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), 
                   rownames(CFMsexCFMaFsig))
genesFG62 <- setdiff(G62.1, rownames(FLXsexCFMsexFsig))
FLXswapFG$g62[which(FLXswapFG$gene.id %in% genesFG62)] <- "Sig"
FLXswapFG$g62 <- factor(FLXswapFG$g62)
FLXswapFG$g62up <- "NotSig"
FLXswapFG$g62down <- "NotSig"
genesFG62up <- intersect(genesFG62, rownames(FLXsexFLXaFsigup))
genesFG62down <- intersect(genesFG62, rownames(FLXsexFLXaFsigdown))
FLXswapFG$g62up[which(FLXswapFG$gene.id %in% genesFG62up)] <- "Sig"
FLXswapFG$g62down[which(FLXswapFG$gene.id %in% genesFG62down)] <- "Sig"
FLXswapFG$g62up <- factor(FLXswapFG$g62up)
FLXswapFG$g62down <- factor(FLXswapFG$g62down)

# Group 63
FLXswapFG$g63 <- "NotSig"
G63.1 <- intersect(intersect(intersect(intersect(rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig)), rownames(FLXaCFMsexFsig)), rownames(FLXaCFMaFsig)), 
                   rownames(CFMsexCFMaFsig))
genesFG63 <- setdiff(G63.1, rownames(FLXsexFLXaFsig))
FLXswapFG$g63[which(FLXswapFG$gene.id %in% genesFG63)] <- "Sig"
FLXswapFG$g63 <- factor(FLXswapFG$g63)
FLXswapFG$g63up <- "NotSig"
FLXswapFG$g63down <- "NotSig"
genesFG63up <- intersect(genesFG63, rownames(FLXsexCFMsexFsigup))
genesFG63down <- intersect(genesFG63, rownames(FLXsexCFMsexFsigdown))
FLXswapFG$g63up[which(FLXswapFG$gene.id %in% genesFG63up)] <- "Sig"
FLXswapFG$g63down[which(FLXswapFG$gene.id %in% genesFG63down)] <- "Sig"
FLXswapFG$g63up <- factor(FLXswapFG$g63up)
FLXswapFG$g63down <- factor(FLXswapFG$g63down)

# Last I'll add the overall significant transcripts for all three classes
FLXswapFG$Overall <- "NotSig"
# To be able to add sig make a list of transcripts
genesFGoverall <- unique(cbind(c(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig), rownames(FLXsexCFMaFsig),
                                 rownames(FLXaCFMsexFsig), rownames(FLXaCFMaFsig), rownames(CFMsexCFMaFsig))))
# Then add Sig for transcripts that are significant
FLXswapFG$Overall[which(FLXswapFG$gene.id %in% genesFGoverall)] <- "Sig"
# Make this a factor
FLXswapFG$Overall <- factor(FLXswapFG$Overall)

#Add chromosome location to the data set as well
#First add the gene names as a new column
FLXswapreadcountUnfilt$id <- rownames(FLXswapreadcountUnfilt)
#Then create a new vector with both chromosome and gene ID 
chrIDswap <- FLXswapreadcountUnfilt[,c(1,98)]
chrIDswap$Chr <- as.character(chrIDswap$Chr)
#Now make a new vector with the chromosomes matching the FLXswapFG data frame
chrdataFswap <- cbind(genesFswap, chrIDswap[match(FLXswapFG$gene.id, rownames(FLXswapreadcount)), c(1,2)])

#Last add the chromosome location to the FLXswapFG data frame
FLXswapFG$chr <- character(16785)
for(i in 1:16785) { FLXswapFG$chr[FLXswapFG$gene.id == chrdataFswap[i, 1]] <- chrdataFswap[i, 2]}

#Clean up the chromosome location column
FLXswapFG$chr <- gsub("\\;.*", "", as.character(FLXswapFG$chr))


#Now save the full new data frame as a csv file
write.csv(FLXswapFG, "FLXswapFemaleGroups.csv", row.names = F)




### FIRST SEVEN FEMALE GROUPS ###

# Time to create another new data set
FLXswapResF <- data.frame(gene.id = genesFswap, gA41 = numeric(16785), gA41up = numeric(16785), gA41down = numeric(16785),
                          gB7 = numeric(16785), gB7up = numeric(16785), gB7down = numeric(16785), 
                          gC6 = numeric(16785), gC6up = numeric(16785), gC6down = numeric(16785),
                          gD22 = numeric(16785), gD22up = numeric(16785), gD22down = numeric(16785),
                          gE4 = numeric(16785), gE4up = numeric(16785), gE4down = numeric(16785),
                          gF18 = numeric(16785), gF18up = numeric(16785), gF18down = numeric(16785),
                          gG19 = numeric(16785), gG19up = numeric(16785), gG19down = numeric(16785),
                          g23 = numeric(16785), g23up = numeric(16785), g23down = numeric(16785))

# Next it's time to add the significant transcripts from the different expressed model
# First I add the transcripts for group A/41 
FLXswapResF$gA41 <- "NotSig"
# To be able to add sig make a list of transcripts
GA41f.1 <- intersect(intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMaFsig)), rownames(CFMsexCFMaFsig))
genesFGA41 <- setdiff(setdiff(setdiff(GA41f.1, rownames(FLXsexFLXaFsig)), rownames(FLXsexCFMsexFsig)), rownames(FLXaCFMsexFsig))
# Then add Sig for transcripts that are significant
FLXswapResF$gA41[which(FLXswapResF$gene.id %in% genesFGA41)] <- "Sig"
# Make this a factor
FLXswapResF$gA41 <- factor(FLXswapResF$gA41)
# Next add the up- and down-regulated transcripts separately
FLXswapResF$gA41up <- "NotSig"
FLXswapResF$gA41down <- "NotSig"
genesFGA41up <- intersect(genesFGA41, rownames(FLXsexCFMaFsigup))
genesFGA41down <- intersect(genesFGA41, rownames(FLXsexCFMaFsigdown))
FLXswapResF$gA41up[which(FLXswapResF$gene.id %in% genesFGA41up)] <- "Sig"
FLXswapResF$gA41down[which(FLXswapResF$gene.id %in% genesFGA41down)] <- "Sig"
FLXswapResF$gA41up <- factor(FLXswapResF$gA41up)
FLXswapResF$gA41down <- factor(FLXswapResF$gA41down)

# Then I add the transcripts from group B/7
FLXswapResF$gB7 <- "NotSig"
genesFGB7 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(CFMsexCFMaFsig), rownames(FLXsexFLXaFsig)), 
                                             rownames(FLXsexCFMsexFsig)), 
                                     rownames(FLXsexCFMaFsig)), 
                             rownames(FLXaCFMsexFsig)), 
                     rownames(FLXaCFMaFsig))
FLXswapResF$gB7[which(FLXswapResF$gene.id %in% genesFGB7)] <- "Sig"
FLXswapResF$gB7 <- factor(FLXswapResF$gB7)
FLXswapResF$gB7up <- "NotSig"
FLXswapResF$gB7down <- "NotSig"
genesFGB7up <- intersect(genesFGB7, rownames(CFMsexCFMaFsigup))
genesFGB7down <- intersect(genesFGB7, rownames(CFMsexCFMaFsigdown))
FLXswapResF$gB7up[which(FLXswapResF$gene.id %in% genesFGB7up)] <- "Sig"
FLXswapResF$gB7down[which(FLXswapResF$gene.id %in% genesFGB7down)] <- "Sig"
FLXswapResF$gB7up <- factor(FLXswapResF$gB7up)
FLXswapResF$gB7down <- factor(FLXswapResF$gB7down)

# Group C/6
FLXswapResF$gC6 <- "NotSig"
genesFGC6 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXaCFMaFsig), rownames(FLXsexFLXaFsig)), 
                                             rownames(FLXsexCFMsexFsig)), 
                                     rownames(FLXsexCFMaFsig)), 
                             rownames(FLXaCFMsexFsig)), 
                     rownames(CFMsexCFMaFsig))
FLXswapResF$gC6[which(FLXswapResF$gene.id %in% genesFGC6)] <- "Sig"
FLXswapResF$gC6 <- factor(FLXswapResF$gC6)
FLXswapResF$gC6up <- "NotSig"
FLXswapResF$gC6down <- "NotSig"
genesFGC6up <- intersect(genesFGC6, rownames(FLXaCFMaFsigup))
genesFGC6down <- intersect(genesFGC6, rownames(FLXaCFMaFsigdown))
FLXswapResF$gC6up[which(FLXswapResF$gene.id %in% genesFGC6up)] <- "Sig"
FLXswapResF$gC6down[which(FLXswapResF$gene.id %in% genesFGC6down)] <- "Sig"
FLXswapResF$gC6up <- factor(FLXswapResF$gC6up)
FLXswapResF$gC6down <- factor(FLXswapResF$gC6down)

# Group D/22
FLXswapResF$gD22 <- "NotSig"
GD22f.1 <- intersect(rownames(FLXaCFMaFsig), rownames(CFMsexCFMaFsig))
genesFGD22 <- setdiff(setdiff(setdiff(setdiff(GD22f.1, rownames(FLXsexFLXaFsig)),
                                      rownames(FLXsexCFMsexFsig)),
                              rownames(FLXsexCFMaFsig)),
                      rownames(FLXaCFMsexFsig))
FLXswapResF$gD22[which(FLXswapResF$gene.id %in% genesFGD22)] <- "Sig"
FLXswapResF$gD22 <- factor(FLXswapResF$gD22)
FLXswapResF$gD22up <- "NotSig"
FLXswapResF$gD22down <- "NotSig"
genesFGD22up <- intersect(genesFGD22, rownames(FLXaCFMaFsigup))
genesFGD22down <- intersect(genesFGD22, rownames(FLXaCFMaFsigdown))
FLXswapResF$gD22up[which(FLXswapResF$gene.id %in% genesFGD22up)] <- "Sig"
FLXswapResF$gD22down[which(FLXswapResF$gene.id %in% genesFGD22down)] <- "Sig"
FLXswapResF$gD22up <- factor(FLXswapResF$gD22up)
FLXswapResF$gD22down <- factor(FLXswapResF$gD22down)

# Group E/4
FLXswapResF$gE4 <- "NotSig"
genesFGE4 <- setdiff(setdiff(setdiff(setdiff(setdiff(rownames(FLXsexCFMaFsig), rownames(FLXsexFLXaFsig)),
                                             rownames(FLXsexCFMsexFsig)),
                                     rownames(FLXaCFMsexFsig)),
                             rownames(FLXaCFMaFsig)),
                     rownames(CFMsexCFMaFsig))
FLXswapResF$gE4[which(FLXswapResF$gene.id %in% genesFGE4)] <- "Sig"
FLXswapResF$gE4 <- factor(FLXswapResF$gE4)
FLXswapResF$gE4up <- "NotSig"
FLXswapResF$gE4down <- "NotSig"
genesFGE4up <- intersect(genesFGE4, rownames(FLXsexCFMaFsigup))
genesFGE4down <- intersect(genesFGE4, rownames(FLXsexCFMaFsigdown))
FLXswapResF$gE4up[which(FLXswapResF$gene.id %in% genesFGE4up)] <- "Sig"
FLXswapResF$gE4down[which(FLXswapResF$gene.id %in% genesFGE4down)] <- "Sig"
FLXswapResF$gE4up <- factor(FLXswapResF$gE4up)
FLXswapResF$gE4down <- factor(FLXswapResF$gE4down)

# Group F/18
FLXswapResF$gF18 <- "NotSig"
GF18f.1 <- intersect(rownames(FLXsexCFMaFsig), rownames(FLXaCFMaFsig))
genesFGF18 <- setdiff(setdiff(setdiff(setdiff(GF18f.1, rownames(FLXsexFLXaFsig)),
                                      rownames(FLXsexCFMsexFsig)),
                              rownames(FLXaCFMsexFsig)),
                      rownames(CFMsexCFMaFsig))
FLXswapResF$gF18[which(FLXswapResF$gene.id %in% genesFGF18)] <- "Sig"
FLXswapResF$gF18 <- factor(FLXswapResF$gF18)
FLXswapResF$gF18up <- "NotSig"
FLXswapResF$gF18down <- "NotSig"
genesFGF18up <- intersect(genesFGF18, rownames(FLXsexCFMaFsigup))
genesFGF18down <- intersect(genesFGF18, rownames(FLXsexCFMaFsigdown))
FLXswapResF$gF18up[which(FLXswapResF$gene.id %in% genesFGF18up)] <- "Sig"
FLXswapResF$gF18down[which(FLXswapResF$gene.id %in% genesFGF18down)] <- "Sig"
FLXswapResF$gF18up <- factor(FLXswapResF$gF18up)
FLXswapResF$gF18down <- factor(FLXswapResF$gF18down)

# Group G/19
FLXswapResF$gG19 <- "NotSig"
GG19f.1 <- intersect(rownames(FLXsexCFMaFsig), rownames(CFMsexCFMaFsig))
genesFGG19 <- setdiff(setdiff(setdiff(setdiff(GG19f.1, rownames(FLXsexFLXaFsig)),
                                      rownames(FLXsexCFMsexFsig)),
                              rownames(FLXaCFMsexFsig)),
                      rownames(FLXaCFMaFsig))
FLXswapResF$gG19[which(FLXswapResF$gene.id %in% genesFGG19)] <- "Sig"
FLXswapResF$gG19 <- factor(FLXswapResF$gG19)
FLXswapResF$gG19up <- "NotSig"
FLXswapResF$gG19down <- "NotSig"
genesFGG19up <- intersect(genesFGG19, rownames(FLXsexCFMaFsigup))
genesFGG19down <- intersect(genesFGG19, rownames(FLXsexCFMaFsigdown))
FLXswapResF$gG19up[which(FLXswapResF$gene.id %in% genesFGG19up)] <- "Sig"
FLXswapResF$gG19down[which(FLXswapResF$gene.id %in% genesFGG19down)] <- "Sig"
FLXswapResF$gG19up <- factor(FLXswapResF$gG19up)
FLXswapResF$gG19down <- factor(FLXswapResF$gG19down)

# Group 23
FLXswapResF$g23 <- "NotSig"
G23f.1 <- intersect(intersect(rownames(FLXsexFLXaFsig), rownames(FLXsexCFMsexFsig)), rownames(FLXsexCFMaFsig))
genesFG23 <- setdiff(setdiff(setdiff(G23f.1, rownames(FLXaCFMsexFsig)),
                                      rownames(FLXaCFMaFsig)),
                              rownames(CFMsexCFMaFsig))
FLXswapResF$g23[which(FLXswapResF$gene.id %in% genesFG23)] <- "Sig"
FLXswapResF$g23 <- factor(FLXswapResF$g23)
FLXswapResF$g23up <- "NotSig"
FLXswapResF$g23down <- "NotSig"
genesFG23up <- intersect(genesFG23, rownames(FLXsexFLXaFsigup))
genesFG23down <- intersect(genesFG23, rownames(FLXsexFLXaFsigdown))
FLXswapResF$g23up[which(FLXswapResF$gene.id %in% genesFG23up)] <- "Sig"
FLXswapResF$g23down[which(FLXswapResF$gene.id %in% genesFG23down)] <- "Sig"
FLXswapResF$g23up <- factor(FLXswapResF$g23up)
FLXswapResF$g23down <- factor(FLXswapResF$g23down)


#Last add the chromosome location to the FLXswapResF data frame
FLXswapResF$chr <- character(16785)
for(i in 1:16785) { FLXswapResF$chr[FLXswapResF$gene.id == chrdataFswap[i, 1]] <- chrdataFswap[i, 2]}

#Clean up the chromosome location column
FLXswapResF$chr <- gsub("\\;.*", "", as.character(FLXswapResF$chr))


#Now save the full new data frame as a csv file
write.csv(FLXswapResF, "FLXswapFres.csv", row.names = F)

