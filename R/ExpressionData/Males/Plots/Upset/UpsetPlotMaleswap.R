################################################################################################
############################# FLX CHROMOSOME SWAP MALES UPSET PLOT #############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Males/Plots/Upset")

# Set up environment
library(DESeq2)
library(ggupset)
library(tidyverse)

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



########## DATA FOR PLOTTING ##########


# Next I have to make a data set for plotting

# Collect transcripts
genesMswap <- rownames(FLXswapM)

# Now create the new data frame
upsetplotM <- data.frame(gene.id = genesMswap,
                         FLXsexFLXauto = numeric(16785), FLXsexCFMsex = numeric(16785), FLXsexCFMauto = numeric(16785),
                         FLXautoCFMsex = numeric(16785), FLXautoCFMauto = numeric(16785), CFMsexCFMauto = numeric(16785))


# FLXsex vs FLXauto
upsetplotM$FLXsexFLXauto <- FALSE
genesFswapFLXsexFLXauto <- rownames(FLXsexFLXaMsig)
upsetplotM$FLXsexFLXauto[which(upsetplotM$gene.id %in% genesMswapFLXsexFLXauto)] <- TRUE
upsetplotM$FLXsexFLXauto <- factor(upsetplotM$FLXsexFLXauto)

# FLXsex vs CFMsex
upsetplotM$FLXsexCFMsex <- FALSE
genesMswapFLXsexCFMsex <- rownames(FLXsexCFMsexMsig)
upsetplotM$FLXsexCFMsex[which(upsetplotM$gene.id %in% genesMswapFLXsexCFMsex)] <- TRUE
upsetplotM$FLXsexCFMsex <- factor(upsetplotM$FLXsexCFMsex)

# FLXsex vs CFMauto
upsetplotM$FLXsexCFMauto <- FALSE
genesMswapFLXsexCFMauto <- rownames(FLXsexCFMaMsig)
upsetplotM$FLXsexCFMauto[which(upsetplotM$gene.id %in% genesMswapFLXsexCFMauto)] <- TRUE
upsetplotM$FLXsexCFMauto <- factor(upsetplotM$FLXsexCFMauto)

# FLXauto vs CFMsex
upsetplotM$FLXautoCFMsex <- FALSE
genesMswapFLXautoCFMsex <- rownames(FLXaCFMsexMsig)
upsetplotM$FLXautoCFMsex[which(upsetplotM$gene.id %in% genesMswapFLXautoCFMsex)] <- TRUE
upsetplotM$FLXautoCFMsex <- factor(upsetplotM$FLXautoCFMsex)

# FLXauto vs CFMauto
upsetplotM$FLXautoCFMauto <- FALSE
genesMswapFLXautoCFMauto <- rownames(FLXaCFMaMsig)
upsetplotM$FLXautoCFMauto[which(upsetplotM$gene.id %in% genesMswapFLXautoCFMauto)] <- TRUE
upsetplotM$FLXautoCFMauto <- factor(upsetplotM$FLXautoCFMauto)

# CFMsex vs CFMauto
upsetplotM$CFMsexCFMauto <- FALSE
genesMswapCFMsexCFMauto <- rownames(CFMsexCFMaMsig)
upsetplotM$CFMsexCFMauto[which(upsetplotM$gene.id %in% genesMswapCFMsexCFMauto)] <- TRUE
upsetplotM$CFMsexCFMauto <- factor(upsetplotM$CFMsexCFMauto)

#Now save the full new data frame as a csv file
write.csv(upsetplotM, "upsetplotM.csv", row.names = F)



#########################################  PLOT DATA  #########################################


# Read in files with data
upsetplotM.data <- read.table(file = "upsetplotM.csv", h = T, sep = ",", stringsAsFactors = T)

# Next add a new column that is a list, I got this code from the internet so not completely sure what the different parts means,
# but it works and makes a list
upsetplotM_com <- upsetplotM.data |>
  mutate(
    Combination = pmap(
      list(FLXsexFLXauto == TRUE, FLXsexCFMsex == TRUE, FLXsexCFMauto == TRUE, FLXautoCFMsex == TRUE, FLXautoCFMauto == TRUE, CFMsexCFMauto == TRUE),
      \(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6) {
        c('FLXsexFLXauto', 'FLXsexCFMsex', 'FLXsexCFMauto', 'FLXautoCFMsex', 'FLXautoCFMauto', 'CFMsexCFMauto')[c(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6)] 
      }
    )
  )


# Next it's time to make the plot, I'm still working on the aesthetics
upsetplotM_com|> 
  ggplot(aes(x = Combination)) +
  geom_bar() +
  ylim(0, 120) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size = 3, family='serif') +
  scale_x_upset(
    sets = c("FLXsexFLXauto", "FLXsexCFMsex", "FLXsexCFMauto", "FLXautoCFMsex", "FLXautoCFMauto", "CFMsexCFMauto")) +
  labs(title = "Male", tag = "B",
       x = "Intersections",
       y = "Number of significant transcripts") +
  theme_classic() +
  theme_combmatrix(
    combmatrix.label.text = element_text(family = 'serif', size=12)
  ) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 18),
    axis.text=element_text(size = 12),
    axis.title=element_text(size = 15))
