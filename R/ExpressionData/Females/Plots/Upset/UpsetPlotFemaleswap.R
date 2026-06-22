################################################################################################
############################ FLX CHROMOSOME SWAP FEMALES UPSET PLOT ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Females/Plots/Upset")

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



########## DATA FOR PLOTTING ##########


# Next I have to make a data set for plotting

# Collect transcripts
genesFswap <- rownames(FLXswapF)

# Now create the new data frame
upsetplotF <- data.frame(gene.id = genesFswap,
                         FLXsexFLXauto = numeric(16785), FLXsexCFMsex = numeric(16785), FLXsexCFMauto = numeric(16785),
                         FLXautoCFMsex = numeric(16785), FLXautoCFMauto = numeric(16785), CFMsexCFMauto = numeric(16785))


# FLXsex vs FLXauto
upsetplotF$FLXsexFLXauto <- FALSE
genesFswapFLXsexFLXauto <- rownames(FLXsexFLXaFsig)
upsetplotF$FLXsexFLXauto[which(upsetplotF$gene.id %in% genesFswapFLXsexFLXauto)] <- TRUE
upsetplotF$FLXsexFLXauto <- factor(upsetplotF$FLXsexFLXauto)

# FLXsex vs CFMsex
upsetplotF$FLXsexCFMsex <- FALSE
genesFswapFLXsexCFMsex <- rownames(FLXsexCFMsexFsig)
upsetplotF$FLXsexCFMsex[which(upsetplotF$gene.id %in% genesFswapFLXsexCFMsex)] <- TRUE
upsetplotF$FLXsexCFMsex <- factor(upsetplotF$FLXsexCFMsex)

# FLXsex vs CFMauto
upsetplotF$FLXsexCFMauto <- FALSE
genesFswapFLXsexCFMauto <- rownames(FLXsexCFMaFsig)
upsetplotF$FLXsexCFMauto[which(upsetplotF$gene.id %in% genesFswapFLXsexCFMauto)] <- TRUE
upsetplotF$FLXsexCFMauto <- factor(upsetplotF$FLXsexCFMauto)

# FLXauto vs CFMsex
upsetplotF$FLXautoCFMsex <- FALSE
genesFswapFLXautoCFMsex <- rownames(FLXaCFMsexFsig)
upsetplotF$FLXautoCFMsex[which(upsetplotF$gene.id %in% genesFswapFLXautoCFMsex)] <- TRUE
upsetplotF$FLXautoCFMsex <- factor(upsetplotF$FLXautoCFMsex)

# FLXauto vs CFMauto
upsetplotF$FLXautoCFMauto <- FALSE
genesFswapFLXautoCFMauto <- rownames(FLXaCFMaFsig)
upsetplotF$FLXautoCFMauto[which(upsetplotF$gene.id %in% genesFswapFLXautoCFMauto)] <- TRUE
upsetplotF$FLXautoCFMauto <- factor(upsetplotF$FLXautoCFMauto)

# CFMsex vs CFMauto
upsetplotF$CFMsexCFMauto <- FALSE
genesFswapCFMsexCFMauto <- rownames(CFMsexCFMaFsig)
upsetplotF$CFMsexCFMauto[which(upsetplotF$gene.id %in% genesFswapCFMsexCFMauto)] <- TRUE
upsetplotF$CFMsexCFMauto <- factor(upsetplotF$CFMsexCFMauto)

#Now save the full new data frame as a csv file
write.csv(upsetplotF, "upsetplotF.csv", row.names = F)



#########################################  PLOT DATA  #########################################


# Read in files with data
upsetplotF.data <- read.table(file = "upsetplotF.csv", h = T, sep = ",", stringsAsFactors = T)

# Next add a new column that is a list, I got this code from the internet so not completely sure what the different parts means,
# but it works and makes a list
upsetplotF_com <- upsetplotF.data |>
  mutate(
    Combination = pmap(
      list(FLXsexFLXauto == TRUE, FLXsexCFMsex == TRUE, FLXsexCFMauto == TRUE, FLXautoCFMsex == TRUE, FLXautoCFMauto == TRUE, CFMsexCFMauto == TRUE),
      \(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6) {
        c('FLXsexFLXauto', 'FLXsexCFMsex', 'FLXsexCFMauto', 'FLXautoCFMsex', 'FLXautoCFMauto', 'CFMsexCFMauto')[c(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6)] 
      }
    )
  )


# Next it's time to make the plot, I'm still working on the aesthetics
upsetplotF_com|> 
  ggplot(aes(x = Combination)) +
  geom_bar() +
  ylim(0, 2500) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size = 3, family='serif') +
  scale_x_upset(
    sets = c("FLXsexFLXauto", "FLXsexCFMsex", "FLXsexCFMauto", "FLXautoCFMsex", "FLXautoCFMauto", "CFMsexCFMauto")) +
  labs(title = "Female", tag = "A",
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


