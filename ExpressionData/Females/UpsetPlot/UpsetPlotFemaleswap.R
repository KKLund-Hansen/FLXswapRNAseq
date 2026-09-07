################################################################################################
############################ FLX CHROMOSOME SWAP FEMALES UPSET PLOT ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/UpsetPlot")

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

#As we want to analyse the data separately per sex, we start by sub-setting the data

#First, separate the data into females
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
# There are 277 transcripts that are different, 193 up-regulated and 84 down-regulated


# FLX-sex chromosomes and CFM-sex chromosomes
FLXsexCFMsexF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "cCFMsCwta"))
# Then select the transcripts that are significantly different
FLXsexCFMsexFsig <- subset(FLXsexCFMsexF, padj < 0.05)
summary(FLXsexCFMsexFsig)
# There are 192 transcripts that are different, 137 up-regulated and 55 down-regulated


# FLX-sex chromosomes and CFM-autosomes
FLXsexCFMaF <- results(FLXswapF, contrast = c("type", "aFLXsCwta", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXsexCFMaFsig <- subset(FLXsexCFMaF, padj < 0.05)
summary(FLXsexCFMaFsig)
# There are 3293 transcripts that are different, 508 up-regulated and 2785 down-regulated


# FLX-autosomes and CFM-sex chromosomes
FLXaCFMsexF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "cCFMsCwta"))
#Then select the transcripts that are significantly different
FLXaCFMsexFsig <- subset(FLXaCFMsexF, padj < 0.05)
summary(FLXaCFMsexFsig)
# There are 411 transcripts that are different, 166 up-regulated and 245 down-regulated


# FLX-autosomes and CFM-autosomes
FLXaCFMaF <- results(FLXswapF, contrast = c("type", "bCwtsFLXa", "dCwtsCFMa"))
# Then select the transcripts that are significantly different
FLXaCFMaFsig <- subset(FLXaCFMaF, padj < 0.05)
summary(FLXaCFMaFsig)
# There are 3519 transcripts that are different, 525 up-regulated and 2994 down-regulated


# CFM-sex chromosomes and CFM-autosomes
CFMsexCFMaF <- results(FLXswapF, contrast = c("type", "cCFMsCwta", "dCwtsCFMa"))
#Then select the transcripts that are significantly different
CFMsexCFMaFsig <- subset(CFMsexCFMaF, padj < 0.05)
summary(CFMsexCFMaFsig)
# There are 4137 transcripts that are different, 1130 up-regulated and 3007 down-regulated



########## DATA FOR PLOTTING ##########


# Next I have to make a data set for plotting

# Collect transcripts
genesFswap <- rownames(FLXswapF)

# Now create the new data frame
upsetplotF <- data.frame(gene.id = genesFswap,
                         FLXsexFLXauto = numeric(15505), FLXsexCFMsex = numeric(15505), FLXsexCFMauto = numeric(15505),
                         FLXautoCFMsex = numeric(15505), FLXautoCFMauto = numeric(15505), CFMsexCFMauto = numeric(15505))


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
  ylim(0, 2100) +
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
