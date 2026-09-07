################################################################################################
############################ TISSUE ENRICHMENT TEST FLX SWAP FEMALES ###########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/TissueTest")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
flyatlas <- read.csv("FlyAtlas2.csv", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First change the name of the first row
colnames(flyatlas)[1] <- "gene.id"

# Next subset the fly atlas data set to the female data only
Femaleflyatlas <- flyatlas[,c(1:4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32:37)]

# Before, I start to do the tests, I first need to rescale the data set to the whole fly baseline expression 
# to ensure an equal signal baseline for all the tissues, by using the using the average expression of unexpressed genes
# (n=599, expression value in the whole fly smaller than 3.4) (Innocenti https://doi.org/10.1371/journal.pbio.1000335)
standard.geneF <- apply(Femaleflyatlas[Femaleflyatlas$Whole.F < 3.4, 5:23], 2, mean)
standard.geneF2 <- standard.geneF - standard.geneF["Whole.F"]
Femaleflyatlas <- as.data.frame(scale(Femaleflyatlas[,5:23], center = standard.geneF2, scale = FALSE))

# Add transcript names as row names
row.names(Femaleflyatlas) <- flyatlas$gene.id


# Next I need to make a new data frame for the actual analysis
Femaleswap.tissue <- stack(Femaleflyatlas)
Femaleswap.tissue$gene.id <- rownames(Femaleflyatlas)
Femaleswap.tissue$gene.id <- gsub(" ", "", as.character(Femaleswap.tissue$gene.id))
Femaleswap.tissue$reference <- Femaleflyatlas$Whole.F
# Transcripts are considered tissue specific if the expression level in the target tissue is 2-fold higher than in the whole fly.
Femaleswap.tissue$over <- (Femaleswap.tissue$reference + 1) < Femaleswap.tissue$values
Femaleswap.tissue$under <- (Femaleswap.tissue$reference + 1) > Femaleswap.tissue$values
# Next specify the level for the tissue, needed for the analysis
levels(Femaleswap.tissue$ind) <- c("Whole fly", "Carcass", "Head", "Eye", "Brain", "Thoracicoabdominal ganglion", "Crop", "Midgut", "Hindgut",
                                   "Malpighian tubule", "Rectal pad", "Fat body", "Heart", "Salivary gland", "Ovary", "Virgin spermatheca",
                                   "Mated spermatheca", "Testes", "Accessory glands")
Femaleswap.tissue$ind <- factor(Femaleswap.tissue$ind, levels = c("Accessory glands", "Brain", "Carcass", "Crop", "Eye", "Fat body", "Head", "Heart", "Hindgut",
                                                                  "Malpighian tubule", "Mated spermatheca", "Midgut", "Ovary", "Rectal pad", "Salivary gland",
                                                                  "Testes", "Thoracicoabdominal ganglion", "Virgin spermatheca", "Whole fly"))

# To test for overabundance of transcripts of interest in a target tissue, I do an one-tailed Fisher's exact test on the observed 
# and expected tissue-specific transcripts of interest compared to the overall number of tissues-specific transcripts in each tissue.
# and collect all the data in a new data frame
Femaleswap.tissuetest <- data.frame("FvCfupsig.oddsratio" = rep(0, 18), row.names = levels(Femaleswap.tissue$ind)[levels(Femaleswap.tissue$ind) != "Whole fly"])
Femaleswap.tissuetest$FvCfupsig.pval <- 0
Femaleswap.tissuetest$FvCfupsig.adjpval <- 0
Femaleswap.tissuetest$FvCfdownsig.oddsratio <- 0
Femaleswap.tissuetest$FvCfdownsig.pval <- 0
Femaleswap.tissuetest$FvCfdownsig.adjpval <- 0
Femaleswap.tissuetest$AvSfupsig.oddsratio <- 0
Femaleswap.tissuetest$AvSfupsig.pval <- 0
Femaleswap.tissuetest$AvSfupsig.adjpval <- 0
Femaleswap.tissuetest$AvSfdownsig.oddsratio <- 0
Femaleswap.tissuetest$AvSfdownsig.pval <- 0
Femaleswap.tissuetest$AvSfdownsig.adjpval <- 0
Femaleswap.tissuetest$FsvCsFupsig.oddsratio <- 0
Femaleswap.tissuetest$FsvCsFupsig.pval <- 0
Femaleswap.tissuetest$FsvCsFupsig.adjpval <- 0
Femaleswap.tissuetest$FsvCsFdownsig.oddsratio <- 0
Femaleswap.tissuetest$FsvCsFdownsig.pval <- 0
Femaleswap.tissuetest$FsvCsFdownsig.adjpval <- 0
Femaleswap.tissuetest$FavCaFupsig.oddsratio <- 0
Femaleswap.tissuetest$FavCaFupsig.pval <- 0
Femaleswap.tissuetest$FavCaFupsig.adjpval <- 0
Femaleswap.tissuetest$FavCaFdownsig.oddsratio <- 0
Femaleswap.tissuetest$FavCaFdownsig.pval <- 0
Femaleswap.tissuetest$FavCaFdownsig.adjpval <- 0



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FvCfupsig <- Femaleswap.tissue$gene.id %in% FvCfupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FvCfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FvCfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & FvCfupsig == T, na.rm = T), sum(under == T & FvCfupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FvCfupsig.adjpval <- Femaleswap.tissuetest$FvCfupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FvCfdownsig <- Femaleswap.tissue$gene.id %in% FvCfdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FvCfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FvCfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FvCfdownsig == T, na.rm = T), sum(under == T & FvCfdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FvCfdownsig.adjpval <- Femaleswap.tissuetest$FvCfdownsig.pval * nrow(Femaleswap.tissuetest)




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$AvSfupsig <- Femaleswap.tissue$gene.id %in% AvSfupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "AvSfupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "AvSfupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & AvSfupsig == T, na.rm = T), sum(under == T & AvSfupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$AvSfupsig.adjpval <- Femaleswap.tissuetest$AvSfupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$AvSfdownsig <- Femaleswap.tissue$gene.id %in% AvSfdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "AvSfdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "AvSfdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & AvSfdownsig == T, na.rm = T), sum(under == T & AvSfdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$AvSfdownsig.adjpval <- Femaleswap.tissuetest$AvSfdownsig.pval * nrow(Femaleswap.tissuetest)




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FsvCsFupsig <- Femaleswap.tissue$gene.id %in% FsvCsFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FsvCsFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FsvCsFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FsvCsFupsig == T, na.rm = T), sum(under == T & FsvCsFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FsvCsFupsig.adjpval <- Femaleswap.tissuetest$FsvCsFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FsvCsFdownsig <- Femaleswap.tissue$gene.id %in% FsvCsFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                             fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                    fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                     fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                    fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                              fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                              fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                   fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                           fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                   fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                        fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FsvCsFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                               fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FsvCsFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FsvCsFdownsig == T, na.rm = T), sum(under == T & FsvCsFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FsvCsFdownsig.adjpval <- Femaleswap.tissuetest$FsvCsFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FavCaFupsig <- Femaleswap.tissue$gene.id %in% FavCaFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FavCaFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FavCaFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FavCaFupsig == T, na.rm = T), sum(under == T & FavCaFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FavCaFupsig.adjpval <- Femaleswap.tissuetest$FavCaFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$FavCaFdownsig <- Femaleswap.tissue$gene.id %in% FavCaFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                             fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                    fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                     fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                    fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                              fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                              fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                           fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                        fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "FavCaFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                               fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "FavCaFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FavCaFdownsig == T, na.rm = T), sum(under == T & FavCaFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$FavCaFdownsig.adjpval <- Femaleswap.tissuetest$FavCaFdownsig.pval * nrow(Femaleswap.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Femaleswap.tissuetest, "FemaleswapTissueTest.csv")

