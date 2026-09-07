################################################################################################
############################# TISSUE ENRICHMENT TEST FLX SWAP MALES ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/TissueTest")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
flyatlas <- read.csv("FlyAtlas2.csv", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First change the name of the first row
colnames(flyatlas)[1] <- "gene.id"

# Next subset the fly atlas data set to the male data only
Maleflyatlas <- flyatlas[,c(1:5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33:37)]

# Before, I start to do the tests, I first need to rescale the data set to the whole fly baseline expression 
# to ensure an equal signal baseline for all the tissues, by using the using the average expression of unexpressed genes
# (n=599, expression value in the whole fly smaller than 3.4) (Innocenti https://doi.org/10.1371/journal.pbio.1000335)
standard.geneM <- apply(Maleflyatlas[Maleflyatlas$Whole.M < 3.4, 5:23], 2, mean)
standard.geneM2 <- standard.geneM - standard.geneM["Whole.M"]
Maleflyatlas <- as.data.frame(scale(Maleflyatlas[,5:23], center = standard.geneM2, scale = FALSE))

# Add transcript names as row names
row.names(Maleflyatlas) <- flyatlas$gene.id


# Next I need to make a new data frame for the actual analysis
Maleswap.tissue <- stack(Maleflyatlas)
Maleswap.tissue$gene.id <- rownames(Maleflyatlas)
Maleswap.tissue$gene.id <- gsub(" ", "", as.character(Maleswap.tissue$gene.id))
Maleswap.tissue$reference <- Maleflyatlas$Whole.M
# Transcripts are considered tissue specific if the expression level in the target tissue is 2-fold higher than in the whole fly.
Maleswap.tissue$over <- (Maleswap.tissue$reference + 1) < Maleswap.tissue$values
Maleswap.tissue$under <- (Maleswap.tissue$reference + 1) > Maleswap.tissue$values
# Next specify the level for the tissue, needed for the analysis
levels(Maleswap.tissue$ind) <- c("Whole fly", "Carcass", "Head", "Eye", "Brain", "Thoracicoabdominal ganglion", "Crop", "Midgut", "Hindgut",
                               "Malpighian tubule", "Rectal pad", "Fat body", "Heart", "Salivary gland", "Ovary", "Virgin spermatheca",
                               "Mated spermatheca", "Testes", "Accessory glands")
Maleswap.tissue$ind <- factor(Maleswap.tissue$ind, levels = c("Accessory glands", "Brain", "Carcass", "Crop", "Eye", "Fat body", "Head", "Heart", "Hindgut",
                                                          "Malpighian tubule", "Mated spermatheca", "Midgut", "Ovary", "Rectal pad",
                                                          "Salivary gland", "Testes", "Thoracicoabdominal ganglion", "Virgin spermatheca", "Whole fly"))

# To test for overabundance of transcripts of interest in a target tissue, I do an one-tailed Fisher's exact test on the observed 
# and expected tissue-specific transcripts of interest compared to the overall number of tissues-specific transcripts in each tissue.
# and collect all the data in a new data frame
Maleswap.tissuetest <- data.frame("FvCmupsig.oddsratio" = rep(0, 18), row.names = levels(Maleswap.tissue$ind)[levels(Maleswap.tissue$ind) != "Whole fly"])
Maleswap.tissuetest$FvCmupsig.pval <- 0
Maleswap.tissuetest$FvCmupsig.adjpval <- 0
Maleswap.tissuetest$FvCmdownsig.oddsratio <- 0
Maleswap.tissuetest$FvCmdownsig.pval <- 0
Maleswap.tissuetest$FvCmdownsig.adjpval <- 0
Maleswap.tissuetest$AvSmupsig.oddsratio <- 0
Maleswap.tissuetest$AvSmupsig.pval <- 0
Maleswap.tissuetest$AvSmupsig.adjpval <- 0
Maleswap.tissuetest$AvSmdownsig.oddsratio <- 0
Maleswap.tissuetest$AvSmdownsig.pval <- 0
Maleswap.tissuetest$AvSmdownsig.adjpval <- 0
Maleswap.tissuetest$FsvCsMupsig.oddsratio <- 0
Maleswap.tissuetest$FsvCsMupsig.pval <- 0
Maleswap.tissuetest$FsvCsMupsig.adjpval <- 0
Maleswap.tissuetest$FsvCsMdownsig.oddsratio <- 0
Maleswap.tissuetest$FsvCsMdownsig.pval <- 0
Maleswap.tissuetest$FsvCsMdownsig.adjpval <- 0
Maleswap.tissuetest$FavCaMupsig.oddsratio <- 0
Maleswap.tissuetest$FavCaMupsig.pval <- 0
Maleswap.tissuetest$FavCaMupsig.adjpval <- 0
Maleswap.tissuetest$FavCaMdownsig.oddsratio <- 0
Maleswap.tissuetest$FavCaMdownsig.pval <- 0
Maleswap.tissuetest$FavCaMdownsig.adjpval <- 0



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FvCmupsig <- Maleswap.tissue$gene.id %in% FvCmupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FvCmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FvCmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & FvCmupsig == T, na.rm = T), sum(under == T & FvCmupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FvCmupsig.adjpval <- Maleswap.tissuetest$FvCmupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FvCmdownsig <- Maleswap.tissue$gene.id %in% FvCmdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FvCmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FvCmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & FvCmdownsig == T, na.rm = T), sum(under == T & FvCmdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FvCmdownsig.adjpval <- Maleswap.tissuetest$FvCmdownsig.pval * nrow(Maleswap.tissuetest)




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$AvSmupsig <- Maleswap.tissue$gene.id %in% AvSmupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "AvSmupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "AvSmupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & AvSmupsig == T, na.rm = T), sum(under == T & AvSmupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$AvSmupsig.adjpval <- Maleswap.tissuetest$AvSmupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$AvSmdownsig <- Maleswap.tissue$gene.id %in% AvSmdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "AvSmdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "AvSmdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & AvSmdownsig == T, na.rm = T), sum(under == T & AvSmdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$AvSmdownsig.adjpval <- Maleswap.tissuetest$AvSmdownsig.pval * nrow(Maleswap.tissuetest)




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FsvCsMupsig <- Maleswap.tissue$gene.id %in% FsvCsMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FsvCsMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FsvCsMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & FsvCsMupsig == T, na.rm = T), sum(under == T & FsvCsMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FsvCsMupsig.adjpval <- Maleswap.tissuetest$FsvCsMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FsvCsMdownsig <- Maleswap.tissue$gene.id %in% FsvCsMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FsvCsMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FsvCsMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FsvCsMdownsig == T, na.rm = T), sum(under == T & FsvCsMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FsvCsMdownsig.adjpval <- Maleswap.tissuetest$FsvCsMdownsig.pval * nrow(Maleswap.tissuetest)




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FavCaMupsig <- Maleswap.tissue$gene.id %in% FavCaMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                         fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                          fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                       fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                    fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FavCaMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                    fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                         fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                        fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                       fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                        fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                         fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                     fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                         fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                               fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FavCaMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & FavCaMupsig == T, na.rm = T), sum(under == T & FavCaMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FavCaMupsig.adjpval <- Maleswap.tissuetest$FavCaMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$FavCaMdownsig <- Maleswap.tissue$gene.id %in% FavCaMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                           fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                   fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                               fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                                  fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                            fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                     fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                         fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                                 fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                      fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "FavCaMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                             fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "FavCaMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & FavCaMdownsig == T, na.rm = T), sum(under == T & FavCaMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$FavCaMdownsig.adjpval <- Maleswap.tissuetest$FavCaMdownsig.pval * nrow(Maleswap.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Maleswap.tissuetest, "MaleswapTissueTest.csv")

