################################################################################################
####################### TISSUE ENRICHMENT TEST FLX SWAP MALES - 7 GROUPS #######################
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
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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

# Add the transcript names as row names
row.names(Maleflyatlas) <- flyatlas$gene.id

# Next I need to make a new data frame for the actual analysis
Mswap7g.tissue <- stack(Maleflyatlas)
Mswap7g.tissue$gene.id <- rownames(Maleflyatlas)
Mswap7g.tissue$gene.id <- gsub(" ", "", as.character(Mswap7g.tissue$gene.id))
Mswap7g.tissue$reference <- Maleflyatlas$Whole.M
# Transcripts are considered tissue specific if the expression level in the target tissue is 2-fold higher than in the whole fly.
Mswap7g.tissue$over <- (Mswap7g.tissue$reference + 1) < Mswap7g.tissue$values
Mswap7g.tissue$under <- (Mswap7g.tissue$reference + 1) > Mswap7g.tissue$values
# Next specify the level for the tissue, needed for the analysis
levels(Mswap7g.tissue$ind) <- c("Whole fly", "Carcass", "Head", "Eye", "Brain", "Thoracicoabdominal ganglion", "Crop", "Midgut", "Hindgut", 
                                "Malpighian tubule", "Rectal pad", "Fat body", "Heart", "Salivary gland", "Ovary", "Virgin spermatheca", 
                                "Mated spermatheca", "Testes", "Accessory glands")
Mswap7g.tissue$ind <- factor(Mswap7g.tissue$ind, levels = c("Accessory glands", "Brain", "Carcass", "Crop", "Eye", "Fat body", "Head", "Heart", "Hindgut", 
                                                            "Malpighian tubule", "Mated spermatheca", "Midgut", "Ovary", "Rectal pad",
                                                            "Salivary gland", "Testes", "Thoracicoabdominal ganglion", "Virgin spermatheca", "Whole fly"))


# To test for overabundance of transcripts of interest in a target tissue, I do an one-tailed Fisher's exact test on the observed 
# and expected tissue-specific transcripts of interest compared to the overall number of tissues-specific transcripts in each tissue.
# And collect all the data in a new data frame
Mswap7g.tissuetest <- data.frame("GAMupsig.oddsratio" = rep(0, 18), row.names = levels(Mswap7g.tissue$ind)[levels(Mswap7g.tissue$ind) != "Whole fly"])
Mswap7g.tissuetest$GAMupsig.pval <- 0
Mswap7g.tissuetest$GAMupsig.adjpval <- 0
Mswap7g.tissuetest$GAMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GAMdownsig.pval <- 0
Mswap7g.tissuetest$GAMdownsig.adjpval <- 0
Mswap7g.tissuetest$GBMupsig.oddsratio <- 0
Mswap7g.tissuetest$GBMupsig.pval <- 0
Mswap7g.tissuetest$GBMupsig.adjpval <- 0
Mswap7g.tissuetest$GBMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GBMdownsig.pval <- 0
Mswap7g.tissuetest$GBMdownsig.adjpval <- 0
Mswap7g.tissuetest$GCMupsig.oddsratio <- 0
Mswap7g.tissuetest$GCMupsig.pval <- 0
Mswap7g.tissuetest$GCMupsig.adjpval <- 0
Mswap7g.tissuetest$GCMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GCMdownsig.pval <- 0
Mswap7g.tissuetest$GCMdownsig.adjpval <- 0
Mswap7g.tissuetest$GDMupsig.oddsratio <- 0
Mswap7g.tissuetest$GDMupsig.pval <- 0
Mswap7g.tissuetest$GDMupsig.adjpval <- 0
Mswap7g.tissuetest$GDMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GDMdownsig.pval <- 0
Mswap7g.tissuetest$GDMdownsig.adjpval <- 0
Mswap7g.tissuetest$GEMupsig.oddsratio <- 0
Mswap7g.tissuetest$GEMupsig.pval <- 0
Mswap7g.tissuetest$GEMupsig.adjpval <- 0
Mswap7g.tissuetest$GEMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GEMdownsig.pval <- 0
Mswap7g.tissuetest$GEMdownsig.adjpval <- 0
Mswap7g.tissuetest$GFMupsig.oddsratio <- 0
Mswap7g.tissuetest$GFMupsig.pval <- 0
Mswap7g.tissuetest$GFMupsig.adjpval <- 0
Mswap7g.tissuetest$GFMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GFMdownsig.pval <- 0
Mswap7g.tissuetest$GFMdownsig.adjpval <- 0
Mswap7g.tissuetest$GGMupsig.oddsratio <- 0
Mswap7g.tissuetest$GGMupsig.pval <- 0
Mswap7g.tissuetest$GGMupsig.adjpval <- 0
Mswap7g.tissuetest$GGMdownsig.oddsratio <- 0
Mswap7g.tissuetest$GGMdownsig.pval <- 0
Mswap7g.tissuetest$GGMdownsig.adjpval <- 0




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GAMupsig <- Mswap7g.tissue$gene.id %in% GAMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GAMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GAMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GAMupsig.adjpval <- Mswap7g.tissuetest$GAMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GAMdownsig <- Mswap7g.tissue$gene.id %in% GAMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GAMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GAMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GAMdownsig.adjpval <- Mswap7g.tissuetest$GAMdownsig.pval * nrow(Mswap7g.tissuetest)




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GBMupsig <- Mswap7g.tissue$gene.id %in% GBMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GBMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GBMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GBMupsig.adjpval <- Mswap7g.tissuetest$GBMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GBMdownsig <- Mswap7g.tissue$gene.id %in% GBMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GBMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GBMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GBMdownsig.adjpval <- Mswap7g.tissuetest$GBMdownsig.pval * nrow(Mswap7g.tissuetest)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GCMupsig <- Mswap7g.tissue$gene.id %in% GCMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GCMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GCMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GCMupsig.adjpval <- Mswap7g.tissuetest$GCMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts iin group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GCMdownsig <- Mswap7g.tissue$gene.id %in% GCMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GCMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GCMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GCMdownsig.adjpval <- Mswap7g.tissuetest$GCMdownsig.pval * nrow(Mswap7g.tissuetest)




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GDMupsig <- Mswap7g.tissue$gene.id %in% GDMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GDMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GDMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GDMupsig.adjpval <- Mswap7g.tissuetest$GDMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GDMdownsig <- Mswap7g.tissue$gene.id %in% GDMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GDMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GDMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GDMdownsig.adjpval <- Mswap7g.tissuetest$GDMdownsig.pval * nrow(Mswap7g.tissuetest)




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GEMupsig <- Mswap7g.tissue$gene.id %in% GEMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GEMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GEMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GEMupsig.adjpval <- Mswap7g.tissuetest$GEMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 20

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GEMdownsig <- Mswap7g.tissue$gene.id %in% GEMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GEMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GEMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GEMdownsig.adjpval <- Mswap7g.tissuetest$GEMdownsig.pval * nrow(Mswap7g.tissuetest)




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GFMupsig <- Mswap7g.tissue$gene.id %in% GFMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GFMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GFMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GFMupsig.adjpval <- Mswap7g.tissuetest$GFMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GFMdownsig <- Mswap7g.tissue$gene.id %in% GFMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GFMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GFMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GFMdownsig.adjpval <- Mswap7g.tissuetest$GFMdownsig.pval * nrow(Mswap7g.tissuetest)





######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GGMupsig <- Mswap7g.tissue$gene.id %in% GGMupsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GGMupsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Mswap7g.tissuetest["Accessory glands", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GGMupsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GGMupsig.adjpval <- Mswap7g.tissuetest$GGMupsig.pval * nrow(Mswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# Add the significant transcripts to the tissue data frame
Mswap7g.tissue$GGMdownsig <- Mswap7g.tissue$gene.id %in% GGMdownsig$gene.id

# Add data to the odd ratio column
Mswap7g.tissuetest["Accessory glands", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Brain", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Carcass", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Crop", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Eye", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Fat body", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Head", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Heart", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Hindgut", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Malpighian tubule", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Mated spermatheca", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Midgut", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Ovary", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Rectal pad", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Salivary gland", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Testes", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Mswap7g.tissuetest["Virgin spermatheca", "GGMdownsig.oddsratio"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Mswap7g.tissuetest["Accessory glands", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Brain", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Carcass", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Crop", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Eye", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Fat body", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Head", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Heart", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Hindgut", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Malpighian tubule", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Mated spermatheca", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Midgut", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Ovary", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Rectal pad", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Salivary gland", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Testes", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Thoracicoabdominal ganglion", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Mswap7g.tissuetest["Virgin spermatheca", "GGMdownsig.pval"] <- with(Mswap7g.tissue[Mswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Mswap7g.tissuetest$GGMdownsig.adjpval <- Mswap7g.tissuetest$GGMdownsig.pval * nrow(Mswap7g.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Mswap7g.tissuetest, "Maleswap7gTissueTest.csv")

