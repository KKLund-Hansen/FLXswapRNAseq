################################################################################################
############################# TISSUE ENRICHMENT TEST FLX SWAP MALES ############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/TissueTest")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
# And collect all the data in a new data frame
Maleswap.tissuetest <- data.frame("GAMupsig.oddsratio" = rep(0, 18), row.names = levels(Maleswap.tissue$ind)[levels(Maleswap.tissue$ind) != "Whole fly"])
Maleswap.tissuetest$GAMupsig.pval <- 0
Maleswap.tissuetest$GAMupsig.adjpval <- 0
Maleswap.tissuetest$GAMdownsig.oddsratio <- 0
Maleswap.tissuetest$GAMdownsig.pval <- 0
Maleswap.tissuetest$GAMdownsig.adjpval <- 0
Maleswap.tissuetest$GBMupsig.oddsratio <- 0
Maleswap.tissuetest$GBMupsig.pval <- 0
Maleswap.tissuetest$GBMupsig.adjpval <- 0
Maleswap.tissuetest$GBMdownsig.oddsratio <- 0
Maleswap.tissuetest$GBMdownsig.pval <- 0
Maleswap.tissuetest$GBMdownsig.adjpval <- 0
Maleswap.tissuetest$GCMupsig.oddsratio <- 0
Maleswap.tissuetest$GCMupsig.pval <- 0
Maleswap.tissuetest$GCMupsig.adjpval <- 0
Maleswap.tissuetest$GCMdownsig.oddsratio <- 0
Maleswap.tissuetest$GCMdownsig.pval <- 0
Maleswap.tissuetest$GCMdownsig.adjpval <- 0
Maleswap.tissuetest$GDMupsig.oddsratio <- 0
Maleswap.tissuetest$GDMupsig.pval <- 0
Maleswap.tissuetest$GDMupsig.adjpval <- 0
Maleswap.tissuetest$GDMdownsig.oddsratio <- 0
Maleswap.tissuetest$GDMdownsig.pval <- 0
Maleswap.tissuetest$GDMdownsig.adjpval <- 0
Maleswap.tissuetest$GEMupsig.oddsratio <- 0
Maleswap.tissuetest$GEMupsig.pval <- 0
Maleswap.tissuetest$GEMupsig.adjpval <- 0
Maleswap.tissuetest$GEMdownsig.oddsratio <- 0
Maleswap.tissuetest$GEMdownsig.pval <- 0
Maleswap.tissuetest$GEMdownsig.adjpval <- 0
Maleswap.tissuetest$GFMupsig.oddsratio <- 0
Maleswap.tissuetest$GFMupsig.pval <- 0
Maleswap.tissuetest$GFMupsig.adjpval <- 0
Maleswap.tissuetest$GFMdownsig.oddsratio <- 0
Maleswap.tissuetest$GFMdownsig.pval <- 0
Maleswap.tissuetest$GFMdownsig.adjpval <- 0
Maleswap.tissuetest$GGMupsig.oddsratio <- 0
Maleswap.tissuetest$GGMupsig.pval <- 0
Maleswap.tissuetest$GGMupsig.adjpval <- 0
Maleswap.tissuetest$GGMdownsig.oddsratio <- 0
Maleswap.tissuetest$GGMdownsig.pval <- 0
Maleswap.tissuetest$GGMdownsig.adjpval <- 0




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group A
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GAMupsig <- Maleswap.tissue$gene.id %in% GAMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GAMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GAMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAMupsig == T, na.rm = T), sum(under == T & GAMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GAMupsig.adjpval <- Maleswap.tissuetest$GAMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GAMdownsig <- Maleswap.tissue$gene.id %in% GAMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GAMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GAMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAMdownsig == T, na.rm = T), sum(under == T & GAMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GAMdownsig.adjpval <- Maleswap.tissuetest$GAMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group B
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GBMupsig <- Maleswap.tissue$gene.id %in% GBMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GBMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GBMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBMupsig == T, na.rm = T), sum(under == T & GBMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GBMupsig.adjpval <- Maleswap.tissuetest$GBMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GBMdownsig <- Maleswap.tissue$gene.id %in% GBMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GBMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GBMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBMdownsig == T, na.rm = T), sum(under == T & GBMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GBMdownsig.adjpval <- Maleswap.tissuetest$GBMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GCMupsig <- Maleswap.tissue$gene.id %in% GCMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GCMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GCMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCMupsig == T, na.rm = T), sum(under == T & GCMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GCMupsig.adjpval <- Maleswap.tissuetest$GCMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts iin group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GCMdownsig <- Maleswap.tissue$gene.id %in% GCMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GCMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GCMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCMdownsig == T, na.rm = T), sum(under == T & GCMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GCMdownsig.adjpval <- Maleswap.tissuetest$GCMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GDMupsig <- Maleswap.tissue$gene.id %in% GDMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GDMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GDMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDMupsig == T, na.rm = T), sum(under == T & GDMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GDMupsig.adjpval <- Maleswap.tissuetest$GDMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GDMdownsig <- Maleswap.tissue$gene.id %in% GDMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GDMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GDMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDMdownsig == T, na.rm = T), sum(under == T & GDMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GDMdownsig.adjpval <- Maleswap.tissuetest$GDMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group E
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GEMupsig <- Maleswap.tissue$gene.id %in% GEMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GEMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GEMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEMupsig == T, na.rm = T), sum(under == T & GEMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GEMupsig.adjpval <- Maleswap.tissuetest$GEMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GEMdownsig <- Maleswap.tissue$gene.id %in% GEMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GEMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GEMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEMdownsig == T, na.rm = T), sum(under == T & GEMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GEMdownsig.adjpval <- Maleswap.tissuetest$GEMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GFMupsig <- Maleswap.tissue$gene.id %in% GFMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GFMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GFMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFMupsig == T, na.rm = T), sum(under == T & GFMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GFMupsig.adjpval <- Maleswap.tissuetest$GFMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GFMdownsig <- Maleswap.tissue$gene.id %in% GFMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GFMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GFMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFMdownsig == T, na.rm = T), sum(under == T & GFMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GFMdownsig.adjpval <- Maleswap.tissuetest$GFMdownsig.pval * nrow(Maleswap.tissuetest)




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group G
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GGMupsig <- Maleswap.tissue$gene.id %in% GGMupsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                 fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GGMupsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Maleswap.tissuetest["Accessory glands", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                    fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                     fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                               fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                            fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GGMupsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGMupsig == T, na.rm = T), sum(under == T & GGMupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GGMupsig.adjpval <- Maleswap.tissuetest$GGMupsig.pval * nrow(Maleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# Add the significant transcripts to the tissue data frame
Maleswap.tissue$GGMdownsig <- Maleswap.tissue$gene.id %in% GGMdownsig$gene.id

# Add data to the odd ratio column
Maleswap.tissuetest["Accessory glands", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Brain", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Carcass", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Crop", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Eye", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Fat body", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Head", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Heart", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Hindgut", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Malpighian tubule", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Mated spermatheca", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Midgut", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Ovary", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Rectal pad", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Salivary gland", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Testes", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Maleswap.tissuetest["Virgin spermatheca", "GGMdownsig.oddsratio"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Maleswap.tissuetest["Accessory glands", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Brain", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Carcass", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Crop", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Eye", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Fat body", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Head", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Heart", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Hindgut", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Malpighian tubule", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Mated spermatheca", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Midgut", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Ovary", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Rectal pad", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Salivary gland", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Testes", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Thoracicoabdominal ganglion", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Maleswap.tissuetest["Virgin spermatheca", "GGMdownsig.pval"] <- with(Maleswap.tissue[Maleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGMdownsig == T, na.rm = T), sum(under == T & GGMdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Maleswap.tissuetest$GGMdownsig.adjpval <- Maleswap.tissuetest$GGMdownsig.pval * nrow(Maleswap.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Maleswap.tissuetest, "MaleswapTissueTest.csv")

