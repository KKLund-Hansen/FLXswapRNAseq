################################################################################################
############################ TISSUE ENRICHMENT TEST FLX SWAP FEMALES ###########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Females/TissueTest")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
                                                                  "Malpighian tubule", "Mated spermatheca", "Midgut", "Ovary", "Rectal pad",
                                                                  "Salivary gland", "Testes", "Thoracicoabdominal ganglion", "Virgin spermatheca", "Whole fly"))


# To test for overabundance of transcripts of interest in a target tissue, I do an one-tailed Fisher's exact test on the observed 
# and expected tissue-specific transcripts of interest compared to the overall number of tissues-specific transcripts in each tissue.
# And collect all the data in a new data frame
Femaleswap.tissuetest <- data.frame("GAFupsig.oddsratio" = rep(0, 18), row.names = levels(Femaleswap.tissue$ind)[levels(Femaleswap.tissue$ind) != "Whole fly"])
Femaleswap.tissuetest$GAFupsig.pval <- 0
Femaleswap.tissuetest$GAFupsig.adjpval <- 0
Femaleswap.tissuetest$GAFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GAFdownsig.pval <- 0
Femaleswap.tissuetest$GAFdownsig.adjpval <- 0
Femaleswap.tissuetest$GBFupsig.oddsratio <- 0
Femaleswap.tissuetest$GBFupsig.pval <- 0
Femaleswap.tissuetest$GBFupsig.adjpval <- 0
Femaleswap.tissuetest$GBFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GBFdownsig.pval <- 0
Femaleswap.tissuetest$GBFdownsig.adjpval <- 0
Femaleswap.tissuetest$GCFupsig.oddsratio <- 0
Femaleswap.tissuetest$GCFupsig.pval <- 0
Femaleswap.tissuetest$GCFupsig.adjpval <- 0
Femaleswap.tissuetest$GCFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GCFdownsig.pval <- 0
Femaleswap.tissuetest$GCFdownsig.adjpval <- 0
Femaleswap.tissuetest$GDFupsig.oddsratio <- 0
Femaleswap.tissuetest$GDFupsig.pval <- 0
Femaleswap.tissuetest$GDFupsig.adjpval <- 0
Femaleswap.tissuetest$GDFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GDFdownsig.pval <- 0
Femaleswap.tissuetest$GDFdownsig.adjpval <- 0
Femaleswap.tissuetest$GEFupsig.oddsratio <- 0
Femaleswap.tissuetest$GEFupsig.pval <- 0
Femaleswap.tissuetest$GEFupsig.adjpval <- 0
Femaleswap.tissuetest$GEFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GEFdownsig.pval <- 0
Femaleswap.tissuetest$GEFdownsig.adjpval <- 0
Femaleswap.tissuetest$GFFupsig.oddsratio <- 0
Femaleswap.tissuetest$GFFupsig.pval <- 0
Femaleswap.tissuetest$GFFupsig.adjpval <- 0
Femaleswap.tissuetest$GFFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GFFdownsig.pval <- 0
Femaleswap.tissuetest$GFFdownsig.adjpval <- 0
Femaleswap.tissuetest$GGFupsig.oddsratio <- 0
Femaleswap.tissuetest$GGFupsig.pval <- 0
Femaleswap.tissuetest$GGFupsig.adjpval <- 0
Femaleswap.tissuetest$GGFdownsig.oddsratio <- 0
Femaleswap.tissuetest$GGFdownsig.pval <- 0
Femaleswap.tissuetest$GGFdownsig.adjpval <- 0




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GAFupsig <- Femaleswap.tissue$gene.id %in% GAFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GAFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GAFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GAFupsig.adjpval <- Femaleswap.tissuetest$GAFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GAFdownsig <- Femaleswap.tissue$gene.id %in% GAFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GAFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GAFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GAFdownsig.adjpval <- Femaleswap.tissuetest$GAFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GBFupsig <- Femaleswap.tissue$gene.id %in% GBFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GBFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GBFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GBFupsig.adjpval <- Femaleswap.tissuetest$GBFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GBFdownsig <- Femaleswap.tissue$gene.id %in% GBFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GBFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GBFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GBFdownsig.adjpval <- Femaleswap.tissuetest$GBFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GCFupsig <- Femaleswap.tissue$gene.id %in% GCFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GCFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GCFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GCFupsig.adjpval <- Femaleswap.tissuetest$GCFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts iin group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GCFdownsig <- Femaleswap.tissue$gene.id %in% GCFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GCFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GCFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GCFdownsig.adjpval <- Femaleswap.tissuetest$GCFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GDFupsig <- Femaleswap.tissue$gene.id %in% GDFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GDFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GDFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GDFupsig.adjpval <- Femaleswap.tissuetest$GDFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GDFdownsig <- Femaleswap.tissue$gene.id %in% GDFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GDFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GDFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GDFdownsig.adjpval <- Femaleswap.tissuetest$GDFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GEFupsig <- Femaleswap.tissue$gene.id %in% GEFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GEFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GEFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GEFupsig.adjpval <- Femaleswap.tissuetest$GEFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GEFdownsig <- Femaleswap.tissue$gene.id %in% GEFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GEFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GEFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GEFdownsig.adjpval <- Femaleswap.tissuetest$GEFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GFFupsig <- Femaleswap.tissue$gene.id %in% GFFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GFFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GFFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GFFupsig.adjpval <- Femaleswap.tissuetest$GFFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GFFdownsig <- Femaleswap.tissue$gene.id %in% GFFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GFFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GFFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GFFdownsig.adjpval <- Femaleswap.tissuetest$GFFdownsig.pval * nrow(Femaleswap.tissuetest)




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GGFupsig <- Femaleswap.tissue$gene.id %in% GGFupsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                   fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GGFupsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Femaleswap.tissuetest["Accessory glands", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                              fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GGFupsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GGFupsig.adjpval <- Femaleswap.tissuetest$GGFupsig.pval * nrow(Femaleswap.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# Add the significant transcripts to the tissue data frame
Femaleswap.tissue$GGFdownsig <- Femaleswap.tissue$gene.id %in% GGFdownsig$gene.id

# Add data to the odd ratio column
Femaleswap.tissuetest["Accessory glands", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Brain", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Carcass", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Crop", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Eye", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Fat body", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Head", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Heart", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Hindgut", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Malpighian tubule", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Mated spermatheca", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Midgut", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                                fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Ovary", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Rectal pad", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Salivary gland", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Testes", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                                fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                     fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Femaleswap.tissuetest["Virgin spermatheca", "GGFdownsig.oddsratio"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Femaleswap.tissuetest["Accessory glands", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Brain", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Carcass", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Crop", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Eye", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Fat body", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Head", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Heart", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Hindgut", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Malpighian tubule", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Mated spermatheca", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Midgut", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Ovary", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Rectal pad", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Salivary gland", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Testes", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Thoracicoabdominal ganglion", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Femaleswap.tissuetest["Virgin spermatheca", "GGFdownsig.pval"] <- with(Femaleswap.tissue[Femaleswap.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Femaleswap.tissuetest$GGFdownsig.adjpval <- Femaleswap.tissuetest$GGFdownsig.pval * nrow(Femaleswap.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Femaleswap.tissuetest, "FemaleSwapTissueTest.csv")

