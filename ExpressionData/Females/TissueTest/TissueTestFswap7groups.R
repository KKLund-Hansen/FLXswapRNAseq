################################################################################################
###################### TISSUE ENRICHMENT TEST FLX SWAP FEMALES - 7 GROUPS ######################
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
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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
Fswap7g.tissue <- stack(Femaleflyatlas)
Fswap7g.tissue$gene.id <- rownames(Femaleflyatlas)
Fswap7g.tissue$gene.id <- gsub(" ", "", as.character(Fswap7g.tissue$gene.id))
Fswap7g.tissue$reference <- Femaleflyatlas$Whole.F
# Transcripts are considered tissue specific if the expression level in the target tissue is 2-fold higher than in the whole fly.
Fswap7g.tissue$over <- (Fswap7g.tissue$reference + 1) < Fswap7g.tissue$values
Fswap7g.tissue$under <- (Fswap7g.tissue$reference + 1) > Fswap7g.tissue$values
# Next specify the level for the tissue, needed for the analysis
levels(Fswap7g.tissue$ind) <- c("Whole fly", "Carcass", "Head", "Eye", "Brain", "Thoracicoabdominal ganglion", "Crop", "Midgut", "Hindgut",
                                "Malpighian tubule", "Rectal pad", "Fat body", "Heart", "Salivary gland", "Ovary", "Virgin spermatheca", 
                                "Mated spermatheca", "Testes", "Accessory glands")
Fswap7g.tissue$ind <- factor(Fswap7g.tissue$ind, levels = c("Accessory glands", "Brain", "Carcass", "Crop", "Eye", "Fat body", "Head", "Heart", "Hindgut", 
                                                            "Malpighian tubule", "Mated spermatheca", "Midgut", "Ovary", "Rectal pad",
                                                            "Salivary gland", "Testes", "Thoracicoabdominal ganglion", "Virgin spermatheca", "Whole fly"))


# To test for overabundance of transcripts of interest in a target tissue, I do an one-tailed Fisher's exact test on the observed 
# and expected tissue-specific transcripts of interest compared to the overall number of tissues-specific transcripts in each tissue.
# And collect all the data in a new data frame
Fswap7g.tissuetest <- data.frame("GAFupsig.oddsratio" = rep(0, 18), row.names = levels(Fswap7g.tissue$ind)[levels(Fswap7g.tissue$ind) != "Whole fly"])
Fswap7g.tissuetest$GAFupsig.pval <- 0
Fswap7g.tissuetest$GAFupsig.adjpval <- 0
Fswap7g.tissuetest$GAFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GAFdownsig.pval <- 0
Fswap7g.tissuetest$GAFdownsig.adjpval <- 0
Fswap7g.tissuetest$GBFupsig.oddsratio <- 0
Fswap7g.tissuetest$GBFupsig.pval <- 0
Fswap7g.tissuetest$GBFupsig.adjpval <- 0
Fswap7g.tissuetest$GBFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GBFdownsig.pval <- 0
Fswap7g.tissuetest$GBFdownsig.adjpval <- 0
Fswap7g.tissuetest$GCFupsig.oddsratio <- 0
Fswap7g.tissuetest$GCFupsig.pval <- 0
Fswap7g.tissuetest$GCFupsig.adjpval <- 0
Fswap7g.tissuetest$GCFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GCFdownsig.pval <- 0
Fswap7g.tissuetest$GCFdownsig.adjpval <- 0
Fswap7g.tissuetest$GDFupsig.oddsratio <- 0
Fswap7g.tissuetest$GDFupsig.pval <- 0
Fswap7g.tissuetest$GDFupsig.adjpval <- 0
Fswap7g.tissuetest$GDFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GDFdownsig.pval <- 0
Fswap7g.tissuetest$GDFdownsig.adjpval <- 0
Fswap7g.tissuetest$GEFupsig.oddsratio <- 0
Fswap7g.tissuetest$GEFupsig.pval <- 0
Fswap7g.tissuetest$GEFupsig.adjpval <- 0
Fswap7g.tissuetest$GEFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GEFdownsig.pval <- 0
Fswap7g.tissuetest$GEFdownsig.adjpval <- 0
Fswap7g.tissuetest$GFFupsig.oddsratio <- 0
Fswap7g.tissuetest$GFFupsig.pval <- 0
Fswap7g.tissuetest$GFFupsig.adjpval <- 0
Fswap7g.tissuetest$GFFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GFFdownsig.pval <- 0
Fswap7g.tissuetest$GFFdownsig.adjpval <- 0
Fswap7g.tissuetest$GGFupsig.oddsratio <- 0
Fswap7g.tissuetest$GGFupsig.pval <- 0
Fswap7g.tissuetest$GGFupsig.adjpval <- 0
Fswap7g.tissuetest$GGFdownsig.oddsratio <- 0
Fswap7g.tissuetest$GGFdownsig.pval <- 0
Fswap7g.tissuetest$GGFdownsig.adjpval <- 0




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GAFupsig <- Fswap7g.tissue$gene.id %in% GAFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GAFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GAFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAFupsig == T, na.rm = T), sum(under == T & GAFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GAFupsig.adjpval <- Fswap7g.tissuetest$GAFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GAFdownsig <- Fswap7g.tissue$gene.id %in% GAFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GAFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GAFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GAFdownsig == T, na.rm = T), sum(under == T & GAFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GAFdownsig.adjpval <- Fswap7g.tissuetest$GAFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GBFupsig <- Fswap7g.tissue$gene.id %in% GBFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GBFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GBFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBFupsig == T, na.rm = T), sum(under == T & GBFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GBFupsig.adjpval <- Fswap7g.tissuetest$GBFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GBFdownsig <- Fswap7g.tissue$gene.id %in% GBFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GBFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GBFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GBFdownsig == T, na.rm = T), sum(under == T & GBFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GBFdownsig.adjpval <- Fswap7g.tissuetest$GBFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GCFupsig <- Fswap7g.tissue$gene.id %in% GCFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GCFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GCFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCFupsig == T, na.rm = T), sum(under == T & GCFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GCFupsig.adjpval <- Fswap7g.tissuetest$GCFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts iin group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GCFdownsig <- Fswap7g.tissue$gene.id %in% GCFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GCFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GCFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GCFdownsig == T, na.rm = T), sum(under == T & GCFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GCFdownsig.adjpval <- Fswap7g.tissuetest$GCFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GDFupsig <- Fswap7g.tissue$gene.id %in% GDFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GDFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GDFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDFupsig == T, na.rm = T), sum(under == T & GDFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GDFupsig.adjpval <- Fswap7g.tissuetest$GDFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GDFdownsig <- Fswap7g.tissue$gene.id %in% GDFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GDFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GDFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GDFdownsig == T, na.rm = T), sum(under == T & GDFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GDFdownsig.adjpval <- Fswap7g.tissuetest$GDFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GEFupsig <- Fswap7g.tissue$gene.id %in% GEFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GEFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GEFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEFupsig == T, na.rm = T), sum(under == T & GEFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GEFupsig.adjpval <- Fswap7g.tissuetest$GEFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GEFdownsig <- Fswap7g.tissue$gene.id %in% GEFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GEFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GEFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GEFdownsig == T, na.rm = T), sum(under == T & GEFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GEFdownsig.adjpval <- Fswap7g.tissuetest$GEFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GFFupsig <- Fswap7g.tissue$gene.id %in% GFFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GFFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GFFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFFupsig == T, na.rm = T), sum(under == T & GFFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GFFupsig.adjpval <- Fswap7g.tissuetest$GFFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GFFdownsig <- Fswap7g.tissue$gene.id %in% GFFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GFFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GFFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GFFdownsig == T, na.rm = T), sum(under == T & GFFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GFFdownsig.adjpval <- Fswap7g.tissuetest$GFFdownsig.pval * nrow(Fswap7g.tissuetest)




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GGFupsig <- Fswap7g.tissue$gene.id %in% GGFupsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GGFupsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P-value column
Fswap7g.tissuetest["Accessory glands", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                     fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                    fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                   fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                    fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                     fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                     fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                      fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                           fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GGFupsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGFupsig == T, na.rm = T), sum(under == T & GGFupsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GGFupsig.adjpval <- Fswap7g.tissuetest$GGFupsig.pval * nrow(Fswap7g.tissuetest)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to look for tissue specificity for down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# Add the significant transcripts to the tissue data frame
Fswap7g.tissue$GGFdownsig <- Fswap7g.tissue$gene.id %in% GGFdownsig$gene.id

# Add data to the odd ratio column
Fswap7g.tissuetest["Accessory glands", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                       fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Brain", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Carcass", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Crop", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Eye", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Fat body", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                               fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                    sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Head", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                           fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Heart", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Hindgut", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                              fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                   sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Malpighian tubule", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Mated spermatheca", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Midgut", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Ovary", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Rectal pad", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                                 fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Salivary gland", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                     fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Testes", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                             fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                                  fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                      sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)
Fswap7g.tissuetest["Virgin spermatheca", "GGFdownsig.oddsratio"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                         fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$estimate)

# Add data to the P value column
Fswap7g.tissuetest["Accessory glands", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Accessory glands", ],
                                                                  fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                       sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Brain", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Brain", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Carcass", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Carcass", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Crop", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Crop", ],
                                                      fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Eye", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Eye", ],
                                                     fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                          sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Fat body", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Fat body", ],
                                                          fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                               sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Head", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Head", ],
                                                      fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                           sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Heart", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Heart", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Hindgut", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Hindgut", ],
                                                         fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                              sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Malpighian tubule", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Malpighian tubule", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Mated spermatheca", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Mated spermatheca", ],
                                                                   fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                        sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Midgut", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Midgut", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Ovary", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Ovary", ],
                                                       fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                            sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Rectal pad", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Rectal pad", ],
                                                            fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                 sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Salivary gland", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Salivary gland", ],
                                                                fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                     sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Testes", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Testes", ],
                                                        fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                             sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Thoracicoabdominal ganglion", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Thoracicoabdominal ganglion", ],
                                                                             fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                                  sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)
Fswap7g.tissuetest["Virgin spermatheca", "GGFdownsig.pval"] <- with(Fswap7g.tissue[Fswap7g.tissue$ind == "Virgin spermatheca", ],
                                                                    fisher.test(matrix(c(sum(over == T & GGFdownsig == T, na.rm = T), sum(under == T & GGFdownsig == T, na.rm = T),
                                                                                         sum(over == T, na.rm = T), sum(under == T, na.rm = T)), ncol = 2), alternative = "greater")$p.val)

# Adjust the P value
Fswap7g.tissuetest$GGFdownsig.adjpval <- Fswap7g.tissuetest$GGFdownsig.pval * nrow(Fswap7g.tissuetest)





# Now save the full new data frame as a csv file
write.csv(Fswap7g.tissuetest, "Femaleswap7gTissueTest.csv")

