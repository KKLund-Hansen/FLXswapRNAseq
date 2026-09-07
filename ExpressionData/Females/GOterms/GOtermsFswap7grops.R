################################################################################################
###############$#### GENE ONTOLOGY TERM ANALYSES FLX SWAP FEMALES - 7 GROUPS ###################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/GOterms")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# First collect the gene IDs in a new vector
genesGAFupsig <- as.vector(GAFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GAFupsigGO.BP <- enrichGO(gene = genesGAFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GAFupsigGO.CC <- enrichGO(gene = genesGAFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFupsigGO.CC)
# 2 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GAFupsigGO.CC <- simplify(GAFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GAFupsigGO.CC)
# 2 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GAFupsigGO.CC, "GAFupsigGO.CC.csv",  row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GAFupsigGO.MF <- enrichGO(gene = genesGAFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# First collect the gene IDs in a new vector
genesGAFdownsig <- as.vector(GAFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GAFdownsigGO.BP <- enrichGO(gene = genesGAFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFdownsigGO.BP)
# 216 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GAFdownsigGO.BP <- simplify(GAFdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GAFdownsigGO.BP)
# 61 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GAFdownsigGO.BP, "GAFdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GAFdownsigGO.CC <- enrichGO(gene = genesGAFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFdownsigGO.CC)
# 57 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GAFdownsigGO.CC <- simplify(GAFdownsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GAFdownsigGO.CC)
# 18 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GAFdownsigGO.CC, "GAFdownsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GAFdownsigGO.MF <- enrichGO(gene = genesGAFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GAFdownsigGO.MF)
# 66 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GAFdownsigGO.MF <- simplify(GAFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GAFdownsigGO.MF)
# 24 MF GO terms

#Save the terms as in a data frame
write.csv(reduced.GAFdownsigGO.MF, "GAFdownsigGO.MF.csv", row.names = F)




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# First collect the gene IDs in a new vector
genesGBFupsig <- as.vector(GBFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GBFupsigGO.BP <- enrichGO(gene = genesGBFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFupsigGO.BP)
# 32 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFupsigGO.BP <- simplify(GBFupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFupsigGO.BP)
# 12 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GBFupsigGO.BP, "GBFupsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GBFupsigGO.CC <- enrichGO(gene = genesGBFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFupsigGO.CC)
# 35 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFupsigGO.CC <- simplify(GBFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFupsigGO.CC)
# 12 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GBFupsigGO.CC, "GBFupsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GBFupsigGO.MF <- enrichGO(gene = genesGBFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFupsigGO.MF)
# 4 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFupsigGO.MF <- simplify(GBFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFupsigGO.MF)
# 4 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GBFupsigGO.MF, "GBFupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# First collect the gene IDs in a new vector
genesGBFdownsig <- as.vector(GBFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GBFdownsigGO.BP <- enrichGO(gene = genesGBFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFdownsigGO.BP)
# 15 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFdownsigGO.BP <- simplify(GBFdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFdownsigGO.BP)
# 2 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GBFdownsigGO.BP, "GBFdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GBFdownsigGO.CC <- enrichGO(gene = genesGBFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFdownsigGO.CC)
# 9 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFdownsigGO.CC <- simplify(GBFdownsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFdownsigGO.CC)
# 3 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GBFdownsigGO.CC, "GBFdownsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GBFdownsigGO.MF <- enrichGO(gene = genesGBFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GBFdownsigGO.MF)
# 1 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GBFdownsigGO.MF <- simplify(GBFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBFdownsigGO.MF)
# 1 MF GO terms

#Save the terms as in a data frame
write.csv(reduced.GBFdownsigGO.MF, "GBFdownsigGO.MF.csv", row.names = F)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# First collect the gene IDs in a new vector
genesGCFupsig <- as.vector(GCFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GCFupsigGO.BP <- enrichGO(gene = genesGCFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GCFupsigGO.CC <- enrichGO(gene = genesGCFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFupsigGO.CC)
# 6 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GCFupsigGO.CC <- simplify(GCFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GCFupsigGO.CC)
# 4 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GCFupsigGO.CC, "GCFupsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GCFupsigGO.MF <- enrichGO(gene = genesGCFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# First collect the gene IDs in a new vector
genesGCFdownsig <- as.vector(GCFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GCFdownsigGO.BP <- enrichGO(gene = genesGCFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFdownsigGO.BP)
# 2 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GCFdownsigGO.BP <- simplify(GCFdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GCFdownsigGO.BP)
# 1 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GCFdownsigGO.BP, "GCFdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GCFdownsigGO.CC <- enrichGO(gene = genesGCFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GCFdownsigGO.MF <- enrichGO(gene = genesGCFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GCFdownsigGO.MF)
# 8 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GCFdownsigGO.MF <- simplify(GCFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GCFdownsigGO.MF)
# 3 MF GO terms

#Save the terms as in a data frame
write.csv(reduced.GCFdownsigGO.MF, "GCFdownsigGO.MF.csv", row.names = F)




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# First collect the gene IDs in a new vector
genesGDFupsig <- as.vector(GDFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GDFupsigGO.BP <- enrichGO(gene = genesGDFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFupsigGO.BP)
# 3 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GDFupsigGO.BP <- simplify(GDFupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GDFupsigGO.BP)
# 2 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GDFupsigGO.BP, "GDFupsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GDFupsigGO.CC <- enrichGO(gene = genesGDFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFupsigGO.CC)
# 14 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GDFupsigGO.CC <- simplify(GDFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GDFupsigGO.CC)
# 4 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GDFupsigGO.CC, "GDFupsigGO.CC.csv",  row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GDFupsigGO.MF <- enrichGO(gene = genesGDFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFupsigGO.MF)
# 1 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GDFupsigGO.MF <- simplify(GDFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GDFupsigGO.MF)
# 1 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GDFupsigGO.MF, "GDFupsigGO.MF.csv",  row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# First collect the gene IDs in a new vector
genesGDFdownsig <- as.vector(GDFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GDFdownsigGO.BP <- enrichGO(gene = genesGDFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GDFdownsigGO.CC <- enrichGO(gene = genesGDFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GDFdownsigGO.MF <- enrichGO(gene = genesGDFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GDFdownsigGO.MF)
# 2 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GDFdownsigGO.MF <- simplify(GDFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GDFdownsigGO.MF)
# 1 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GDFdownsigGO.MF, "GDFdownsigGO.MF.csv", row.names = F)




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# First collect the gene IDs in a new vector
genesGEFupsig <- as.vector(GEFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GEFupsigGO.BP <- enrichGO(gene = genesGEFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFupsigGO.BP)
# 4 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GEFupsigGO.BP <- simplify(GEFupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GEFupsigGO.BP)
# 4 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GEFupsigGO.BP, "GEFupsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GEFupsigGO.CC <- enrichGO(gene = genesGEFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFupsigGO.CC)
# 1 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GEFupsigGO.CC <- simplify(GEFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GEFupsigGO.CC)
# 1 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GEFupsigGO.CC, "GEFupsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GEFupsigGO.MF <- enrichGO(gene = genesGEFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFupsigGO.MF)
# 2 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GEFupsigGO.MF <- simplify(GEFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GEFupsigGO.MF)
# 1 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GEFupsigGO.MF, "GEFupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# First collect the gene IDs in a new vector
genesGEFdownsig <- as.vector(GEFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GEFdownsigGO.BP <- enrichGO(gene = genesGEFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GEFdownsigGO.CC <- enrichGO(gene = genesGEFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GEFdownsigGO.MF <- enrichGO(gene = genesGEFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GEFdownsigGO.MF)
# 0 MF GO terms




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# First collect the gene IDs in a new vector
genesGFFupsig <- as.vector(GFFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GFFupsigGO.BP <- enrichGO(gene = genesGFFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GFFupsigGO.CC <- enrichGO(gene = genesGFFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFupsigGO.CC)
# 12 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GFFupsigGO.CC <- simplify(GFFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFFupsigGO.CC)
# 4 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GFFupsigGO.CC, "GFFupsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GFFupsigGO.MF <- enrichGO(gene = genesGFFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFupsigGO.MF)
# 6 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GFFupsigGO.MF <- simplify(GFFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFFupsigGO.MF)
# 3 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GFFupsigGO.MF, "GFFupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# First collect the gene IDs in a new vector
genesGFFdownsig <- as.vector(GFFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GFFdownsigGO.BP <- enrichGO(gene = genesGFFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFdownsigGO.BP)
# 14 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GFFdownsigGO.BP <- simplify(GFFdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFFdownsigGO.BP)
# 6 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GFFdownsigGO.BP, "GFFdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GFFdownsigGO.CC <- enrichGO(gene = genesGFFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GFFdownsigGO.MF <- enrichGO(gene = genesGFFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GFFdownsigGO.MF)
# 4 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GFFdownsigGO.MF <- simplify(GFFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFFdownsigGO.MF)
# 3 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GFFdownsigGO.MF, "GFFdownsigGO.MF.csv", row.names = F)




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# First collect the gene IDs in a new vector
genesGGFupsig <- as.vector(GGFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GGFupsigGO.BP <- enrichGO(gene = genesGGFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GGFupsigGO.CC <- enrichGO(gene = genesGGFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFupsigGO.CC)
# 3 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GGFupsigGO.CC <- simplify(GGFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GGFupsigGO.CC)
# 2 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GGFupsigGO.CC, "GGFupsigGO.CC.csv",  row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GGFupsigGO.MF <- enrichGO(gene = genesGGFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# First collect the gene IDs in a new vector
genesGGFdownsig <- as.vector(GGFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GGFdownsigGO.BP <- enrichGO(gene = genesGGFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GGFdownsigGO.CC <- enrichGO(gene = genesGGFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFdownsigGO.CC)
# 5 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GGFdownsigGO.CC <- simplify(GGFdownsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GGFdownsigGO.CC)
# 3 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GGFdownsigGO.CC, "GGFdownsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GGFdownsigGO.MF <- enrichGO(gene = genesGGFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(GGFdownsigGO.MF)
# 0 MF GO terms

