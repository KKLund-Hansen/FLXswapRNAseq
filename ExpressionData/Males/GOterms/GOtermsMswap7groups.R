################################################################################################
##################### GENE ONTOLOGY TERM ANALYSES FLX SWAP MALES - 7 GROUPS ####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/GOterms")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(clusterProfiler)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# First collect the gene IDs in a new vector
genesGAMupsig <- as.vector(GAMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GAMupsigGO.BP <- enrichGO(gene = genesGAMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GAMupsigGO.CC <- enrichGO(gene = genesGAMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GAMupsigGO.MF <- enrichGO(gene = genesGAMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMupsigGO.MF)
# 2 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GAMupsigGO.MF <- simplify(GAMupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GAMupsigGO.MF)
# 1 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GAMupsigGO.MF, "GAMupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# First collect the gene IDs in a new vector
genesGAMdownsig <- as.vector(GAMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GAMdownsigGO.BP <- enrichGO(gene = genesGAMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GAMdownsigGO.CC <- enrichGO(gene = genesGAMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GAMdownsigGO.MF <- enrichGO(gene = genesGAMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GAMdownsigGO.MF)
# 0 MF GO terms




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# First collect the gene IDs in a new vector
genesGBMupsig <- as.vector(GBMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GBMupsigGO.BP <- enrichGO(gene = genesGBMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMupsigGO.BP)
# 6 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GBMupsigGO.BP <- simplify(GBMupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBMupsigGO.BP)
# 4 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GBMupsigGO.BP, "GBMupsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GBMupsigGO.CC <- enrichGO(gene = genesGBMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GBMupsigGO.MF <- enrichGO(gene = genesGBMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMupsigGO.MF)
# 5 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GBMupsigGO.MF <- simplify(GBMupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBMupsigGO.MF)
# 3 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GBMupsigGO.MF, "GBMupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# First collect the gene IDs in a new vector
genesGBMdownsig <- as.vector(GBMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GBMdownsigGO.BP <- enrichGO(gene = genesGBMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMdownsigGO.BP)
# 4 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GBMdownsigGO.BP <- simplify(GBMdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBMdownsigGO.BP)
# 3 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GBMdownsigGO.BP, "GBMdownsigGO.BP.csv",  row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GBMdownsigGO.CC <- enrichGO(gene = genesGBMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GBMdownsigGO.MF <- enrichGO(gene = genesGBMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GBMdownsigGO.MF)
# 1 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GBMdownsigGO.MF <- simplify(GBMdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GBMdownsigGO.MF)
# 1 MF GO terms

#Save the terms as in a data frame
write.csv(reduced.GBMdownsigGO.MF, "GBMdownsigGO.MF.csv",  row.names = F)




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# First collect the gene IDs in a new vector
genesGCMupsig <- as.vector(GCMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GCMupsigGO.BP <- enrichGO(gene = genesGCMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GCMupsigGO.CC <- enrichGO(gene = genesGCMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GCMupsigGO.MF <- enrichGO(gene = genesGCMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMupsigGO.MF)
# 0 MF GO terms




##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# First collect the gene IDs in a new vector
genesGCMdownsig <- as.vector(GCMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GCMdownsigGO.BP <- enrichGO(gene = genesGCMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GCMdownsigGO.CC <- enrichGO(gene = genesGCMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GCMdownsigGO.MF <- enrichGO(gene = genesGCMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GCMdownsigGO.MF)
# 0 MF GO terms




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# First collect the gene IDs in a new vector
genesGDMupsig <- as.vector(GDMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GDMupsigGO.BP <- enrichGO(gene = genesGDMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GDMupsigGO.CC <- enrichGO(gene = genesGDMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GDMupsigGO.MF <- enrichGO(gene = genesGDMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# First collect the gene IDs in a new vector
genesGDMdownsig <- as.vector(GDMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GDMdownsigGO.BP <- enrichGO(gene = genesGDMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GDMdownsigGO.CC <- enrichGO(gene = genesGDMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GDMdownsigGO.MF <- enrichGO(gene = genesGDMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GDMdownsigGO.MF)
# 3 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GDMdownsigGO.MF <- simplify(GDMdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GDMdownsigGO.MF)
# 1 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GDMdownsigGO.MF, "GDMdownsigGO.MF.csv", row.names = F)



######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# First collect the gene IDs in a new vector
genesGEMupsig <- as.vector(GEMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GEMupsigGO.BP <- enrichGO(gene = genesGEMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GEMupsigGO.CC <- enrichGO(gene = genesGEMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GEMupsigGO.MF <- enrichGO(gene = genesGEMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# First collect the gene IDs in a new vector
genesGEMdownsig <- as.vector(GEMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GEMdownsigGO.BP <- enrichGO(gene = genesGEMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GEMdownsigGO.CC <- enrichGO(gene = genesGEMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GEMdownsigGO.MF <- enrichGO(gene = genesGEMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GEMdownsigGO.MF)
# 0 MF GO terms




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# First collect the gene IDs in a new vector
genesGFMupsig <- as.vector(GFMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GFMupsigGO.BP <- enrichGO(gene = genesGFMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GFMupsigGO.CC <- enrichGO(gene = genesGFMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GFMupsigGO.MF <- enrichGO(gene = genesGFMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# First collect the gene IDs in a new vector
genesGFMdownsig <- as.vector(GFMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GFMdownsigGO.BP <- enrichGO(gene = genesGFMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMdownsigGO.BP)
# 9 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.GFMdownsigGO.BP <- simplify(GFMdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFMdownsigGO.BP)
# 6 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.GFMdownsigGO.BP, "GFMdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GFMdownsigGO.CC <- enrichGO(gene = genesGFMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMdownsigGO.CC)
# 1 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.GFMdownsigGO.CC <- simplify(GFMdownsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFMdownsigGO.CC)
# 1 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.GFMdownsigGO.CC, "GFMdownsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GFMdownsigGO.MF <- enrichGO(gene = genesGFMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GFMdownsigGO.MF)
# 8 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.GFMdownsigGO.MF <- simplify(GFMdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.GFMdownsigGO.MF)
# 4 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.GFMdownsigGO.MF, "GFMdownsigGO.MF.csv",  row.names = F)




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# First collect the gene IDs in a new vector
genesGGMupsig <- as.vector(GGMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GGMupsigGO.BP <- enrichGO(gene = genesGGMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GGMupsigGO.CC <- enrichGO(gene = genesGGMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GGMupsigGO.MF <- enrichGO(gene = genesGGMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# First collect the gene IDs in a new vector
genesGGMdownsig <- as.vector(GGMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
GGMdownsigGO.BP <- enrichGO(gene = genesGGMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
GGMdownsigGO.CC <- enrichGO(gene = genesGGMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
GGMdownsigGO.MF <- enrichGO(gene = genesGGMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                            ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(GGMdownsigGO.MF)
# 0 MF GO terms

