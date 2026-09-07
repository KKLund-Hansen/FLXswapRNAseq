################################################################################################
########################## GENE ONTOLOGY TERM ANALYSES FLX SWAP MALES ##########################
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
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# First collect the gene IDs in a new vector
genesFvCmupsig <- as.vector(FvCmupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FvCmupsigGO.BP <- enrichGO(gene = genesFvCmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCmupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FvCmupsigGO.CC <- enrichGO(gene = genesFvCmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCmupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FvCmupsigGO.MF <- enrichGO(gene = genesFvCmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCmupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# First collect the gene IDs in a new vector
genesFvCmdownsig <- as.vector(FvCmdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FvCmdownsigGO.BP <- enrichGO(gene = genesFvCmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FvCmdownsigGO.BP)
# 1 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FvCmdownsigGO.BP <- simplify(FvCmdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FvCmdownsigGO.BP)
# 1 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FvCmdownsigGO.BP, "FvCmdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FvCmdownsigGO.CC <- enrichGO(gene = genesFvCmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FvCmdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FvCmdownsigGO.MF <- enrichGO(gene = genesFvCmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FvCmdownsigGO.MF)
# 0 MF GO terms




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# First collect the gene IDs in a new vector
genesAvSmupsig <- as.vector(AvSmupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
AvSmupsigGO.BP <- enrichGO(gene = genesAvSmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmupsigGO.BP)
# 11 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.AvSmupsigGO.BP <- simplify(AvSmupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.AvSmupsigGO.BP)
# 5 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.AvSmupsigGO.BP, "AvSmupsigGO.BP.csv",  row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
AvSmupsigGO.CC <- enrichGO(gene = genesAvSmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
AvSmupsigGO.MF <- enrichGO(gene = genesAvSmupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# First collect the gene IDs in a new vector
genesAvSmdownsig <- as.vector(AvSmdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
AvSmdownsigGO.BP <- enrichGO(gene = genesAvSmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
AvSmdownsigGO.CC <- enrichGO(gene = genesAvSmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
AvSmdownsigGO.MF <- enrichGO(gene = genesAvSmdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(AvSmdownsigGO.MF)
# 0 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.AvSmdownsigGO.MF <- simplify(AvSmdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.AvSmdownsigGO.MF)
# 1 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.AvSmdownsigGO.MF, "AvSmdownsigGO.MF.csv",  row.names = F)




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# First collect the gene IDs in a new vector
genesFsvCsMupsig <- as.vector(FsvCsMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FsvCsMupsigGO.BP <- enrichGO(gene = genesFsvCsMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMupsigGO.BP)
# 3 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FsvCsMupsigGO.BP<- simplify(FsvCsMupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FsvCsMupsigGO.BP)
# 2 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FsvCsMupsigGO.BP, "FsvCsMupsigGO.BP.csv",  row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FsvCsMupsigGO.CC <- enrichGO(gene = genesFsvCsMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FsvCsMupsigGO.MF <- enrichGO(gene = genesFsvCsMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMupsigGO.MF)
# 6 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.FsvCsMupsigGO.MF <- simplify(FsvCsMupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FsvCsMupsigGO.MF)
# 3 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.FsvCsMupsigGO.MF, "FsvCsMupsigGO.MF.csv",  row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# First collect the gene IDs in a new vector
genesFsvCsMdownsig <- as.vector(FsvCsMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FsvCsMdownsigGO.BP <- enrichGO(gene = genesFsvCsMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMdownsigGO.BP)
# 2 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FsvCsMdownsigGO.BP <- simplify(FsvCsMdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FsvCsMdownsigGO.BP)
# 2 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FsvCsMdownsigGO.BP, "FsvCsMdownsigGO.BP.csv",  row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FsvCsMdownsigGO.CC <- enrichGO(gene = genesFsvCsMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FsvCsMdownsigGO.MF <- enrichGO(gene = genesFsvCsMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FsvCsMdownsigGO.MF)
# 0 MF GO terms




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65


# First collect the gene IDs in a new vector
genesFavCaMupsig <- as.vector(FavCaMupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FavCaMupsigGO.BP <- enrichGO(gene = genesFavCaMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FavCaMupsigGO.CC <- enrichGO(gene = genesFavCaMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FavCaMupsigGO.MF <- enrichGO(gene = genesFavCaMupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# First collect the gene IDs in a new vector
genesFavCaMdownsig <- as.vector(FavCaMdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FavCaMdownsigGO.BP <- enrichGO(gene = genesFavCaMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMdownsigGO.BP)
# 1 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaMdownsigGO.BP <- simplify(FavCaMdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaMdownsigGO.BP)
# 1 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaMdownsigGO.BP, "FavCaMdownsigGO.BP.csv",  row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FavCaMdownsigGO.CC <- enrichGO(gene = genesFavCaMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FavCaMdownsigGO.MF <- enrichGO(gene = genesFavCaMdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
# View the data
dim(FavCaMdownsigGO.MF)
# 12 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaMdownsigGO.MF <- simplify(FavCaMdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaMdownsigGO.MF)
# 2 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaMdownsigGO.MF, "FavCaMdownsigGO.MF.csv",  row.names = F)


