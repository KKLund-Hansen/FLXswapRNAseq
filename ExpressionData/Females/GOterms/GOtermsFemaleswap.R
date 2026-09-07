################################################################################################
######################### GENE ONTOLOGY TERM ANALYSES FLX SWAP FEMALES #########################
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
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# First collect the gene IDs in a new vector
genesFvCfupsig <- as.vector(FvCfupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FvCfupsigGO.BP <- enrichGO(gene = genesFvCfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FvCfupsigGO.CC <- enrichGO(gene = genesFvCfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FvCfupsigGO.MF <- enrichGO(gene = genesFvCfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfupsigGO.MF)
# 0 MF GO terms



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# First collect the gene IDs in a new vector
genesFvCfdownsig <- as.vector(FvCfdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FvCfdownsigGO.BP <- enrichGO(gene = genesFvCfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FvCfdownsigGO.CC <- enrichGO(gene = genesFvCfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FvCfdownsigGO.MF <- enrichGO(gene = genesFvCfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FvCfdownsigGO.MF)
# 0 MF GO terms




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# First collect the gene IDs in a new vector
genesAvSfupsig <- as.vector(AvSfupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
AvSfupsigGO.BP <- enrichGO(gene = genesAvSfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfupsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
AvSfupsigGO.CC <- enrichGO(gene = genesAvSfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
AvSfupsigGO.MF <- enrichGO(gene = genesAvSfupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                           ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfupsigGO.MF)
# 11 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.AvSfupsigGO.MF <- simplify(AvSfupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.AvSfupsigGO.MF)
# 7 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.AvSfupsigGO.MF, "AvSfupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# First collect the gene IDs in a new vector
genesAvSfdownsig <- as.vector(AvSfdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
AvSfdownsigGO.BP <- enrichGO(gene = genesAvSfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
AvSfdownsigGO.CC <- enrichGO(gene = genesAvSfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
AvSfdownsigGO.MF <- enrichGO(gene = genesAvSfdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(AvSfdownsigGO.MF)
# 1 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.AvSfdownsigGO.MF <- simplify(AvSfdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.AvSfdownsigGO.MF)
# 1 MF GO terms

#Save the terms as in a data frame
write.csv(reduced.AvSfdownsigGO.MF, "AvSfdownsigGO.MF.csv", row.names = F)




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# First collect the gene IDs in a new vector
genesFsvCsFupsig <- as.vector(FsvCsFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FsvCsFupsigGO.BP <- enrichGO(gene = genesFsvCsFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE", 
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFupsigGO.BP)
# 14 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FsvCsFupsigGO.BP<- simplify(FsvCsFupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FsvCsFupsigGO.BP)
# 7 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FsvCsFupsigGO.BP, "FsvCsFupsigGO.BP.csv", row.names = F)



### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FsvCsFupsigGO.CC <- enrichGO(gene = genesFsvCsFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFupsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FsvCsFupsigGO.MF <- enrichGO(gene = genesFsvCsFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFupsigGO.MF)
# 5 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.FsvCsFupsigGO.MF <- simplify(FsvCsFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FsvCsFupsigGO.MF)
# 3 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.FsvCsFupsigGO.MF, "FsvCsFupsigGO.MF.csv",  row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# First collect the gene IDs in a new vector
genesFsvCsFdownsig <- as.vector(FsvCsFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FsvCsFdownsigGO.BP <- enrichGO(gene = genesFsvCsFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFdownsigGO.BP)
# 0 BP GO terms


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FsvCsFdownsigGO.CC <- enrichGO(gene = genesFsvCsFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFdownsigGO.CC)
# 0 CC GO terms


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FsvCsFdownsigGO.MF <- enrichGO(gene = genesFsvCsFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FsvCsFdownsigGO.MF)
# 0 MF GO terms




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# First collect the gene IDs in a new vector
genesFavCaFupsig <- as.vector(FavCaFupsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FavCaFupsigGO.BP <- enrichGO(gene = genesFavCaFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFupsigGO.BP)
# 8 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFupsigGO.BP <- simplify(FavCaFupsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFupsigGO.BP)
# 6 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFupsigGO.BP, "FavCaFupsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FavCaFupsigGO.CC <- enrichGO(gene = genesFavCaFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFupsigGO.CC)
# 26 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFupsigGO.CC <- simplify(FavCaFupsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFupsigGO.CC)
# 13 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFupsigGO.CC, "FavCaFupsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FavCaFupsigGO.MF <- enrichGO(gene = genesFavCaFupsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFupsigGO.MF)
# 5 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFupsigGO.MF <- simplify(FavCaFupsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFupsigGO.MF)
# 5 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFupsigGO.MF, "FavCaFupsigGO.MF.csv", row.names = F)



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify GO terms for down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# First collect the gene IDs in a new vector
genesFavCaFdownsig <- as.vector(FavCaFdownsig$gene.id)


### Biological Process ###

# Significantly enriched Gene Ontology terms for Biological Process
FavCaFdownsigGO.BP <- enrichGO(gene = genesFavCaFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFdownsigGO.BP)
# 252 BP GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFdownsigGO.BP <- simplify(FavCaFdownsigGO.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFdownsigGO.BP)
# 77 BP GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFdownsigGO.BP, "FavCaFdownsigGO.BP.csv", row.names = F)


### Cellular Component ###

# Significantly enriched Gene Ontology terms for Cellular Component
FavCaFdownsigGO.CC <- enrichGO(gene = genesFavCaFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFdownsigGO.CC)
# 60 CC GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFdownsigGO.CC <- simplify(FavCaFdownsigGO.CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFdownsigGO.CC)
# 21 CC GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFdownsigGO.CC, "FavCaFdownsigGO.CC.csv", row.names = F)


### Molecular Function ###

# Significantly enriched Gene Ontology terms for Molecular Function
FavCaFdownsigGO.MF <- enrichGO(gene = genesFavCaFdownsig, OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
                               ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)
# View the data
dim(FavCaFdownsigGO.MF)
# 102 MF GO terms

# Reduce redundancy of enriched GO terms
reduced.FavCaFdownsigGO.MF <- simplify(FavCaFdownsigGO.MF, cutoff = 0.7, by = "p.adjust", select_fun = min)

# View the data
dim(reduced.FavCaFdownsigGO.MF)
# 40 MF GO terms

# Save the terms as in a data frame
write.csv(reduced.FavCaFdownsigGO.MF, "FavCaFdownsigGO.MF.csv",  row.names = F)


