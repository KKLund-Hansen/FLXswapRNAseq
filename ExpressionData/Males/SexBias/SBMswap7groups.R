################################################################################################
########## SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP MALES - 7 GROUPS ##########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexBias")

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)
sebida.database <- read.table(file = "SEBIDA.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# Subset the database to only include transcripts with an identified biased class
sebida.data <- subset(sebida.database, MetaClass != "Unclassified")
# 12688

# To be able to identify the different biased classes first subset into the three biased classes:
# Female
sebida.febi <- subset(sebida.data, MetaClass == "Female")
# 4744

# Male
sebida.mabi <- subset(sebida.data, MetaClass == "Male")
# 3684

# Unbiased
sebida.unbi <- subset(sebida.data, MetaClass == "Unbiased")
# 4260 


# Next calculate the proportion of transcripts in either class
# Female
prop.female <- nrow(sebida.febi)/nrow(sebida.data)
# 0.3738966

# Male
prop.male <- nrow(sebida.mabi)/nrow(sebida.data)
# 0.2903531

# Unbiased
prop.unbi <- nrow(sebida.unbi)/nrow(sebida.data)
# 0.3357503



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group A
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# Subset the data to only include transcripts that overlap between group A up-regulated and sebida database
GAMupsigBias <- subset(GAMupsig, GAMupsig$gene.id %in% sebida.data$gene.id)
# 2

# Get the sum for the transcripts which are significant for three bias classes
GAMupsigBiasF <- sum(GAMupsigBias$gene.id %in% sebida.febi$gene.id)
# 0
GAMupsigBiasM <- sum(GAMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GAMupsigBiasU <- sum(GAMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 2


# Last create a new data frame with the observed and expected numbers
sexbi.GAMup <- rbind( c(GAMupsigBiasF, GAMupsigBiasM, GAMupsigBiasU), 
                      c(prop.female * nrow(GAMupsigBias), prop.male * nrow(GAMupsigBias), prop.unbi * nrow(GAMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GAMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GAMup["Obs",], p = sexbi.GAMup["Exp",], rescale.p = T)
# P = 0.1383

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GAMupF <- rbind( c(GAMupsigBiasF, GAMupsigBiasM, GAMupsigBiasU),
                       c(nrow(sebida.febi)-GAMupsigBiasF, nrow(sebida.mabi)-GAMupsigBiasM, nrow(sebida.unbi)-GAMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GAMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GAMupF)
# P = 0.197

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
# 16

# Subset the data to only include transcripts that overlap between group A down-regulated and sebida database
GAMdownsigBias <- subset(GAMdownsig, GAMdownsig$gene.id %in% sebida.data$gene.id)
# 15

# Get the sum for the transcripts which are significant for three bias classes
GAMdownsigBiasF <- sum(GAMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 6
GAMdownsigBiasM <- sum(GAMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 3
GAMdownsigBiasU <- sum(GAMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 6


# Last create a new data frame with the observed and expected numbers
sexbi.GAMdown <- rbind( c(GAMdownsigBiasF, GAMdownsigBiasM, GAMdownsigBiasU),
                        c(prop.female * nrow(GAMdownsigBias), prop.male * nrow(GAMdownsigBias), prop.unbi * nrow(GAMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GAMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GAMdown["Obs",], p = sexbi.GAMdown["Exp",], rescale.p = T)
# P = 0.7285

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GAMdownF <- rbind( c(GAMdownsigBiasF, GAMdownsigBiasM, GAMdownsigBiasU),
                         c(nrow(sebida.febi)-GAMdownsigBiasF, nrow(sebida.mabi)-GAMdownsigBiasM, nrow(sebida.unbi)-GAMdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GAMdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GAMdownF)
# P = 0.7652

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group B
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# Subset the data to only include transcripts that overlap between group B up-regulated and sebida database
GBMupsigBias <- subset(GBMupsig, GBMupsig$gene.id %in% sebida.data$gene.id)
# 32

# Get the sum for the transcripts which are significant for three bias classes
GBMupsigBiasF <- sum(GBMupsigBias$gene.id %in% sebida.febi$gene.id)
# 6
GBMupsigBiasM <- sum(GBMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 7
GBMupsigBiasU <- sum(GBMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 19


# Last create a new data frame with the observed and expected numbers
sexbi.GBMup <- rbind( c(GBMupsigBiasF, GBMupsigBiasM, GBMupsigBiasU), 
                      c(prop.female * nrow(GBMupsigBias), prop.male * nrow(GBMupsigBias), prop.unbi * nrow(GBMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GBMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GBMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GBMup["Obs",], p = sexbi.GBMup["Exp",], rescale.p = T)
# P = 0.007145



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# Subset the data to only include transcripts that overlap between group B down-regulated and sebida database
GBMdownsigBias <- subset(GBMdownsig, GBMdownsig$gene.id %in% sebida.data$gene.id)
# 51

# Get the sum for the transcripts which are significant for three bias classes
GBMdownsigBiasF <- sum(GBMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 11
GBMdownsigBiasM <- sum(GBMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 2
GBMdownsigBiasU <- sum(GBMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 38


# Last create a new data frame with the observed and expected numbers
sexbi.GBMdown <- rbind( c(GBMdownsigBiasF, GBMdownsigBiasM, GBMdownsigBiasU),
                        c(prop.female * nrow(GBMdownsigBias), prop.male * nrow(GBMdownsigBias), prop.unbi * nrow(GBMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GBMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GBMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GBMdown["Obs",], p = sexbi.GBMdown["Exp",], rescale.p = T)
# P = 2.118e-09




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
# 27

# Subset the data to only include transcripts that overlap between group C up-regulated and sebida database
GCMupsigBias <- subset(GCMupsig, GCMupsig$gene.id %in% sebida.data$gene.id)
# 20

# Get the sum for the transcripts which are significant for three bias classes
GCMupsigBiasF <- sum(GCMupsigBias$gene.id %in% sebida.febi$gene.id)
# 8
GCMupsigBiasM <- sum(GCMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 3
GCMupsigBiasU <- sum(GCMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 9


# Last create a new data frame with the observed and expected numbers
sexbi.GCMup <- rbind( c(GCMupsigBiasF, GCMupsigBiasM, GCMupsigBiasU), 
                      c(prop.female * nrow(GCMupsigBias), prop.male * nrow(GCMupsigBias), prop.unbi * nrow(GCMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GCMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GCMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GCMup["Obs",], p = sexbi.GCMup["Exp",], rescale.p = T)
# P = 0.3378



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# Subset the data to only include transcripts that overlap between group C down-regulated and sebida database
GCMdownsigBias <- subset(GCMdownsig, GCMdownsig$gene.id %in% sebida.data$gene.id)
# 35

# Get the sum for the transcripts which are significant for three bias classes
GCMdownsigBiasF <- sum(GCMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 10
GCMdownsigBiasM <- sum(GCMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 4
GCMdownsigBiasU <- sum(GCMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 21


# Last create a new data frame with the observed and expected numbers
sexbi.GCMdown <- rbind( c(GCMdownsigBiasF, GCMdownsigBiasM, GCMdownsigBiasU),
                        c(prop.female * nrow(GCMdownsigBias), prop.male * nrow(GCMdownsigBias), prop.unbi * nrow(GCMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GCMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GCMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GCMdown["Obs",], p = sexbi.GCMdown["Exp",], rescale.p = T)
# P = 0.002817




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# Subset the data to only include transcripts that overlap between group D up-regulated and sebida database
GDMupsigBias <- subset(GDMupsig, GDMupsig$gene.id %in% sebida.data$gene.id)
# 10

# Get the sum for the transcripts which are significant for three bias classes
GDMupsigBiasF <- sum(GDMupsigBias$gene.id %in% sebida.febi$gene.id)
# 3
GDMupsigBiasM <- sum(GDMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GDMupsigBiasU <- sum(GDMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 7


# Last create a new data frame with the observed and expected numbers
sexbi.GDMup <- rbind( c(GDMupsigBiasF, GDMupsigBiasM, GDMupsigBiasU), 
                      c(prop.female * nrow(GDMupsigBias), prop.male * nrow(GDMupsigBias), prop.unbi * nrow(GDMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GDMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GDMup["Obs",], p = sexbi.GDMup["Exp",], rescale.p = T)
# P = 0.03018

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GDMupF <- rbind( c(GDMupsigBiasF, GDMupsigBiasM, GDMupsigBiasU),
                         c(nrow(sebida.febi)-GDMupsigBiasF, nrow(sebida.mabi)-GDMupsigBiasM, nrow(sebida.unbi)-GDMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GDMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GDMupF)
# P = 0.0324

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
# 16

# Subset the data to only include transcripts that overlap between group D down-regulated and sebida database
GDMdownsigBias <- subset(GDMdownsig, GDMdownsig$gene.id %in% sebida.data$gene.id)
# 16

# Get the sum for the transcripts which are significant for three bias classes
GDMdownsigBiasF <- sum(GDMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
GDMdownsigBiasM <- sum(GDMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GDMdownsigBiasU <- sum(GDMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 11


# Last create a new data frame with the observed and expected numbers
sexbi.GDMdown <- rbind( c(GDMdownsigBiasF, GDMdownsigBiasM, GDMdownsigBiasU),
                        c(prop.female * nrow(GDMdownsigBias), prop.male * nrow(GDMdownsigBias), prop.unbi * nrow(GDMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GDMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GDMdown["Obs",], p = sexbi.GDMdown["Exp",], rescale.p = T)
# P = 0.009032

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GDMdownF <- rbind( c(GDMdownsigBiasF, GDMdownsigBiasM, GDMdownsigBiasU),
                         c(nrow(sebida.febi)-GDMdownsigBiasF, nrow(sebida.mabi)-GDMdownsigBiasM, nrow(sebida.unbi)-GDMdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GDMdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GDMdownF)
# P = 0.01088

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group E
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# Subset the data to only include transcripts that overlap between group E up-regulated and sebida database
GEMupsigBias <- subset(GEMupsig, GEMupsig$gene.id %in% sebida.data$gene.id)
# 12

# Get the sum for the transcripts which are significant for three bias classes
GEMupsigBiasF <- sum(GEMupsigBias$gene.id %in% sebida.febi$gene.id)
# 1
GEMupsigBiasM <- sum(GEMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 4
GEMupsigBiasU <- sum(GEMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 7


# Last create a new data frame with the observed and expected numbers
sexbi.GEMup <- rbind( c(GEMupsigBiasF, GEMupsigBiasM, GEMupsigBiasU), 
                      c(prop.female * nrow(GEMupsigBias), prop.male * nrow(GEMupsigBias), prop.unbi * nrow(GEMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEMup["Obs",], p = sexbi.GEMup["Exp",], rescale.p = T)
# P = 0.08304

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GEMupF <- rbind( c(GEMupsigBiasF, GEMupsigBiasM, GEMupsigBiasU),
                       c(nrow(sebida.febi)-GEMupsigBiasF, nrow(sebida.mabi)-GEMupsigBiasM, nrow(sebida.unbi)-GEMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GEMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GEMupF)
# P = 0.06354

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# Subset the data to only include transcripts that overlap between group E down-regulated and sebida database
GEMdownsigBias <- subset(GEMdownsig, GEMdownsig$gene.id %in% sebida.data$gene.id)
# 22

# Get the sum for the transcripts which are significant for three bias classes
GEMdownsigBiasF <- sum(GEMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 8
GEMdownsigBiasM <- sum(GEMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GEMdownsigBiasU <- sum(GEMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 13


# Last create a new data frame with the observed and expected numbers
sexbi.GEMdown <- rbind( c(GEMdownsigBiasF, GEMdownsigBiasM, GEMdownsigBiasU),
                        c(prop.female * nrow(GEMdownsigBias), prop.male * nrow(GEMdownsigBias), prop.unbi * nrow(GEMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEMdown["Obs",], p = sexbi.GEMdown["Exp",], rescale.p = T)
# P = 0.01218




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
# 2

# Subset the data to only include transcripts that overlap between group F up-regulated and sebida database
GFMupsigBias <- subset(GFMupsig, GFMupsig$gene.id %in% sebida.data$gene.id)
# 1

# Get the sum for the transcripts which are significant for three bias classes
GFMupsigBiasF <- sum(GFMupsigBias$gene.id %in% sebida.febi$gene.id)
# 0
GFMupsigBiasM <- sum(GFMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GFMupsigBiasU <- sum(GFMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 1


# Last create a new data frame with the observed and expected numbers
sexbi.GFMup <- rbind( c(GFMupsigBiasF, GFMupsigBiasM, GFMupsigBiasU), 
                      c(prop.female * nrow(GFMupsigBias), prop.male * nrow(GFMupsigBias), prop.unbi * nrow(GFMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GFMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GFMup["Obs",], p = sexbi.GFMup["Exp",], rescale.p = T)
# P = 0.3719

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GFMupF <- rbind( c(GFMupsigBiasF, GFMupsigBiasM, GFMupsigBiasU),
                       c(nrow(sebida.febi)-GFMupsigBiasF, nrow(sebida.mabi)-GFMupsigBiasM, nrow(sebida.unbi)-GFMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GFMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GFMupF)
# P = 0.6261

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group F
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
# 3

# Subset the data to only include transcripts that overlap between group F down-regulated and sebida database
GFMdownsigBias <- subset(GFMdownsig, GFMdownsig$gene.id %in% sebida.data$gene.id)
# 3

# Get the sum for the transcripts which are significant for three bias classes
GFMdownsigBiasF <- sum(GFMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 0
GFMdownsigBiasM <- sum(GFMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GFMdownsigBiasU <- sum(GFMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 3


# Last create a new data frame with the observed and expected numbers
sexbi.GFMdown <- rbind( c(GFMdownsigBiasF, GFMdownsigBiasM, GFMdownsigBiasU),
                        c(prop.female * nrow(GFMdownsigBias), prop.male * nrow(GFMdownsigBias), prop.unbi * nrow(GFMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GFMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GFMdown["Obs",], p = sexbi.GFMdown["Exp",], rescale.p = T)
# P = 0.05143

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GFMdownF <- rbind( c(GFMdownsigBiasF, GFMdownsigBiasM, GFMdownsigBiasU),
                         c(nrow(sebida.febi)-GFMdownsigBiasF, nrow(sebida.mabi)-GFMdownsigBiasM, nrow(sebida.unbi)-GFMdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GFMdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GFMdownF)
# P = 0.06229

# We get the same results, so that is good



######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group G
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
# 10

# Subset the data to only include transcripts that overlap between group G up-regulated and sebida database
GGMupsigBias <- subset(GGMupsig, GGMupsig$gene.id %in% sebida.data$gene.id)
# 5

# Get the sum for the transcripts which are significant for three bias classes
GGMupsigBiasF <- sum(GGMupsigBias$gene.id %in% sebida.febi$gene.id)
# 2
GGMupsigBiasM <- sum(GGMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GGMupsigBiasU <- sum(GGMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 3


# Last create a new data frame with the observed and expected numbers
sexbi.GGMup <- rbind( c(GGMupsigBiasF, GGMupsigBiasM, GGMupsigBiasU), 
                      c(prop.female * nrow(GGMupsigBias), prop.male * nrow(GGMupsigBias), prop.unbi * nrow(GGMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GGMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GGMup["Obs",], p = sexbi.GGMup["Exp",], rescale.p = T)
# P = 0.2864

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GGMupF <- rbind( c(GGMupsigBiasF, GGMupsigBiasM, GGMupsigBiasU),
                       c(nrow(sebida.febi)-GGMupsigBiasF, nrow(sebida.mabi)-GGMupsigBiasM, nrow(sebida.unbi)-GGMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GGMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GGMupF)
# P = 0.3328

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group G
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# Subset the data to only include transcripts that overlap between group G down-regulated and sebida database
GGMdownsigBias <- subset(GGMdownsig, GGMdownsig$gene.id %in% sebida.data$gene.id)
# 20

# Get the sum for the transcripts which are significant for three bias classes
GGMdownsigBiasF <- sum(GGMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 5
GGMdownsigBiasM <- sum(GGMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 2
GGMdownsigBiasU <- sum(GGMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 13


# Last create a new data frame with the observed and expected numbers
sexbi.GGMdown <- rbind( c(GGMdownsigBiasF, GGMdownsigBiasM, GGMdownsigBiasU),
                        c(prop.female * nrow(GGMdownsigBias), prop.male * nrow(GGMdownsigBias), prop.unbi * nrow(GGMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GGMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GGMdown["Obs",], p = sexbi.GGMdown["Exp",], rescale.p = T)
# P = 0.01005

