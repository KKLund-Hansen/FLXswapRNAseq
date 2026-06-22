################################################################################################
################ SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP MALES ###############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Males/SexBias")

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 4

# Subset the data to only include transcripts that overlap between group A up-regulated and sebida database
GAMupsigBias <- subset(GAMupsig, GAMupsig$gene.id %in% sebida.data$gene.id)
# 1

# Get the sum for the transcripts which are significant for three bias classes
GAMupsigBiasF <- sum(GAMupsigBias$gene.id %in% sebida.febi$gene.id)
# 0
GAMupsigBiasM <- sum(GAMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GAMupsigBiasU <- sum(GAMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 1


# Last create a new data frame with the observed and expected numbers
sexbi.GAMup <- rbind( c(GAMupsigBiasF, GAMupsigBiasM, GAMupsigBiasU), 
                      c(prop.female * nrow(GAMupsigBias), prop.male * nrow(GAMupsigBias), prop.unbi * nrow(GAMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GAMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GAMup["Obs",], p = sexbi.GAMup["Exp",], rescale.p = T)
# P = 0.3719


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GAMupF <- rbind( c(GAMupsigBiasF, GAMupsigBiasM, GAMupsigBiasU),
                       c(nrow(sebida.febi)-GAMupsigBiasF, nrow(sebida.mabi)-GAMupsigBiasM, nrow(sebida.unbi)-GAMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GAMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GAMupF)
# P = 0.6261

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
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
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# Subset the data to only include transcripts that overlap between group B up-regulated and sebida database
GBMupsigBias <- subset(GBMupsig, GBMupsig$gene.id %in% sebida.data$gene.id)
# 34

# Get the sum for the transcripts which are significant for three bias classes
GBMupsigBiasF <- sum(GBMupsigBias$gene.id %in% sebida.febi$gene.id)
# 7
GBMupsigBiasM <- sum(GBMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 7
GBMupsigBiasU <- sum(GBMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 20


# Last create a new data frame with the observed and expected numbers
sexbi.GBMup <- rbind( c(GBMupsigBiasF, GBMupsigBiasM, GBMupsigBiasU), 
                      c(prop.female * nrow(GBMupsigBias), prop.male * nrow(GBMupsigBias), prop.unbi * nrow(GBMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GBMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GBMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GBMup["Obs",], p = sexbi.GBMup["Exp",], rescale.p = T)
# P = 0.007234



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# Subset the data to only include transcripts that overlap between group B down-regulated and sebida database
GBMdownsigBias <- subset(GBMdownsig, GBMdownsig$gene.id %in% sebida.data$gene.id)
# 53

# Get the sum for the transcripts which are significant for three bias classes
GBMdownsigBiasF <- sum(GBMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 12
GBMdownsigBiasM <- sum(GBMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 3
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
# P = 1.524e-08




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
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
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

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
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# Subset the data to only include transcripts that overlap between group D up-regulated and sebida database
GDMupsigBias <- subset(GDMupsig, GDMupsig$gene.id %in% sebida.data$gene.id)
# 12

# Get the sum for the transcripts which are significant for three bias classes
GDMupsigBiasF <- sum(GDMupsigBias$gene.id %in% sebida.febi$gene.id)
# 3
GDMupsigBiasM <- sum(GDMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GDMupsigBiasU <- sum(GDMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 9


# Last create a new data frame with the observed and expected numbers
sexbi.GDMup <- rbind( c(GDMupsigBiasF, GDMupsigBiasM, GDMupsigBiasU), 
                      c(prop.female * nrow(GDMupsigBias), prop.male * nrow(GDMupsigBias), prop.unbi * nrow(GDMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GDMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GDMup["Obs",], p = sexbi.GDMup["Exp",], rescale.p = T)
# P = 0.006377


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GDMupF <- rbind( c(GDMupsigBiasF, GDMupsigBiasM, GDMupsigBiasU),
                         c(nrow(sebida.febi)-GDMupsigBiasF, nrow(sebida.mabi)-GDMupsigBiasM, nrow(sebida.unbi)-GDMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GDMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GDMupF)
# P = 0.005592

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
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
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# Subset the data to only include transcripts that overlap between group E up-regulated and sebida database
GEMupsigBias <- subset(GEMupsig, GEMupsig$gene.id %in% sebida.data$gene.id)
# 11

# Get the sum for the transcripts which are significant for three bias classes
GEMupsigBiasF <- sum(GEMupsigBias$gene.id %in% sebida.febi$gene.id)
# 1
GEMupsigBiasM <- sum(GEMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 4
GEMupsigBiasU <- sum(GEMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 6


# Last create a new data frame with the observed and expected numbers
sexbi.GEMup <- rbind( c(GEMupsigBiasF, GEMupsigBiasM, GEMupsigBiasU), 
                      c(prop.female * nrow(GEMupsigBias), prop.male * nrow(GEMupsigBias), prop.unbi * nrow(GEMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEMup["Obs",], p = sexbi.GEMup["Exp",], rescale.p = T)
# P = 0.1353


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GEMupF <- rbind( c(GEMupsigBiasF, GEMupsigBiasM, GEMupsigBiasU),
                       c(nrow(sebida.febi)-GEMupsigBiasF, nrow(sebida.mabi)-GEMupsigBiasM, nrow(sebida.unbi)-GEMupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GEMupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GEMupF)
# P = 0.09549

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# Subset the data to only include transcripts that overlap between group E down-regulated and sebida database
GEMdownsigBias <- subset(GEMdownsig, GEMdownsig$gene.id %in% sebida.data$gene.id)
# 23

# Get the sum for the transcripts which are significant for three bias classes
GEMdownsigBiasF <- sum(GEMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 8
GEMdownsigBiasM <- sum(GEMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GEMdownsigBiasU <- sum(GEMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 14


# Last create a new data frame with the observed and expected numbers
sexbi.GEMdown <- rbind( c(GEMdownsigBiasF, GEMdownsigBiasM, GEMdownsigBiasU),
                        c(prop.female * nrow(GEMdownsigBias), prop.male * nrow(GEMdownsigBias), prop.unbi * nrow(GEMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEMdown["Obs",], p = sexbi.GEMdown["Exp",], rescale.p = T)
# P = 0.006829




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
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
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
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
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
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
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# Subset the data to only include transcripts that overlap between group G down-regulated and sebida database
GGMdownsigBias <- subset(GGMdownsig, GGMdownsig$gene.id %in% sebida.data$gene.id)
# 17

# Get the sum for the transcripts which are significant for three bias classes
GGMdownsigBiasF <- sum(GGMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
GGMdownsigBiasM <- sum(GGMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GGMdownsigBiasU <- sum(GGMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 12


# Last create a new data frame with the observed and expected numbers
sexbi.GGMdown <- rbind( c(GGMdownsigBiasF, GGMdownsigBiasM, GGMdownsigBiasU),
                        c(prop.female * nrow(GGMdownsigBias), prop.male * nrow(GGMdownsigBias), prop.unbi * nrow(GGMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GGMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GGMdown["Obs",], p = sexbi.GGMdown["Exp",], rescale.p = T)
# P = 0.004193


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.GGMdownF <- rbind( c(GGMdownsigBiasF, GGMdownsigBiasM, GGMdownsigBiasU),
                         c(nrow(sebida.febi)-GGMdownsigBiasF, nrow(sebida.mabi)-GGMdownsigBiasM, nrow(sebida.unbi)-GGMdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.GGMdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.GGMdownF)
# P = 0.00524

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# Subset the data to only include transcripts that overlap between group 2 up-regulated and sebida database
G2MupsigBias <- subset(G2Mupsig, G2Mupsig$gene.id %in% sebida.data$gene.id)
# 24

# Get the sum for the transcripts which are significant for three bias classes
G2MupsigBiasF <- sum(G2MupsigBias$gene.id %in% sebida.febi$gene.id)
# 3
G2MupsigBiasM <- sum(G2MupsigBias$gene.id %in% sebida.mabi$gene.id)
# 12
G2MupsigBiasU <- sum(G2MupsigBias$gene.id %in% sebida.unbi$gene.id)
# 99


# Last create a new data frame with the observed and expected numbers
sexbi.G2Mup <- rbind( c(G2MupsigBiasF, G2MupsigBiasM, G2MupsigBiasU), 
                      c(prop.female * nrow(G2MupsigBias), prop.male * nrow(G2MupsigBias), prop.unbi * nrow(G2MupsigBias)))
# Add names to the columns and rows
colnames(sexbi.G2Mup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G2Mup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G2Mup["Obs",], p = sexbi.G2Mup["Exp",], rescale.p = T)
# P = 0.02107



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# Subset the data to only include transcripts that overlap between group 2 down-regulated and sebida database
G2MdownsigBias <- subset(G2Mdownsig, G2Mdownsig$gene.id %in% sebida.data$gene.id)
# 29

# Get the sum for the transcripts which are significant for three bias classes
G2MdownsigBiasF <- sum(G2MdownsigBias$gene.id %in% sebida.febi$gene.id)
# 7
G2MdownsigBiasM <- sum(G2MdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 2
G2MdownsigBiasU <- sum(G2MdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 20


# Last create a new data frame with the observed and expected numbers
sexbi.G2Mdown <- rbind( c(G2MdownsigBiasF, G2MdownsigBiasM, G2MdownsigBiasU),
                        c(prop.female * nrow(G2MdownsigBias), prop.male * nrow(G2MdownsigBias), prop.unbi * nrow(G2MdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.G2Mdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G2Mdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G2Mdown["Obs",], p = sexbi.G2Mdown["Exp",], rescale.p = T)
# P = 0.0001959




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# Subset the data to only include transcripts that overlap between group 11 up-regulated and sebida database
G11MupsigBias <- subset(G11Mupsig, G11Mupsig$gene.id %in% sebida.data$gene.id)
# 21

# Get the sum for the transcripts which are significant for three bias classes
G11MupsigBiasF <- sum(G11MupsigBias$gene.id %in% sebida.febi$gene.id)
# 2
G11MupsigBiasM <- sum(G11MupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
G11MupsigBiasU <- sum(G11MupsigBias$gene.id %in% sebida.unbi$gene.id)
# 19


# Last create a new data frame with the observed and expected numbers
sexbi.G11Mup <- rbind( c(G11MupsigBiasF, G11MupsigBiasM, G11MupsigBiasU),
                       c(prop.female * nrow(G11MupsigBias), prop.male * nrow(G11MupsigBias), prop.unbi * nrow(G11MupsigBias)))
# Add names to the columns and rows
colnames(sexbi.G11Mup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G11Mup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G11Mup["Obs",], p = sexbi.G11Mup["Exp",], rescale.p = T)
# P = 2.145e-07



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# Subset the data to only include transcripts that overlap between group 11 down-regulated and sebida database
G11MdownsigBias <- subset(G11Mdownsig, G11Mdownsig$gene.id %in% sebida.data$gene.id)
# 7

# Get the sum for the transcripts which are significant for three bias classes
G11MdownsigBiasF <- sum(G11MdownsigBias$gene.id %in% sebida.febi$gene.id)
# 0
G11MdownsigBiasM <- sum(G11MdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
G11MdownsigBiasU <- sum(G11MdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 6


# Last create a new data frame with the observed and expected numbers
sexbi.G11Mdown <- rbind( c(G11MdownsigBiasF, G11MdownsigBiasM, G11MdownsigBiasU),
                         c(prop.female * nrow(G11MdownsigBias), prop.male * nrow(G11MdownsigBias), prop.unbi * nrow(G11MdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.G11Mdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G11Mdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G11Mdown["Obs",], p = sexbi.G11Mdown["Exp",], rescale.p = T)
# P = 0.01222

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.G11MdownF <- rbind( c(G11MdownsigBiasF, G11MdownsigBiasM, G11MdownsigBiasU),
                          c(nrow(sebida.febi)-G11MdownsigBiasF, nrow(sebida.mabi)-G11MdownsigBiasM, nrow(sebida.unbi)-G11MdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.G11MdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.G11MdownF)
# P = 0.007547

# We get the same results, so that is good

