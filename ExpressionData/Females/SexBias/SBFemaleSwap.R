################################################################################################
############### SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP FEMALES ##############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexBias")

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# Subset the data to only include transcripts that overlap between up-regulated group and sebida database
FvCfupsigBias <- subset(FvCfupsig, FvCfupsig$gene.id %in% sebida.data$gene.id)
# 28

# Get the sum for the transcripts which are significant for three bias classes
FvCfupsigBiasF <- sum(FvCfupsigBias$gene.id %in% sebida.febi$gene.id)
# 9
FvCfupsigBiasM <- sum(FvCfupsigBias$gene.id %in% sebida.mabi$gene.id)
# 4
FvCfupsigBiasU <- sum(FvCfupsigBias$gene.id %in% sebida.unbi$gene.id)
# 15


# Last create a new data frame with the observed and expected numbers
sexbi.FvCfup <- rbind( c(FvCfupsigBiasF, FvCfupsigBiasM, FvCfupsigBiasU), 
                       c(prop.female * nrow(FvCfupsigBias), prop.male * nrow(FvCfupsigBias), prop.unbi * nrow(FvCfupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FvCfup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FvCfup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FvCfup["Obs",], p = sexbi.FvCfup["Exp",], rescale.p = T)
# P = 0.05964



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# Subset the data to only include transcripts that overlap between down-regulated group and sebida database
FvCfdownsigBias <- subset(FvCfdownsig, FvCfdownsig$gene.id %in% sebida.data$gene.id)
# 24

# Get the sum for the transcripts which are significant for three bias classes
FvCfdownsigBiasF <- sum(FvCfdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
FvCfdownsigBiasM <- sum(FvCfdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 5
FvCfdownsigBiasU <- sum(FvCfdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 15

# Last create a new data frame with the observed and expected numbers
sexbi.FvCfdown <- rbind( c(FvCfdownsigBiasF, FvCfdownsigBiasM, FvCfdownsigBiasU),
                         c(prop.female * nrow(FvCfdownsigBias), prop.male * nrow(FvCfdownsigBias), prop.unbi * nrow(FvCfdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FvCfdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FvCfdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FvCfdown["Obs",], p = sexbi.FvCfdown["Exp",], rescale.p = T)
# P = 0.009594




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the autovsex group 
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# Subset the data to only include transcripts that overlap between up-regulated group and sebida database
AvSfupsigBias <- subset(AvSfupsig, AvSfupsig$gene.id %in% sebida.data$gene.id)
# 187

# Get the sum for the transcripts which are significant for three bias classes
AvSfupsigBiasF <- sum(AvSfupsigBias$gene.id %in% sebida.febi$gene.id)
# 14
AvSfupsigBiasM <- sum(AvSfupsigBias$gene.id %in% sebida.mabi$gene.id)
# 63
AvSfupsigBiasU <- sum(AvSfupsigBias$gene.id %in% sebida.unbi$gene.id)
# 110

# Last create a new data frame with the observed and expected numbers
sexbi.AvSfup <- rbind( c(AvSfupsigBiasF, AvSfupsigBiasM, AvSfupsigBiasU), 
                       c(prop.female * nrow(AvSfupsigBias), prop.male * nrow(AvSfupsigBias), prop.unbi * nrow(AvSfupsigBias)))
# Add names to the columns and rows
colnames(sexbi.AvSfup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.AvSfup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.AvSfup["Obs",], p = sexbi.AvSfup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the autovsex group 
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# Subset the data to only include transcripts that overlap between down-regulated group and sebida database
AvSfdownsigBias <- subset(AvSfdownsig, AvSfdownsig$gene.id %in% sebida.data$gene.id)
# 129

# Get the sum for the transcripts which are significant for three bias classes
AvSfdownsigBiasF <- sum(AvSfdownsigBias$gene.id %in% sebida.febi$gene.id)
# 41
AvSfdownsigBiasM <- sum(AvSfdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 17
AvSfdownsigBiasU <- sum(AvSfdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 71

7 <- rbind( c(AvSfdownsigBiasF, AvSfdownsigBiasM, AvSfdownsigBiasU),
                         c(prop.female * nrow(AvSfdownsigBias), prop.male * nrow(AvSfdownsigBias), prop.unbi * nrow(AvSfdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.AvSfdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.AvSfdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.AvSfdown["Obs",], p = sexbi.AvSfdown["Exp",], rescale.p = T)
# P = 3.127e-07




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# Subset the data to only include transcripts that overlap between up-regulated group and sebida database
FsvCsFupsigBias <- subset(FsvCsFupsig, FsvCsFupsig$gene.id %in% sebida.data$gene.id)
# 127

# Get the sum for the transcripts which are significant for three bias classes
FsvCsFupsigBiasF <- sum(FsvCsFupsigBias$gene.id %in% sebida.febi$gene.id)
# 9
FsvCsFupsigBiasM <- sum(FsvCsFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 12
FsvCsFupsigBiasU <- sum(FsvCsFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 106

# Last create a new data frame with the observed and expected numbers
sexbi.FsvCsFup <- rbind( c(FsvCsFupsigBiasF, FsvCsFupsigBiasM, FsvCsFupsigBiasU), 
                         c(prop.female * nrow(FsvCsFupsigBias), prop.male * nrow(FsvCsFupsigBias), prop.unbi * nrow(FsvCsFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FsvCsFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FsvCsFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FsvCsFup["Obs",], p = sexbi.FsvCsFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# Subset the data to only include transcripts that overlap between down-regulated and sebida database
FsvCsFdownsigBias <- subset(FsvCsFdownsig, FsvCsFdownsig$gene.id %in% sebida.data$gene.id)
# 37

# Get the sum for the transcripts which are significant for three bias classes
FsvCsFdownsigBiasF <- sum(FsvCsFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 9
FsvCsFdownsigBiasM <- sum(FsvCsFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
FsvCsFdownsigBiasU <- sum(FsvCsFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 27

# Last create a new data frame with the observed and expected numbers
sexbi.FsvCsFdown <- rbind( c(FsvCsFdownsigBiasF, FsvCsFdownsigBiasM, FsvCsFdownsigBiasU),
                           c(prop.female * nrow(FsvCsFdownsigBias), prop.male * nrow(FsvCsFdownsigBias), prop.unbi * nrow(FsvCsFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FsvCsFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FsvCsFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FsvCsFdown["Obs",], p = sexbi.FsvCsFdown["Exp",], rescale.p = T)
# P = 1e-06




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# Subset the data to only include transcripts that overlap between up-regulated group and sebida database
FavCaFupsigBias <- subset(FavCaFupsig, FavCaFupsig$gene.id %in% sebida.data$gene.id)
# 478

# Get the sum for the transcripts which are significant for three bias classes
FavCaFupsigBiasF <- sum(FavCaFupsigBias$gene.id %in% sebida.febi$gene.id)
# 421
FavCaFupsigBiasM <- sum(FavCaFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
FavCaFupsigBiasU <- sum(FavCaFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 56

# Last create a new data frame with the observed and expected numbers
sexbi.FavCaFup <- rbind( c(FavCaFupsigBiasF, FavCaFupsigBiasM, FavCaFupsigBiasU), 
                         c(prop.female * nrow(FavCaFupsigBias), prop.male * nrow(FavCaFupsigBias), prop.unbi * nrow(FavCaFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FavCaFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FavCaFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FavCaFup["Obs",], p = sexbi.FavCaFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# Subset the data to only include transcripts that overlap between down-regulated group and sebida database
FavCaFdownsigBias <- subset(FavCaFdownsig, FavCaFdownsig$gene.id %in% sebida.data$gene.id)
# 2592

# Get the sum for the transcripts which are significant for three bias classes
FavCaFdownsigBiasF <- sum(FavCaFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 19
FavCaFdownsigBiasM <- sum(FavCaFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1604
FavCaFdownsigBiasU <- sum(FavCaFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 969

# Last create a new data frame with the observed and expected numbers
sexbi.FavCaFdown <- rbind( c(FavCaFdownsigBiasF, FavCaFdownsigBiasM, FavCaFdownsigBiasU),
                           c(prop.female * nrow(FavCaFdownsigBias), prop.male * nrow(FavCaFdownsigBias), prop.unbi * nrow(FavCaFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FavCaFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FavCaFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FavCaFdown["Obs",], p = sexbi.FavCaFdown["Exp",], rescale.p = T)
# P < 2.2e-16

