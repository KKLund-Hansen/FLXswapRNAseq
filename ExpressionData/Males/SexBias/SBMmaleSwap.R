################################################################################################
################ SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP MALES ###############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/SexBias")

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)
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
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# Subset the data to only include transcripts that overlap between up-regulated and sebida database
FvCmupsigBias <- subset(FvCmupsig, FvCmupsig$gene.id %in% sebida.data$gene.id)
# 42

# Get the sum for the transcripts which are significant for three bias classes
FvCmupsigBiasF <- sum(FvCmupsigBias$gene.id %in% sebida.febi$gene.id)
# 14
FvCmupsigBiasM <- sum(FvCmupsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
FvCmupsigBiasU <- sum(FvCmupsigBias$gene.id %in% sebida.unbi$gene.id)
# 27


# Last create a new data frame with the observed and expected numbers
sexbi.FvCmup <- rbind( c(FvCmupsigBiasF, FvCmupsigBiasM, FvCmupsigBiasU), 
                       c(prop.female * nrow(FvCmupsigBias), prop.male * nrow(FvCmupsigBias), prop.unbi * nrow(FvCmupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FvCmup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FvCmup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FvCmup["Obs",], p = sexbi.FvCmup["Exp",], rescale.p = T)
# P = 1.467e-05



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXvCFM group
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

# Subset the data to only include transcripts that overlap between down-regulated and sebida database
FvCmdownsigBias <- subset(FvCmdownsig, FvCmdownsig$gene.id %in% sebida.data$gene.id)
# 32

# Get the sum for the transcripts which are significant for three bias classes
FvCmdownsigBiasF <- sum(FvCmdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
FvCmdownsigBiasM <- sum(FvCmdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 8
FvCmdownsigBiasU <- sum(FvCmdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 20

# Last create a new data frame with the observed and expected numbers
sexbi.FvCmdown <- rbind( c(FvCmdownsigBiasF, FvCmdownsigBiasM, FvCmdownsigBiasU),
                         c(prop.female * nrow(FvCmdownsigBias), prop.male * nrow(FvCmdownsigBias), prop.unbi * nrow(FvCmdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FvCmdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FvCmdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FvCmdown["Obs",], p = sexbi.FvCmdown["Exp",], rescale.p = T)
# P = 0.001197




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# Subset the data to only include transcripts that overlap between up-regulated and sebida database
AvSmupsigBias <- subset(AvSmupsig, AvSmupsig$gene.id %in% sebida.data$gene.id)
# 171

# Get the sum for the transcripts which are significant for three bias classes
AvSmupsigBiasF <- sum(AvSmupsigBias$gene.id %in% sebida.febi$gene.id)
# 61
AvSmupsigBiasM <- sum(AvSmupsigBias$gene.id %in% sebida.mabi$gene.id)
# 20
AvSmupsigBiasU <- sum(AvSmupsigBias$gene.id %in% sebida.unbi$gene.id)
# 90

# Last create a new data frame with the observed and expected numbers
sexbi.AvSmup <- rbind( c(AvSmupsigBiasF, AvSmupsigBiasM, AvSmupsigBiasU), 
                       c(prop.female * nrow(AvSmupsigBias), prop.male * nrow(AvSmupsigBias), prop.unbi * nrow(AvSmupsigBias)))
# Add names to the columns and rows
colnames(sexbi.AvSmup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.AvSmup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.AvSmup["Obs",], p = sexbi.AvSmup["Exp",], rescale.p = T)
# P = 1.287e-08



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated n the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# Subset the data to only include transcripts that overlap between down-regulated and sebida database
AvSmdownsigBias <- subset(AvSmdownsig, AvSmdownsig$gene.id %in% sebida.data$gene.id)
# 162

# Get the sum for the transcripts which are significant for three bias classes
AvSmdownsigBiasF <- sum(AvSmdownsigBias$gene.id %in% sebida.febi$gene.id)
# 16
AvSmdownsigBiasM <- sum(AvSmdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 48
AvSmdownsigBiasU <- sum(AvSmdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 98

# Last create a new data frame with the observed and expected numbers
sexbi.AvSmdown <- rbind( c(AvSmdownsigBiasF, AvSmdownsigBiasM, AvSmdownsigBiasU),
                         c(prop.female * nrow(AvSmdownsigBias), prop.male * nrow(AvSmdownsigBias), prop.unbi * nrow(AvSmdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.AvSmdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.AvSmdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.AvSmdown["Obs",], p = sexbi.AvSmdown["Exp",], rescale.p = T)
# P = 1.913e-15




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# Subset the data to only include transcripts that overlap between up-regulated and sebida database
FsvCsMupsigBias <- subset(FsvCsMupsig, FsvCsMupsig$gene.id %in% sebida.data$gene.id)
# 47

# Get the sum for the transcripts which are significant for three bias classes
FsvCsMupsigBiasF <- sum(FsvCsMupsigBias$gene.id %in% sebida.febi$gene.id)
# 6
FsvCsMupsigBiasM <- sum(FsvCsMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 5
FsvCsMupsigBiasU <- sum(FsvCsMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 36

# Last create a new data frame with the observed and expected numbers
sexbi.FsvCsMup <- rbind( c(FsvCsMupsigBiasF, FsvCsMupsigBiasM, FsvCsMupsigBiasU), 
                         c(prop.female * nrow(FsvCsMupsigBias), prop.male * nrow(FsvCsMupsigBias), prop.unbi * nrow(FsvCsMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FsvCsMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FsvCsMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FsvCsMup["Obs",], p = sexbi.FsvCsMup["Exp",], rescale.p = T)
# P = 3.384e-09



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# Subset the data to only include transcripts that overlap between down-regulated and sebida database
FsvCsMdownsigBias <- subset(FsvCsMdownsig, FsvCsMdownsig$gene.id %in% sebida.data$gene.id)
# 29

# Get the sum for the transcripts which are significant for three bias classes
FsvCsMdownsigBiasF <- sum(FsvCsMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 5
FsvCsMdownsigBiasM <- sum(FsvCsMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 2
FsvCsMdownsigBiasU <- sum(FsvCsMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 22

# Last create a new data frame with the observed and expected numbers
sexbi.FsvCsMdown <- rbind( c(FsvCsMdownsigBiasF, FsvCsMdownsigBiasM, FsvCsMdownsigBiasU),
                           c(prop.female * nrow(FsvCsMdownsigBias), prop.male * nrow(FsvCsMdownsigBias), prop.unbi * nrow(FsvCsMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FsvCsMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FsvCsMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FsvCsMdown["Obs",], p = sexbi.FsvCsMdown["Exp",], rescale.p = T)
# P = 7.932e-06




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# Subset the data to only include transcripts that overlap between up-regulated and sebida database
FavCaMupsigBias <- subset(FavCaMupsig, FavCaMupsig$gene.id %in% sebida.data$gene.id)
# 51

# Get the sum for the transcripts which are significant for three bias classes
FavCaMupsigBiasF <- sum(FavCaMupsigBias$gene.id %in% sebida.febi$gene.id)
# 13
FavCaMupsigBiasM <- sum(FavCaMupsigBias$gene.id %in% sebida.mabi$gene.id)
# 5
FavCaMupsigBiasU <- sum(FavCaMupsigBias$gene.id %in% sebida.unbi$gene.id)
# 33

# Last create a new data frame with the observed and expected numbers
sexbi.FavCaMup <- rbind( c(FavCaMupsigBiasF, FavCaMupsigBiasM, FavCaMupsigBiasU), 
                         c(prop.female * nrow(FavCaMupsigBias), prop.male * nrow(FavCaMupsigBias), prop.unbi * nrow(FavCaMupsigBias)))
# Add names to the columns and rows
colnames(sexbi.FavCaMup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FavCaMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FavCaMup["Obs",], p = sexbi.FavCaMup["Exp",], rescale.p = T)
# P = 9.405e-06



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# Subset the data to only include transcripts that overlap between down-regulated and sebida database
FavCaMdownsigBias <- subset(FavCaMdownsig, FavCaMdownsig$gene.id %in% sebida.data$gene.id)
# 110

# Get the sum for the transcripts which are significant for three bias classes
FavCaMdownsigBiasF <- sum(FavCaMdownsigBias$gene.id %in% sebida.febi$gene.id)
# 25
FavCaMdownsigBiasM <- sum(FavCaMdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 10
FavCaMdownsigBiasU <- sum(FavCaMdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 75

# Last create a new data frame with the observed and expected numbers
sexbi.FavCaMdown <- rbind( c(FavCaMdownsigBiasF, FavCaMdownsigBiasM, FavCaMdownsigBiasU),
                           c(prop.female * nrow(FavCaMdownsigBias), prop.male * nrow(FavCaMdownsigBias), prop.unbi * nrow(FavCaMdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.FavCaMdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FavCaMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FavCaMdown["Obs",], p = sexbi.FavCaMdown["Exp",], rescale.p = T)
# P = 6.822e-14

