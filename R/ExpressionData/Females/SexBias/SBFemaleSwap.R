################################################################################################
############### SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP FEMALES ##############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Females/SexBias")

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
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
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# Subset the data to only include transcripts that overlap between group A up-regulated and sebida database
GAFupsigBias <- subset(GAFupsig, GAFupsig$gene.id %in% sebida.data$gene.id)
# 110

# Get the sum for the transcripts which are significant for three bias classes
GAFupsigBiasF <- sum(GAFupsigBias$gene.id %in% sebida.febi$gene.id)
# 98
GAFupsigBiasM <- sum(GAFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GAFupsigBiasU <- sum(GAFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 12


# Last create a new data frame with the observed and expected numbers
sexbi.GAFup <- rbind( c(GAFupsigBiasF, GAFupsigBiasM, GAFupsigBiasU), 
                      c(prop.female * nrow(GAFupsigBias), prop.male * nrow(GAFupsigBias), prop.unbi * nrow(GAFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GAFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GAFup["Obs",], p = sexbi.GAFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

# Subset the data to only include transcripts that overlap between group A down-regulated and sebida database
GAFdownsigBias <- subset(GAFdownsig, GAFdownsig$gene.id %in% sebida.data$gene.id)
# 1776

# Get the sum for the transcripts which are significant for three bias classes
GAFdownsigBiasF <- sum(GAFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
GAFdownsigBiasM <- sum(GAFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1220
GAFdownsigBiasU <- sum(GAFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 552


# Last create a new data frame with the observed and expected numbers
sexbi.GAFdown <- rbind( c(GAFdownsigBiasF, GAFdownsigBiasM, GAFdownsigBiasU),
                        c(prop.female * nrow(GAFdownsigBias), prop.male * nrow(GAFdownsigBias), prop.unbi * nrow(GAFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GAFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GAFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GAFdown["Obs",], p = sexbi.GAFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group B
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# Subset the data to only include transcripts that overlap between group B up-regulated and sebida database
GBFupsigBias <- subset(GBFupsig, GBFupsig$gene.id %in% sebida.data$gene.id)
# 517

# Get the sum for the transcripts which are significant for three bias classes
GBFupsigBiasF <- sum(GBFupsigBias$gene.id %in% sebida.febi$gene.id)
# 469
GBFupsigBiasM <- sum(GBFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GBFupsigBiasU <- sum(GBFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 48


# Last create a new data frame with the observed and expected numbers
sexbi.GBFup <- rbind( c(GBFupsigBiasF, GBFupsigBiasM, GBFupsigBiasU), 
                      c(prop.female * nrow(GBFupsigBias), prop.male * nrow(GBFupsigBias), prop.unbi * nrow(GBFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GBFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GBFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GBFup["Obs",], p = sexbi.GBFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# Subset the data to only include transcripts that overlap between group B down-regulated and sebida database
GBFdownsigBias <- subset(GBFdownsig, GBFdownsig$gene.id %in% sebida.data$gene.id)
# 314

# Get the sum for the transcripts which are significant for three bias classes
GBFdownsigBiasF <- sum(GBFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 10
GBFdownsigBiasM <- sum(GBFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 98
GBFdownsigBiasU <- sum(GBFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 206


# Last create a new data frame with the observed and expected numbers
sexbi.GBFdown <- rbind( c(GBFdownsigBiasF, GBFdownsigBiasM, GBFdownsigBiasU),
                        c(prop.female * nrow(GBFdownsigBias), prop.male * nrow(GBFdownsigBias), prop.unbi * nrow(GBFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GBFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GBFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GBFdown["Obs",], p = sexbi.GBFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# Subset the data to only include transcripts that overlap between group C up-regulated and sebida database
GCFupsigBias <- subset(GCFupsig, GCFupsig$gene.id %in% sebida.data$gene.id)
# 159

# Get the sum for the transcripts which are significant for three bias classes
GCFupsigBiasF <- sum(GCFupsigBias$gene.id %in% sebida.febi$gene.id)
# 146
GCFupsigBiasM <- sum(GCFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GCFupsigBiasU <- sum(GCFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 13


# Last create a new data frame with the observed and expected numbers
sexbi.GCFup <- rbind( c(GCFupsigBiasF, GCFupsigBiasM, GCFupsigBiasU), 
                      c(prop.female * nrow(GCFupsigBias), prop.male * nrow(GCFupsigBias), prop.unbi * nrow(GCFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GCFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GCFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GCFup["Obs",], p = sexbi.GCFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# Subset the data to only include transcripts that overlap between group C down-regulated and sebida database
GCFdownsigBias <- subset(GCFdownsig, GCFdownsig$gene.id %in% sebida.data$gene.id)
# 302

# Get the sum for the transcripts which are significant for three bias classes
GCFdownsigBiasF <- sum(GCFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 9
GCFdownsigBiasM <- sum(GCFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 137
GCFdownsigBiasU <- sum(GCFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 156


# Last create a new data frame with the observed and expected numbers
sexbi.GCFdown <- rbind( c(GCFdownsigBiasF, GCFdownsigBiasM, GCFdownsigBiasU),
                        c(prop.female * nrow(GCFdownsigBias), prop.male * nrow(GCFdownsigBias), prop.unbi * nrow(GCFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GCFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GCFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GCFdown["Obs",], p = sexbi.GCFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# Subset the data to only include transcripts that overlap between group D up-regulated and sebida database
GDFupsigBias <- subset(GDFupsig, GDFupsig$gene.id %in% sebida.data$gene.id)
# 181

# Get the sum for the transcripts which are significant for three bias classes
GDFupsigBiasF <- sum(GDFupsigBias$gene.id %in% sebida.febi$gene.id)
# 170
GDFupsigBiasM <- sum(GDFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GDFupsigBiasU <- sum(GDFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 10


# Last create a new data frame with the observed and expected numbers
sexbi.GDFup <- rbind( c(GDFupsigBiasF, GDFupsigBiasM, GDFupsigBiasU), 
                      c(prop.female * nrow(GDFupsigBias), prop.male * nrow(GDFupsigBias), prop.unbi * nrow(GDFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GDFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GDFup["Obs",], p = sexbi.GDFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# Subset the data to only include transcripts that overlap between group D down-regulated and sebida database
GDFdownsigBias <- subset(GDFdownsig, GDFdownsig$gene.id %in% sebida.data$gene.id)
# 237

# Get the sum for the transcripts which are significant for three bias classes
GDFdownsigBiasF <- sum(GDFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 2
GDFdownsigBiasM <- sum(GDFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 112
GDFdownsigBiasU <- sum(GDFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 123


# Last create a new data frame with the observed and expected numbers
sexbi.GDFdown <- rbind( c(GDFdownsigBiasF, GDFdownsigBiasM, GDFdownsigBiasU),
                        c(prop.female * nrow(GDFdownsigBias), prop.male * nrow(GDFdownsigBias), prop.unbi * nrow(GDFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GDFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GDFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GDFdown["Obs",], p = sexbi.GDFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# Subset the data to only include transcripts that overlap between group E up-regulated and sebida database
GEFupsigBias <- subset(GEFupsig, GEFupsig$gene.id %in% sebida.data$gene.id)
# 87

# Get the sum for the transcripts which are significant for three bias classes
GEFupsigBiasF <- sum(GEFupsigBias$gene.id %in% sebida.febi$gene.id)
# 53
GEFupsigBiasM <- sum(GEFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 3
GEFupsigBiasU <- sum(GEFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 31


# Last create a new data frame with the observed and expected numbers
sexbi.GEFup <- rbind( c(GEFupsigBiasF, GEFupsigBiasM, GEFupsigBiasU), 
                      c(prop.female * nrow(GEFupsigBias), prop.male * nrow(GEFupsigBias), prop.unbi * nrow(GEFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEFup["Obs",], p = sexbi.GEFup["Exp",], rescale.p = T)
# P = 8.298e-08



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# Subset the data to only include transcripts that overlap between group E down-regulated and sebida database
GEFdownsigBias <- subset(GEFdownsig, GEFdownsig$gene.id %in% sebida.data$gene.id)
# 221

# Get the sum for the transcripts which are significant for three bias classes
GEFdownsigBiasF <- sum(GEFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 23
GEFdownsigBiasM <- sum(GEFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 47
GEFdownsigBiasU <- sum(GEFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 151


# Last create a new data frame with the observed and expected numbers
sexbi.GEFdown <- rbind( c(GEFdownsigBiasF, GEFdownsigBiasM, GEFdownsigBiasU),
                        c(prop.female * nrow(GEFdownsigBias), prop.male * nrow(GEFdownsigBias), prop.unbi * nrow(GEFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEFdown["Obs",], p = sexbi.GEFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group F
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# Subset the data to only include transcripts that overlap between group F up-regulated and sebida database
GFFupsigBias <- subset(GFFupsig, GFFupsig$gene.id %in% sebida.data$gene.id)
# 35

# Get the sum for the transcripts which are significant for three bias classes
GFFupsigBiasF <- sum(GFFupsigBias$gene.id %in% sebida.febi$gene.id)
# 24
GFFupsigBiasM <- sum(GFFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GFFupsigBiasU <- sum(GFFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 11


# Last create a new data frame with the observed and expected numbers
sexbi.GFFup <- rbind( c(GFFupsigBiasF, GFFupsigBiasM, GFFupsigBiasU), 
                      c(prop.female * nrow(GFFupsigBias), prop.male * nrow(GFFupsigBias), prop.unbi * nrow(GFFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GFFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GFFup["Obs",], p = sexbi.GFFup["Exp",], rescale.p = T)
# P = 6.404e-05



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# Subset the data to only include transcripts that overlap between group F down-regulated and sebida database
GFFdownsigBias <- subset(GFFdownsig, GFFdownsig$gene.id %in% sebida.data$gene.id)
# 232

# Get the sum for the transcripts which are significant for three bias classes
GFFdownsigBiasF <- sum(GFFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 2
GFFdownsigBiasM <- sum(GFFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 118
GFFdownsigBiasU <- sum(GFFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 112


# Last create a new data frame with the observed and expected numbers
sexbi.GFFdown <- rbind( c(GFFdownsigBiasF, GFFdownsigBiasM, GFFdownsigBiasU),
                        c(prop.female * nrow(GFFdownsigBias), prop.male * nrow(GFFdownsigBias), prop.unbi * nrow(GFFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GFFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GFFdown["Obs",], p = sexbi.GFFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# Subset the data to only include transcripts that overlap between group G up-regulated and sebida database
GGFupsigBias <- subset(GGFupsig, GGFupsig$gene.id %in% sebida.data$gene.id)
# 81

# Get the sum for the transcripts which are significant for three bias classes
GGFupsigBiasF <- sum(GGFupsigBias$gene.id %in% sebida.febi$gene.id)
# 67
GGFupsigBiasM <- sum(GGFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GGFupsigBiasU <- sum(GGFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 14


# Last create a new data frame with the observed and expected numbers
sexbi.GGFup <- rbind( c(GGFupsigBiasF, GGFupsigBiasM, GGFupsigBiasU), 
                      c(prop.female * nrow(GGFupsigBias), prop.male * nrow(GGFupsigBias), prop.unbi * nrow(GGFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GGFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GGFup["Obs",], p = sexbi.GGFup["Exp",], rescale.p = T)
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# Subset the data to only include transcripts that overlap between group G down-regulated and sebida database
GGFdownsigBias <- subset(GGFdownsig, GGFdownsig$gene.id %in% sebida.data$gene.id)
# 163

# Get the sum for the transcripts which are significant for three bias classes
GGFdownsigBiasF <- sum(GGFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 9
GGFdownsigBiasM <- sum(GGFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 61
GGFdownsigBiasU <- sum(GGFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 93


# Last create a new data frame with the observed and expected numbers
sexbi.GGFdown <- rbind( c(GGFdownsigBiasF, GGFdownsigBiasM, GGFdownsigBiasU),
                        c(prop.female * nrow(GGFdownsigBias), prop.male * nrow(GGFdownsigBias), prop.unbi * nrow(GGFdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.GGFdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GGFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GGFdown["Obs",], p = sexbi.GGFdown["Exp",], rescale.p = T)
# P < 2.2e-16




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# Subset the data to only include transcripts that overlap between group 23 up-regulated and sebida database
G23FupsigBias <- subset(G23Fupsig, G23Fupsig$gene.id %in% sebida.data$gene.id)
# 14

# Get the sum for the transcripts which are significant for three bias classes
G23FupsigBiasF <- sum(G23FupsigBias$gene.id %in% sebida.febi$gene.id)
# 2
G23FupsigBiasM <- sum(G23FupsigBias$gene.id %in% sebida.mabi$gene.id)
# 2
G23FupsigBiasU <- sum(G23FupsigBias$gene.id %in% sebida.unbi$gene.id)
# 10


# Last create a new data frame with the observed and expected numbers
sexbi.G23Fup <- rbind( c(G23FupsigBiasF, G23FupsigBiasM, G23FupsigBiasU), 
                       c(prop.female * nrow(G23FupsigBias), prop.male * nrow(G23FupsigBias), prop.unbi * nrow(G23FupsigBias)))
# Add names to the columns and rows
colnames(sexbi.G23Fup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G23Fup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G23Fup["Obs",], p = sexbi.G23Fup["Exp",], rescale.p = T)
# P = 0.01098

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.G23FupF <- rbind( c(G23FupsigBiasF, G23FupsigBiasM, G23FupsigBiasU),
                        c(nrow(sebida.febi)-G23FupsigBiasF, nrow(sebida.mabi)-G23FupsigBiasM, nrow(sebida.unbi)-G23FupsigBiasU))
# Add names to the columns and rows
colnames(sexbi.G23FupF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G23FupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.G23FupF)
# P = 0.01595

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# Subset the data to only include transcripts that overlap between group 23 down-regulated and sebida database
G23FdownsigBias <- subset(G23Fdownsig, G23Fdownsig$gene.id %in% sebida.data$gene.id)
# 4

# Get the sum for the transcripts which are significant for three bias classes
G23FdownsigBiasF <- sum(G23FdownsigBias$gene.id %in% sebida.febi$gene.id)
# 1
G23FdownsigBiasM <- sum(G23FdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
G23FdownsigBiasU <- sum(G23FdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 2


# Last create a new data frame with the observed and expected numbers
sexbi.G23Fdown <- rbind( c(G23FdownsigBiasF, G23FdownsigBiasM, G23FdownsigBiasU),
                         c(prop.female * nrow(G23FdownsigBias), prop.male * nrow(G23FdownsigBias), prop.unbi * nrow(G23FdownsigBias)))
# Add names to the columns and rows
colnames(sexbi.G23Fdown) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G23Fdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.G23Fdown["Obs",], p = sexbi.G23Fdown["Exp",], rescale.p = T)
# P = 0.7757

# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.G23FdownF <- rbind( c(G23FdownsigBiasF, G23FdownsigBiasM, G23FdownsigBiasU),
                          c(nrow(sebida.febi)-G23FdownsigBiasF, nrow(sebida.mabi)-G23FdownsigBiasM, nrow(sebida.unbi)-G23FdownsigBiasU))
# Add names to the columns and rows
colnames(sexbi.G23FdownF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.G23FdownF)
# P = 0.8364

# We get the same results, so that is good


