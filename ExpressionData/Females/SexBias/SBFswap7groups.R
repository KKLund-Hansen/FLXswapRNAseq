################################################################################################
######### SEXUAL BIAS EXPRESSION USING THE SEBIDA DATABASE FLX SWAP FEMALES - 7 GROUPS #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexBias")

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
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
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# Subset the data to only include transcripts that overlap between group A up-regulated and sebida database
GAFupsigBias <- subset(GAFupsig, GAFupsig$gene.id %in% sebida.data$gene.id)
# 104

# Get the sum for the transcripts which are significant for three bias classes
GAFupsigBiasF <- sum(GAFupsigBias$gene.id %in% sebida.febi$gene.id)
# 92
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
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# Subset the data to only include transcripts that overlap between group A down-regulated and sebida database
GAFdownsigBias <- subset(GAFdownsig, GAFdownsig$gene.id %in% sebida.data$gene.id)
# 1737

# Get the sum for the transcripts which are significant for three bias classes
GAFdownsigBiasF <- sum(GAFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 4
GAFdownsigBiasM <- sum(GAFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 1200
GAFdownsigBiasU <- sum(GAFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 533


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
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# Subset the data to only include transcripts that overlap between group B up-regulated and sebida database
GBFupsigBias <- subset(GBFupsig, GBFupsig$gene.id %in% sebida.data$gene.id)
# 517

# Get the sum for the transcripts which are significant for three bias classes
GBFupsigBiasF <- sum(GBFupsigBias$gene.id %in% sebida.febi$gene.id)
# 468
GBFupsigBiasM <- sum(GBFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GBFupsigBiasU <- sum(GBFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 49


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
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# Subset the data to only include transcripts that overlap between group B down-regulated and sebida database
GBFdownsigBias <- subset(GBFdownsig, GBFdownsig$gene.id %in% sebida.data$gene.id)
# 315

# Get the sum for the transcripts which are significant for three bias classes
GBFdownsigBiasF <- sum(GBFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 11
GBFdownsigBiasM <- sum(GBFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 96
GBFdownsigBiasU <- sum(GBFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 208


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
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# Subset the data to only include transcripts that overlap between group C up-regulated and sebida database
GCFupsigBias <- subset(GCFupsig, GCFupsig$gene.id %in% sebida.data$gene.id)
# 145

# Get the sum for the transcripts which are significant for three bias classes
GCFupsigBiasF <- sum(GCFupsigBias$gene.id %in% sebida.febi$gene.id)
# 134
GCFupsigBiasM <- sum(GCFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GCFupsigBiasU <- sum(GCFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 11


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
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# Subset the data to only include transcripts that overlap between group C down-regulated and sebida database
GCFdownsigBias <- subset(GCFdownsig, GCFdownsig$gene.id %in% sebida.data$gene.id)
# 297

# Get the sum for the transcripts which are significant for three bias classes
GCFdownsigBiasF <- sum(GCFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 9
GCFdownsigBiasM <- sum(GCFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 131
GCFdownsigBiasU <- sum(GCFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 157


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
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# Subset the data to only include transcripts that overlap between group D up-regulated and sebida database
GDFupsigBias <- subset(GDFupsig, GDFupsig$gene.id %in% sebida.data$gene.id)
# 167

# Get the sum for the transcripts which are significant for three bias classes
GDFupsigBiasF <- sum(GDFupsigBias$gene.id %in% sebida.febi$gene.id)
# 155
GDFupsigBiasM <- sum(GDFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 1
GDFupsigBiasU <- sum(GDFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 11


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
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# Subset the data to only include transcripts that overlap between group D down-regulated and sebida database
GDFdownsigBias <- subset(GDFdownsig, GDFdownsig$gene.id %in% sebida.data$gene.id)
# 237

# Get the sum for the transcripts which are significant for three bias classes
GDFdownsigBiasF <- sum(GDFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 2
GDFdownsigBiasM <- sum(GDFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 119
GDFdownsigBiasU <- sum(GDFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 116


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
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# Subset the data to only include transcripts that overlap between group E up-regulated and sebida database
GEFupsigBias <- subset(GEFupsig, GEFupsig$gene.id %in% sebida.data$gene.id)
# 84

# Get the sum for the transcripts which are significant for three bias classes
GEFupsigBiasF <- sum(GEFupsigBias$gene.id %in% sebida.febi$gene.id)
# 51
GEFupsigBiasM <- sum(GEFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 3
GEFupsigBiasU <- sum(GEFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 30


# Last create a new data frame with the observed and expected numbers
sexbi.GEFup <- rbind( c(GEFupsigBiasF, GEFupsigBiasM, GEFupsigBiasU), 
                      c(prop.female * nrow(GEFupsigBias), prop.male * nrow(GEFupsigBias), prop.unbi * nrow(GEFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GEFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GEFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GEFup["Obs",], p = sexbi.GEFup["Exp",], rescale.p = T)
# P = 1.769e-07



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# Subset the data to only include transcripts that overlap between group E down-regulated and sebida database
GEFdownsigBias <- subset(GEFdownsig, GEFdownsig$gene.id %in% sebida.data$gene.id)
# 223

# Get the sum for the transcripts which are significant for three bias classes
GEFdownsigBiasF <- sum(GEFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 23
GEFdownsigBiasM <- sum(GEFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 46
GEFdownsigBiasU <- sum(GEFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 154


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
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# Subset the data to only include transcripts that overlap between group F up-regulated and sebida database
GFFupsigBias <- subset(GFFupsig, GFFupsig$gene.id %in% sebida.data$gene.id)
# 34

# Get the sum for the transcripts which are significant for three bias classes
GFFupsigBiasF <- sum(GFFupsigBias$gene.id %in% sebida.febi$gene.id)
# 24
GFFupsigBiasM <- sum(GFFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GFFupsigBiasU <- sum(GFFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 10


# Last create a new data frame with the observed and expected numbers
sexbi.GFFup <- rbind( c(GFFupsigBiasF, GFFupsigBiasM, GFFupsigBiasU), 
                      c(prop.female * nrow(GFFupsigBias), prop.male * nrow(GFFupsigBias), prop.unbi * nrow(GFFupsigBias)))
# Add names to the columns and rows
colnames(sexbi.GFFup) <- c("Female", "Male", "Unbiased")
rownames(sexbi.GFFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.GFFup["Obs",], p = sexbi.GFFup["Exp",], rescale.p = T)
# P = 4.384e-05



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# Subset the data to only include transcripts that overlap between group F down-regulated and sebida database
GFFdownsigBias <- subset(GFFdownsig, GFFdownsig$gene.id %in% sebida.data$gene.id)
# 230

# Get the sum for the transcripts which are significant for three bias classes
GFFdownsigBiasF <- sum(GFFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 2
GFFdownsigBiasM <- sum(GFFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 118
GFFdownsigBiasU <- sum(GFFdownsigBias$gene.id %in% sebida.unbi$gene.id)
# 110


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
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# Subset the data to only include transcripts that overlap between group G up-regulated and sebida database
GGFupsigBias <- subset(GGFupsig, GGFupsig$gene.id %in% sebida.data$gene.id)
# 83

# Get the sum for the transcripts which are significant for three bias classes
GGFupsigBiasF <- sum(GGFupsigBias$gene.id %in% sebida.febi$gene.id)
# 70
GGFupsigBiasM <- sum(GGFupsigBias$gene.id %in% sebida.mabi$gene.id)
# 0
GGFupsigBiasU <- sum(GGFupsigBias$gene.id %in% sebida.unbi$gene.id)
# 13


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
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# Subset the data to only include transcripts that overlap between group G down-regulated and sebida database
GGFdownsigBias <- subset(GGFdownsig, GGFdownsig$gene.id %in% sebida.data$gene.id)
# 164

# Get the sum for the transcripts which are significant for three bias classes
GGFdownsigBiasF <- sum(GGFdownsigBias$gene.id %in% sebida.febi$gene.id)
# 8
GGFdownsigBiasM <- sum(GGFdownsigBias$gene.id %in% sebida.mabi$gene.id)
# 63
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

