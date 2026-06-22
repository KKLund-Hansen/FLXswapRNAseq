################################################################################################
######### SEXUAL BIASED EXPRESSION USING THE SEBIDA DATABASE FLX SWAP MALE GROUPS ONLY #########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/SexBias")

# Read in files with data
SwapMaleOnly.data <- read.table(file = "FLXswapMaleGroups.csv", h = T, sep = ",", stringsAsFactors = T)
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



######################## SEX BIAS IN GENERAL ########################

# Subset the overall data to only include transcripts that are significant and overlap with the sebida data base
MaleSig <- subset(SwapMaleOnly.data, Overall == "Sig")
# 588
MaleSigBias <- subset(MaleSig, MaleSig$gene.id %in% sebida.data$gene.id)
# 467

# How many are proportional sex biased?
MaleSigSexBias <- sum(MaleSigBias$gene.id %in% sebida.febi$gene.id) + sum(MaleSigBias$gene.id %in% sebida.mabi$gene.id)
# 165
# Proportion
prop.MaleSigSexBias <- MaleSigSexBias/length(MaleSigBias$gene.id)
# 0.3533191
# How many are not?
MaleSigNotSexBias <-  sum(MaleSigBias$gene.id %in% sebida.unbi$gene.id)
# 302
# Proportion
prop.MaleSigNotSexBias <- MaleSigNotSexBias/length(MaleSigBias$gene.id)
# 0.6466809


# Subset to significant transcripts in the female only groups
MaOnSig <- subset(MaleSigBias, g26 == "Sig" | g33 == "Sig" | g36 == "Sig" | g55 == "Sig")
# 3

# How many are sex biased?
MaOnSigSexBias <- sum(MaOnSig$gene.id %in% sebida.febi$gene.id) + sum(MaOnSig$gene.id %in% sebida.mabi$gene.id)
# 1
# The number of expected transcripts which are sex biased
length(MaOnSig$gene.id)*prop.MaleSigSexBias
# 1.059957
# How many are not?
MaOnSigNotSexBias <- sum(MaOnSig$gene.id %in% sebida.unbi$gene.id)
# 2


# Chi-squared Test
# Do the test, with the number of sex biased transcripts and not sex biased transcripts and the proportion of overall sex biased transcripts 
chisq.test(c(MaOnSigSexBias, MaOnSigNotSexBias), p = c(prop.MaleSigSexBias, prop.MaleSigNotSexBias))
# P = 0.9423


# Fisher test
SexBias.MaleOnly <- rbind( c(MaOnSigSexBias, MaOnSigNotSexBias),
                           c(MaleSigSexBias-MaOnSigSexBias, MaleSigNotSexBias-MaOnSigNotSexBias ))
# Add names to the columns and rows
colnames(SexBias.MaleOnly) <- c("SexBias", "NotSexBias")
rownames(SexBias.MaleOnly) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SexBias.MaleOnly)
# P = 1



######################## MALE ONLY GROUPS ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts
MaOnUpSig <- subset(SwapMaleOnly.data, g26up == "Sig" | g33up == "Sig" | g36up == "Sig" | g55up == "Sig")
# 2


# Subset the data to only include transcripts that overlap between male only groups and the SEBIDA data base
MaOnUpSigBias <- subset(MaOnUpSig, MaOnUpSig$gene.id %in% sebida.data$gene.id)
# 2

#Get the sum for the transcripts which are significant for three bias classes
MaOnUpSigBiasF <- sum(MaOnUpSigBias$gene.id %in% sebida.febi$gene.id)
# 0
MaOnUpSigBiasM <- sum(MaOnUpSigBias$gene.id %in% sebida.mabi$gene.id)
# 0
MaOnUpSigBiasU <- sum(MaOnUpSigBias$gene.id %in% sebida.unbi$gene.id)
# 2


# Last create a new data frame with the observed and expected numbers
sexbi.MaOnUpSig <- rbind( c(MaOnUpSigBiasF, MaOnUpSigBiasM, MaOnUpSigBiasU), 
                          c(prop.female * nrow(MaOnUpSigBias), prop.male * nrow(MaOnUpSigBias), prop.unbi * nrow(MaOnUpSigBias)))
#Add names to the columns and rows
colnames(sexbi.MaOnUpSig) <- c("Female", "Male", "Unbiased")
rownames(sexbi.MaOnUpSig) <- c("Obs", "Exp")

#It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.MaOnUpSig["Obs",], p = sexbi.MaOnUpSig["Exp",], rescale.p = T)
# P = 0.1383


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.MaOnUpSigF <- rbind( c(MaOnUpSigBiasF, MaOnUpSigBiasM , MaOnUpSigBiasU),
                           c(nrow(sebida.febi)-MaOnUpSigBiasF, nrow(sebida.mabi)-MaOnUpSigBiasM, nrow(sebida.unbi)-MaOnUpSigBiasU))
#Add names to the columns and rows
colnames(sexbi.MaOnUpSigF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.MaOnUpSigF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.MaOnUpSigF)
# P = 0.197

# We get the same results, so that is good




##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify the observed sex-biased classes for down-regulated transcripts
MaOnDownSig <- subset(SwapMaleOnly.data, g26down == "Sig" | g33down == "Sig" | g36down == "Sig" | g55down == "Sig")
# 3


# Subset the data to only include transcripts that overlap
MaOnDownSigBias <- subset(MaOnDownSig, MaOnDownSig$gene.id %in% sebida.data$gene.id)
# 1

# Get the sum for the transcripts which are significant for three bias classes
MaOnDownSigBiasF <- sum(MaOnDownSigBias$gene.id %in% sebida.febi$gene.id)
# 0
MaOnDownSigBiasM <- sum(MaOnDownSigBias$gene.id %in% sebida.mabi$gene.id)
# 1
MaOnDownSigBiasU <- sum(MaOnDownSigBias$gene.id %in% sebida.unbi$gene.id)
# 0


#Last create a new data frame with the observed and expected numbers
sexbi.MaOnDownSig <- rbind( c(MaOnDownSigBiasF, MaOnDownSigBiasM, MaOnDownSigBiasU), 
                            c(prop.female * nrow(MaOnDownSigBias), prop.male * nrow(MaOnDownSigBias), prop.unbi * nrow(MaOnDownSigBias)))
#Add names to the columns and rows
colnames(sexbi.MaOnDownSig) <- c("Female", "Male", "Unbiased")
rownames(sexbi.MaOnDownSig) <- c("Obs", "Exp")

#It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.MaOnDownSig["Obs",], p = sexbi.MaOnDownSig["Exp",], rescale.p = T)
# P = 0.2946


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
sexbi.MaOnDownSigF <- rbind( c(MaOnDownSigBiasF, MaOnDownSigBiasM , MaOnDownSigBiasU),
                             c(nrow(sebida.febi)-MaOnDownSigBiasF, nrow(sebida.mabi)-MaOnDownSigBiasM, nrow(sebida.unbi)-MaOnDownSigBiasU))
#Add names to the columns and rows
colnames(sexbi.MaOnDownSigF) <- c("Female", "Male", "Unbiased")
rownames(sexbi.MaOnDownSigF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(sexbi.MaOnDownSigF)
# P = 0.2904

# We get the same results, so that is good



#########################################  PLOT DATA  #########################################
