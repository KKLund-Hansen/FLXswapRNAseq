################################################################################################
######## SEXUAL BIASED EXPRESSION USING THE SEBIDA DATABASE FLX SWAP FEMALE GROUPS ONLY ########
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Flies/FLX/5.FLXchromosomeSwap/R/Females/SexBias")

# Read in files with data
SwapFemaleOnly.data <- read.table(file = "FLXswapFemaleGroups.csv", h = T, sep = ",", stringsAsFactors = T)
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



######################## SEX BIAS IN COMPARISON ########################


# Subset the overall data to only include transcripts that are significant and overlap with the sebida data base
FemaleSig <- subset(SwapFemaleOnly.data, Overall == "Sig")
# 5725
FemaleSigBias <- subset(FemaleSig, FemaleSig$gene.id %in% sebida.data$gene.id)
# 4960

# How many are proportional sex biased?
FemaleSigSexBias <- sum(FemaleSigBias$gene.id %in% sebida.febi$gene.id) + sum(FemaleSigBias$gene.id %in% sebida.mabi$gene.id)
# 3027
# Proportion
prop.FemaleSigSexBias <- FemaleSigSexBias/length(FemaleSigBias$gene.id)
# 0.6102823
# How many are not?
FemaleSigNotSexBias <-  sum(FemaleSigBias$gene.id %in% sebida.unbi$gene.id)
# 1933
# Proportion
prop.FemaleSigNotSexBias <- FemaleSigNotSexBias/length(FemaleSigBias$gene.id)
# 0.3897177


# Subset to significant transcripts in the female only groups
FeOnSig <- subset(FemaleSigBias, g35 == "Sig" | g39 == "Sig" | g40 == "Sig" | g42 == "Sig" | g44 == "Sig" | g49 == "Sig" | 
                                 g56 == "Sig" | g57 == "Sig" | g59 == "Sig" | g62 == "Sig" | g63 == "Sig")
# 71

# How many are sex biased?
FeOnSigSexBias <- sum(FeOnSig$gene.id %in% sebida.febi$gene.id) + sum(FeOnSig$gene.id %in% sebida.mabi$gene.id)
# 23
# The number of expected transcripts which are sex biased
length(FeOnSig$gene.id)*prop.FemaleSigSexBias
# 43.33004
# How many are not?
FeOnSigNotSexBias <-  sum(FeOnSig$gene.id %in% sebida.unbi$gene.id)
# 48


# Chi-squared Test
# Do the test, with the number of sex biased transcripts and not sex biased transcripts and the proportion of overall sex biased transcripts 
chisq.test(c(FeOnSigSexBias, FeOnSigNotSexBias), p = c(prop.FemaleSigSexBias, prop.FemaleSigNotSexBias))
# P = 7.525e-07


# Fisher test
SexBias.FemaleOnly <- rbind( c(FeOnSigSexBias, FeOnSigNotSexBias),
                             c(FemaleSigSexBias-FeOnSigSexBias, FemaleSigNotSexBias-FeOnSigNotSexBias))
# Add names to the columns and rows
colnames(SexBias.FemaleOnly) <- c("SexBias", "NotSexBias")
rownames(SexBias.FemaleOnly) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(SexBias.FemaleOnly)
# P = 1.015e-06




######################## FEMALE ONLY GROUPS ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation of sex biased transcripts

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLX effect class
FeOnUpSig <- subset(SwapFemaleOnly.data, g35up == "Sig" | g39up == "Sig" | g40up == "Sig" | g42up == "Sig" | g44up == "Sig" | g49up == "Sig" |
                                         g56up == "Sig" | g57up == "Sig" | g59up == "Sig" | g62up == "Sig" | g63up == "Sig")
# 46


# Subset the data to only include transcripts that overlap between the FLX up-regulated and sebida database
FeOnUpSigBias <- subset(FeOnUpSig, FeOnUpSig$gene.id %in% sebida.data$gene.id)
# 43

# Get the sum for the transcripts which are significant for three bias classes
FeOnUpSigBiasF <- sum(FeOnUpSigBias$gene.id %in% sebida.febi$gene.id)
# 7
FeOnUpSigBiasM <- sum(FeOnUpSigBias$gene.id %in% sebida.mabi$gene.id)
# 3
FeOnUpSigBiasU <- sum(FeOnUpSigBias$gene.id %in% sebida.unbi$gene.id)
# 33


# Last create a new data frame with the observed and expected numbers
sexbi.FeOnUpSig <- rbind( c(FeOnUpSigBiasF, FeOnUpSigBiasM, FeOnUpSigBiasU), 
                          c(prop.female * nrow(FeOnUpSigBias), prop.male * nrow(FeOnUpSigBias), prop.unbi * nrow(FeOnUpSigBias)))
# Add names to the columns and rows
colnames(sexbi.FeOnUpSig) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FeOnUpSig) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FeOnUpSig["Obs",], p = sexbi.FeOnUpSig["Exp",], rescale.p = T)
# P = 1.379e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify the observed sex-biased classes for up-regulated transcripts in the FLX effect class
FeOnDownSig <- subset(SwapFemaleOnly.data, g35down == "Sig" | g39down == "Sig" | g40down == "Sig" | g42down == "Sig" | g44down == "Sig" | g49down == "Sig" |
                                           g56down == "Sig" | g57down == "Sig" | g59down == "Sig" | g62down == "Sig" | g63down == "Sig")
# 40


# Subset the data to only include transcripts that overlap between the FLX up-regulated and sebida database
FeOnDownSigBias <- subset(FeOnDownSig, FeOnDownSig$gene.id %in% sebida.data$gene.id)
# 28

# Get the sum for the transcripts which are significant for three bias classes
FeOnDownSigBiasF <- sum(FeOnDownSigBias$gene.id %in% sebida.febi$gene.id)
# 1
FeOnDownSigBiasM <- sum(FeOnDownSigBias$gene.id %in% sebida.mabi$gene.id)
# 12
FeOnDownSigBiasU <- sum(FeOnDownSigBias$gene.id %in% sebida.unbi$gene.id)
# 15


# Last create a new data frame with the observed and expected numbers
sexbi.FeOnDownSig <- rbind( c(FeOnDownSigBiasF, FeOnDownSigBiasM, FeOnDownSigBiasU), 
                            c(prop.female * nrow(FeOnDownSigBias), prop.male * nrow(FeOnDownSigBias), prop.unbi * nrow(FeOnDownSigBias)))
# Add names to the columns and rows
colnames(sexbi.FeOnDownSig) <- c("Female", "Male", "Unbiased")
rownames(sexbi.FeOnDownSig) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(sexbi.FeOnDownSig["Obs",], p = sexbi.FeOnDownSig["Exp",], rescale.p = T)
# P = 0.001038




#########################################  PLOT DATA  #########################################


# First collect the data in a matrix
FemaleOnly.plot <- as.matrix(cbind(t(sexbi.FeOnUpSig), t(sexbi.FeOnDownSig)))

# Next plot
pdf("FemaleOnly.pdf", family = "serif", width = 9, height = 7)

par(mar = c(4, 5, 2, 2))
par(lwd = 2)
barplot(FemaleOnly.plot[,c(1, 3)], col = "#FFADAD", border = NA, beside = TRUE, space = c(0.1, 1.2), 
        ylim = c(0, 40), las = 1, cex.axis = 1.2, xaxt = "n", main = "Exclusive female groups", cex.main = 1.8)
barplot(FemaleOnly.plot[,c(2, 4)], beside = TRUE, col = rgb(0, 0, 0, 0), space = c(0.1, 1.2), 
        xaxt = "n", yaxt = "n", add = T)
legend("topright", legend = c("Observed", "Expected"), pch = c(15, 22), col = c("#FFADAD", "black"))
mtext(c("Female bias", "Male bias", "Unbiased"), side = 1, line = 0.3, at = c(1.7, 2.8, 3.9, 6.1, 7.2, 8.3), cex = 1.1)
mtext(c("Up regulated", "Down regulated"), side = 1, line = 1.9, at = c(2.8, 7.2), cex = 1.5)

dev.off()

mtext("A", side = 3, line = 0.8, at = 0, cex = 1.6)
mtext("*", side = 3, line = -15.5, at = 6.1, cex = 1.2)
mtext("*", side = 3, line = -4.8, at = 7.2, cex = 1.2)


