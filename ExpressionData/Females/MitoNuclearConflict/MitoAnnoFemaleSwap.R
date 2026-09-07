################################################################################################
########################## MITO-ANNOTATED TRANSCRIPTS FLX SWAP FEMALES #########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/MitoNuclearConflict")

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)
MA.data <- read.table(file = "mitoannotated.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First I'll subset the mito-annotated data set to match the female data set
MAfemaleswap <- subset(MA.data, MA.data$From %in% FLXswapF.data$gene.id)
# 463 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-annotated transcripts out of all my transcripts
MApropF <- length(MAfemaleswap$From)/length(FLXswapF.data$gene.id)
# 0.03

# And for the chisq.test I also have to calculate
NotMApropF <- 1-MApropF
# 0.97



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# How many mito-annotated transcripts are there?
FvCfupsigMA <- sum(FvCfupsig$gene.id %in% MAfemaleswap$From)
# 0
# How many are not?
FvCfupsigNotMA <- length(FvCfupsig$gene.id)-FvCfupsigMA
# 35

# The number of expected mito-annotated transcripts
length(FvCfupsig$gene.id)*MApropF
# 1.045147

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfupsigMA, FvCfupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.2993

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FvCfupF <- rbind(c(FvCfupsigMA, FvCfupsigNotMA), 
                    c(length(MAfemaleswap$From)-FvCfupsigMA, (length(FLXswapF.data$gene.id)-length(MAfemaleswap$From))-FvCfupsigNotMA))
#Add names to the columns and rows
colnames(MA.FvCfupF) <- c("MA", "NotMA")
rownames(MA.FvCfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FvCfupF)
# P = 0.6268

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the FLXvCFM group 
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

# How many mito-annotated transcripts are there?
FvCfdownsigMA <- sum(FvCfdownsig$gene.id %in% MAfemaleswap$From)
# 0
# How many are not?
FvCfdownsigNotMA <- length(FvCfdownsig$gene.id)-FvCfdownsigMA
# 40

# The number of expected mito-annotated transcripts
length(FvCfdownsig$gene.id)*MApropF
# 1.194453

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FvCfdownsigMA, FvCfdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.2672

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FvCfdownF <- rbind(c(FvCfdownsigMA, FvCfdownsigNotMA), 
                      c(length(MAfemaleswap$From)-FvCfdownsigMA, (length(FLXswapF.data$gene.id)-length(MAfemaleswap$From))-FvCfdownsigNotMA))
#Add names to the columns and rows
colnames(MA.FvCfdownF) <- c("MA", "NotMA")
rownames(MA.FvCfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FvCfdownF)
# P = 0.6334

# We get the same results, so that is good



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# How many mito-annotated transcripts are there?
AvSfupsigMA <- sum(AvSfupsig$gene.id %in% MAfemaleswap$From)
# 7
# How many are not?
AvSfupsigNotMA <- length(AvSfupsig$gene.id)-AvSfupsigMA
# 237

# The number of expected mito-annotated transcripts
length(AvSfupsig$gene.id)*MApropF
# 7.286166

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfupsigMA, AvSfupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.9143



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# How many mito-annotated transcripts are there?
AvSfdownsigMA <- sum(AvSfdownsig$gene.id %in% MAfemaleswap$From)
# 3
# How many are not?
AvSfdownsigNotMA <- length(AvSfdownsig$gene.id)-AvSfdownsigMA
# 158

# The number of expected mito-annotated transcripts
length(AvSfdownsig$gene.id)*MApropF
# 4.807675

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(AvSfdownsigMA, AvSfdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.4026

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.AvSfdownF <- rbind(c(AvSfdownsigMA, AvSfdownsigNotMA), 
                    c(length(MAfemaleswap$From)-AvSfdownsigMA, (length(FLXswapF.data$gene.id)-length(MAfemaleswap$From))-AvSfdownsigNotMA))
#Add names to the columns and rows
colnames(MA.AvSfdownF) <- c("MA", "NotMA")
rownames(MA.AvSfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.AvSfdownF)
# P = 0.6371

# We get the same results, so that is good




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# How many mito-annotated transcripts are there?
FsvCsFupsigMA <- sum(FsvCsFupsig$gene.id %in% MAfemaleswap$From)
# 2
# How many are not?
FsvCsFupsigNotMA <- length(FsvCsFupsig$gene.id)-FsvCsFupsigMA
# 135

# The number of expected mito-annotated transcripts
length(FsvCsFupsig$gene.id)*MApropF
# 4.091003

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFupsigMA, FsvCsFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.2939

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FsvCsFupF <- rbind(c(FsvCsFupsigMA, FsvCsFupsigNotMA), 
                      c(length(MAfemaleswap$From)-FsvCsFupsigMA, (length(FLXswapF.data$gene.id)-length(MAfemaleswap$From))-FsvCsFupsigNotMA))
#Add names to the columns and rows
colnames(MA.FsvCsFupF) <- c("MA", "NotMA")
rownames(MA.FsvCsFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FsvCsFupF)
# P = 0.4457

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# How many mito-annotated transcripts are there?
FsvCsFdownsigMA <- sum(FsvCsFdownsig$gene.id %in% MAfemaleswap$From)
# 2
# How many are not?
FsvCsFdownsigNotMA <- length(FsvCsFdownsig$gene.id)-FsvCsFdownsigMA
# 53

# The number of expected mito-annotated transcripts
length(FsvCsFdownsig$gene.id)*MApropF
# 1.642373

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FsvCsFdownsigMA, FsvCsFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 0.7769

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MA.FsvCsFdownF <- rbind(c(FsvCsFdownsigMA, FsvCsFdownsigNotMA), 
                        c(length(MAfemaleswap$From)-FsvCsFdownsigMA, (length(FLXswapF.data$gene.id)-length(MAfemaleswap$From))-FsvCsFdownsigNotMA))
#Add names to the columns and rows
colnames(MA.FsvCsFdownF) <- c("MA", "NotMA")
rownames(MA.FsvCsFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MA.FsvCsFdownF)
# P = 0.6803

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# How many mito-annotated transcripts are there?
FavCaFupsigMA <- sum(FavCaFupsig$gene.id %in% MAfemaleswap$From)
# 38
# How many are not?
FavCaFupsigNotMA <- length(FavCaFupsig$gene.id)-FavCaFupsigMA
# 487

# The number of expected mito-annotated transcripts
length(FavCaFupsig$gene.id)*MApropF
# 15.6772

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFupsigMA,FavCaFupsigNotMA), p = c(MApropF, NotMApropF))
# P = 1.041e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-annotated transcripts among the down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# How many mito-annotated transcripts are there?
FavCaFdownsigMA <- sum(FavCaFdownsig$gene.id %in% MAfemaleswap$From)
# 53
# How many are not?
FavCaFdownsigNotMA <- length(FavCaFdownsig$gene.id)-FavCaFdownsigMA
# 2941

# The number of expected mito-annotated transcripts
length(FavCaFdownsig$gene.id)*MApropF
# 89.40484

# Do the test, with the number of mito-annotated transcripts and not and the proportion of overall transcripts
chisq.test(c(FavCaFdownsigMA, FavCaFdownsigNotMA), p = c(MApropF, NotMApropF))
# P = 9.269e-05


