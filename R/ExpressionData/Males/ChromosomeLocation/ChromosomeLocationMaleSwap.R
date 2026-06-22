################################################################################################
############################## CHROMOSOME LOCTAION FLX SWAP MALES ##############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Males/ChromosomeLocation")

# Read in files with data
RNAswapmale.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


## Under/over representation on any chromosome

# To be able to calculate the expected number I first have collect the number of overall transcripts located on each chromosome

# Each chromosome
chrXM <- sum(grepl("X", RNAswapmale.data$chr))
# 2553
chr2M <- sum(grepl("2L|2R", RNAswapmale.data$chr))
# 6608
chr3M <- sum(grepl("3L|3R", RNAswapmale.data$chr))
# 7350
chr4M <-sum(grepl("\\<4\\>", RNAswapmale.data$chr))
# 111

# And then collect
chrM <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", RNAswapmale.data$chr))
# 16622



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group A 
GAMupsig <- subset(RNAswapmale.data, gA41up == "Sig")
# 129

# Next collect the observed number for each chromosome specifically
GAMupsigX <- sum(grepl("X", GAMupsig$chr))
# 1
GAMupsig2 <- sum(grepl("2L|2R", GAMupsig$chr))
# 1
GAMupsig3 <- sum(grepl("3L|3R", GAMupsig$chr))
# 2
GAMupsig4 <- sum(grepl("\\<4\\>", GAMupsig$chr))
# 0

# And for all chromosomes collectedly 
GAMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GAMupsig$chr))
# 4


# Last create a new data frame with the observed and expected numbers 
chrloc.GAMup <- rbind( c(GAMupsigX, GAMupsig2, GAMupsig3, GAMupsig4),
                       c(chrXM/chrM*GAMupsigChr, chr2M/chrM*GAMupsigChr, chr3M/chrM*GAMupsigChr, chr4M/chrM*GAMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GAMup) <- c("X", "2", "3", "4")
rownames(chrloc.GAMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GAMup["Obs",], p = chrloc.GAMup["Exp",], rescale.p = T)
# P = 0.9149


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GAMupF <- rbind( c(GAMupsigX, GAMupsig2, GAMupsig3, GAMupsig4),
                        c(chrXM-GAMupsigX, chr2M-GAMupsig2, chr3M-GAMupsig3, chr4M-GAMupsig4))
# Add names to the columns and rows
colnames(chrloc.GAMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GAMupF)
# P = 0.8146

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group A
GAMdownsig <- subset(RNAswapmale.data, gA41down == "Sig")
# 16

#Next collect the observed number for each chromosome specifically
GAMdownsigX <- sum(grepl("X", GAMdownsig$chr))
# 3
GAMdownsig2 <- sum(grepl("2L|2R", GAMdownsig$chr))
# 5
GAMdownsig3 <- sum(grepl("3L|3R", GAMdownsig$chr))
# 8
GAMdownsig4 <- sum(grepl("\\<4\\>", GAMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GAMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GAMdownsig$chr))
# 16


# Last create a new data frame with the observed and expected numbers 
chrloc.GAMdown <- rbind( c(GAMdownsigX, GAMdownsig2, GAMdownsig3, GAMdownsig4),
                         c(chrXM/chrM*GAMdownsigChr, chr2M/chrM*GAMdownsigChr, chr3M/chrM*GAMdownsigChr, chr4M/chrM*GAMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GAMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GAMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GAMdown["Obs",], p = chrloc.GAMdown["Exp",], rescale.p = T)
# P = 0.8875


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GAMdownF <- rbind( c(GAMdownsigX, GAMdownsig2, GAMdownsig3, GAMdownsig4),
                          c(chrXM-GAMdownsigX, chr2M-GAMdownsig2, chr3M-GAMdownsig3, chr4M-GAMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GAMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GAMdownF)
# P = 0.7745

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group B 
GBMupsig <- subset(RNAswapmale.data, gB7up == "Sig")
# 54

# Next collect the observed number for each chromosome specifically
GBMupsigX <- sum(grepl("X", GBMupsig$chr))
# 5
GBMupsig2 <- sum(grepl("2L|2R", GBMupsig$chr))
# 24
GBMupsig3 <- sum(grepl("3L|3R", GBMupsig$chr))
# 23
GBMupsig4 <- sum(grepl("\\<4\\>", GBMupsig$chr))
# 2

# And for all chromosomes collectedly 
GBMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBMupsig$chr))
# 54


# Last create a new data frame with the observed and expected numbers 
chrloc.GBMup <- rbind( c(GBMupsigX, GBMupsig2, GBMupsig3, GBMupsig4),
                       c(chrXM/chrM*GBMupsigChr, chr2M/chrM*GBMupsigChr, chr3M/chrM*GBMupsigChr, chr4M/chrM*GBMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GBMup) <- c("X", "2", "3", "4")
rownames(chrloc.GBMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBMup["Obs",], p = chrloc.GBMup["Exp",], rescale.p = T)
# P = 0.02809


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBMupF <- rbind( c(GBMupsigX, GBMupsig2, GBMupsig3, GBMupsig4),
                        c(chrXM-GBMupsigX, chr2M-GBMupsig2, chr3M-GBMupsig3, chr4M-GBMupsig4))
# Add names to the columns and rows
colnames(chrloc.GBMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBMupF)
# P = 0.06703

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group B
GBMdownsig <- subset(RNAswapmale.data, gB7down == "Sig")
# 58

# Next collect the observed number for each chromosome specifically
GBMdownsigX <- sum(grepl("X", GBMdownsig$chr))
# 7
GBMdownsig2 <- sum(grepl("2L|2R", GBMdownsig$chr))
# 24
GBMdownsig3 <- sum(grepl("3L|3R", GBMdownsig$chr))
# 27
GBMdownsig4 <- sum(grepl("\\<4\\>", GBMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GBMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBMdownsig$chr))
# 58


# Last create a new data frame with the observed and expected numbers 
chrloc.GBMdown <- rbind( c(GBMdownsigX, GBMdownsig2, GBMdownsig3, GBMdownsig4),
                         c(chrXM/chrM*GBMdownsigChr, chr2M/chrM*GBMdownsigChr, chr3M/chrM*GBMdownsigChr, chr4M/chrM*GBMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GBMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GBMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBMdown["Obs",], p = chrloc.GBMdown["Exp",], rescale.p = T)
# P = 0.824


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBMdownF <- rbind( c(GBMdownsigX, GBMdownsig2, GBMdownsig3, GBMdownsig4),
                          c(chrXM-GBMdownsigX, chr2M-GBMdownsig2, chr3M-GBMdownsig3, chr4M-GBMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GBMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBMdownF)
# P = 0.8787

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group C
GCMupsig <- subset(RNAswapmale.data, gC6up == "Sig")
# 27

# Next collect the observed number for each chromosome specifically
GCMupsigX <- sum(grepl("X", GCMupsig$chr))
# 2
GCMupsig2 <- sum(grepl("2L|2R", GCMupsig$chr))
# 12
GCMupsig3 <- sum(grepl("3L|3R", GCMupsig$chr))
# 13
GCMupsig4 <- sum(grepl("\\<4\\>", GCMupsig$chr))
# 0

# And for all chromosomes collectedly 
GCMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GCMupsig$chr))
# 27


# Last create a new data frame with the observed and expected numbers 
chrloc.GCMup <- rbind( c(GCMupsigX, GCMupsig2, GCMupsig3, GCMupsig4),
                       c(chrXM/chrM*GCMupsigChr, chr2M/chrM*GCMupsigChr, chr3M/chrM*GCMupsigChr, chr4M/chrM*GCMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GCMup) <- c("X", "2", "3", "4")
rownames(chrloc.GCMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GCMup["Obs",], p = chrloc.GCMup["Exp",], rescale.p = T)
# P = 0.6741


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCMupF <- rbind( c(GCMupsigX, GCMupsig2, GCMupsig3, GCMupsig4),
                        c(chrXM-GCMupsigX, chr2M-GCMupsig2, chr3M-GCMupsig3, chr4M-GCMupsig4))
# Add names to the columns and rows
colnames(chrloc.GCMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCMupF)
# P = 0.6591

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group C
GCMdownsig <- subset(RNAswapmale.data, gC6down == "Sig")
# 37

# Next collect the observed number for each chromosome specifically
GCMdownsigX <- sum(grepl("X", GCMdownsig$chr))
# 7
GCMdownsig2 <- sum(grepl("2L|2R", GCMdownsig$chr))
# 14
GCMdownsig3 <- sum(grepl("3L|3R", GCMdownsig$chr))
# 7
GCMdownsig4 <- sum(grepl("\\<4\\>", GCMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GCMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GCMdownsig$chr))
# 37


# Last create a new data frame with the observed and expected numbers 
chrloc.GCMdown <- rbind( c(GCMdownsigX, GCMdownsig2, GCMdownsig3, GCMdownsig4),
                         c(chrXM/chrM*GCMdownsigChr, chr2M/chrM*GCMdownsigChr, chr3M/chrM*GCMdownsigChr, chr4M/chrM*GCMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GCMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GCMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GCMdown["Obs",], p = chrloc.GCMdown["Exp",], rescale.p = T)
# P = 0.8977


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCMdownF <- rbind( c(GCMdownsigX, GCMdownsig2, GCMdownsig3, GCMdownsig4),
                          c(chrXM-GCMdownsigX, chr2M-GCMdownsig2, chr3M-GCMdownsig3, chr4M-GCMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GCMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCMdownF)
# P = 0.872

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group D
GDMupsig <- subset(RNAswapmale.data, gD22up == "Sig")
# 14

# Next collect the observed number for each chromosome specifically
GDMupsigX <- sum(grepl("X", GDMupsig$chr))
# 1
GDMupsig2 <- sum(grepl("2L|2R", GDMupsig$chr))
# 5
GDMupsig3 <- sum(grepl("3L|3R", GDMupsig$chr))
# 7
GDMupsig4 <- sum(grepl("\\<4\\>", GDMupsig$chr))
# 1

# And for all chromosomes collectedly 
GDMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GDMupsig$chr))
# 14


# Last create a new data frame with the observed and expected numbers 
chrloc.GDMup <- rbind( c(GDMupsigX, GDMupsig2, GDMupsig3, GDMupsig4),
                       c(chrXM/chrM*GDMupsigChr, chr2M/chrM*GDMupsigChr, chr3M/chrM*GDMupsigChr, chr4M/chrM*GDMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GDMup) <- c("X", "2", "3", "4")
rownames(chrloc.GDMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GDMup["Obs",], p = chrloc.GDMup["Exp",], rescale.p = T)
# P = 0.02261


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDMupF <- rbind( c(GDMupsigX, GDMupsig2, GDMupsig3, GDMupsig4),
                        c(chrXM-GDMupsigX, chr2M-GDMupsig2, chr3M-GDMupsig3, chr4M-GDMupsig4))
# Add names to the columns and rows
colnames(chrloc.GDMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GDMupF)
# P = 0.1225

# It went from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group D
GDMdownsig <- subset(RNAswapmale.data, gD22down == "Sig")
# 16

# Next collect the observed number for each chromosome specifically
GDMdownsigX <- sum(grepl("X", GDMdownsig$chr))
# 2
GDMdownsig2 <- sum(grepl("2L|2R", GDMdownsig$chr))
# 6
GDMdownsig3 <- sum(grepl("3L|3R", GDMdownsig$chr))
# 8
GDMdownsig4 <- sum(grepl("\\<4\\>", GDMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GDMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GDMdownsig$chr))
# 16


# Last create a new data frame with the observed and expected numbers 
chrloc.GDMdown <- rbind( c(GDMdownsigX, GDMdownsig2, GDMdownsig3, GDMdownsig4),
                         c(chrXM/chrM*GDMdownsigChr, chr2M/chrM*GDMdownsigChr, chr3M/chrM*GDMdownsigChr, chr4M/chrM*GDMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GDMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GDMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GDMdown["Obs",], p = chrloc.GDMdown["Exp",], rescale.p = T)
# P = 0.9536


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDMdownF <- rbind( c(GDMdownsigX, GDMdownsig2, GDMdownsig3, GDMdownsig4),
                          c(chrXM-GDMdownsigX, chr2M-GDMdownsig2, chr3M-GDMdownsig3, chr4M-GDMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GDMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Misher's exact test
fisher.test(chrloc.GDMdownF)
# P = 0.9496

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group E 
GEMupsig <- subset(RNAswapmale.data, gE4up == "Sig")
# 17

# Next collect the observed number for each chromosome specifically
GEMupsigX <- sum(grepl("X", GEMupsig$chr))
# 0
GEMupsig2 <- sum(grepl("2L|2R", GEMupsig$chr))
# 12
GEMupsig3 <- sum(grepl("3L|3R", GEMupsig$chr))
# 4
GEMupsig4 <- sum(grepl("\\<4\\>", GEMupsig$chr))
# 0

# And for all chromosomes collectedly 
GEMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GEMupsig$chr))
# 16


# Last create a new data frame with the observed and expected numbers 
chrloc.GEMup <- rbind( c(GEMupsigX, GEMupsig2, GEMupsig3, GEMupsig4),
                       c(chrXM/chrM*GEMupsigChr, chr2M/chrM*GEMupsigChr, chr3M/chrM*GEMupsigChr, chr4M/chrM*GEMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GEMup) <- c("X", "2", "3", "4")
rownames(chrloc.GEMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GEMup["Obs",], p = chrloc.GEMup["Exp",], rescale.p = T)
# P = 0.03064


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEMupF <- rbind( c(GEMupsigX, GEMupsig2, GEMupsig3, GEMupsig4),
                        c(chrXM-GEMupsigX, chr2M-GEMupsig2, chr3M-GEMupsig3, chr4M-GEMupsig4))
# Add names to the columns and rows
colnames(chrloc.GEMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEMupF)
# P = 0.03418

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group E
GEMdownsig <- subset(RNAswapmale.data, gE4down == "Sig")
# 26

# Next collect the observed number for each chromosome specifically
GEMdownsigX <- sum(grepl("X", GEMdownsig$chr))
# 4
GEMdownsig2 <- sum(grepl("2L|2R", GEMdownsig$chr))
# 12
GEMdownsig3 <- sum(grepl("3L|3R", GEMdownsig$chr))
# 10
GEMdownsig4 <- sum(grepl("\\<4\\>", GEMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GEMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GEMdownsig$chr))
# 26


# Last create a new data frame with the observed and expected numbers 
chrloc.GEMdown <- rbind( c(GEMdownsigX, GEMdownsig2, GEMdownsig3, GEMdownsig4),
                         c(chrXM/chrM*GEMdownsigChr, chr2M/chrM*GEMdownsigChr, chr3M/chrM*GEMdownsigChr, chr4M/chrM*GEMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GEMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GEMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GEMdown["Obs",], p = chrloc.GEMdown["Exp",], rescale.p = T)
# P = 0.8881


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEMdownF <- rbind( c(GEMdownsigX, GEMdownsig2, GEMdownsig3, GEMdownsig4),
                          c(chrXM-GEMdownsigX, chr2M-GEMdownsig2, chr3M-GEMdownsig3, chr4M-GEMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GEMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEMdownF)
# P = 0.8118

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group F
GFMupsig <- subset(RNAswapmale.data, gF18up == "Sig")
# 2

# Next collect the observed number for each chromosome specifically
GFMupsigX <- sum(grepl("X", GFMupsig$chr))
# 0
GFMupsig2 <- sum(grepl("2L|2R", GFMupsig$chr))
# 1
GFMupsig3 <- sum(grepl("3L|3R", GFMupsig$chr))
# 1
GFMupsig4 <- sum(grepl("\\<4\\>", GFMupsig$chr))
# 0

# And for all chromosomes collectedly 
GFMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GFMupsig$chr))
# 2


# Last create a new data frame with the observed and expected numbers 
chrloc.GFMup <- rbind( c(GFMupsigX, GFMupsig2, GFMupsig3, GFMupsig4),
                       c(chrXM/chrM*GFMupsigChr, chr2M/chrM*GFMupsigChr, chr3M/chrM*GFMupsigChr, chr4M/chrM*GFMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GFMup) <- c("X", "2", "3", "4")
rownames(chrloc.GFMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GFMup["Obs",], p = chrloc.GFMup["Exp",], rescale.p = T)
# P = 0.9426


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GFMupF <- rbind( c(GFMupsigX, GFMupsig2, GFMupsig3, GFMupsig4),
                        c(chrXM-GFMupsigX, chr2M-GFMupsig2, chr3M-GFMupsig3, chr4M-GFMupsig4))
# Add names to the columns and rows
colnames(chrloc.GFMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GFMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GFMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group F
GFMdownsig <- subset(RNAswapmale.data, gF18down == "Sig")
# 3

# Next collect the observed number for each chromosome specifically
GFMdownsigX <- sum(grepl("X", GFMdownsig$chr))
# 1
GFMdownsig2 <- sum(grepl("2L|2R", GFMdownsig$chr))
# 1
GFMdownsig3 <- sum(grepl("3L|3R", GFMdownsig$chr))
# 1
GFMdownsig4 <- sum(grepl("\\<4\\>", GFMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GFMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GFMdownsig$chr))
# 3


# Last create a new data frame with the observed and expected numbers 
chrloc.GFMdown <- rbind( c(GFMdownsigX, GFMdownsig2, GFMdownsig3, GFMdownsig4),
                         c(chrXM/chrM*GFMdownsigChr, chr2M/chrM*GFMdownsigChr, chr3M/chrM*GFMdownsigChr, chr4M/chrM*GFMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GFMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GFMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GFMdown["Obs",], p = chrloc.GFMdown["Exp",], rescale.p = T)
# P = 0.8584


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GFMdownF <- rbind( c(GFMdownsigX, GFMdownsig2, GFMdownsig3, GFMdownsig4),
                          c(chrXM-GFMdownsigX, chr2M-GFMdownsig2, chr3M-GFMdownsig3, chr4M-GFMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GFMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GFMdownF)
# P = 0.5571

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group G 
GGMupsig <- subset(RNAswapmale.data, gG19up == "Sig")
# 10

# Next collect the observed number for each chromosome specifically
GGMupsigX <- sum(grepl("X", GGMupsig$chr))
# 1
GGMupsig2 <- sum(grepl("2L|2R", GGMupsig$chr))
# 5
GGMupsig3 <- sum(grepl("3L|3R", GGMupsig$chr))
# 4
GGMupsig4 <- sum(grepl("\\<4\\>", GGMupsig$chr))
# 0

# And for all chromosomes collectedly 
GGMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GGMupsig$chr))
# 10


# Last create a new data frame with the observed and expected numbers 
chrloc.GGMup <- rbind( c(GGMupsigX, GGMupsig2, GGMupsig3, GGMupsig4),
                       c(chrXM/chrM*GGMupsigChr, chr2M/chrM*GGMupsigChr, chr3M/chrM*GGMupsigChr, chr4M/chrM*GGMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GGMup) <- c("X", "2", "3", "4")
rownames(chrloc.GGMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GGMup["Obs",], p = chrloc.GGMup["Exp",], rescale.p = T)
# P = 0.906


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GGMupF <- rbind( c(GGMupsigX, GGMupsig2, GGMupsig3, GGMupsig4),
                        c(chrXM-GGMupsigX, chr2M-GGMupsig2, chr3M-GGMupsig3, chr4M-GGMupsig4))
# Add names to the columns and rows
colnames(chrloc.GGMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GGMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GGMupF)
# P = 0.9182

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group G
GGMdownsig <- subset(RNAswapmale.data, gG19down == "Sig")
# 19

# Next collect the observed number for each chromosome specifically
GGMdownsigX <- sum(grepl("X", GGMdownsig$chr))
# 3
GGMdownsig2 <- sum(grepl("2L|2R", GGMdownsig$chr))
# 9
GGMdownsig3 <- sum(grepl("3L|3R", GGMdownsig$chr))
# 7
GGMdownsig4 <- sum(grepl("\\<4\\>", GGMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GGMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GGMdownsig$chr))
# 19


# Last create a new data frame with the observed and expected numbers 
chrloc.GGMdown <- rbind( c(GGMdownsigX, GGMdownsig2, GGMdownsig3, GGMdownsig4),
                         c(chrXM/chrM*GGMdownsigChr, chr2M/chrM*GGMdownsigChr, chr3M/chrM*GGMdownsigChr, chr4M/chrM*GGMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GGMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GGMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GGMdown["Obs",], p = chrloc.GGMdown["Exp",], rescale.p = T)
# P = 0.8872


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GGMdownF <- rbind( c(GGMdownsigX, GGMdownsig2, GGMdownsig3, GGMdownsig4),
                          c(chrXM-GGMdownsigX, chr2M-GGMdownsig2, chr3M-GGMdownsig3, chr4M-GGMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GGMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GGMdownF)
# P = 0.8076

# We get the same results, so that is good




######################## GROUP 2 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group 2
G2Mupsig <- subset(RNAswapmale.data, g2up == "Sig")
# 29

# Next collect the observed number for each chromosome specifically
G2MupsigX <- sum(grepl("X", G2Mupsig$chr))
# 3
G2Mupsig2 <- sum(grepl("2L|2R", G2Mupsig$chr))
# 12
G2Mupsig3 <- sum(grepl("3L|3R", G2Mupsig$chr))
# 14
G2Mupsig4 <- sum(grepl("\\<4\\>", G2Mupsig$chr))
# 0

# And for all chromosomes collectedly 
G2MupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G2Mupsig$chr))
# 29


# Last create a new data frame with the observed and expected numbers 
chrloc.G2Mup <- rbind( c(G2MupsigX, G2Mupsig2, G2Mupsig3, G2Mupsig4),
                       c(chrXM/chrM*G2MupsigChr, chr2M/chrM*G2MupsigChr, chr3M/chrM*G2MupsigChr, chr4M/chrM*G2MupsigChr))
# Add names to the columns and rows
colnames(chrloc.G2Mup) <- c("X", "2", "3", "4")
rownames(chrloc.G2Mup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G2Mup["Obs",], p = chrloc.G2Mup["Exp",], rescale.p = T)
# P = 0.8505


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G2MupF <- rbind( c(G2MupsigX, G2Mupsig2, G2Mupsig3, G2Mupsig4),
                        c(chrXM-G2MupsigX, chr2M-G2Mupsig2, chr3M-G2Mupsig3, chr4M-G2Mupsig4))
# Add names to the columns and rows
colnames(chrloc.G2MupF) <- c("X", "2", "3", "4")
rownames(chrloc.G2MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G2MupF)
# P = 0.8759

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group 2
G2Mdownsig <- subset(RNAswapmale.data, g2down == "Sig")
# 36

# Next collect the observed number for each chromosome specifically
G2MdownsigX <- sum(grepl("X", G2Mdownsig$chr))
# 3
G2Mdownsig2 <- sum(grepl("2L|2R", G2Mdownsig$chr))
# 14
G2Mdownsig3 <- sum(grepl("3L|3R", G2Mdownsig$chr))
# 19
G2Mdownsig4 <- sum(grepl("\\<4\\>", G2Mdownsig$chr))
# 0

# And for all chromosomes collectedly 
G2MdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G2Mdownsig$chr))
# 36


# Last create a new data frame with the observed and expected numbers 
chrloc.G2Mdown <- rbind( c(G2MdownsigX, G2Mdownsig2, G2Mdownsig3, G2Mdownsig4),
                         c(chrXM/chrM*G2MdownsigChr, chr2M/chrM*G2MdownsigChr, chr3M/chrM*G2MdownsigChr, chr4M/chrM*G2MdownsigChr))
# Add names to the columns and rows
colnames(chrloc.G2Mdown) <- c("X", "2", "3", "4")
rownames(chrloc.G2Mdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G2Mdown["Obs",], p = chrloc.G2Mdown["Exp",], rescale.p = T)
# P = 0.5723


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G2MdownF <- rbind( c(G2MdownsigX, G2Mdownsig2, G2Mdownsig3, G2Mdownsig4),
                          c(chrXM-G2MdownsigX, chr2M-G2Mdownsig2, chr3M-G2Mdownsig3, chr4M-G2Mdownsig4))
# Add names to the columns and rows
colnames(chrloc.G2MdownF) <- c("X", "2", "3", "4")
rownames(chrloc.G2MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G2MdownF)
# P = 0.5702

# We get the same results, so that is good




######################## GROUP 11 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group 11
G11Mupsig <- subset(RNAswapmale.data, g11up == "Sig")
# 21

# Next collect the observed number for each chromosome specifically
G11MupsigX <- sum(grepl("X", G11Mupsig$chr))
# 1
G11Mupsig2 <- sum(grepl("2L|2R", G11Mupsig$chr))
# 8
G11Mupsig3 <- sum(grepl("3L|3R", G11Mupsig$chr))
# 12
G11Mupsig4 <- sum(grepl("\\<4\\>", G11Mupsig$chr))
# 0

# And for all chromosomes collectedly 
G11MupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G11Mupsig$chr))
# 21


# Last create a new data frame with the observed and expected numbers 
chrloc.G11Mup <- rbind( c(G11MupsigX, G11Mupsig2, G11Mupsig3, G11Mupsig4),
                        c(chrXM/chrM*G11MupsigChr, chr2M/chrM*G11MupsigChr, chr3M/chrM*G11MupsigChr, chr4M/chrM*G11MupsigChr))
# Add names to the columns and rows
colnames(chrloc.G11Mup) <- c("X", "2", "3", "4")
rownames(chrloc.G11Mup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G11Mup["Obs",], p = chrloc.G11Mup["Exp",], rescale.p = T)
# P = 0.4783


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G11MupF <- rbind( c(G11MupsigX, G11Mupsig2, G11Mupsig3, G11Mupsig4),
                         c(chrXM-G11MupsigX, chr2M-G11Mupsig2, chr3M-G11Mupsig3, chr4M-G11Mupsig4))
# Add names to the columns and rows
colnames(chrloc.G11MupF) <- c("X", "2", "3", "4")
rownames(chrloc.G11MupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G11MupF)
# P = 0.4505

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group 11
G11Mdownsig <- subset(RNAswapmale.data, g11down == "Sig")
# 7

# Next collect the observed number for each chromosome specifically
G11MdownsigX <- sum(grepl("X", G11Mdownsig$chr))
# 1
G11Mdownsig2 <- sum(grepl("2L|2R", G11Mdownsig$chr))
# 3
G11Mdownsig3 <- sum(grepl("3L|3R", G11Mdownsig$chr))
# 3
G11Mdownsig4 <- sum(grepl("\\<4\\>", G11Mdownsig$chr))
# 0

# And for all chromosomes collectedly 
G11MdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G11Mdownsig$chr))
# 7


# Last create a new data frame with the observed and expected numbers 
chrloc.G11Mdown <- rbind( c(G11MdownsigX, G11Mdownsig2, G11Mdownsig3, G11Mdownsig4),
                          c(chrXM/chrM*G11MdownsigChr, chr2M/chrM*G11MdownsigChr, chr3M/chrM*G11MdownsigChr, chr4M/chrM*G11MdownsigChr))
# Add names to the columns and rows
colnames(chrloc.G11Mdown) <- c("X", "2", "3", "4")
rownames(chrloc.G11Mdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G11Mdown["Obs",], p = chrloc.G11Mdown["Exp",], rescale.p = T)
# P = 0.995


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G11MdownF <- rbind( c(G11MdownsigX, G11Mdownsig2, G11Mdownsig3, G11Mdownsig4),
                           c(chrXM-G11MdownsigX, chr2M-G11Mdownsig2, chr3M-G11Mdownsig3, chr4M-G11Mdownsig4))
# Add names to the columns and rows
colnames(chrloc.G11MdownF) <- c("X", "2", "3", "4")
rownames(chrloc.G11MdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G11MdownF)
# P = 1

# We get the same results, so that is good

