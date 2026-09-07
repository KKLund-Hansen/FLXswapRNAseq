################################################################################################
######################### CHROMOSOME LOCATION FLX SWAP MALES - 7 GROUPS ########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/ChromosomeLocation")

# Read in files with data
FLXswap7gM.data <- read.table(file = "FLXswapMres7G.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


## Under/over representation on any chromosome

# To be able to calculate the expected number I first have collect the number of overall transcripts located on each chromosome

# Each chromosome
chrXM <- sum(grepl("X", FLXswap7gM.data$chr))
# 2372
chr2M <- sum(grepl("2L|2R", FLXswap7gM.data$chr))
# 6091
chr3M <- sum(grepl("3L|3R", FLXswap7gM.data$chr))
# 6865
chr4M <-sum(grepl("\\<4\\>", FLXswap7gM.data$chr))
# 107

# And then collect
chrM <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FLXswap7gM.data$chr))
# 15435



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group A 
GAMupsig <- subset(FLXswap7gM.data, gAup == "Sig")
# 5

# Next collect the observed number for each chromosome specifically
GAMupsigX <- sum(grepl("X", GAMupsig$chr))
# 1
GAMupsig2 <- sum(grepl("2L|2R", GAMupsig$chr))
# 2
GAMupsig3 <- sum(grepl("3L|3R", GAMupsig$chr))
# 2
GAMupsig4 <- sum(grepl("\\<4\\>", GAMupsig$chr))
# 0

# And for all chromosomes collectedly 
GAMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GAMupsig$chr))
# 5


# Last create a new data frame with the observed and expected numbers 
chrloc.GAMup <- rbind( c(GAMupsigX, GAMupsig2, GAMupsig3, GAMupsig4),
                       c(chrXM/chrM*GAMupsigChr, chr2M/chrM*GAMupsigChr, chr3M/chrM*GAMupsigChr, chr4M/chrM*GAMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GAMup) <- c("X", "2", "3", "4")
rownames(chrloc.GAMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GAMup["Obs",], p = chrloc.GAMup["Exp",], rescale.p = T)
# P = 0.9884

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GAMupF <- rbind( c(GAMupsigX, GAMupsig2, GAMupsig3, GAMupsig4),
                        c(chrXM-GAMupsigX, chr2M-GAMupsig2, chr3M-GAMupsig3, chr4M-GAMupsig4))
# Add names to the columns and rows
colnames(chrloc.GAMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GAMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GAMupF)
# P = 1

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group A
GAMdownsig <- subset(FLXswap7gM.data, gAdown == "Sig")
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
# P = 0.8934

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GAMdownF <- rbind( c(GAMdownsigX, GAMdownsig2, GAMdownsig3, GAMdownsig4),
                          c(chrXM-GAMdownsigX, chr2M-GAMdownsig2, chr3M-GAMdownsig3, chr4M-GAMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GAMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GAMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GAMdownF)
# P = 0.7764

# We get the same results, so that is good




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group B 
GBMupsig <- subset(FLXswap7gM.data, gBup == "Sig")
# 51

# Next collect the observed number for each chromosome specifically
GBMupsigX <- sum(grepl("X", GBMupsig$chr))
# 4
GBMupsig2 <- sum(grepl("2L|2R", GBMupsig$chr))
# 24
GBMupsig3 <- sum(grepl("3L|3R", GBMupsig$chr))
# 21
GBMupsig4 <- sum(grepl("\\<4\\>", GBMupsig$chr))
# 2

# And for all chromosomes collectedly 
GBMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBMupsig$chr))
# 51


# Last create a new data frame with the observed and expected numbers 
chrloc.GBMup <- rbind( c(GBMupsigX, GBMupsig2, GBMupsig3, GBMupsig4),
                       c(chrXM/chrM*GBMupsigChr, chr2M/chrM*GBMupsigChr, chr3M/chrM*GBMupsigChr, chr4M/chrM*GBMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GBMup) <- c("X", "2", "3", "4")
rownames(chrloc.GBMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBMup["Obs",], p = chrloc.GBMup["Exp",], rescale.p = T)
# P = 0.01533

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBMupF <- rbind( c(GBMupsigX, GBMupsig2, GBMupsig3, GBMupsig4),
                        c(chrXM-GBMupsigX, chr2M-GBMupsig2, chr3M-GBMupsig3, chr4M-GBMupsig4))
# Add names to the columns and rows
colnames(chrloc.GBMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GBMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBMupF)
# P = 0.03926

# It goes from significant to non-significant



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group B
GBMdownsig <- subset(FLXswap7gM.data, gBdown == "Sig")
# 56

# Next collect the observed number for each chromosome specifically
GBMdownsigX <- sum(grepl("X", GBMdownsig$chr))
# 7
GBMdownsig2 <- sum(grepl("2L|2R", GBMdownsig$chr))
# 24
GBMdownsig3 <- sum(grepl("3L|3R", GBMdownsig$chr))
# 25
GBMdownsig4 <- sum(grepl("\\<4\\>", GBMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GBMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBMdownsig$chr))
# 56


# Last create a new data frame with the observed and expected numbers 
chrloc.GBMdown <- rbind( c(GBMdownsigX, GBMdownsig2, GBMdownsig3, GBMdownsig4),
                         c(chrXM/chrM*GBMdownsigChr, chr2M/chrM*GBMdownsigChr, chr3M/chrM*GBMdownsigChr, chr4M/chrM*GBMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GBMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GBMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBMdown["Obs",], p = chrloc.GBMdown["Exp",], rescale.p = T)
# P = 0.837

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBMdownF <- rbind( c(GBMdownsigX, GBMdownsig2, GBMdownsig3, GBMdownsig4),
                          c(chrXM-GBMdownsigX, chr2M-GBMdownsig2, chr3M-GBMdownsig3, chr4M-GBMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GBMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GBMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBMdownF)
# P = 0.9041

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group C
GCMupsig <- subset(FLXswap7gM.data, gCup == "Sig")
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
# P = 0.6703

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCMupF <- rbind( c(GCMupsigX, GCMupsig2, GCMupsig3, GCMupsig4),
                        c(chrXM-GCMupsigX, chr2M-GCMupsig2, chr3M-GCMupsig3, chr4M-GCMupsig4))
# Add names to the columns and rows
colnames(chrloc.GCMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GCMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCMupF)
# P = 0.6614

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group C
GCMdownsig <- subset(FLXswap7gM.data, gCdown == "Sig")
# 38

# Next collect the observed number for each chromosome specifically
GCMdownsigX <- sum(grepl("X", GCMdownsig$chr))
# 7
GCMdownsig2 <- sum(grepl("2L|2R", GCMdownsig$chr))
# 15
GCMdownsig3 <- sum(grepl("3L|3R", GCMdownsig$chr))
# 7
GCMdownsig4 <- sum(grepl("\\<4\\>", GCMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GCMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GCMdownsig$chr))
# 38


# Last create a new data frame with the observed and expected numbers 
chrloc.GCMdown <- rbind( c(GCMdownsigX, GCMdownsig2, GCMdownsig3, GCMdownsig4),
                         c(chrXM/chrM*GCMdownsigChr, chr2M/chrM*GCMdownsigChr, chr3M/chrM*GCMdownsigChr, chr4M/chrM*GCMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GCMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GCMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GCMdown["Obs",], p = chrloc.GCMdown["Exp",], rescale.p = T)
# P = 0.9096

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCMdownF <- rbind( c(GCMdownsigX, GCMdownsig2, GCMdownsig3, GCMdownsig4),
                          c(chrXM-GCMdownsigX, chr2M-GCMdownsig2, chr3M-GCMdownsig3, chr4M-GCMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GCMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GCMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCMdownF)
# P = 0.8771

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group D
GDMupsig <- subset(FLXswap7gM.data, gDup == "Sig")
# 12

# Next collect the observed number for each chromosome specifically
GDMupsigX <- sum(grepl("X", GDMupsig$chr))
# 1
GDMupsig2 <- sum(grepl("2L|2R", GDMupsig$chr))
# 3
GDMupsig3 <- sum(grepl("3L|3R", GDMupsig$chr))
# 7
GDMupsig4 <- sum(grepl("\\<4\\>", GDMupsig$chr))
# 1

# And for all chromosomes collectedly 
GDMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GDMupsig$chr))
# 12


# Last create a new data frame with the observed and expected numbers 
chrloc.GDMup <- rbind( c(GDMupsigX, GDMupsig2, GDMupsig3, GDMupsig4),
                       c(chrXM/chrM*GDMupsigChr, chr2M/chrM*GDMupsigChr, chr3M/chrM*GDMupsigChr, chr4M/chrM*GDMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GDMup) <- c("X", "2", "3", "4")
rownames(chrloc.GDMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GDMup["Obs",], p = chrloc.GDMup["Exp",], rescale.p = T)
# P = 0.008705

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDMupF <- rbind( c(GDMupsigX, GDMupsig2, GDMupsig3, GDMupsig4),
                        c(chrXM-GDMupsigX, chr2M-GDMupsig2, chr3M-GDMupsig3, chr4M-GDMupsig4))
# Add names to the columns and rows
colnames(chrloc.GDMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GDMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GDMupF)
# P = 0.07799

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group D
GDMdownsig <- subset(FLXswap7gM.data, gDdown == "Sig")
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
# P = 0.9559

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDMdownF <- rbind( c(GDMdownsigX, GDMdownsig2, GDMdownsig3, GDMdownsig4),
                          c(chrXM-GDMdownsigX, chr2M-GDMdownsig2, chr3M-GDMdownsig3, chr4M-GDMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GDMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GDMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Misher's exact test
fisher.test(chrloc.GDMdownF)
# P = 0.9501

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group E 
GEMupsig <- subset(FLXswap7gM.data, gEup == "Sig")
# 20

# Next collect the observed number for each chromosome specifically
GEMupsigX <- sum(grepl("X", GEMupsig$chr))
# 0
GEMupsig2 <- sum(grepl("2L|2R", GEMupsig$chr))
# 13
GEMupsig3 <- sum(grepl("3L|3R", GEMupsig$chr))
# 6
GEMupsig4 <- sum(grepl("\\<4\\>", GEMupsig$chr))
# 0

# And for all chromosomes collectedly 
GEMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GEMupsig$chr))
# 19


# Last create a new data frame with the observed and expected numbers 
chrloc.GEMup <- rbind( c(GEMupsigX, GEMupsig2, GEMupsig3, GEMupsig4),
                       c(chrXM/chrM*GEMupsigChr, chr2M/chrM*GEMupsigChr, chr3M/chrM*GEMupsigChr, chr4M/chrM*GEMupsigChr))
# Add names to the columns and rows
colnames(chrloc.GEMup) <- c("X", "2", "3", "4")
rownames(chrloc.GEMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GEMup["Obs",], p = chrloc.GEMup["Exp",], rescale.p = T)
# P = 0.05033

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEMupF <- rbind( c(GEMupsigX, GEMupsig2, GEMupsig3, GEMupsig4),
                        c(chrXM-GEMupsigX, chr2M-GEMupsig2, chr3M-GEMupsig3, chr4M-GEMupsig4))
# Add names to the columns and rows
colnames(chrloc.GEMupF) <- c("X", "2", "3", "4")
rownames(chrloc.GEMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEMupF)
# P = 0.04883

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group E
GEMdownsig <- subset(FLXswap7gM.data, gEdown == "Sig")
# 26

# Next collect the observed number for each chromosome specifically
GEMdownsigX <- sum(grepl("X", GEMdownsig$chr))
# 4
GEMdownsig2 <- sum(grepl("2L|2R", GEMdownsig$chr))
# 13
GEMdownsig3 <- sum(grepl("3L|3R", GEMdownsig$chr))
# 9
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
# P = 0.6868

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEMdownF <- rbind( c(GEMdownsigX, GEMdownsig2, GEMdownsig3, GEMdownsig4),
                          c(chrXM-GEMdownsigX, chr2M-GEMdownsig2, chr3M-GEMdownsig3, chr4M-GEMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GEMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GEMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEMdownF)
# P = 0.5838

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group F
GFMupsig <- subset(FLXswap7gM.data, gFup == "Sig")
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
# P = 0.9421

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
GFMdownsig <- subset(FLXswap7gM.data, gFdown == "Sig")
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
# P = 0.8582

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GFMdownF <- rbind( c(GFMdownsigX, GFMdownsig2, GFMdownsig3, GFMdownsig4),
                          c(chrXM-GFMdownsigX, chr2M-GFMdownsig2, chr3M-GFMdownsig3, chr4M-GFMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GFMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GFMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GFMdownF)
# P = 0.558

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group G 
GGMupsig <- subset(FLXswap7gM.data, gGup == "Sig")
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
# P = 0.9003

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
GGMdownsig <- subset(FLXswap7gM.data, gGdown == "Sig")
# 22

# Next collect the observed number for each chromosome specifically
GGMdownsigX <- sum(grepl("X", GGMdownsig$chr))
# 3
GGMdownsig2 <- sum(grepl("2L|2R", GGMdownsig$chr))
# 10
GGMdownsig3 <- sum(grepl("3L|3R", GGMdownsig$chr))
# 9
GGMdownsig4 <- sum(grepl("\\<4\\>", GGMdownsig$chr))
# 0

# And for all chromosomes collectedly 
GGMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GGMdownsig$chr))
# 22


# Last create a new data frame with the observed and expected numbers 
chrloc.GGMdown <- rbind( c(GGMdownsigX, GGMdownsig2, GGMdownsig3, GGMdownsig4),
                         c(chrXM/chrM*GGMdownsigChr, chr2M/chrM*GGMdownsigChr, chr3M/chrM*GGMdownsigChr, chr4M/chrM*GGMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GGMdown) <- c("X", "2", "3", "4")
rownames(chrloc.GGMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GGMdown["Obs",], p = chrloc.GGMdown["Exp",], rescale.p = T)
# P = 0.9279

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GGMdownF <- rbind( c(GGMdownsigX, GGMdownsig2, GGMdownsig3, GGMdownsig4),
                          c(chrXM-GGMdownsigX, chr2M-GGMdownsig2, chr3M-GGMdownsig3, chr4M-GGMdownsig4))
# Add names to the columns and rows
colnames(chrloc.GGMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GGMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GGMdownF)
# P = 0.9301

# We get the same results, so that is good

