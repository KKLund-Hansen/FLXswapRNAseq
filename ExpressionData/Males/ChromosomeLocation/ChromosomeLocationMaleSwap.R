################################################################################################
############################## CHROMOSOME LOCATION FLX SWAP MALES ##############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Males/ChromosomeLocation")

# Read in files with data
FLXswapM.data <- read.table(file = "FLXswapMres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


## Under/over representation on any chromosome

# To be able to calculate the expected number I first have collect the number of overall transcripts located on each chromosome

# Each chromosome
chrXM <- sum(grepl("X", FLXswapM.data$chr))
# 2372
chr2M <- sum(grepl("2L|2R", FLXswapM.data$chr))
# 6091
chr3M <- sum(grepl("3L|3R", FLXswapM.data$chr))
# 6865
chr4M <-sum(grepl("\\<4\\>", FLXswapM.data$chr))
# 107

# And then collect
chrM <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FLXswapM.data$chr))
# 15435



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in FLX vs CFM
FvCmupsig <- subset(FLXswapM.data, FLX.CFMup == "Sig")
# 57

# Next collect the observed number for each chromosome specifically
FvCmupsigX <- sum(grepl("X", FvCmupsig$chr))
# 13
FvCmupsig2 <- sum(grepl("2L|2R", FvCmupsig$chr))
# 25
FvCmupsig3 <- sum(grepl("3L|3R", FvCmupsig$chr))
# 17
FvCmupsig4 <- sum(grepl("\\<4\\>", FvCmupsig$chr))
# 0

# And for all chromosomes collectedly 
FvCmupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FvCmupsig$chr))
# 55


# Last create a new data frame with the observed and expected numbers 
chrloc.FvCmup <- rbind( c(FvCmupsigX, FvCmupsig2, FvCmupsig3, FvCmupsig4),
                        c(chrXM/chrM*FvCmupsigChr, chr2M/chrM*FvCmupsigChr, chr3M/chrM*FvCmupsigChr, chr4M/chrM*FvCmupsigChr))
# Add names to the columns and rows
colnames(chrloc.FvCmup) <- c("X", "2", "3", "4")
rownames(chrloc.FvCmup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FvCmup["Obs",], p = chrloc.FvCmup["Exp",], rescale.p = T)
# P = 0.1325

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FvCmupF <- rbind( c(FvCmupsigX, FvCmupsig2, FvCmupsig3, FvCmupsig4),
                         c(chrXM-FvCmupsigX, chr2M-FvCmupsig2, chr3M-FvCmupsig3, chr4M-FvCmupsig4))
# Add names to the columns and rows
colnames(chrloc.FvCmupF) <- c("X", "2", "3", "4")
rownames(chrloc.FvCmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FvCmupF)
# P = 0.1398

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in FLX vs CFM
FvCmdownsig <- subset(FLXswapM.data, FLX.CFMdown == "Sig")
# 57

#Next collect the observed number for each chromosome specifically
FvCmdownsigX <- sum(grepl("X", FvCmdownsig$chr))
# 5
FvCmdownsig2 <- sum(grepl("2L|2R", FvCmdownsig$chr))
# 26
FvCmdownsig3 <- sum(grepl("3L|3R", FvCmdownsig$chr))
# 26
FvCmdownsig4 <- sum(grepl("\\<4\\>", FvCmdownsig$chr))
# 0

# And for all chromosomes collectedly 
FvCmdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FvCmdownsig$chr))
# 57


# Last create a new data frame with the observed and expected numbers 
chrloc.FvCmdown <- rbind( c(FvCmdownsigX, FvCmdownsig2, FvCmdownsig3, FvCmdownsig4),
                          c(chrXM/chrM*FvCmdownsigChr, chr2M/chrM*FvCmdownsigChr, chr3M/chrM*FvCmdownsigChr, chr4M/chrM*FvCmdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FvCmdown) <- c("X", "2", "3", "4")
rownames(chrloc.FvCmdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FvCmdown["Obs",], p = chrloc.FvCmdown["Exp",], rescale.p = T)
# P = 0.4624

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FvCmdownF <- rbind( c(FvCmdownsigX, FvCmdownsig2, FvCmdownsig3, FvCmdownsig4),
                           c(chrXM-FvCmdownsigX, chr2M-FvCmdownsig2, chr3M-FvCmdownsig3, chr4M-FvCmdownsig4))
# Add names to the columns and rows
colnames(chrloc.FvCmdownF) <- c("X", "2", "3", "4")
rownames(chrloc.FvCmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FvCmdownF)
# P = 0.5341

# We get the same results, so that is good



######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the autovsex group 
AvSmupsig <- subset(FLXswapM.data, auto.sexup == "Sig")
# 215

# Next collect the observed number for each chromosome specifically
AvSmupsigX <- sum(grepl("X", AvSmupsig$chr))
# 20
AvSmupsig2 <- sum(grepl("2L|2R", AvSmupsig$chr))
# 97
AvSmupsig3 <- sum(grepl("3L|3R", AvSmupsig$chr))
# 98
AvSmupsig4 <- sum(grepl("\\<4\\>", AvSmupsig$chr))
# 0

# And for all chromosomes collectedly 
AvSmupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", AvSmupsig$chr))
# 215


# Last create a new data frame with the observed and expected numbers 
chrloc.AvSmup <- rbind( c(AvSmupsigX, AvSmupsig2, AvSmupsig3, AvSmupsig4),
                        c(chrXM/chrM*AvSmupsigChr, chr2M/chrM*AvSmupsigChr, chr3M/chrM*AvSmupsigChr, chr4M/chrM*AvSmupsigChr))
# Add names to the columns and rows
colnames(chrloc.AvSmup) <- c("X", "2", "3", "4")
rownames(chrloc.AvSmup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.AvSmup["Obs",], p = chrloc.AvSmup["Exp",], rescale.p = T)
# P = 0.03778

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.AvSmupF <- rbind( c(AvSmupsigX, AvSmupsig2, AvSmupsig3, AvSmupsig4),
                         c(chrXM-AvSmupsigX, chr2M-AvSmupsig2, chr3M-AvSmupsig3, chr4M-AvSmupsig4))
# Add names to the columns and rows
colnames(chrloc.AvSmupF) <- c("X", "2", "3", "4")
rownames(chrloc.AvSmupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.AvSmupF)
# P = 0.03525

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the autovsex group 
AvSmdownsig <- subset(FLXswapM.data, auto.sexdown == "Sig")
# 262

# Next collect the observed number for each chromosome specifically
AvSmdownsigX <- sum(grepl("X", AvSmdownsig$chr))
# 16
AvSmdownsig2 <- sum(grepl("2L|2R", AvSmdownsig$chr))
# 142
AvSmdownsig3 <- sum(grepl("3L|3R", AvSmdownsig$chr))
# 101
AvSmdownsig4 <- sum(grepl("\\<4\\>", AvSmdownsig$chr))
# 2

# And for all chromosomes collectedly 
AvSmdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", AvSmdownsig$chr))
# 261


# Last create a new data frame with the observed and expected numbers 
chrloc.AvSmdown <- rbind( c(AvSmdownsigX, AvSmdownsig2, AvSmdownsig3, AvSmdownsig4),
                          c(chrXM/chrM*AvSmdownsigChr, chr2M/chrM*AvSmdownsigChr, chr3M/chrM*AvSmdownsigChr, chr4M/chrM*AvSmdownsigChr))
# Add names to the columns and rows
colnames(chrloc.AvSmdown) <- c("X", "2", "3", "4")
rownames(chrloc.AvSmdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.AvSmdown["Obs",], p = chrloc.AvSmdown["Exp",], rescale.p = T)
# P = 7.557e-07

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.AvSmdownF <- rbind( c(AvSmdownsigX, AvSmdownsig2, AvSmdownsig3, AvSmdownsig4),
                           c(chrXM-AvSmdownsigX, chr2M-AvSmdownsig2, chr3M-AvSmdownsig3, chr4M-AvSmdownsig4))
# Add names to the columns and rows
colnames(chrloc.AvSmdownF) <- c("X", "2", "3", "4")
rownames(chrloc.AvSmdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.AvSmdownF)
# P = 1.616e-07

# We get the same results, so that is good




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsMupsig <- subset(FLXswapM.data, FLXs.CFMsup == "Sig")
# 56

# Next collect the observed number for each chromosome specifically
FsvCsMupsigX <- sum(grepl("X", FsvCsMupsig$chr))
# 6
FsvCsMupsig2 <- sum(grepl("2L|2R", FsvCsMupsig$chr))
# 22
FsvCsMupsig3 <- sum(grepl("3L|3R", FsvCsMupsig$chr))
# 28
FsvCsMupsig4 <- sum(grepl("\\<4\\>", FsvCsMupsig$chr))
# 0

# And for all chromosomes collectedly 
FsvCsMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FsvCsMupsig$chr))
# 56


# Last create a new data frame with the observed and expected numbers 
chrloc.FsvCsMup <- rbind( c(FsvCsMupsigX, FsvCsMupsig2, FsvCsMupsig3, FsvCsMupsig4),
                          c(chrXM/chrM*FsvCsMupsigChr, chr2M/chrM*FsvCsMupsigChr, chr3M/chrM*FsvCsMupsigChr, chr4M/chrM*FsvCsMupsigChr))
# Add names to the columns and rows
colnames(chrloc.FsvCsMup) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FsvCsMup["Obs",], p = chrloc.FsvCsMup["Exp",], rescale.p = T)
# P = 0.6681

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FsvCsMupF <- rbind( c(FsvCsMupsigX, FsvCsMupsig2, FsvCsMupsig3, FsvCsMupsig4),
                           c(chrXM-FsvCsMupsigX, chr2M-FsvCsMupsig2, chr3M-FsvCsMupsig3, chr4M-FsvCsMupsig4))
# Add names to the columns and rows
colnames(chrloc.FsvCsMupF) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FsvCsMupF)
# P = 0.7293

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsMdownsig <- subset(FLXswapM.data, FLXs.CFMsdown == "Sig")
# 49

# Next collect the observed number for each chromosome specifically
FsvCsMdownsigX <- sum(grepl("X", FsvCsMdownsig$chr))
# 8
FsvCsMdownsig2 <- sum(grepl("2L|2R", FsvCsMdownsig$chr))
# 17
FsvCsMdownsig3 <- sum(grepl("3L|3R", FsvCsMdownsig$chr))
# 24
FsvCsMdownsig4 <- sum(grepl("\\<4\\>", FsvCsMdownsig$chr))
# 0

# And for all chromosomes collectedly 
FsvCsMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FsvCsMdownsig$chr))
# 49


# Last create a new data frame with the observed and expected numbers 
chrloc.FsvCsMdown <- rbind( c(FsvCsMdownsigX, FsvCsMdownsig2, FsvCsMdownsig3, FsvCsMdownsig4),
                            c(chrXM/chrM*FsvCsMdownsigChr, chr2M/chrM*FsvCsMdownsigChr, chr3M/chrM*FsvCsMdownsigChr, chr4M/chrM*FsvCsMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FsvCsMdown) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FsvCsMdown["Obs",], p = chrloc.FsvCsMdown["Exp",], rescale.p = T)
# P = 0.8315


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FsvCsMdownF <- rbind( c(FsvCsMdownsigX, FsvCsMdownsig2, FsvCsMdownsig3, FsvCsMdownsig4),
                             c(chrXM-FsvCsMdownsigX, chr2M-FsvCsMdownsig2, chr3M-FsvCsMdownsig3, chr4M-FsvCsMdownsig4))
# Add names to the columns and rows
colnames(chrloc.FsvCsMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FsvCsMdownF)
# P = 0.8201

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the FLXautovCFMauto group
FavCaMupsig <- subset(FLXswapM.data, FLXa.CFMaup == "Sig")
# 65

# Next collect the observed number for each chromosome specifically
FavCaMupsigX <- sum(grepl("X", FavCaMupsig$chr))
# 5
FavCaMupsig2 <- sum(grepl("2L|2R", FavCaMupsig$chr))
# 28
FavCaMupsig3 <- sum(grepl("3L|3R", FavCaMupsig$chr))
# 31
FavCaMupsig4 <- sum(grepl("\\<4\\>", FavCaMupsig$chr))
# 1

# And for all chromosomes collectedly 
FavCaMupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FavCaMupsig$chr))
# 65


# Last create a new data frame with the observed and expected numbers 
chrloc.FavCaMup <- rbind( c(FavCaMupsigX, FavCaMupsig2, FavCaMupsig3, FavCaMupsig4),
                          c(chrXM/chrM*FavCaMupsigChr, chr2M/chrM*FavCaMupsigChr, chr3M/chrM*FavCaMupsigChr, chr4M/chrM*FavCaMupsigChr))
# Add names to the columns and rows
colnames(chrloc.FavCaMup) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaMup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FavCaMup["Obs",], p = chrloc.FavCaMup["Exp",], rescale.p = T)
# P = 0.3172

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FavCaMupF <- rbind( c(FavCaMupsigX, FavCaMupsig2, FavCaMupsig3, FavCaMupsig4),
                           c(chrXM-FavCaMupsigX, chr2M-FavCaMupsig2, chr3M-FavCaMupsig3, chr4M-FavCaMupsig4))
# Add names to the columns and rows
colnames(chrloc.FavCaMupF) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaMupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FavCaMupF)
# P = 0.2107

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the FLXautovCFMauto group
FavCaMdownsig <- subset(FLXswapM.data, FLXa.CFMadown == "Sig")
# 121

# Next collect the observed number for each chromosome specifically
FavCaMdownsigX <- sum(grepl("X", FavCaMdownsig$chr))
# 16
FavCaMdownsig2 <- sum(grepl("2L|2R", FavCaMdownsig$chr))
# 47
FavCaMdownsig3 <- sum(grepl("3L|3R", FavCaMdownsig$chr))
# 58
FavCaMdownsig4 <- sum(grepl("\\<4\\>", FavCaMdownsig$chr))
# 0

# And for all chromosomes collectedly 
FavCaMdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FavCaMdownsig$chr))
# 121


# Last create a new data frame with the observed and expected numbers 
chrloc.FavCaMdown <- rbind( c(FavCaMdownsigX, FavCaMdownsig2, FavCaMdownsig3, FavCaMdownsig4),
                            c(chrXM/chrM*FavCaMdownsigChr, chr2M/chrM*FavCaMdownsigChr, chr3M/chrM*FavCaMdownsigChr, chr4M/chrM*FavCaMdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FavCaMdown) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaMdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FavCaMdown["Obs",], p = chrloc.FavCaMdown["Exp",], rescale.p = T)
# P = 0.6736

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FavCaMdownF <- rbind( c(FavCaMdownsigX, FavCaMdownsig2, FavCaMdownsig3, FavCaMdownsig4),
                             c(chrXM-FavCaMdownsigX, chr2M-FavCaMdownsig2, chr3M-FavCaMdownsig3, chr4M-FavCaMdownsig4))
# Add names to the columns and rows
colnames(chrloc.FavCaMdownF) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaMdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FavCaMdownF)
# P = 0.8411

# We get the same results, so that is good

