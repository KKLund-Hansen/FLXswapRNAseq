################################################################################################
############################# CHROMOSOME LOCATION FLX SWAP FEMALES #############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/ChromosomeLocation")

# Read in files with data
FLXswapF.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


## Under/over representation on any chromosome

# To be able to calculate the expected number I first have collect the number of overall transcripts located on each chromosome

# Each chromosome
chrXF <- sum(grepl("X", FLXswapF.data$chr))
# 2372
chr2F <- sum(grepl("2L|2R", FLXswapF.data$chr))
# 6091
chr3F <- sum(grepl("3L|3R", FLXswapF.data$chr))
# 6865
chr4F <-sum(grepl("\\<4\\>", FLXswapF.data$chr))
# 107

# And then collect
chrF <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FLXswapF.data$chr))
# 15435



######################## FLX vs CFM ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the FLXvCFM group
FvCfupsig <- subset(FLXswapF.data, FLX.CFMup == "Sig")
# 35

# Next collect the observed number for each chromosome specifically
FvCfupsigX <- sum(grepl("X", FvCfupsig$chr))
# 9
FvCfupsig2 <- sum(grepl("2L|2R", FvCfupsig$chr))
# 13
FvCfupsig3 <- sum(grepl("3L|3R", FvCfupsig$chr))
# 13
FvCfupsig4 <- sum(grepl("\\<4\\>", FvCfupsig$chr))
# 0

# And for all chromosomes collectedly 
FvCfupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FvCfupsig$chr))
# 35


# Last create a new data frame with the observed and expected numbers 
chrloc.FvCfup <- rbind( c(FvCfupsigX, FvCfupsig2, FvCfupsig3, FvCfupsig4),
                        c(chrXF/chrF*FvCfupsigChr, chr2F/chrF*FvCfupsigChr, chr3F/chrF*FvCfupsigChr, chr4F/chrF*FvCfupsigChr))
# Add names to the columns and rows
colnames(chrloc.FvCfup) <- c("X", "2", "3", "4")
rownames(chrloc.FvCfup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FvCfup["Obs",], p = chrloc.FvCfup["Exp",], rescale.p = T)
# P = 0.3688

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FvCfupF <- rbind( c(FvCfupsigX, FvCfupsig2, FvCfupsig3, FvCfupsig4),
                         c(chrXF-FvCfupsigX, chr2F-FvCfupsig2, chr3F-FvCfupsig3, chr4F-FvCfupsig4))
# Add names to the columns and rows
colnames(chrloc.FvCfupF) <- c("X", "2", "3", "4")
rownames(chrloc.FvCfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FvCfupF)
# P = 0.4051

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the FLXvCFM group
FvCfdownsig <- subset(FLXswapF.data, FLX.CFMdown == "Sig")
# 40

#Next collect the observed number for each chromosome specifically
FvCfdownsigX <- sum(grepl("X", FvCfdownsig$chr))
# 9
FvCfdownsig2 <- sum(grepl("2L|2R", FvCfdownsig$chr))
# 20
FvCfdownsig3 <- sum(grepl("3L|3R", FvCfdownsig$chr))
# 9
FvCfdownsig4 <- sum(grepl("\\<4\\>", FvCfdownsig$chr))
# 0

# And for all chromosomes collectedly 
FvCfdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FvCfdownsig$chr))
# 40


# Last create a new data frame with the observed and expected numbers 
chrloc.FvCfdown <- rbind( c(FvCfdownsigX, FvCfdownsig2, FvCfdownsig3, FvCfdownsig4),
                          c(chrXF/chrF*FvCfdownsigChr, chr2F/chrF*FvCfdownsigChr, chr3F/chrF*FvCfdownsigChr, chr4F/chrF*FvCfdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FvCfdown) <- c("X", "2", "3", "4")
rownames(chrloc.FvCfdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FvCfdown["Obs",], p = chrloc.FvCfdown["Exp",], rescale.p = T)
# P = 0.1499

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FvCfdownF <- rbind( c(FvCfdownsigX, FvCfdownsig2, FvCfdownsig3, FvCfdownsig4),
                          c(chrXF-FvCfdownsigX, chr2F-FvCfdownsig2, chr3F-FvCfdownsig3, chr4F-FvCfdownsig4))
# Add names to the columns and rows
colnames(chrloc.FvCfdownF) <- c("X", "2", "3", "4")
rownames(chrloc.FvCfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FvCfdownF)
# P = 0.1327

# We get the same results, so that is good




######################## AUTOSOMES VS SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcriptsin the autovsex group
AvSfupsig <- subset(FLXswapF.data, auto.sexup == "Sig")
# 244

# Next collect the observed number for each chromosome specifically
AvSfupsigX <- sum(grepl("X", AvSfupsig$chr))
# 29
AvSfupsig2 <- sum(grepl("2L|2R", AvSfupsig$chr))
# 110
AvSfupsig3 <- sum(grepl("3L|3R", AvSfupsig$chr))
# 29
AvSfupsig4 <- sum(grepl("\\<4\\>", AvSfupsig$chr))
# 1

# And for all chromosomes collectedly 
AvSfupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", AvSfupsig$chr))
# 244


# Last create a new data frame with the observed and expected numbers 
chrloc.AvSfup <- rbind( c(AvSfupsigX, AvSfupsig2, AvSfupsig3, AvSfupsig4),
                        c(chrXF/chrF*AvSfupsigChr, chr2F/chrF*AvSfupsigChr, chr3F/chrF*AvSfupsigChr, chr4F/chrF*AvSfupsigChr))
# Add names to the columns and rows
colnames(chrloc.AvSfup) <- c("X", "2", "3", "4")
rownames(chrloc.AvSfup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.AvSfup["Obs",], p = chrloc.AvSfup["Exp",], rescale.p = T)
# P = 0.2261

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.AvSfupF <- rbind(c(AvSfupsigX, AvSfupsig2, AvSfupsig3, AvSfupsig4),
                        c(chrXF-AvSfupsigX, chr2F-AvSfupsig2, chr3F-AvSfupsig3, chr4F-AvSfupsig4))
# Add names to the columns and rows
colnames(chrloc.AvSfupF) <- c("X", "2", "3", "4")
rownames(chrloc.AvSfupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.AvSfupF)
# P = 0.2311

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcriptsin the autovsex group
AvSfdownsig <- subset(FLXswapF.data, auto.sexdown == "Sig")
# 161

# Next collect the observed number for each chromosome specifically
AvSfdownsigX <- sum(grepl("X", AvSfdownsig$chr))
# 14
AvSfdownsig2 <- sum(grepl("2L|2R", AvSfdownsig$chr))
# 76
AvSfdownsig3 <- sum(grepl("3L|3R", AvSfdownsig$chr))
# 69
AvSfdownsig4 <- sum(grepl("\\<4\\>", AvSfdownsig$chr))
# 2

# And for all chromosomes collectedly 
AvSfdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", AvSfdownsig$chr))
# 161


# Last create a new data frame with the observed and expected numbers 
chrloc.AvSfdown <- rbind(c(AvSfdownsigX, AvSfdownsig2, AvSfdownsig3, AvSfdownsig4),
                        c(chrXF/chrF*AvSfdownsigChr, chr2F/chrF*AvSfdownsigChr, chr3F/chrF*AvSfdownsigChr, chr4F/chrF*AvSfdownsigChr))
# Add names to the columns and rows
colnames(chrloc.AvSfdown) <- c("X", "2", "3", "4")
rownames(chrloc.AvSfdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.AvSfdown["Obs",], p = chrloc.AvSfdown["Exp",], rescale.p = T)
# P = 0.04803

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.AvSfdownF <- rbind( c(AvSfdownsigX, AvSfdownsig2, AvSfdownsig3, AvSfdownsig4),
                           c(chrXF-AvSfdownsigX, chr2F-AvSfdownsig2, chr3F-AvSfdownsig3, chr4F-AvSfdownsig4))
# Add names to the columns and rows
colnames(chrloc.AvSfdownF) <- c("X", "2", "3", "4")
rownames(chrloc.AvSfdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.AvSfdownF)
# P = 0.02784

# We get the same results, so that is good




######################## FLX SEX CHROMSOMES VS CFM SEX CHROMOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the FLXsexvCFMsex group
FsvCsFupsig <- subset(FLXswapF.data, FLXs.CFMsup == "Sig")
# 137

# Next collect the observed number for each chromosome specifically
FsvCsFupsigX <- sum(grepl("X", FsvCsFupsig$chr))
# 17
FsvCsFupsig2 <- sum(grepl("2L|2R", FsvCsFupsig$chr))
# 44
FsvCsFupsig3 <- sum(grepl("3L|3R", FsvCsFupsig$chr))
# 76
FsvCsFupsig4 <- sum(grepl("\\<4\\>", FsvCsFupsig$chr))
# 0

# And for all chromosomes collectedly 
FsvCsFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FsvCsFupsig$chr))
# 137


# Last create a new data frame with the observed and expected numbers 
chrloc.FsvCsFup <- rbind( c(FsvCsFupsigX, FsvCsFupsig2, FsvCsFupsig3, FsvCsFupsig4),
                          c(chrXF/chrF*FsvCsFupsigChr, chr2F/chrF*FsvCsFupsigChr, chr3F/chrF*FsvCsFupsigChr, chr4F/chrF*FsvCsFupsigChr))
# Add names to the columns and rows
colnames(chrloc.FsvCsFup) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FsvCsFup["Obs",], p = chrloc.FsvCsFup["Exp",], rescale.p = T)
# P = 0.06212

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FsvCsFupF <- rbind( c(FsvCsFupsigX, FsvCsFupsig2, FsvCsFupsig3, FsvCsFupsig4),
                           c(chrXF-FsvCsFupsigX, chr2F-FsvCsFupsig2, chr3F-FsvCsFupsig3, chr4F-FsvCsFupsig4))
# Add names to the columns and rows
colnames(chrloc.FsvCsFupF) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FsvCsFupF)
# P = 0.08575

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the FLXsexvCFMsex group
FsvCsFdownsig <- subset(FLXswapF.data, FLXs.CFMsdown == "Sig")
# 55

# Next collect the observed number for each chromosome specifically
FsvCsFdownsigX <- sum(grepl("X", FsvCsFdownsig$chr))
# 8
FsvCsFdownsig2 <- sum(grepl("2L|2R", FsvCsFdownsig$chr))
# 22
FsvCsFdownsig3 <- sum(grepl("3L|3R", FsvCsFdownsig$chr))
# 24
FsvCsFdownsig4 <- sum(grepl("\\<4\\>", FsvCsFdownsig$chr))
# 1

# And for all chromosomes collectedly 
FsvCsFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FsvCsFdownsig$chr))
# 55


# Last create a new data frame with the observed and expected numbers 
chrloc.FsvCsFdown <- rbind( c(FsvCsFdownsigX, FsvCsFdownsig2, FsvCsFdownsig3, FsvCsFdownsig4),
                            c(chrXF/chrF*FsvCsFdownsigChr, chr2F/chrF*FsvCsFdownsigChr, chr3F/chrF*FsvCsFdownsigChr, chr4F/chrF*FsvCsFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FsvCsFdown) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FsvCsFdown["Obs",], p = chrloc.FsvCsFdown["Exp",], rescale.p = T)
# P = 0.7913

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FsvCsFdownF <- rbind( c(FsvCsFdownsigX, FsvCsFdownsig2, FsvCsFdownsig3, FsvCsFdownsig4),
                             c(chrXF-FsvCsFdownsigX, chr2F-FsvCsFdownsig2, chr3F-FsvCsFdownsig3, chr4F-FsvCsFdownsig4))
# Add names to the columns and rows
colnames(chrloc.FsvCsFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.FsvCsFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FsvCsFdownF)
# P = 0.5801

# We get the same results, so that is good




######################## FLX AUTOSOMES VS CFM AUTOSOMES ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in the FLXautovCFMauto group
FavCaFupsig <- subset(FLXswapF.data, FLXa.CFMaup == "Sig")
# 525

# Next collect the observed number for each chromosome specifically
FavCaFupsigX <- sum(grepl("X", FavCaFupsig$chr))
# 85
FavCaFupsig2 <- sum(grepl("2L|2R", FavCaFupsig$chr))
# 200
FavCaFupsig3 <- sum(grepl("3L|3R", FavCaFupsig$chr))
# 237
FavCaFupsig4 <- sum(grepl("\\<4\\>", FavCaFupsig$chr))
# 2

# And for all chromosomes collectedly 
FavCaFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FavCaFupsig$chr))
# 524


# Last create a new data frame with the observed and expected numbers 
chrloc.FavCaFup <- rbind( c(FavCaFupsigX, FavCaFupsig2, FavCaFupsig3, FavCaFupsig4),
                          c(chrXF/chrF*FavCaFupsigChr, chr2F/chrF*FavCaFupsigChr, chr3F/chrF*FavCaFupsigChr, chr4F/chrF*FavCaFupsigChr))
# Add names to the columns and rows
colnames(chrloc.FavCaFup) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FavCaFup["Obs",], p = chrloc.FavCaFup["Exp",], rescale.p = T)
# P = 0.736

# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.FavCaFupF <- rbind( c(FavCaFupsigX, FavCaFupsig2, FavCaFupsig3, FavCaFupsig4),
                           c(chrXF-FavCaFupsigX, chr2F-FavCaFupsig2, chr3F-FavCaFupsig3, chr4F-FavCaFupsig4))
# Add names to the columns and rows
colnames(chrloc.FavCaFupF) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.FavCaFupF)
# P = 0.793

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in the FLXautovCFMauto group
FavCaFdownsig <- subset(FLXswapF.data, FLXa.CFMadown == "Sig")
# 2994

# Next collect the observed number for each chromosome specifically
FavCaFdownsigX <- sum(grepl("X", FavCaFdownsig$chr))
# 399
FavCaFdownsig2 <- sum(grepl("2L|2R", FavCaFdownsig$chr))
# 1196
FavCaFdownsig3 <- sum(grepl("3L|3R", FavCaFdownsig$chr))
# 1365
FavCaFdownsig4 <- sum(grepl("\\<4\\>", FavCaFdownsig$chr))
# 32

# And for all chromosomes collectedly 
FavCaFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", FavCaFdownsig$chr))
# 2992


# Last create a new data frame with the observed and expected numbers 
chrloc.FavCaFdown <- rbind( c(FavCaFdownsigX, FavCaFdownsig2, FavCaFdownsig3, FavCaFdownsig4),
                            c(chrXF/chrF*FavCaFdownsigChr, chr2F/chrF*FavCaFdownsigChr, chr3F/chrF*FavCaFdownsigChr, chr4F/chrF*FavCaFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.FavCaFdown) <- c("X", "2", "3", "4")
rownames(chrloc.FavCaFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.FavCaFdown["Obs",], p = chrloc.FavCaFdown["Exp",], rescale.p = T)
# P = 0.00163


