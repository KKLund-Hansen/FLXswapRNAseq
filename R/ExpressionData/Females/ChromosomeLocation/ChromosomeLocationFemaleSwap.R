################################################################################################
############################# CHROMOSOME LOCTAION FLX SWAP FEMALES #############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiment/FLX/5.FLXchromosomeSwap/R/Females/ChromosomeLocation")

# Read in files with data
RNAswapfemale.data <- read.table(file = "FLXswapFres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


## Under/over representation on any chromosome

# To be able to calculate the expected number I first have collect the number of overall transcripts located on each chromosome

# Each chromosome
chrXF <- sum(grepl("X", RNAswapfemale.data$chr))
# 2553
chr2F <- sum(grepl("2L|2R", RNAswapfemale.data$chr))
# 6608
chr3F <- sum(grepl("3L|3R", RNAswapfemale.data$chr))
# 7350
chr4F <-sum(grepl("\\<4\\>", RNAswapfemale.data$chr))
# 111

# And then collect
chrF <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", RNAswapfemale.data$chr))
# 16622



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in in group A
GAFupsig <- subset(RNAswapfemale.data, gA41up == "Sig")
# 129

# Next collect the observed number for each chromosome specifically
GAFupsigX <- sum(grepl("X", GAFupsig$chr))
# 35
GAFupsig2 <- sum(grepl("2L|2R", GAFupsig$chr))
# 43
GAFupsig3 <- sum(grepl("3L|3R", GAFupsig$chr))
# 49
GAFupsig4 <- sum(grepl("\\<4\\>", GAFupsig$chr))
# 2

# And for all chromosomes collectedly 
GAFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GAFupsig$chr))
# 129


# Last create a new data frame with the observed and expected numbers 
chrloc.GAFup <- rbind( c(GAFupsigX, GAFupsig2, GAFupsig3, GAFupsig4),
                       c(chrXF/chrF*GAFupsigChr, chr2F/chrF*GAFupsigChr, chr3F/chrF*GAFupsigChr, chr4F/chrF*GAFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GAFup) <- c("X", "2", "3", "4")
rownames(chrloc.GAFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GAFup["Obs",], p = chrloc.GAFup["Exp",], rescale.p = T)
# P = 0.001359


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GAFupF <- rbind( c(GAFupsigX, GAFupsig2, GAFupsig3, GAFupsig4),
                        c(chrXF-GAFupsigX, chr2F-GAFupsig2, chr3F-GAFupsig3, chr4F-GAFupsig4))
# Add names to the columns and rows
colnames(chrloc.GAFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GAFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GAFupF)
# P = 0.00207

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group A
GAFdownsig <- subset(RNAswapfemale.data, gA41down == "Sig")
# 2021

#Next collect the observed number for each chromosome specifically
GAFdownsigX <- sum(grepl("X", GAFdownsig$chr))
# 272
GAFdownsig2 <- sum(grepl("2L|2R", GAFdownsig$chr))
# 813
GAFdownsig3 <- sum(grepl("3L|3R", GAFdownsig$chr))
# 914
GAFdownsig4 <- sum(grepl("\\<4\\>", GAFdownsig$chr))
# 22

# And for all chromosomes collectedly 
GAFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GAFdownsig$chr))
# 2021


# Last create a new data frame with the observed and expected numbers 
chrloc.GAFdown <- rbind( c(GAFdownsigX, GAFdownsig2, GAFdownsig3, GAFdownsig4),
                             c(chrXF/chrF*GAFdownsigChr, chr2F/chrF*GAFdownsigChr, chr3F/chrF*GAFdownsigChr, chr4F/chrF*GAFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GAFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GAFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GAFdown["Obs",], p = chrloc.GAFdown["Exp",], rescale.p = T)
# P = 0.01354




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group B 
GBFupsig <- subset(RNAswapfemale.data, gB7up == "Sig")
# 570

# Next collect the observed number for each chromosome specifically
GBFupsigX <- sum(grepl("X", GBFupsig$chr))
# 76
GBFupsig2 <- sum(grepl("2L|2R", GBFupsig$chr))
# 241
GBFupsig3 <- sum(grepl("3L|3R", GBFupsig$chr))
# 247
GBFupsig4 <- sum(grepl("\\<4\\>", GBFupsig$chr))
# 6

# And for all chromosomes collectedly 
GBFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBFupsig$chr))
# 570


# Last create a new data frame with the observed and expected numbers 
chrloc.GBFup <- rbind( c(GBFupsigX, GBFupsig2, GBFupsig3, GBFupsig4),
                       c(chrXF/chrF*GBFupsigChr, chr2F/chrF*GBFupsigChr, chr3F/chrF*GBFupsigChr, chr4F/chrF*GBFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GBFup) <- c("X", "2", "3", "4")
rownames(chrloc.GBFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBFup["Obs",], p = chrloc.GBFup["Exp",], rescale.p = T)
# P = 0.2835


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBFupF <- rbind( c(GBFupsigX, GBFupsig2, GBFupsig3, GBFupsig4),
                       c(chrXF-GBFupsigX, chr2F-GBFupsig2, chr3F-GBFupsig3, chr4F-GBFupsig4))
# Add names to the columns and rows
colnames(chrloc.GBFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GBFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBFupF)
# P = 0.2386

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group B
GBFdownsig <- subset(RNAswapfemale.data, gB7down == "Sig")
# 382

# Next collect the observed number for each chromosome specifically
GBFdownsigX <- sum(grepl("X", GBFdownsig$chr))
# 57
GBFdownsig2 <- sum(grepl("2L|2R", GBFdownsig$chr))
# 154
GBFdownsig3 <- sum(grepl("3L|3R", GBFdownsig$chr))
# 169
GBFdownsig4 <- sum(grepl("\\<4\\>", GBFdownsig$chr))
# 2

# And for all chromosomes collectedly 
GBFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GBFdownsig$chr))
# 382


# Last create a new data frame with the observed and expected numbers 
chrloc.GBFdown <- rbind( c(GBFdownsigX, GBFdownsig2, GBFdownsig3, GBFdownsig4),
                         c(chrXF/chrF*GBFdownsigChr, chr2F/chrF*GBFdownsigChr, chr3F/chrF*GBFdownsigChr, chr4F/chrF*GBFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GBFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GBFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GBFdown["Obs",], p = chrloc.GBFdown["Exp",], rescale.p = T)
# P = 0.9781


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GBFdownF <- rbind( c(GBFdownsigX, GBFdownsig2, GBFdownsig3, GBFdownsig4),
                          c(chrXF-GBFdownsigX, chr2F-GBFdownsig2, chr3F-GBFdownsig3, chr4F-GBFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GBFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GBFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GBFdownF)
# P = 0.9915

# We get the same results, so that is good




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group C
GCFupsig <- subset(RNAswapfemale.data, gC6up == "Sig")
# 177

# Next collect the observed number for each chromosome specifically
GCFupsigX <- sum(grepl("X", GCFupsig$chr))
# 35
GCFupsig2 <- sum(grepl("2L|2R", GCFupsig$chr))
# 70
GCFupsig3 <- sum(grepl("3L|3R", GCFupsig$chr))
# 71
GCFupsig4 <- sum(grepl("\\<4\\>", GCFupsig$chr))
# 0

# And for all chromosomes collectedly 
GCFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GCFupsig$chr))
# 176


# Last create a new data frame with the observed and expected numbers 
chrloc.GCFup <- rbind( c(GCFupsigX, GCFupsig2, GCFupsig3, GCFupsig4),
                       c(chrXF/chrF*GCFupsigChr, chr2F/chrF*GCFupsigChr, chr3F/chrF*GCFupsigChr, chr4F/chrF*GCFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GCFup) <- c("X", "2", "3", "4")
rownames(chrloc.GCFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GCFup["Obs",], p = chrloc.GCFup["Exp",], rescale.p = T)
# P = 0.2486


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCFupF <- rbind( c(GCFupsigX, GCFupsig2, GCFupsig3, GCFupsig4),
                        c(chrXF-GCFupsigX, chr2F-GCFupsig2, chr3F-GCFupsig3, chr4F-GCFupsig4))
# Add names to the columns and rows
colnames(chrloc.GCFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GCFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCFupF)
# P = 0.3094

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group C
GCFdownsig <- subset(RNAswapfemale.data, gC6down == "Sig")
# 372

# Next collect the observed number for each chromosome specifically
GCFdownsigX <- sum(grepl("X", GCFdownsig$chr))
# 39
GCFdownsig2 <- sum(grepl("2L|2R", GCFdownsig$chr))
# 151
GCFdownsig3 <- sum(grepl("3L|3R", GCFdownsig$chr))
# 176
GCFdownsig4 <- sum(grepl("\\<4\\>", GCFdownsig$chr))
# 5

# And for all chromosomes collectedly 
GCFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GCFdownsig$chr))
# 371


# Last create a new data frame with the observed and expected numbers 
chrloc.GCFdown <- rbind( c(GCFdownsigX, GCFdownsig2, GCFdownsig3, GCFdownsig4),
                         c(chrXF/chrF*GCFdownsigChr, chr2F/chrF*GCFdownsigChr, chr3F/chrF*GCFdownsigChr, chr4F/chrF*GCFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GCFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GCFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GCFdown["Obs",], p = chrloc.GCFdown["Exp",], rescale.p = T)
# P = 0.02678


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GCFdownF <- rbind( c(GCFdownsigX, GCFdownsig2, GCFdownsig3, GCFdownsig4),
                          c(chrXF-GCFdownsigX, chr2F-GCFdownsig2, chr3F-GCFdownsig3, chr4F-GCFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GCFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GCFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GCFdownF)
# P = 0.01669

# We get the same results, so that is good




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group D
GDFupsig <- subset(RNAswapfemale.data, gD22up == "Sig")
# 191

# Next collect the observed number for each chromosome specifically
GDFupsigX <- sum(grepl("X", GDFupsig$chr))
# 10
GDFupsig2 <- sum(grepl("2L|2R", GDFupsig$chr))
# 76
GDFupsig3 <- sum(grepl("3L|3R", GDFupsig$chr))
# 105
GDFupsig4 <- sum(grepl("\\<4\\>", GDFupsig$chr))
# 0

# And for all chromosomes collectedly 
GDFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GDFupsig$chr))
# 191


# Last create a new data frame with the observed and expected numbers 
chrloc.GDFup <- rbind( c(GDFupsigX, GDFupsig2, GDFupsig3, GDFupsig4),
                       c(chrXF/chrF*GDFupsigChr, chr2F/chrF*GDFupsigChr, chr3F/chrF*GDFupsigChr, chr4F/chrF*GDFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GDFup) <- c("X", "2", "3", "4")
rownames(chrloc.GDFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GDFup["Obs",], p = chrloc.GDFup["Exp",], rescale.p = T)
# P = 0.0002712


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDFupF <- rbind( c(GDFupsigX, GDFupsig2, GDFupsig3, GDFupsig4),
                        c(chrXF-GDFupsigX, chr2F-GDFupsig2, chr3F-GDFupsig3, chr4F-GDFupsig4))
# Add names to the columns and rows
colnames(chrloc.GDFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GDFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GDFupF)
# P = 6.555e-05

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group D
GDFdownsig <- subset(RNAswapfemale.data, gD22down == "Sig")
# 284

# Next collect the observed number for each chromosome specifically
GDFdownsigX <- sum(grepl("X", GDFdownsig$chr))
# 37
GDFdownsig2 <- sum(grepl("2L|2R", GDFdownsig$chr))
# 113
GDFdownsig3 <- sum(grepl("3L|3R", GDFdownsig$chr))
# 133
GDFdownsig4 <- sum(grepl("\\<4\\>", GDFdownsig$chr))
# 1

# And for all chromosomes collectedly 
GDFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GDFdownsig$chr))
# 284


# Last create a new data frame with the observed and expected numbers 
chrloc.GDFdown <- rbind( c(GDFdownsigX, GDFdownsig2, GDFdownsig3, GDFdownsig4),
                         c(chrXF/chrF*GDFdownsigChr, chr2F/chrF*GDFdownsigChr, chr3F/chrF*GDFdownsigChr, chr4F/chrF*GDFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GDFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GDFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GDFdown["Obs",], p = chrloc.GDFdown["Exp",], rescale.p = T)
# P = 0.6005


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GDFdownF <- rbind( c(GDFdownsigX, GDFdownsig2, GDFdownsig3, GDFdownsig4),
                          c(chrXF-GDFdownsigX, chr2F-GDFdownsig2, chr3F-GDFdownsig3, chr4F-GDFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GDFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GDFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GDFdownF)
# P = 0.6872

# We get the same results, so that is good




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group E
GEFupsig <- subset(RNAswapfemale.data, gE4up == "Sig")
# 98

# Next collect the observed number for each chromosome specifically
GEFupsigX <- sum(grepl("X", GEFupsig$chr))
# 18
GEFupsig2 <- sum(grepl("2L|2R", GEFupsig$chr))
# 31
GEFupsig3 <- sum(grepl("3L|3R", GEFupsig$chr))
# 48
GEFupsig4 <- sum(grepl("\\<4\\>", GEFupsig$chr))
# 1

# And for all chromosomes collectedly 
GEFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GEFupsig$chr))
# 98


# Last create a new data frame with the observed and expected numbers 
chrloc.GEFup <- rbind( c(GEFupsigX, GEFupsig2, GEFupsig3, GEFupsig4),
                       c(chrXF/chrF*GEFupsigChr, chr2F/chrF*GEFupsigChr, chr3F/chrF*GEFupsigChr, chr4F/chrF*GEFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GEFup) <- c("X", "2", "3", "4")
rownames(chrloc.GEFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GEFup["Obs",], p = chrloc.GEFup["Exp",], rescale.p = T)
# P = 0.4092


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEFupF <- rbind( c(GEFupsigX, GEFupsig2, GEFupsig3, GEFupsig4),
                        c(chrXF-GEFupsigX, chr2F-GEFupsig2, chr3F-GEFupsig3, chr4F-GEFupsig4))
# Add names to the columns and rows
colnames(chrloc.GEFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GEFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEFupF)
# P = 0.2725

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group E
GEFdownsig <- subset(RNAswapfemale.data, gE4down == "Sig")
# 264

# Next collect the observed number for each chromosome specifically
GEFdownsigX <- sum(grepl("X", GEFdownsig$chr))
# 40
GEFdownsig2 <- sum(grepl("2L|2R", GEFdownsig$chr))
# 116
GEFdownsig3 <- sum(grepl("3L|3R", GEFdownsig$chr))
# 105
GEFdownsig4 <- sum(grepl("\\<4\\>", GEFdownsig$chr))
# 2

# And for all chromosomes collectedly 
GEFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GEFdownsig$chr))
# 263


# Last create a new data frame with the observed and expected numbers 
chrloc.GEFdown <- rbind( c(GEFdownsigX, GEFdownsig2, GEFdownsig3, GEFdownsig4),
                         c(chrXF/chrF*GEFdownsigChr, chr2F/chrF*GEFdownsigChr, chr3F/chrF*GEFdownsigChr, chr4F/chrF*GEFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GEFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GEFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GEFdown["Obs",], p = chrloc.GEFdown["Exp",], rescale.p = T)
# P = 0.496


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GEFdownF <- rbind( c(GEFdownsigX, GEFdownsig2, GEFdownsig3, GEFdownsig4),
                          c(chrXF-GEFdownsigX, chr2F-GEFdownsig2, chr3F-GEFdownsig3, chr4F-GEFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GEFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GEFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GEFdownF)
# P = 0.4286

# We get the same results, so that is good




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group F 
GFFupsig <- subset(RNAswapfemale.data, gF18up == "Sig")
# 40

# Next collect the observed number for each chromosome specifically
GFFupsigX <- sum(grepl("X", GFFupsig$chr))
# 8
GFFupsig2 <- sum(grepl("2L|2R", GFFupsig$chr))
# 15
GFFupsig3 <- sum(grepl("3L|3R", GFFupsig$chr))
# 17
GFFupsig4 <- sum(grepl("\\<4\\>", GFFupsig$chr))
# 0

# And for all chromosomes collectedly 
GFFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GFFupsig$chr))
# 40


# Last create a new data frame with the observed and expected numbers 
chrloc.GFFup <- rbind( c(GFFupsigX, GFFupsig2, GFFupsig3, GFFupsig4),
                       c(chrXF/chrF*GFFupsigChr, chr2F/chrF*GFFupsigChr, chr3F/chrF*GFFupsigChr, chr4F/chrF*GFFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GFFup) <- c("X", "2", "3", "4")
rownames(chrloc.GFFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GFFup["Obs",], p = chrloc.GFFup["Exp",], rescale.p = T)
# P = 0.824


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GFFupF <- rbind( c(GFFupsigX, GFFupsig2, GFFupsig3, GFFupsig4),
                        c(chrXF-GFFupsigX, chr2F-GFFupsig2, chr3F-GFFupsig3, chr4F-GFFupsig4))
# Add names to the columns and rows
colnames(chrloc.GFFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GFFupF)
# P = 0.7834

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group F
GFFdownsig <- subset(RNAswapfemale.data, gF18down == "Sig")
# 264

# Next collect the observed number for each chromosome specifically
GFFdownsigX <- sum(grepl("X", GFFdownsig$chr))
# 41
GFFdownsig2 <- sum(grepl("2L|2R", GFFdownsig$chr))
# 103
GFFdownsig3 <- sum(grepl("3L|3R", GFFdownsig$chr))
# 117
GFFdownsig4 <- sum(grepl("\\<4\\>", GFFdownsig$chr))
# 3

# And for all chromosomes collectedly 
GFFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GFFdownsig$chr))
# 263


# Last create a new data frame with the observed and expected numbers 
chrloc.GFFdown <- rbind( c(GFFdownsigX, GFFdownsig2, GFFdownsig3, GFFdownsig4),
                         c(chrXF/chrF*GFFdownsigChr, chr2F/chrF*GFFdownsigChr, chr3F/chrF*GFFdownsigChr, chr4F/chrF*GFFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GFFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GFFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GFFdown["Obs",], p = chrloc.GFFdown["Exp",], rescale.p = T)
# P = 0.8118


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GFFdownF <- rbind( c(GFFdownsigX, GFFdownsig2, GFFdownsig3, GFFdownsig4),
                          c(chrXF-GFFdownsigX, chr2F-GFFdownsig2, chr3F-GFFdownsig3, chr4F-GFFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GFFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GFFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GFFdownF)
# P = 0.7041

# We get the same results, so that is good




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group G
GGFupsig <- subset(RNAswapfemale.data, gG19up == "Sig")
# 92

# Next collect the observed number for each chromosome specifically
GGFupsigX <- sum(grepl("X", GGFupsig$chr))
# 16
GGFupsig2 <- sum(grepl("2L|2R", GGFupsig$chr))
# 32
GGFupsig3 <- sum(grepl("3L|3R", GGFupsig$chr))
# 42
GGFupsig4 <- sum(grepl("\\<4\\>", GGFupsig$chr))
# 2

# And for all chromosomes collectedly 
GGFupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GGFupsig$chr))
# 92


# Last create a new data frame with the observed and expected numbers 
chrloc.GGFup <- rbind( c(GGFupsigX, GGFupsig2, GGFupsig3, GGFupsig4),
                       c(chrXF/chrF*GGFupsigChr, chr2F/chrF*GGFupsigChr, chr3F/chrF*GGFupsigChr, chr4F/chrF*GGFupsigChr))
# Add names to the columns and rows
colnames(chrloc.GGFup) <- c("X", "2", "3", "4")
rownames(chrloc.GGFup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GGFup["Obs",], p = chrloc.GGFup["Exp",], rescale.p = T)
# P = 0.2628


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GGFupF <- rbind( c(GGFupsigX, GGFupsig2, GGFupsig3, GGFupsig4),
                        c(chrXF-GGFupsigX, chr2F-GGFupsig2, chr3F-GGFupsig3, chr4F-GGFupsig4))
# Add names to the columns and rows
colnames(chrloc.GGFupF) <- c("X", "2", "3", "4")
rownames(chrloc.GGFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GGFupF)
# P = 0.2144

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group G
GGFdownsig <- subset(RNAswapfemale.data, gG19down == "Sig")
# 199

# Next collect the observed number for each chromosome specifically
GGFdownsigX <- sum(grepl("X", GGFdownsig$chr))
# 32
GGFdownsig2 <- sum(grepl("2L|2R", GGFdownsig$chr))
# 67
GGFdownsig3 <- sum(grepl("3L|3R", GGFdownsig$chr))
# 100
GGFdownsig4 <- sum(grepl("\\<4\\>", GGFdownsig$chr))
# 0

# And for all chromosomes collectedly 
GGFdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", GGFdownsig$chr))
# 199


# Last create a new data frame with the observed and expected numbers 
chrloc.GGFdown <- rbind( c(GGFdownsigX, GGFdownsig2, GGFdownsig3, GGFdownsig4),
                         c(chrXF/chrF*GGFdownsigChr, chr2F/chrF*GGFdownsigChr, chr3F/chrF*GGFdownsigChr, chr4F/chrF*GGFdownsigChr))
# Add names to the columns and rows
colnames(chrloc.GGFdown) <- c("X", "2", "3", "4")
rownames(chrloc.GGFdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.GGFdown["Obs",], p = chrloc.GGFdown["Exp",], rescale.p = T)
# P = 0.1802


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.GGFdownF <- rbind( c(GGFdownsigX, GGFdownsig2, GGFdownsig3, GGFdownsig4),
                          c(chrXF-GGFdownsigX, chr2F-GGFdownsig2, chr3F-GGFdownsig3, chr4F-GGFdownsig4))
# Add names to the columns and rows
colnames(chrloc.GGFdownF) <- c("X", "2", "3", "4")
rownames(chrloc.GGFdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.GGFdownF)
# P = 0.2145

# We get the same results, so that is good




######################## GROUP 23 ########################


##### UP-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for up-regulated transcripts in group 23
G23Fupsig <- subset(RNAswapfemale.data, g23up == "Sig")
# 16

# Next collect the observed number for each chromosome specifically
G23FupsigX <- sum(grepl("X", G23Fupsig$chr))
# 1
G23Fupsig2 <- sum(grepl("2L|2R", G23Fupsig$chr))
# 8
G23Fupsig3 <- sum(grepl("3L|3R", G23Fupsig$chr))
# 7
G23Fupsig4 <- sum(grepl("\\<4\\>", G23Fupsig$chr))
# 0

# And for all chromosomes collectedly 
G23FupsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G23Fupsig$chr))
# 16


# Last create a new data frame with the observed and expected numbers 
chrloc.G23Fup <- rbind( c(G23FupsigX, G23Fupsig2, G23Fupsig3, G23Fupsig4),
                        c(chrXF/chrF*G23FupsigChr, chr2F/chrF*G23FupsigChr, chr3F/chrF*G23FupsigChr, chr4F/chrF*G23FupsigChr))
# Add names to the columns and rows
colnames(chrloc.G23Fup) <- c("X", "2", "3", "4")
rownames(chrloc.G23Fup) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G23Fup["Obs",], p = chrloc.G23Fup["Exp",], rescale.p = T)
# P = 0.7068


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G23FupF <- rbind( c(G23FupsigX, G23Fupsig2, G23Fupsig3, G23Fupsig4),
                         c(chrXF-G23FupsigX, chr2F-G23Fupsig2, chr3F-G23Fupsig3, chr4F-G23Fupsig4))
# Add names to the columns and rows
colnames(chrloc.G23FupF) <- c("X", "2", "3", "4")
rownames(chrloc.G23FupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G23FupF)
# P = 0.6639

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Under/over representation on any chromosome

# First subset the data to be able to identify the observed chromosome location for down-regulated transcripts in group 23
G23Fdownsig <- subset(RNAswapfemale.data, g23down == "Sig")
# 6

# Next collect the observed number for each chromosome specifically
G23FdownsigX <- sum(grepl("X", G23Fdownsig$chr))
# 2
G23Fdownsig2 <- sum(grepl("2L|2R", G23Fdownsig$chr))
# 2
G23Fdownsig3 <- sum(grepl("3L|3R", G23Fdownsig$chr))
# 2
G23Fdownsig4 <- sum(grepl("\\<4\\>", G23Fdownsig$chr))
# 0

# And for all chromosomes collectedly 
G23FdownsigChr <- sum(grepl("X|2L|2R|3L|3R|\\<4\\>", G23Fdownsig$chr))
# 6


# Last create a new data frame with the observed and expected numbers 
chrloc.G23Fdown <- rbind( c(G23FdownsigX, G23Fdownsig2, G23Fdownsig3, G23Fdownsig4),
                          c(chrXF/chrF*G23FdownsigChr, chr2F/chrF*G23FdownsigChr, chr3F/chrF*G23FdownsigChr, chr4F/chrF*G23FdownsigChr))
# Add names to the columns and rows
colnames(chrloc.G23Fdown) <- c("X", "2", "3", "4")
rownames(chrloc.G23Fdown) <- c("Obs", "Exp")

# It is now possible to use the new data frame to do a chi-squared test
chisq.test(chrloc.G23Fdown["Obs",], p = chrloc.G23Fdown["Exp",], rescale.p = T)
# P = 0.6765


# But I get a warning message that I might not be able to use Chi-squared because the numbers are too small.
# So I have to use a Fisher's exact test instead and to do so I have to reformat the data
chrloc.G23FdownF <- rbind( c(G23FdownsigX, G23Fdownsig2, G23Fdownsig3, G23Fdownsig4),
                           c(chrXF-G23FdownsigX, chr2F-G23Fdownsig2, chr3F-G23Fdownsig3, chr4F-G23Fdownsig4))
# Add names to the columns and rows
colnames(chrloc.G23FdownF) <- c("X", "2", "3", "4")
rownames(chrloc.G23FdownF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(chrloc.G23FdownF)
# P = 0.4182

# We get the same results, so that is good

