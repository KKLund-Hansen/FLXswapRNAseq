################################################################################################
##################### MITO-PROTEOME TRANSCRIPTS FLX SWAP FEMALES - 7 GROUPS ####################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/MitoNuclearConflict")

# First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
MP.data <- read.table(file = "mitoproteome.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# There are some gene IDs that are multiple in the MP.data so first I'll remove the doubles
str(MP.data)
MP.data <- unique(MP.data)

# Next add the CG names as row names
rownames(MP.data) <- MP.data[,1]

# Next add a new column with the flybase names
MP.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MP.data), keytype = "FLYBASECG", column = "FLYBASE")

# Now I can finally subset the mito-proteome data set to match the female data set
MPfswap7g <- subset(MP.data, MP.data$gene.id %in% FLXswap7gF.data$gene.id)
# 762 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-proteome transcriptsout of all my transcripts
MPpropF <- length(MPfswap7g$gene.id)/length(FLXswap7gF.data$gene.id)
# 0.05

# And for the chisq.test I also have to calculate
NotMPpropF <- 1-MPpropF
# 0.95




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many mito-proteome transcripts are there?
GAFupsigMP <- sum(GAFupsig$gene.id %in% MPfswap7g$gene.id)
# 7
# How many are not?
GAFupsigNotMP <- length(GAFupsig$gene.id)-GAFupsigMP
# 115

# The number of expected mito-proteome transcripts
length(GAFupsig$gene.id)*MPpropF
# 5.995743


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMP, GAFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.674



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many mito-proteome transcripts are there?
GAFdownsigMP <- sum(GAFdownsig$gene.id %in% MPfswap7g$gene.id)
# 90
# How many are not?
GAFdownsigNotMP <- length(GAFdownsig$gene.id)-GAFdownsigMP
# 1881

# The number of expected mito-proteome transcripts
length(GAFdownsig$gene.id)*MPpropF
# 96.86566


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMP, GAFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.4744




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many mito-proteome transcripts are there?
GBFupsigMP <- sum(GBFupsig$gene.id %in% MPfswap7g$gene.id)
# 57
# How many are not?
GBFupsigNotMP <- length(GBFupsig$gene.id)-GBFupsigMP
# 506

# The number of expected mito-proteome transcripts
length(GBFupsig$gene.id)*MPpropF
# 27.66888


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMP, GBFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 1.075e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many mito-proteome transcripts are there?
GBFdownsigMP <- sum(GBFdownsig$gene.id %in% MPfswap7g$gene.id)
# 34
# How many are not?
GBFdownsigNotMP <- length(GBFdownsig$gene.id)-GBFdownsigMP
# 348

# The number of expected mito-proteome transcripts
length(GBFdownsig$gene.id)*MPpropF
# 18.77356


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMP, GBFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.0003135




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many mito-proteome transcripts are there?
GCFupsigMP <- sum(GCFupsig$gene.id %in% MPfswap7g$gene.id)
# 23
# How many are not?
GCFupsigNotMP <- length(GCFupsig$gene.id)-GCFupsigMP
# 139

# The number of expected mito-proteome transcripts
length(GCFupsig$gene.id)*MPpropF
# 7.961561


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMP, GCFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 4.611e-08



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many mito-proteome transcripts are there?
GCFdownsigMP <- sum(GCFdownsig$gene.id %in% MPfswap7g$gene.id)
# 14
# How many are not?
GCFdownsigNotMP <- length(GCFdownsig$gene.id)-GCFdownsigMP
# 351

# The number of expected mito-proteome transcripts
length(GCFdownsig$gene.id)*MPpropF
# 17.93808


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMP, GCFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.3403




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many mito-proteome transcripts are there?
GDFupsigMP <- sum(GDFupsig$gene.id %in% MPfswap7g$gene.id)
# 33
# How many are not?
GDFupsigNotMP <- length(GDFupsig$gene.id)-GDFupsigMP
# 143

# The number of expected mito-proteome transcripts
length(GDFupsig$gene.id)*MPpropF
# 8.649597


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMP, GDFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P < 2.2e-16



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many mito-proteome transcripts are there?
GDFdownsigMP <- sum(GDFdownsig$gene.id %in% MPfswap7g$gene.id)
# 10
# How many are not?
GDFdownsigNotMP <- length(GDFdownsig$gene.id)-GDFdownsigMP
# 273

# The number of expected mito-proteome transcripts
length(GDFdownsig$gene.id)*MPpropF
# 13.90816


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMP, GDFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.2825




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many mito-proteome transcripts are there?
GEFupsigMP <- sum(GEFupsig$gene.id %in% MPfswap7g$gene.id)
# 5
# How many are not?
GEFupsigNotMP <- length(GEFupsig$gene.id)-GEFupsigMP
# 86

# The number of expected mito-proteome transcripts
length(GEFupsig$gene.id)*MPpropF
# 4.472235


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMP, GEFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.798

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GEFupF <- rbind(c(GEFupsigMP, GEFupsigNotMP), 
                   c(length(MPfswap7g$gene.id)-GEFupsigMP, (length(FLXswap7gF.data$gene.id)-length(MPfswap7g$gene.id))-GEFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GEFupF) <- c("MP", "NotMP")
rownames(MP.GEFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GEFupF)
# P = 0.8051

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many mito-proteome transcripts are there?
GEFdownsigMP <- sum(GEFdownsig$gene.id %in% MPfswap7g$gene.id)
# 20
# How many are not?
GEFdownsigNotMP <- length(GEFdownsig$gene.id)-GEFdownsigMP
# 243

# The number of expected mito-proteome transcripts
length(GEFdownsig$gene.id)*MPpropF
# 12.92525


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMP, GEFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.04358




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many mito-proteome transcripts are there?
GFFupsigMP <- sum(GFFupsig$gene.id %in% MPfswap7g$gene.id)
# 4
# How many are not?
GFFupsigNotMP <- length(GFFupsig$gene.id)-GFFupsigMP
# 33

# The number of expected mito-proteome transcripts
length(GFFupsig$gene.id)*MPpropF
# 1.818381


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMP, GFFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.09709

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GFFupF <- rbind(c(GFFupsigMP, GFFupsigNotMP), 
                   c(length(MPfswap7g$gene.id)-GFFupsigMP, (length(FLXswap7gF.data$gene.id)-length(MPfswap7g$gene.id))-GFFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GFFupF) <- c("MP", "NotMP")
rownames(MP.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GFFupF)
# P = 0.1065

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many mito-proteome transcripts are there?
GFFdownsigMP <- sum(GFFdownsig$gene.id %in% MPfswap7g$gene.id)
# 8
# How many are not?
GFFdownsigNotMP <- length(GFFdownsig$gene.id)-GFFdownsigMP
# 255

# The number of expected mito-proteome transcripts
length(GFFdownsig$gene.id)*MPpropF
# 12.92525


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMP, GFFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.5431




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many mito-proteome transcripts are there?
GGFupsigMP <- sum(GGFupsig$gene.id %in% MPfswap7g$gene.id)
# 7
# How many are not?
GGFupsigNotMP <- length(GGFupsig$gene.id)-GGFupsigMP
# 88

# The number of expected mito-proteome transcripts
length(GGFupsig$gene.id)*MPpropF
# 4.668817


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMP, GGFupsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.2685

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MP.GGFupF <- rbind(c(GGFupsigMP, GGFupsigNotMP), 
                   c(length(MPfswap7g$gene.id)-GGFupsigMP, (length(FLXswap7gF.data$gene.id)-length(MPfswap7g$gene.id))-GGFupsigNotMP))
#Add names to the columns and rows
colnames(MP.GGFupF) <- c("MP", "NotMP")
rownames(MP.GGFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MP.GGFupF)
# P = 0.2348

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-proteome transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many mito-proteome transcripts are there?
GGFdownsigMP <- sum(GGFdownsig$gene.id %in% MPfswap7g$gene.id)
# 21
# How many are not?
GGFdownsigNotMP <- length(GGFdownsig$gene.id)-GGFdownsigMP
# 173

# The number of expected mito-proteome transcripts
length(GGFdownsig$gene.id)*MPpropF
# 9.534215


# Do the test, with the number of mito-proteome transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMP, GGFdownsigNotMP), p = c(MPpropF, NotMPpropF))
# P = 0.0001401

