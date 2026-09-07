################################################################################################
#################### MITO-SENSITIVE TRANSCRIPTS FLX SWAP FEMALES - 7 GROUPS ####################
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
MS.data <- read.table(file = "mitosensitive.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(MS.data) <- MS.data[,1]

# Next add a new column with the flybase names
MS.data$gene.id <- mapIds(drosophila2.db, keys = row.names(MS.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can subset the mito-sensitive data set to match the female data set
MSfswap7g<- subset(MS.data, MS.data$gene.id %in% FLXswap7gF.data$gene.id)
# 1210 gene IDs overlap

# Next I'll calculate the proportion of the number of mito-sensitive transcripts out of all my transcripts
MSpropF <- length(MSfswap7g$gene.id)/length(FLXswap7gF.data$gene.id)
# 0.08

# And for the chisq.test I also have to calculate
NotMSpropF <- 1-MSpropF
# 0.92




######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many mito-sensitive transcripts are there?
GAFupsigMS <- sum(GAFupsig$gene.id %in% MSfswap7g$gene.id)
# 6
# How many are not?
GAFupsigNotMS <- length(GAFupsig$gene.id)-GAFupsigMS
# 116

# The number of expected mito-sensitive transcripts
length(GAFupsig$gene.id)*MSpropF
# 9.383914


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFupsigMS, GAFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2347



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many mito-sensitive transcripts are there?
GAFdownsigMS <- sum(GAFdownsig$gene.id %in% MSfswap7g$gene.id)
# 213
# How many are not?
GAFdownsigNotMS <- length(GAFdownsig$gene.id)-GAFdownsigMS
# 1758

# The number of expected mito-sensitive transcripts
length(GAFdownsig$gene.id)*MSpropF
# 153.8155


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigMS, GAFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 2.046e-09




######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many mito-sensitive transcripts are there?
GBFupsigMS <- sum(GBFupsig$gene.id %in% MSfswap7g$gene.id)
# 53
# How many are not?
GBFupsigNotMS <- length(GBFupsig$gene.id)-GBFupsigMS
# 510

# The number of expected mito-sensitive transcripts
length(GBFupsig$gene.id)*MSpropF
# 43.93615


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFupsigMS, GBFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.1544



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many mito-sensitive transcripts are there?
GBFdownsigMS <- sum(GBFdownsig$gene.id %in% MSfswap7g$gene.id)
# 18
# How many are not?
GBFdownsigNotMS <- length(GBFdownsig$gene.id)-GBFdownsigMS
# 364

# The number of expected mito-sensitive transcripts
length(GBFdownsig$gene.id)*MSpropF
# 29.81103


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigMS, GBFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.02427




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many mito-sensitive transcripts are there?
GCFupsigMS <- sum(GCFupsig$gene.id %in% MSfswap7g$gene.id)
# 15
# How many are not?
GCFupsigNotMS <- length(GCFupsig$gene.id)-GCFupsigMS
# 147

# The number of expected mito-sensitive transcripts
length(GCFupsig$gene.id)*MSpropF
# 12.64237


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFupsigMS, GCFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.4898



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many mito-sensitive transcripts are there?
GCFdownsigMS <- sum(GCFdownsig$gene.id %in% MSfswap7g$gene.id)
# 15
# How many are not?
GCFdownsigNotMS <- length(GCFdownsig$gene.id)-GCFdownsigMS
# 350

# The number of expected mito-sensitive transcripts
length(GCFdownsig$gene.id)*MSpropF
# 28.48436


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigMS, GCFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.008506




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many mito-sensitive transcripts are there?
GDFupsigMS <- sum(GDFupsig$gene.id %in% MSfswap7g$gene.id)
# 15
# How many are not?
GDFupsigNotMS <- length(GDFupsig$gene.id)-GDFupsigMS
# 161

# The number of expected mito-sensitive transcripts
length(GDFupsig$gene.id)*MSpropF
# 13.73492


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFupsigMS, GDFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.7222



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many mito-sensitive transcripts are there?
GDFdownsigMS <- sum(GDFdownsig$gene.id %in% MSfswap7g$gene.id)
# 10
# How many are not?
GDFdownsigNotMS <- length(GDFdownsig$gene.id)-GDFdownsigMS
# 273

# The number of expected mito-sensitive transcripts
length(GDFdownsig$gene.id)*MSpropF
# 22.08513


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigMS, GDFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.007402




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many mito-sensitive transcripts are there?
GEFupsigMS <- sum(GEFupsig$gene.id %in% MSfswap7g$gene.id)
# 4
# How many are not?
GEFupsigNotMS <- length(GEFupsig$gene.id)-GEFupsigMS
# 87

# The number of expected mito-sensitive transcripts
length(GEFupsig$gene.id)*MSpropF
# 7.10158


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFupsigMS, GEFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2255



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many mito-sensitive transcripts are there?
GEFdownsigMS <- sum(GEFdownsig$gene.id %in% MSfswap7g$gene.id)
# 16
# How many are not?
GEFdownsigNotMS <- length(GEFdownsig$gene.id)-GEFdownsigMS
# 247

# The number of expected mito-sensitive transcripts
length(GEFdownsig$gene.id)*MSpropF
# 20.52435


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigMS, GEFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.2983




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many mito-sensitive transcripts are there?
GFFupsigMS <- sum(GFFupsig$gene.id %in% MSfswap7g$gene.id)
# 3
# How many are not?
GFFupsigNotMS <- length(GFFupsig$gene.id)-GFFupsigMS
# 34

# The number of expected mito-sensitive transcripts
length(GFFupsig$gene.id)*MSpropF
# 2.887456


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFupsigMS, GFFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.945

# However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
MS.GFFupF <- rbind(c(GFFupsigMS, GFFupsigNotMS), 
                   c(length(MSfswap7g$gene.id)-GFFupsigMS, (length(FLXswap7gF.data$gene.id)-length(MSfswap7g$gene.id))-GFFupsigNotMS))
#Add names to the columns and rows
colnames(MS.GFFupF) <- c("MS", "NotMS")
rownames(MS.GFFupF) <- c("Sig", "NotSig")

# It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(MS.GFFupF)
# P = 0.7638

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many mito-sensitive transcripts are there?
GFFdownsigMS <- sum(GFFdownsig$gene.id %in% MSfswap7g$gene.id)
# 20
# How many are not?
GFFdownsigNotMS <- length(GFFdownsig$gene.id)-GFFdownsigMS
# 243

# The number of expected mito-sensitive transcripts
length(GFFdownsig$gene.id)*MSpropF
# 20.52435


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigMS, GFFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.9041




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many mito-sensitive transcripts are there?
GGFupsigMS <- sum(GGFupsig$gene.id %in% MSfswap7g$gene.id)
# 4
# How many are not?
GGFupsigNotMS <- length(GGFupsig$gene.id)-GGFupsigMS
# 91

# The number of expected mito-sensitive transcripts
length(GGFupsig$gene.id)*MSpropF
# 7.413738


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFupsigMS, GGFupsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.1916



##### DOWN-REGULATED TRANSCRIPTS #####

# First subset the data to be able to identify mito-sensitive transcripts among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many mito-sensitive transcripts are there?
GGFdownsigMS <- sum(GGFdownsig$gene.id %in% MSfswap7g$gene.id)
# 19
# How many are not?
GGFdownsigNotMS <- length(GGFdownsig$gene.id)-GGFdownsigMS
# 175

# The number of expected mito-sensitive transcripts
length(GGFdownsig$gene.id)*MSpropF
# 15.13963


# Do the test, with the number of mito-sensitive transcripts and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigMS, GGFdownsigNotMS), p = c(MSpropF, NotMSpropF))
# P = 0.3015

