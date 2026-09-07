################################################################################################
############## CHANGE IN INTENCITY OF SEXUAL CONFLICT FLX SWAP FEMALES - 7 GROUPS ##############
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Females/SexualConflict")

#First install the software from Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("drosophila2.db")
BiocManager::install("org.Dm.eg.db")

# Set up environment
library(drosophila2.db)
library(org.Dm.eg.db)

# Read in files with data
FLXswap7gF.data <- read.table(file = "FLXswapFres7G.csv", h = T, sep = ",", stringsAsFactors = T)
InnocentiSC.data <- read.table(file = "InnocentietalSexualConflict.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### PREPARE THE DATA #####

# First add the probe names as row names
rownames(InnocentiSC.data) <- InnocentiSC.data[,1]

# Next add a new column with the flybase names
InnocentiSC.data$gene.id <- mapIds(drosophila2.db, keys = row.names(InnocentiSC.data), keytype = "PROBEID", column = "FLYBASE")

# Now I can finally subset the Innocenti et al. data set to match the female data set
InnocentiSCF <- subset(InnocentiSC.data, InnocentiSC.data$gene.id %in% FLXswap7gF.data$gene.id)
# 1106 transcripts


# Next I'll calculate the proportion of the number of SC transcripts out of all my transcripts
ISCpropF <- length(InnocentiSCF$gene.id)/length(FLXswap7gF.data$gene.id)
# 0.07
# And for the chisq.test I also have to calculate
NotISCpropF <- 1-ISCpropF
# 0.93



######################## GROUP A ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group A
GAFupsig <- subset(FLXswap7gF.data, gAup == "Sig")
# 122

# How many transcripts are related to sexual conflict?
GAFupsigISC <- sum(GAFupsig$gene.id %in% InnocentiSCF$gene.id)
# 19
# How many are not?
GAFupsigNotISC <- length(GAFupsig$gene.id)-GAFupsigISC
# 103

# The number of expected transcripts which are related to sexual conflict
length(GAFupsig$gene.id)*ISCpropF
# 8.702483

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFupsigISC, GAFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.000292



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group A
GAFdownsig <- subset(FLXswap7gF.data, gAdown == "Sig")
# 1971

# How many transcripts are related to sexual conflict?
GAFdownsigISC <- sum(GAFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 177
# How many are not?
GAFdownsigNotISC <- length(GAFdownsig$gene.id)-GAFdownsigISC
# 1794

# The number of expected transcripts which are related to sexual conflict
length(GAFdownsig$gene.id)*ISCpropF
# 140.595

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GAFdownsigISC, GAFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.001443


######################## GROUP B ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group B
GBFupsig <- subset(FLXswap7gF.data, gBup == "Sig")
# 563

# How many transcripts are related to sexual conflict?
GBFupsigISC <- sum(GBFupsig$gene.id %in% InnocentiSCF$gene.id)
# 52
# How many are not?
GBFupsigNotISC <- length(GBFupsig$gene.id)-GBFupsigISC
# 511

# The number of expected transcripts which are related to sexual conflict
length(GBFupsig$gene.id)*ISCpropF
# 40.15982

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFupsigISC, GBFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.05253



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group B
GBFdownsig <- subset(FLXswap7gF.data, gBdown == "Sig")
# 382

# How many transcripts are related to sexual conflict?
GBFdownsigISC <- sum(GBFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 49
# How many are not?
GBFdownsigNotISC <- length(GBFdownsig$gene.id)-GBFdownsigISC
# 333

# The number of expected transcripts which are related to sexual conflict
length(GBFdownsig$gene.id)*ISCpropF
# 27.24876

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GBFdownsigISC, GBFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 1.171e-06




######################## GROUP C ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group C
GCFupsig <- subset(FLXswap7gF.data, gCup == "Sig")
# 162

# How many transcripts are related to sexual conflict?
GCFupsigISC <- sum(GCFupsig$gene.id %in% InnocentiSCF$gene.id)
# 9
# How many are not?
GCFupsigNotISC <- length(GCFupsig$gene.id)-GCFupsigISC
# 153

# The number of expected transcripts which are related to sexual conflict
length(GCFupsig$gene.id)*ISCpropF
# 11.55576

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFupsigISC, GCFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.4353



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group C
GCFdownsig <- subset(FLXswap7gF.data, gCdown == "Sig")
# 365

# How many transcripts are related to sexual conflict?
GCFdownsigISC <- sum(GCFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 41
# How many are not?
GCFdownsigNotISC <- length(GCFdownsig$gene.id)-GCFdownsigISC
# 324

# The number of expected transcripts which are related to sexual conflict
length(GCFdownsig$gene.id)*ISCpropF
# 26.03612

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GCFdownsigISC, GCFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.002341




######################## GROUP D ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group D
GDFupsig <- subset(FLXswap7gF.data, gDup == "Sig")
# 176

# How many transcripts are related to sexual conflict?
GDFupsigISC <- sum(GDFupsig$gene.id %in% InnocentiSCF$gene.id)
# 10
# How many are not?
GDFupsigNotISC <- length(GDFupsig$gene.id)-GDFupsigISC
# 166

# The number of expected transcripts which are related to sexual conflict
length(GDFupsig$gene.id)*ISCpropF
# 12.5544

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFupsigISC, GDFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.4544



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group D
GDFdownsig <- subset(FLXswap7gF.data, gDdown == "Sig")
# 283

# How many transcripts are related to sexual conflict?
GDFdownsigISC <- sum(GDFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 29
# How many are not?
GDFdownsigNotISC <- length(GDFdownsig$gene.id)-GDFdownsigISC
# 254

# The number of expected transcripts which are related to sexual conflict
length(GDFdownsig$gene.id)*ISCpropF
# 20.18691

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GDFdownsigISC, GDFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.0418




######################## GROUP E ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group E
GEFupsig <- subset(FLXswap7gF.data, gEup == "Sig")
# 91

# How many transcripts are related to sexual conflict?
GEFupsigISC <- sum(GEFupsig$gene.id %in% InnocentiSCF$gene.id)
# 13
# How many are not?
GEFupsigNotISC <- length(GEFupsig$gene.id)-GEFupsigISC
# 78

# The number of expected transcripts which are related to sexual conflict
length(GEFupsig$gene.id)*ISCpropF
# 6.491196

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFupsigISC, GEFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.008026



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group E
GEFdownsig <- subset(FLXswap7gF.data, gEdown == "Sig")
# 263

# How many transcripts are related to sexual conflict?
GEFdownsigISC <- sum(GEFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 18
# How many are not?
GEFdownsigNotISC <- length(GEFdownsig$gene.id)-GEFdownsigISC
# 245

# The number of expected transcripts which are related to sexual conflict
length(GEFdownsig$gene.id)*ISCpropF
# 18.76027

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GEFdownsigISC, GEFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.8555




######################## GROUP F ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group F
GFFupsig <- subset(FLXswap7gF.data, gFup == "Sig")
# 37

# How many transcripts are related to sexual conflict?
GFFupsigISC <- sum(GFFupsig$gene.id %in% InnocentiSCF$gene.id)
# 5
# How many are not?
GFFupsigNotISC <- length(GFFupsig$gene.id)-GFFupsigISC
# 32

# The number of expected transcripts which are related to sexual conflict
length(GFFupsig$gene.id)*ISCpropF
# 2.639278

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFupsigISC, GFFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.1316

#However I get a warning message, so I'll use Fisher's exact test instead because of the low numbers
ISC.GFFupF <- rbind(c(GFFupsigISC, GFFupsigNotISC),
                    c(length(InnocentiSCF$gene.id)-GFFupsigISC, (length(FLXswap7gF.data$gene.id)-length(InnocentiSCF$gene.id))-GFFupsigNotISC))
#Add names to the columns and rows
colnames(ISC.GFFupF) <- c("SA", "NotSA")
rownames(ISC.GFFupF) <- c("Sig", "NotSig")

#It is now possible to use the new data frame to do a Fisher's exact test
fisher.test(ISC.GFFupF)
# P = 0.185

# We get the same results, so that is good



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group F
GFFdownsig <- subset(FLXswap7gF.data, gFdown == "Sig")
# 263

# How many transcripts are related to sexual conflict?
GFFdownsigISC <- sum(GFFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 31
# How many are not?
GFFdownsigNotISC <- length(GFFdownsig$gene.id)-GFFdownsigISC
# 232

# The number of expected transcripts which are related to sexual conflict
length(GFFdownsig$gene.id)*ISCpropF
# 18.76027

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GFFdownsigISC, GFFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.003364




######################## GROUP G ########################


##### UP-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify genes related to sexual conflict among the up-regulated transcripts in group G
GGFupsig <- subset(FLXswap7gF.data, gGup == "Sig")
# 95

# How many transcripts are related to sexual conflict?
GGFupsigISC <- sum(GGFupsig$gene.id %in% InnocentiSCF$gene.id)
# 8
# How many are not?
GGFupsigNotISC <- length(GGFupsig$gene.id)-GGFupsigISC
# 87

# The number of expected transcripts which are related to sexual conflict
length(GGFupsig$gene.id)*ISCpropF
# 6.776524

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFupsigISC, GGFupsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.6258



##### DOWN-REGULATED TRANSCRIPTS #####

## Transcripts related to sexual conflict

# First subset the data to be able to identify transcripts related to sexual conflict among the down-regulated transcripts in group G
GGFdownsig <- subset(FLXswap7gF.data, gGdown == "Sig")
# 194

# How many transcripts are related to sexual conflict?
GGFdownsigISC <- sum(GGFdownsig$gene.id %in% InnocentiSCF$gene.id)
# 20
# How many are not?
GGFdownsigNotISC <- length(GGFdownsig$gene.id)-GGFdownsigISC
# 174

# The number of expected transcripts which are related to sexual conflict
length(GGFdownsig$gene.id)*ISCpropF
# 13.83837

# Do the test, with the number of transcripts related to sexual conflict and not and the proportion of overall transcripts
chisq.test(c(GGFdownsigISC, GGFdownsigNotISC), p = c(ISCpropF, NotISCpropF))
# P = 0.08565

