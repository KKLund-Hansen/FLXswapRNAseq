################################################################################################
########################### DIFFERENCES IN LOG2FOLD BETWEEN THE SEXES ##########################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Both/LogFold")

# Read in files with data
FLXswapL2F.data <- read.table(file = "FLXswapL2Fres.csv", h = T, sep = ",", stringsAsFactors = T)
FLXswapL2Fabs.data <- read.table(file = "FLXswapL2Fabsres.csv", h = T, sep = ",", stringsAsFactors = T)
Fsig.data <- read.table(file = "AllsigFemaleres.csv", h = T, sep = ",", stringsAsFactors = T)
Msig.data <- read.table(file = "AllsigMaleres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


##### MAX LOG2FOLD IN FEMALES #####

# First subset the data to match the female significant data set
SigFl2f <- subset(FLXswapL2F.data, FLXswapL2F.data$gene.id %in% Fsig.data$V1)

# Next subset the data to find the max in only female and male values 
SigFl2fmaxf <- SigFl2f[,c(2:7)]
SigFl2fmaxm <- SigFl2f[,c(8:13)]

# Find the max value for each gene among the female contrasts
maxcolumnF <- colnames(SigFl2fmaxf)[max.col(abs(SigFl2fmaxf))]

# Next make a new data set
maxvaluesF <- data.frame(colName = maxcolumnF)
maxvaluesF$gene.id <- SigFl2f$gene.id
maxvaluesF$valueF <- character(5573)
for(i in 1:5573) { maxvaluesF$valueF[i] <- SigFl2fmaxf[i, names(SigFl2fmaxf[i,]) %in% maxvaluesF$colName[i] ] }
maxvaluesF$valueF <- as.numeric(maxvaluesF$valueF)
# Change female to male in the names so they can match with the male data
maxvaluesF$colName <- gsub("f", "m", maxvaluesF$colName)
# Now add the male
maxvaluesF$valueM <- character(5573)
for(i in 1:5573) { maxvaluesF$valueM[i] <- as.numeric(SigFl2fmaxm[i, names(SigFl2fmaxm[i,]) %in% maxvaluesF$colName[i] ]) }
maxvaluesF$valueM <- as.numeric(maxvaluesF$valueM)

# Next test for normality
# Female
qqnorm(maxvaluesF$valueF)
qqline(maxvaluesF$valueF)
# No

# Male
qqnorm(maxvaluesF$valueM)
qqline(maxvaluesF$valueM)
# No


# Because the data is not normally distributed, I'll use spearman for the method
# Do a correlation test between the number of transcripts between groups in females and males
cor.test(maxvaluesF$valueF, maxvaluesF$valueM, method = "spearman")
# P < 2.2e-16

# Is the correlation negative or positive?
cor.test(maxvaluesF$valueF, maxvaluesF$valueM, method = "spearman", alternative = "less")
# Yes, negative, P < 2.2e-16

# We also do a Wilcoxon rank sum test
wilcox.test(abs(maxvaluesF$valueF), abs(maxvaluesF$valueM))
# P < 2.2e-16

median(abs(maxvaluesF$valueF))
median(abs(maxvaluesF$valueM))


##### MAX LOG2FOLD IN MALES #####

# First subset the data to match the female significant data set
SigMl2f <- subset(FLXswapL2F.data, FLXswapL2F.data$gene.id %in% Msig.data$V1)

# Next subset the data to find the max in only male and female values 
SigMl2fmaxm <- SigMl2f[,c(8:13)]
SigMl2fmaxf <- SigMl2f[,c(2:7)]

# Find the max value for each gene among the male contrasts
maxcolumnM <- colnames(SigMl2fmaxm )[max.col(abs(SigMl2fmaxm ))]

# Next make a new data set
maxvaluesM <- data.frame(colName = maxcolumnM)
maxvaluesM$gene.id <- SigMl2f$gene.id
maxvaluesM$valueM <- character(596)
for(i in 1:596) { maxvaluesM$valueM[i] <- SigMl2fmaxm[i, names(SigMl2fmaxm[i,]) %in% maxvaluesM$colName[i] ] }
maxvaluesM$valueM <- as.numeric(maxvaluesM$valueM)

# Change male to female in the names so they can match with the female data
maxvaluesM$colName <- gsub("m", "f", maxvaluesM$colName)

# Now add the male
maxvaluesM$valueF <- character(596)
for(i in 1:596) { maxvaluesM$valueF[i] <- as.numeric(SigMl2fmaxf[i, names(SigMl2fmaxf[i,]) %in% maxvaluesM$colName[i] ]) }
maxvaluesM$valueF <- as.numeric(maxvaluesM$valueF)


# Next test for normality
# Female
shapiro.test(maxvaluesM$valueF)
# P < 2.2e-16
qqnorm(maxvaluesM$valueF)
qqline(maxvaluesM$valueF)
# No

# Male
shapiro.test(maxvaluesM$valueM)
# P = 7.329e-11
qqnorm(maxvaluesM$valueM)
qqline(maxvaluesM$valueM)
# No

# Because the data is not normally distributed, I'll use spearman for the method
# Do a correlation test between the number of transcripts between groups in females and males
cor.test(maxvaluesM$valueF, maxvaluesM$valueM, method = "spearman")
# P < 2.2e-16

# Is the correlation negative or positive?
cor.test(maxvaluesM$valueF, maxvaluesM$valueM, method = "spearman", alternative = "greater")
# Yes, positive, P < 2.2e-16


# We also do a Wilcoxon rank sum test
wilcox.test(abs(maxvaluesM$valueF), abs(maxvaluesM$valueM))
# P < 2.2e-16

median(abs(maxvaluesM$valueF))
median(abs(maxvaluesM$valueM))



#########################################  PLOT DATA  #########################################


##### PLOT FEMALE DATA #####

# Add an extra column for plotting
maxvaluesF$fx <- as.factor(rep("A", 5573))
maxvaluesF$mx <- as.factor(rep("B", 5573))

# The plot
pdf("L2FchangeF.pdf", family = "serif", width = 7, height = 7)

par(mar = c(3, 5, 2, 2))
plot(NULL, xlim = c(0.5, 2.5), xlab = "", xaxt = "n", ylim = c(-60, 10), ylab = "Log2 fold change", cex.axis = 1.25, cex.lab = 1.4, las = 1,
     main = "Significant female transcripts", cex.main = 1.5, font.main = 1)
# Add points
abline(h=0)
segments(1, maxvaluesF$valueF, 2, maxvaluesF$valueM, cex = 1.3, col = "gray15")
stripchart(maxvaluesF$valueF ~ maxvaluesF$fx, add = T, vertical = T, pch = 16, cex = 1.5, col = "black", at = 1)
stripchart(maxvaluesF$valueM ~ maxvaluesF$mx, add = T, vertical = T, pch = 16, cex = 1.5, col = "black", at = 2)
axis(1, at = c(1:2), cex.axis = 1.4, labels = c("Female", "Male"))
mtext("A", side = 3, line = 0.3, at = 0.15, cex = 1.4)

dev.off()



##### PLOT MALE DATA #####

# Add an extra column for plotting
maxvaluesM$fx <- as.factor(rep("A", 596))
maxvaluesM$mx <- as.factor(rep("B", 596))

# The plot
pdf("L2FchangeM.pdf", family = "serif", width = 7, height = 7)

par(mar = c(3, 5, 2, 2))
plot(NULL, xlim = c(0.5, 2.5), xlab = "", xaxt = "n", ylim = c(-60, 10), ylab = "Log2 fold change", cex.axis = 1.25, cex.lab = 1.4, las = 1, 
     main = "Significant male transcripts", cex.main = 1.4, font.main = 1)
# Add points and lines
abline(h=0)
segments(1, maxvaluesM$valueF, 2, maxvaluesM$valueM, cex = 1.3, col = "gray15")
stripchart(maxvaluesM$valueF ~ maxvaluesM$fx, add = T, vertical = T, pch = 16, cex = 1.5, col = "black", at = 1)
stripchart(maxvaluesM$valueM ~ maxvaluesM$mx, add = T, vertical = T, pch = 16, cex = 1.5, col = "black", at = 2)
axis(1, at = c(1:2), cex.axis = 1.4, labels = c("Female", "Male"))
mtext("B", side = 3, line = 0.3, at = 0.15, cex = 1.4)

dev.off()


##### COMBINED PLOT #####

# The plot
pdf("L2FchangeBoth.pdf", family = "serif", width = 14, height = 7)

layout(matrix(c(1, 2), 2, 2, byrow = TRUE))
par(mar = c(3, 5, 2, 2))
# First plot, females
plot(NULL, xlim = c(0.5, 2.5), xlab = "", xaxt = "n", ylim = c(-60, 10), ylab = "Log2 fold change", cex.axis = 1.35, cex.lab = 1.6, las = 1,
     main = "Significant female transcripts", cex.main = 1.8, font.main = 1)
# Add points
abline(h=0)
segments(1, maxvaluesF$valueF, 2, maxvaluesF$valueM, cex = 1.3, col = "gray15")
stripchart(maxvaluesF$valueF ~ maxvaluesF$fx, add = T, vertical = T, pch = 16, cex = 1.7, col = "black", at = 1)
stripchart(maxvaluesF$valueM ~ maxvaluesF$mx, add = T, vertical = T, pch = 16, cex = 1.7, col = "black", at = 2)
axis(1, at = c(1:2), cex.axis = 1.6, labels = c("Female", "Male"))
mtext("A", side = 3, line = 0.3, at = 0.2, cex = 1.4)
# Second plot, males
plot(NULL, xlim = c(0.5, 2.5), xlab = "", xaxt = "n", ylim = c(-60, 10), ylab = "Log2 fold change", cex.axis = 1.35, cex.lab = 1.6, las = 1, 
     main = "Significant male transcripts", cex.main = 1.8, font.main = 1)
# Add points and lines
abline(h=0)
segments(1, maxvaluesM$valueF, 2, maxvaluesM$valueM, cex = 1.3, col = "gray15")
stripchart(maxvaluesM$valueF ~ maxvaluesM$fx, add = T, vertical = T, pch = 16, cex = 1.7, col = "black", at = 1)
stripchart(maxvaluesM$valueM ~ maxvaluesM$mx, add = T, vertical = T, pch = 16, cex = 1.7, col = "black", at = 2)
axis(1, at = c(1:2), cex.axis = 1.6, labels = c("Female", "Male"))
mtext("B", side = 3, line = 0.3, at = 0.2, cex = 1.4)

dev.off()


