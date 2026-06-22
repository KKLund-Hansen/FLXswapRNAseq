################################################################################################
################################ FLX CHROMOSOME SWAP CORRELATION ###############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/Both/Correlation")

# Read in files with data
corr.data <- read.table(file = "correlation.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


# First look at the data, is it normally distributed?

# Female
shapiro.test(corr.data$female)
# P = 1.242e-10
qqnorm(corr.data$female)
qqline(corr.data$female)
# No

# Male
shapiro.test(corr.data$male)
# P = 7.389e-08
qqnorm(corr.data$male)
qqline(corr.data$male)
# No


# Because the data is not normally distributed, I'll use spearman for the method
# Do a correlation test between the number of transcripts between groups in females and males
cor.test(corr.data$female, corr.data$male, method = "spearman")
# P = 1.255e-06


#########################################  PLOT DATA  #########################################


## The plot

pdf("Correlation.pdf", width = 7, height = 7)

par(mar = c(5, 5, 2, 2))
#
plot(corr.data$rank_female, corr.data$rank_male, xlim = c(0, 40), xlab = "Female", ylim = c(0, 40), 
     ylab = "Male", cex.axis = 1.2, cex.lab = 1.5, las = 1, pch = 16, col = "black", cex = 1.5)
abline(lm(corr.data$rank_male ~ corr.data$rank_female), col = "black", lwd = 1.5)

dev.off()
