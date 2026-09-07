################################################################################################
################################ FLX CHROMOSOME SWAP CORRELATION ###############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Both/Correlation")

# Read in files with data
corr.data <- read.table(file = "correlation.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################

# First look at the data, is it normally distributed?

# Female
shapiro.test(corr.data$female)
# P = 7.583e-11
qqnorm(corr.data$female)
qqline(corr.data$female)
# No

# Male
shapiro.test(corr.data$male)
# P = 6.439594e-08
qqnorm(corr.data$male)
qqline(corr.data$male)
# No


# Because the data is not normally distributed, I'll use spearman for the method
# Do a correlation test between the number of transcripts between groups in females and males
cor.test(corr.data$female, corr.data$male, method = "spearman")
# P = 8.434e-07


#########################################  PLOT DATA  #########################################


### PLOT - RANK ###

# The plot
pdf("Correlation.pdf", family = "serif", width = 7, height = 7)

par(mar = c(5, 5, 2, 2))
plot(corr.data$rank_female, corr.data$rank_male, xlim = c(0, 40), xlab = "Female", ylim = c(0, 40), 
     ylab = "Male", cex.axis = 1.2, cex.lab = 1.5, las = 1, pch = 16, col = "black", cex = 1.5)
abline(lm(corr.data$rank_male ~ corr.data$rank_female), col = "black", lwd = 1.5)

dev.off()



### PLOT - LOG ###

# First add log to the data set
corr.data$log_female <- log10(corr.data$female)
corr.data$log_male <- log10(corr.data$male)


# The plot
pdf("Correlation_log.pdf", family = "serif", width = 7, height = 7)

par(mar = c(5, 5, 2, 2))
plot(corr.data$log_female, corr.data$log_male, xlim = c(0, 4), xlab = "log(Female)", ylim = c(0, 2.5), 
     ylab = "log(Male)", cex.axis = 1.25, cex.lab = 1.4, las = 1, pch = 16, col = "black", cex = 1.5)
abline(lm(corr.data$log_male ~ corr.data$log_female), col = "black", lwd = 1.5)

dev.off()


