################################################################################################
######################## DIFFERENCES IN ERROR VARIANCE BETWEEN THE SEXES #######################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/ExpressionData/Both/StandardError")

# Read in files with data
FLXswapSE.data <- read.table(file = "FLXswapSEres.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################

# First average across the sexes
FLXswapSE.data$femaleM <- rowMeans(FLXswapSE.data[,c(2:7)])
FLXswapSE.data$maleM <- rowMeans(FLXswapSE.data[,c(8:13)])


# Next test for normality
# Female
qqnorm(FLXswapSE.data$femaleM)
qqline(FLXswapSE.data$femaleM)
# No

# Male
qqnorm(FLXswapSE.data$maleM)
qqline(FLXswapSE.data$maleM)
# No


# Because the data is not normally distributed we do a Wilcoxon Rank Sum test 
# First test if the standard error are significant different 
wilcox.test(FLXswapSE.data$femaleM, FLXswapSE.data$maleM)
# Yes, P < 2.2e-16

# Is there more variance in males?
wilcox.test(FLXswapSE.data$femaleM, FLXswapSE.data$maleM, alternative = "less")
# No, P = 1

# Or females?
wilcox.test(FLXswapSE.data$femaleM, FLXswapSE.data$maleM, alternative = "greater")
# Yes, P < 2.2e-16


