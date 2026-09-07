################################################################################################
################################### FLX SWAP LOCOMOTION ASSAY ##################################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/PhenotypicData/Locomotion")

# Set up environment
library(car)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(Hmisc)
library(lme4)

# Read in csv file with data
Locoswap.data <- read.csv(file = "locomotionswap.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


### GENERAL LINEAR MIXED MODELS ###

# General linear mixed-effect model testing if genotype and sex have a significant effect on locomotion
model.loco <- glmer(cbind(active, inactive) ~ type * sex + day + (1|type:rep_pop) + (1|time), data = Locoswap.data, family = binomial)
# ANOVA
Anova(model.loco, Type = "III")
# Genotype is not significant, P = 0.4742668
# Sex is significant, P = 4.972e-06 
# Day is not significant, P = 0.0902279
# The interaction is significant, P = 0.0002698

# Test for overdispersion
testDispersion(model.loco)
simulateResiduals(fittedModel = model.loco, plot = T)
# Model is overdispersed and residuals don't look good.

# Adjust for zeros with a zero-inflated model
model.loco.zi <- glmmTMB(cbind(active, inactive) ~ type * sex + day + (1|type:rep_pop) + (1|time), 
                         ziformula = ~1, data = Locoswap.data, family = binomial)
# ANOVA
Anova(model.loco.zi, Type = "III")
# Genotype is not significant, P = 0.19981
# Sex is significant, P = 0.01574 
# Day is not significant, P = 0.33773
# The interaction is significant, P = 6.45e-05

# Test for overdispersion
testDispersion(model.loco.zi)
simulateResiduals(fittedModel = model.loco.zi, plot = T)
# Fit is a bit better but still not ideal

# Last check if the zero-inflation is non-randomly distributed across genotypes
model.loco.zi2 <- glmmTMB(cbind(active, inactive) ~ type * sex + day + (1|type:rep_pop) + (1|time),
                          ziformula = ~ type * sex + day + (1|type:rep_pop) + (1|time), data = Locoswap.data, family = binomial)
# ANOVA
Anova(model.loco.zi2, Type = "III")
# Genotype is not significant, P = 0.19497
# Sex is significant, P = 0.01692
# Day is not significant, P = 0.34508 
# The interaction is significant, P = 7.581e-05

# Test for overdispersion
testDispersion(model.loco.zi2)
simulateResiduals(fittedModel = model.loco.zi2, plot = T)
# Fit is a bit better but still not ideal


# Last we compare the two zero-inflated models
AIC(model.loco.zi)
# 3238.217
AIC(model.loco.zi2)
# 3237.934
# The AIC values are quite similar, so we'll use the simpler model


# Use emmeans to test if  are significant different from each other
emmeans(model.loco.zi, list(pairwise ~ type | sex), adjust = "tukey")
# FLXsex is significant different from CFMsex and CFM auto in females



##########################################  PLOT DATA  #########################################


# Create new data frame for plotting
Locoswap.plot <- Locoswap.data[c(4, 10)]
# Remove NAs
Locoswap.plot <- Locoswap.plot[which(is.na(Locoswap.plot$locoactive) == F),]
# Add fitted values
Locoswap.plot$fitted <- fitted(model.loco.zi)

# Use emmeans to get SE
emmeans(model.loco.zi, list(pairwise ~ type | sex), adjust = "tukey", type = "response")

# Calculate mean
LocoSwap.M <- tapply(Locoswap.plot$fitted, Locoswap.plot$typesex, mean)
# SE
LocoSwap.SE <- c(0.0634, 0.0649, 0.0614, 0.0612, 0.0628, 0.0592, 0.0620, 0.0616) 


# The plot
pdf("Locomotion.pdf", family = "serif", width = 13, height = 11)

par(mar = c(5, 5, 2, 2))
plot(NULL, xlim = c(0.5, 9), xlab = "", xaxt = "n", ylim = c(0, 1.03), ylab = "Locomotion Activity",
     cex.lab = 1.8, cex.axis = 1.5, las = 1)
# Add points
stripchart(Locoswap.plot$locoactive ~ Locoswap.plot$typesex, add = T, vertical = T, 
           method = "jitter", pch = 16, cex = 1.8, col = "#f3f3f3", at = c(1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5))
# Add errorbar
X <- c(1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5)
errbar(X, LocoSwap.M, LocoSwap.M + LocoSwap.SE, LocoSwap.M - LocoSwap.SE, pch = 16, cex = 3, lwd = 3, add = T)
# Add axis
axis(1, at = c(1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5), cex.axis = 1.5, labels = c("FLXsex", "FLXauto", "CFMsex", "CFMauto",
                                                                           "FLXsex", "FLXauto", "CFMsex", "CFMauto"))
mtext("Female", side = 1, line = 3, at = 2.5, cex = 1.8)
mtext("Male", side = 1, line = 3, at = 7, cex = 1.9)
#SIG
mtext("b", side = 3, line = -1.6, at = 1, cex = 1.5)
mtext("ab", side = 3, line = -1.6, at = 2, cex = 1.5)
mtext("a", side = 3, line = -1.6, at = 3, cex = 1.5)
mtext("a", side = 3, line = -1.6, at = 4, cex = 1.5)
mtext("ab", side = 3, line = -1.6, at = 5.5, cex = 1.5)
mtext("ab", side = 3, line = -1.6, at = 6.5, cex = 1.5)
mtext("ab", side = 3, line = -1.6, at = 7.5, cex = 1.5)
mtext("ab", side = 3, line = -1.6, at = 8.5, cex = 1.5)

dev.off()

