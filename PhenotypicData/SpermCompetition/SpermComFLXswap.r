################################################################################################
############################### FLX SWAP SPERM COMPETITION ASSAY ###############################
################################################################################################


# Set working directory
setwd("~/Library/CloudStorage/Box-Box/Work/Experiments/Drosophila/FLX/5.FLXchromosomeSwap/R/PhenotypicData/SpermCompetition")

# Set up environment
library(car)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(Hmisc)
library(lme4)

# Read in csv file with data
defswap.data <- read.csv(file = "defenceswap.csv", h = T, sep = ",", stringsAsFactors = T)
offswap.data <- read.csv(file = "offenceswap.csv", h = T, sep = ",", stringsAsFactors = T)


##########################################  STATISTIC  #########################################


######################## SPERM DEFENCE ########################

### GENERAL LINEAR MIXED MODELS ###

# General linear mixed-effect model testing if genotype and sex have a significant effect on sperm defence
model.def <- glmer(cbind(red, brown) ~ type + block + (1|type:rep_pop), data = defswap.data, family = binomial)
# ANOVA
Anova(model.def, Type = "III")
# Genotype is not significant, P = 0.4191
# Block is significant, P < 2e-16 

# Test for overdispersion
testDispersion(model.def)
simulateResiduals(fittedModel = model.def, plot = T)
# Clearly over dispersed


# Adjust for zeros with a zero-inflated model
model.def.zi <- glmmTMB(cbind(red, brown) ~ type + block + (1|type:rep_pop), ziformula = ~1, data = defswap.data, family = binomial)
# ANOVA
Anova(model.def.zi, Type = "III")
# Genotype is not significant, P = 0.1077
# Block is significant, P = 3.879e-14

# Test for overdispersion
testDispersion(model.def.zi)
simulateResiduals(fittedModel = model.def.zi, plot = T)
# Assumptions are good


# Last check if the zero-inflation is non-randomly distributed across genotypes
model.def.zi2 <- glmmTMB(cbind(red, brown) ~ type + block + (1|type:rep_pop), 
                         ziformula = ~ type + block + (1|type:rep_pop), data = defswap.data, family = binomial)
# ANOVA
Anova(model.def.zi2, Type = "III")
# Genotype is not significant, P = 0.107
# Block is significant, P = 1.484e-13

# Test for overdispersion
testDispersion(model.def.zi2)
simulateResiduals(fittedModel = model.def.zi2, plot = T)
# Assumptions are good

# Last we compare the two zero-inflated models
AIC(model.def.zi)
# 2140.287
AIC(model.def.zi2)
# 2149.283
# The simpler model is better


# Use emmeans to test if  are significant different from each other
emmeans(model.def.zi, list(pairwise ~ type), adjust = "tukey")
# Nothing is significant different from each other




######################## SPERM OFFENCE ########################

### GENERAL LINEAR MIXED MODELS ###

# General linear mixed-effect model testing if genotype and sex have a significant effect on sperm offence
model.off <- glmer(cbind(red, brown) ~ type + block + (1|type:rep_pop), data = offswap.data, family = binomial)
# ANOVA
Anova(model.off, Type = "III")
# Genotype is not significant, P = 0.9798
# Block is not significant, P = 0.3674 

# Test for overdispersion
testDispersion(model.off)
simulateResiduals(fittedModel = model.off, plot = T)
# Model is overdispersed and residuals don't look good


# Adjust for zeros with a zero-inflated model
model.off.zi <- glmmTMB(cbind(red, brown) ~ type + block + (1|type:rep_pop), ziformula = ~1, data = offswap.data, family = binomial)
# ANOVA
Anova(model.off.zi, Type = "III")
# Genotype is not significant, P = 0.8335
# Block is not significant, P = 0.6035

# Test for overdispersion
testDispersion(model.off.zi)
simulateResiduals(fittedModel = model.off.zi, plot = T)
# Fit is a bit better but still not ideal


# Last check if the zero-inflation is non-randomly distributed across genotypes
model.off.zi2 <- glmmTMB(cbind(red, brown) ~ type + block + (1|type:rep_pop), 
                         ziformula = ~ type + block + (1|type:rep_pop), data = offswap.data, family = binomial)
# ANOVA
Anova(model.off.zi2, Type = "III")
# Genotype is not significant, P = 0.8325
# Block is not significant, P = 0.6035

# Test for overdispersion
testDispersion(model.off.zi2)
simulateResiduals(fittedModel = model.off.zi2, plot = T)
# Assumptions acceptable

# Last we compare the two zero-inflated models
AIC(model.off.zi)
# 1402.367
AIC(model.off.zi2)
# 1407.998
# The simpler model is better


# Use emmeans to test if  are significant different from each other
emmeans(model.off.zi, list(pairwise ~ type), adjust = "tukey")
# Nothing is significant different from each other




##########################################  PLOT DATA  #########################################


######################## SPERM DEFENCE ########################


# Create new data frame for plotting
Defswap.plot <- defswap.data[c(2, 8)]
# Remove NAs
Defswap.plot <- Defswap.plot[which(is.na(Defswap.plot$prop_red) == F),]
# Add fitted values
Defswap.plot$fitted <- fitted(model.def.zi)

# Use emmeans to get SE
emmeans(model.def.zi, list(pairwise ~ type), adjust = "tukey", type = "response")

# Calculate mean
DefSwap.M <- tapply(Defswap.plot$fitted, Defswap.plot$type, mean)
# SE
DefSwap.SE <- c(0.0855, 0.0370, 0.0221, 0.0361) 


# The plot
pdf("Defence.pdf", family = "serif", width = 7, height = 7)

par(mar = c(3, 5, 2, 2))
plot(NULL, xlim = c(0.5, 4.5), xlab = "", xaxt = "n", ylim = c(0, 1.03), ylab = "Proportion of offspring",
     cex.lab = 1.8, cex.axis = 1.5, las = 1)
#Add points
stripchart(Defswap.plot$prop_red ~ Defswap.plot$type, add = T, vertical = T, 
           method = "jitter", pch = 16, cex = 1.8, col = "#f3f3f3", at = c(1, 2, 3, 4))
#Add errorbar
Xdef <- c(1, 2, 3, 4)
errbar(Xdef, DefSwap.M, DefSwap.M + DefSwap.SE, DefSwap.M - DefSwap.SE,
       pch = 16, cex = 3, lwd = 3, add = T)
axis(1, at = c(1, 2, 3, 4), cex.axis = 1.5, labels = c("FLXsex", "FLXauto", "CFMsex", "CFMauto"))
#ADD LETTERS
mtext("A", side = 3, line = 0.3, at = -0.2, cex = 1.5)

dev.off()



######################## SPERM OFFENCE ########################


# Create new data frame for plotting
Offswap.plot <- offswap.data[c(2, 8)]
# Remove NAs
Offswap.plot <- Offswap.plot[which(is.na(Offswap.plot$prop_red) == F),]
# Add fitted values
Offswap.plot$fitted <- fitted(model.off.zi)

# Use emmeans to get SE
emmeans(model.off.zi, list(pairwise ~ type), adjust = "tukey", type = "response")

# Calculate mean
OffSwap.M <- tapply(Offswap.plot$fitted, Offswap.plot$type, mean)
# SE
OffSwap.SE <- c(0.0126, 0.0123, 0.0211, 0.0113) 


# The plot
pdf("Offence.pdf", family = "serif", width = 7, height = 7)

par(mar = c(3, 5, 2, 2))
plot(NULL, xlim = c(0.5, 4.5), xlab = "", xaxt = "n", ylim = c(0, 1.03), ylab = "Proportion of offspring",
     cex.lab = 1.8, cex.axis = 1.5, las = 1)
#Add points
stripchart(Offswap.plot$prop_red ~ Offswap.plot$type, add = T, vertical = T, 
           method = "jitter", pch = 16, cex = 1.8, col = "#f3f3f3", at = c(1, 2, 3, 4))
#Add errorbar
Xoff <- c(1, 2, 3, 4)
errbar(Xoff, OffSwap.M, OffSwap.M + OffSwap.SE, OffSwap.M - OffSwap.SE,
       pch = 16, cex = 3, lwd = 3, add = T)
axis(1, at = c(1, 2, 3, 4), cex.axis = 1.5, labels = c("FLXsex", "FLXauto", "CFMsex", "CFMauto"))
#ADD LETTERS
mtext("B", side = 3, line = 0.3, at = -0.2, cex = 1.5)

dev.off()


