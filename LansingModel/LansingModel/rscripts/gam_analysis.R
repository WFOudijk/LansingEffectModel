################################################################################
# GAM ANALYSIS
################################################################################
library(tidyverse)
library(mgcv)

# combined data 

# quality-only. MutationProbAgeSpecificGenes, meanMutBias and sdMutEffectSize varied. 
d <- sampled_data_tot 
# damage-only. MutationProbGametes and mutationProbStemCells varied together. 
d <- sampled_data
# investment-only. MutationProbInvestmentGenes varied. 
d <- sampled_data_investment_tot
# weight investment varied
d <- sampled_data_weight_investment

# check number of parents.
length(unique(d$ID))  

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(ageOfParent)))
table(z$na) 

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation (but this selects ...)
d <- d %>% filter(na>1)

# change column ID into a factor 
d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If expectedAgeAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath)) 
# scale the expected age at deaths to the maximum age. 
d <- d %>% mutate(y2 = y1/40)
## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

# GAM analysis
# remove some data otherwise it takes too long 
d2 <- d %>% filter(na>15)

##########################
#for damage-only scenario
##########################

## Model without interaction
m1z <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k= 5) +
             s(mutationProbGametes, k = 5),
           #family=betar(link="logit"),
           data=d2,
           method = "REML")
summary(m1z)

## Model with interaction (ti term)
m1zb <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) +
              s(mutationProbGametes, k = 5) +
              ti(ageOfParent, mutationProbGametes, k = c(5,5)) + 
              ti(ageOfParent, mutationProbSC, k = c(7,7)),
            #family=betar(link="logit"),
            data=d2,
            method = "REML")
summary(m1zb)
## ti effect significant (barely)
# ti between Sc and age of parent significant, but ti between gamete and age of parent 
# not ... Why? They mutate at the same rate. 

# fit linear model to see if it is a better fit 
# choose one mutation prob to perform both model analysis on 
dat <- d2[d2$mutationProbGametes == 0.0004,]
m1zd <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5),
            data = dat,
            method = "REML") # bam analysis
# use bam for linear model 
# random slope model 
mlin <- bam(y3 ~ ageOfParent + s(ageOfParent, ID, bs = "re"), data = dat, method = "REML")
AIC(mlin, m1zd)
# non - linear is better.  

AIC(m1z,m1zb)
## Appears to be an interaction since m1zb better
## Surprising that m1zb has smaller df ...

gratia::draw(m1zb)

## Removing the random effect to be able to use vis.gam for 3D plot
m1zc <- bam(y3 ~ s(ageOfParent, k = 5) + 
              s(mutationProbGametes, k = 5) +
              ti(ageOfParent, mutationProbGametes, k = c(5,5)),
            #family=betar(link="logit"),
            data=d2,
            method = "REML")

vis.gam(m1zc,theta=130,phi=20,ticktype="detailed",n.grid=40)
## Not very clear evidence interaction

# transform the predictor parameters into factors 
d2 <- d2 %>% mutate(ID = factor(ID)) 
d2 <- d2 %>% mutate(mutationProbGametes = factor(mutationProbGametes)) 
# new predicted data for determining the percentage decrease
pred_data = expand.grid(ageOfParent = c(0,40), ID = levels(d2$ID)[1], 
                        mutationProbGametes = levels(d2$mutationProbGametes)[1])

# GENERATE NEW DATA FOR PLOT
len <- 20
length(unique(d2$mutationProbGametes)) #10
length(unique(d2$mutationProbSC)) #10

# I'm using 20 age values, and 4 for each of the 2 other predictors.
# Could use all of course, but then the plot becomes a bit messy.
# A single ID value because we will ignore the term involving ID anyway.
d.pred <- expand.grid(ageOfParent=seq(0,40,length=len),
                      ID=levels(d2$ID)[1],
                      mutationProbGametes=unique(d2$mutationProbGametes),
                      mutationProbSC=unique(d2$mutationProbSC))

# Compute all terms, excluding the one with ID
z <- predict(m1zb,newdata = d.pred,
             type = "terms",
             exclude = c("s(ageOfParent,ID)"))
str(z) # contains 9 terms (the s and ti terms, not including the intercept)
attr(z,"constant") # here is the intercept

# Add intercept and all terms, transform back
d.pred$y3 <- attr(z,"constant") + apply(z,1,sum) 
logist <- function(x) 1/(1+exp(-x))
d.pred$y4 <- 40*logist(d.pred$y3)

# Use predictors as factors for use in ggplot
d.pred <- d.pred %>% mutate(mutationProbGametes=factor(mutationProbGametes),
                            mutationProbSC=factor(mutationProbSC))

ggplot(d.pred,aes(x=ageOfParent,y=y4,color=mutationProbGametes)) +
  theme_cowplot() +
  geom_line() +
  facet_wrap(d.pred$mutationProbSC) +
  background_grid(major = "xy")


##########################
#for quality-only scenario
##########################
# model without interaction
m1z <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) 
           + s(mutationProbAgeGenes, k = 5) 
           + s(meanMutBias, k = 5) 
           + s(sdMutEffectSize, k = 5),
           data=d2, 
           method="REML") 
summary(m1z)

# model with interaction
m1zb <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) +
             s(mutationProbAgeGenes, k = 5) +
             s(meanMutBias, k = 5) +
             s(sdMutEffectSize, k = 5) +
             ti(ageOfParent,mutationProbAgeGenes,k=c(5,5)) +
             ti(ageOfParent,meanMutBias,k=c(5,5)) +
             ti(ageOfParent,sdMutEffectSize,k=c(5,5)) +
             ti(mutationProbAgeGenes,meanMutBias, k=c(5,5)) +
             ti(mutationProbAgeGenes, sdMutEffectSize, k=c(5,5)) +
             ti(meanMutBias, sdMutEffectSize, k=c(9,9)),
             #ti(meanMutBias, mutationProbAgeGenes, k=c(5,5)) + 
             #ti(sdMutEffectSize, mutationProbAgeGenes, k=c(5,5)),
             #ti(sdMutEffectSize, meanMutBias, k=c(5,5)), 
           data=d2, 
           method="REML") 
summary(m1zb)
# For every smooth and interaction term p-val < 2e-16

AIC(m1z,m1zb)
# m1zb smaller dfs. m1z higher AIC (barely). 

gratia::draw(m1zb)

# removing random effect for 3D plot 
m1zc <- bam(y3 ~ s(ageOfParent, k = 5) +
              s(mutationProbAgeGenes, k = 5) +
              s(meanMutBias, k = 5) +
              s(sdMutEffectSize, k = 5) +
              ti(ageOfParent,mutationProbAgeGenes,k=c(5,5)) +
              ti(ageOfParent,meanMutBias,k=c(5,5)) +
              ti(ageOfParent,sdMutEffectSize,k=c(5,5)), # could add more ti terms here
            data=d2, 
            method="REML") 

vis.gam(m1zc,theta=130,phi=20,ticktype="detailed",n.grid=40)


# transforms the parameter model predictors to factors 
d2 <- d2 %>% mutate(ID = factor(ID))
d2 <- d2 %>% mutate(mutationProbAgeGenes = factor(mutationProbAgeGenes))
d2 <- d2 %>% mutate(meanMutBias = factor(meanMutBias))
d2 <- d2 %>% mutate(sdMutEffectSize = factor(sdMutEffectSize))
# new predicted data for determining the percentage decrease 
pred_data = expand.grid(ageOfParent = c(0,40), ID = levels(d2$ID)[1], 
                        mutationProbAgeGenes = levels(d2$mutationProbAgeGenes)[1],
                        meanMutBias = levels(d2$meanMutBias)[1],
                        sdMutEffectSize = levels(d2$sdMutEffectSize)[1])

### GENERATE NEW DATA FOR PLOT
len <- 20
length(unique(d2$mutationProbAgeGenes)) #11
length(unique(d2$meanMutBias)) #10
length(unique(d2$sdMutEffectSize)) #10

# I'm using 20 age values, and 4 for each of the 3 other predictors.
# Could use all of course, but then the plot becomes a bit messy.
# A single ID value because we will ignore the term involving ID anyway.
d.pred <- expand.grid(ageOfParent=seq(0,40,length=len),
                      ID=levels(d2$ID)[1],
                      mutationProbAgeGenes=unique(d2$mutationProbAgeGenes),
                      meanMutBias=unique(d2$meanMutBias)[c(1,4,5,6)],
                      sdMutEffectSize=unique(d2$sdMutEffectSize)[c(1,4,7,10)])
d.pred <- expand.grid(ageOfParent=seq(0,40,length=len),
                      ID=levels(d2$ID)[1],
                      mutationProbAgeGenes=unique(d2$mutationProbAgeGenes),
                      meanMutBias=unique(d2$meanMutBias),
                      sdMutEffectSize=unique(d2$sdMutEffectSize))

# Compute all terms, excluding the one with ID
z <- predict(m1zb,newdata = d.pred,
             type = "terms",
             exclude = c("s(ageOfParent,ID)"))
str(z) # contains 9 terms (the s and ti terms, not including the intercept)
attr(z,"constant") # here is the intercept

# Add intercept and all terms, transform back
d.pred$y3 <- attr(z,"constant") + apply(z,1,sum) 
logist <- function(x) 1/(1+exp(-x))
d.pred$y4 <- 40*logist(d.pred$y3)

# Use predictors as factors for use in ggplot
d.pred <- d.pred %>% mutate(mutationProbAgeGenes=factor(mutationProbAgeGenes),
                            meanMutBias=factor(meanMutBias),
                            sdMutEffectSize=factor(sdMutEffectSize))

#library(cowplot)
ggplot(d.pred,aes(x=ageOfParent,y=y4,color=mutationProbAgeGenes)) +
  theme_cowplot() +
  geom_line() +
  facet_wrap(meanMutBias~sdMutEffectSize) +
  background_grid(major = "xy")
# the 'bump' is peculiar. However, the combination of parameters which show the bump is not present 
# in original data. Might be explanation that model fits incorrectly. 
# running the model with combination which shows bump. 
# mean = -0.015; sd = 0.054; mutation rate age genes = 0.003
dataTest <- read.table("/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/CombiningDataForGam/outputLETrackedIndividuals.txt")
dataTest <- read.table("/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/CombiningDataForGam/outputLETrackedIndividualsTest.txt")
colnames(dataTest) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "mean", "sd", "mutationProbAgeGenes")

# try smooth function from ggplot 
ggplot(dataTest, aes(ageOfParent, expectedAgeAtDeath)) +
  geom_smooth() 
# no bump 

d <- dataTest
# perform d transformations as mentioned above 
m1f <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5),
           data = d2, 
           method = "REML")
d.pred.f <- expand.grid(ageOfParent = seq(0,40),
                        ID = levels(d2$ID)[1])
x <- predict(m1f, d.pred.f, type = "terms", 
             exclude = c("s(ageOfParent,ID)"))
d.pred.f$y3 <- attr(x, "constant") + apply(x, 1, sum)
d.pred.f$y4 <- 40*logist(d.pred.f$y3)

ggplot(d.pred.f, aes(ageOfParent, y4)) + geom_point() + geom_line(alpha = 0.5)
# No bump. 

# per parameter. 
t <- subset(d.pred, d.pred$meanMutBias == -0.022 & d.pred$sdMutEffectSize == 0.024)
ggplot(t, aes(x = ageOfParent, y = y4, color = mutationProbAgeGenes)) +
  geom_line()
t <- subset(d.pred, d.pred$mutationProbAgeGenes == 0.003 & d.pred$sdMutEffectSize == 0.024)
ggplot(t, aes(x = ageOfParent, y = y4, color = meanMutBias)) +
  geom_line()
t <- subset(d.pred, d.pred$mutationProbAgeGenes == 0.003 & d.pred$meanMutBias == -0.022)
ggplot(t, aes(x = ageOfParent, y = y4, color = sdMutEffectSize)) +
  geom_line()
# they resemble the ggplot geom_smooth a lot. 

# number of parents per age class
ggplot(d, aes(x = ageOfParent)) +
  geom_density()

# get averages
averages <- aggregate(d$expectedAgeAtDeath, by = list(d$ageOfParent, d$mutationProbAgeGenes, d$meanMutBias, d$sdMutEffectSize), mean)
mean(d[d$mutationProbAgeGenes == 0.0004 & d$meanMutBias == -0.015 & d$sdMutEffectSize == 0.054,]$expectedAgeAtDeath)

#######################
## Predict, only using the first term (i.e. without the "random" effect)
z <- predict(m1zb, newdata = pred_data, type="terms",terms = c("s(ageOfParent)", "ti(ageOfParent,mutationProbAgeGenes)"))
z <- predict(m1zb, newdata = d2, type="terms",terms = c("s(ageOfParent)", "ti(ageOfParent,mutationProbAgeGenes)"))
z
str(z)
## Use logist transformed:
logist <- function(x) 40/(1+exp(-x))
zl <- logist(z)

## Percentage decrease:
perc_decr <- 100*(zl[1]-zl[2])/zl[1]
perc_decr

summary(m1z)
ggplot(data = data.frame(zl), aes("s(ageOfParent)", "ti(ageOfParent,mutationProbAgeGenes)")) +
  geom_point()

#gratia::draw(m1z, fun = logist, main = "test")
plot(m1z, trans = logist, 
     #xlab = "Parental age",
     #ylab = "Expected age at death of offspring",
     #main = paste("% decrease = ", perc_decr)
)

################################
# combining both model mechanisms
################################
# first one more look at quality-only 
# mean = - 0.022; sd = 0.024; mutation rate = 0.003
quality_only <- read.table("/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/CombiningDataForGam/outputLETrackedIndividualsOptParams.txt")
colnames(quality_only) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "mean", "sd", "mutationProbAgeGenes")

d <- quality_only

# check number of parents.
length(unique(d$ID))  

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(ageOfParent)))
table(z$na) 

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation (but this selects ...)
d <- d %>% filter(na>1)

# change column ID into a factor 
d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If expectedAgeAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath)) 
# scale the expected age at deaths to the maximum age. 
d <- d %>% mutate(y2 = y1/40)
## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

# GAM analysis
# remove some data otherwise it takes too long 
d2 <- d %>% filter(na>11)

m1a <- bam(y3 ~ s(ageOfParent, k = 5) +
             s(ageOfParent, ID, bs = "fs", k = 5 ),
           data = d2, 
           method = "REML")
summary(m1a)
logist <- function(x) 40/(1+exp(-x))

d2 <- d2 %>% mutate(ID = factor(ID)) 

pred_data = expand.grid(ageOfParent = c(0,40), ID = levels(d2$ID)[1])

## Predict, only using the first term (i.e. without the "random" effect)
z <- predict(m1a, newdata = pred_data, type="terms",terms = c("s(ageOfParent)"))
z
str(z)
## Use logist transformed:
zl <- logist(z)

## Percentage decrease:
perc_decr <- 100*(zl[1]-zl[2])/zl[1]
perc_decr

plot(m1a, trans = logist,
     #main = paste("Quality-only, % decrease =", perc_decr),
     main = paste("combining damage and quality. % decrease =", perc_decr),
     subtitle = "Mutation rate = 0.003; mean = -0.022; sd = 0.024",
     xlab = "Age of parent",
     ylab = "Expected age at death")

# look at damage-only with gametes = stem cells = 0.0024
damage_only <- read.table("/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/CombiningDataForGam/outputLETrackedIndividualsOptParamsDamage.txt")
colnames(damage_only) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "mean", "sd", "mutationProbAgeGenes")

d <- damage_only

# check number of parents.
length(unique(d$ID))  

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(ageOfParent)))
table(z$na) 

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation (but this selects ...)
d <- d %>% filter(na>1)

# change column ID into a factor 
d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If expectedAgeAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath)) 
# scale the expected age at deaths to the maximum age. 
d <- d %>% mutate(y2 = y1/40)
## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

# GAM analysis
# remove some data otherwise it takes too long 
d2 <- d %>% filter(na>11)

m1a <- bam(y3 ~ s(ageOfParent, k = 5) +
             s(ageOfParent, ID, bs = "fs", k = 5 ),
           data = d2, 
           method = "REML")
summary(m1a)
logist <- function(x) 40/(1+exp(-x))

d2 <- d2 %>% mutate(ID = factor(ID)) 

pred_data = expand.grid(ageOfParent = c(0,40), ID = levels(d2$ID)[1])

## Predict, only using the first term (i.e. without the "random" effect)
z <- predict(m1a, newdata = pred_data, type="terms",terms = c("s(ageOfParent)"))
z
str(z)
## Use logist transformed:
zl <- logist(z)

## Percentage decrease:
perc_decr <- 100*(zl[1]-zl[2])/zl[1]
perc_decr

plot(m1a, trans = logist,
     main = paste("Quality-only, % decrease =", perc_decr),
     subtitle = "Mutation rate = 0.003; mean = -0.022; sd = 0.024",
     xlab = "Age of parent",
     ylab = "Expected age at death")


##########################
#for investment-only scenario
##########################
# model without interaction
m1z <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) 
           + s(mutationProbInvestmentGenes, k = 5) +
             s(sdInvestmentGenes, k = 5),
           data=d2, 
           method="REML") 
summary(m1z)

# model with interaction
m1zb <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) +
              s(mutationProbInvestmentGenes, k = 5) +
              s(sdInvestmentGenes, k = 5) +
              ti(ageOfParent,mutationProbInvestmentGenes, k=c(5,5)) +
              ti(ageOfParent, sdInvestmentGenes, k= c(5,5)),
            data=d2, 
            method="REML") 
summary(m1zb)

AIC(m1z,m1zb)
# m1zb has lower value so better fit 

### GENERATE NEW DATA FOR PLOT
len <- 20
length(unique(d2$mutationProbInvestmentGenes)) #10
length(unique(d2$sdInvestmentGenes))#10

# I'm using 20 age values, and 4 for each of the 2 other predictors.
# Could use all of course, but then the plot becomes a bit messy.
# A single ID value because we will ignore the term involving ID anyway.
d.pred <- expand.grid(ageOfParent=seq(0,40,length=len),
                      ID=levels(d2$ID)[1],
                      mutationProbInvestmentGenes=unique(d2$mutationProbInvestmentGenes),
                      sdInvestmentGenes = unique(d2$sdInvestmentGenes)[c(1,2,3,5,7,9)])

# Compute all terms, excluding the one with ID
z <- predict(m1zb,newdata = d.pred,
             type = "terms",
             exclude = c("s(ageOfParent,ID)"))
str(z) # contains 9 terms (the s and ti terms, not including the intercept)
attr(z,"constant") # here is the intercept

# Add intercept and all terms, transform back
d.pred$y3 <- attr(z,"constant") + apply(z,1,sum) 
logist <- function(x) 1/(1+exp(-x))
d.pred$y4 <- 40*logist(d.pred$y3)

# Use predictors as factors for use in ggplot
d.pred <- d.pred %>% mutate(mutationProbInvestmentGenes=factor(mutationProbInvestmentGenes))
d.pred <- d.pred %>% mutate(sdInvestmentGenes=factor(sdInvestmentGenes))

#library(cowplot)
ggplot(d.pred,aes(x=ageOfParent,y=y4,color=mutationProbInvestmentGenes)) +
  labs(title = "resource-only scenario",
       subtitle = "where every panel shows a different sd",
       y = "expected age at death offspring") +
  theme_cowplot() +
  geom_line() +
  facet_wrap(d.pred$sdInvestmentGenes) +
  background_grid(major = "xy")


