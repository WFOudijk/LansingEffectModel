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
d2 <- d %>% filter(na>12)

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
              ti(ageOfParent, mutationProbGametes, k = c(5,5)),
            #family=betar(link="logit"),
            data=d2,
            method = "REML")
summary(m1zb)
## ti effect significant (barely)

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

# test 
d2 <- d2 %>% mutate(mutationProbGametes = factor(mutationProbGametes))
d2 <- d2 %>% mutate(ID = factor(ID))
d2 <- d2 %>% mutate(ageOfParent = factor(ageOfParent))

pred_data = expand.grid(ageOfParent = levels(d2$ageOfParent), ID = levels(d2$ID),
                        mutationProbGametes = levels(d2$mutationProbGametes))

pred_data$z <- predict(m1zb, newdata = pred_data, type="terms",terms = c("s(ageOfParent)", "ti(ageOfParent,mutationProbGametes)"))
pred_data$z <- logist(pred_data$z)

d2$z <- predict(m1zb)
d2$z <- logist(d2$z)
ggplot(d2, aes(ageOfParent, z, group = factor(mutationProbGametes), color = factor(mutationProbGametes))) + geom_smooth() 

ggplot(pred_data, aes(ageOfParent, z[,1], group = factor(mutationProbGametes), color = factor(mutationProbGametes))) +
  geom_point()
d3 <- d2
d3$z <- predict(m1zb, newdata = pred_data, type="terms",terms = c("s(ageOfParent)", "ti(ageOfParent,mutationProbGametes)"))
ggplot(data = d3, aes(x = ageOfParent, z, group = factor(mutationProbGametes), color = factor(mutationProbGametes))) +
  geom_smooth()

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
             ti(mutationProbAgeGenes, sdMutEffectSize, k=c(5,5)),
             #ti(meanMutBias, sdMutEffectSize, k=c(5,5)),
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
pred_data = expand.grid(ageOfParent = seq(0,40), ID = levels(d2$ID),
                        expectedAgeAtDeath = d2$expectedAgeAtDeath,
                        mutationProbAgeGenes = levels(d2$mutationProbAgeGenes),
                        meanMutBias = levels(d2$meanMutBias)[5],
                        sdMutEffectSize = levels(d2$sdMutEffectSize)[5])
t = expand.grid(ageOfParent = seq(0, 40), ID = levels(d2$ID)[1], 
                 mutationProbAgeGenes = levels(d2$mutationProbAgeGenes))
t$z <- predict(m1zb, t)
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
# try bayesian 
################################

library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

bf2 <- bf(y2 ~ s(ageOfParent) + s(ageOfParent, ID, bs = "fs", k = 5) +
  s(mutationProbGametes, k = 5) +
  t2(ageOfParent, mutationProbGametes, k = c(5,5))) +
  zero_one_inflated_beta()

my_priors <- prior(normal(0,1), class = Intercept) +
  prior(normal(0,1), class = b)

m2 <- brm(bf2,
          data = d,
          prior = my_priors,
          warmup = 1000,
          iter = 2000,
          chains = 4,
          cores = 4,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99,
                         max_treedepth = 10),
          seed = 667,
          file = "model1")
### FATAL ERROR: a precompiled header has been changed. 
# please rebuild precompiled header 'stan/src/stan/model/model_header.hpp.gch'



