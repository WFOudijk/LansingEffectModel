library(tidyverse)
library(mgcv)

#d <- read.csv("maternal_LE.csv")
d <- sub_maternal_LE
d <- sub_paternal_LE
d <- sub_both_LE
d <- myLongitudinalData
d <- subSexSpec

ggplot(d,aes(x=expectedAgeAtDeath)) + geom_histogram(bins=50)

## Let's cut off at 40
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath))

ggplot(d,aes(x=y1)) + geom_histogram(bins=50)
## Large peak at 40
## Not exactly looking like a beta distribution because of that.
## Better would be a zero one inflated beta. WO: [0, 1]
## This is only possible with brms as far as I know, WO: bayesian regression models 
## at least in combination with smoothers. WO: estimating a smooth trend by weighing averages of observations 

## Transform to [0,1]
d <- d %>% mutate(y2 = y1/40)

ggplot(d,aes(x=y2)) + geom_histogram(bins=50)


## Beta regression with spline "smoother"
m1 <- gam(y2 ~ s(ageOfMother), 
          family=betar(link="logit"), 
          data=d, 
          eps=0.001,
          method = "ML")
m1 <- gam(y2 ~ s(ageOfFather), 
          family=betar(link="logit"), 
          data=d, 
          eps=0.001,
          method = "ML")
m1 <- gam(y2 ~ s(ageOfFather) + s(ageOfMother), 
          family=betar(link="logit"), 
          data=d, 
          eps=0.001,
          method = "ML")

# add id to the dataframe 
d$identification <- NA
for (i in seq(1, 4000, 2)){
  d$identification[i] <- as.integer(i)
  d$identification[(i + 1)] <- as.integer(i)
}

m1 <- gam(y2 ~ sexOfParent + s(ageOfParent, by = as.factor(sexOfParent)) + s(identification, bs="re"), # the "re" is random effect. To correct for the individuals being present duplicates 
          family=betar(link="logit"), 
          data=d, 
          eps=0.001,
          method = "ML")

m1 <- gam(y2 ~ s(ageOfParent, by = as.factor(ID)),
          family=betar(link="logit"), 
          data=d, 
          eps=0.001,
          method = "ML")
# WO: by adding a smoothing spline, you can minimize the residual 
# GAM = generalized additive model; gam is a linear model which allows to learn non-linear features. 
# beta coefficients from linear model are replaced by a flexible, non-linear, function.
# this flexible function = spline. Sum of many splines form a gam. 
# smooth function = spline; can include a basis expansion, i.e, x^0, x^1, x^2 etc. In the smooth function
# there can be several weights and functions per variable. > more flexible. 
# gam then consists of several of these splines. 
# each spline also has a weight to prevent overfitting. 

# a link function can be seen as a transformation of the model's predictions 
# for example: a link function can result in predictions falling between 0 and 1 for analysis. 
# The underlying model allows for bigger or smaller values. 

# QUESTION: why the link logit? We already transformed the data to be [0, 1]

## Some slightly worrying warning. 
## Maybe disappears for better-looking distribution

summary(m1)
## significant effect age, appears to be linear (edf=1)
# WO: edf = effective degrees of freedom: 1 means straight line = linear 

## Transformation to plot fitted curve at original scale:
logist <- function(x) 40/(1+exp(-x))

## Plot the fit plus confidence band
plot(m1,all.terms = T,trans = logist)

## Seems a bit more extreme decline than suggested by the box plots


##########################################################################

## brms version; bayesian regression models using 'Stan'

library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

bf2 <- bf(y2 ~ s(ageOfMother)) + zero_one_inflated_beta()
bf2 <- bf(y2 ~ s(ageOfFather)) + zero_one_inflated_beta()
bf2 <- bf(y2 ~ s(ageOfParent, by = ID)) + zero_one_inflated_beta()


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
          file = "m2.n")
# j is with sex specific effects 
# k is longitudinal over ageOfParent 
# kl is bf2 with 'by = ID'
# m is males, by = ID. 

summary(m2)

pp_check(m2) 

## Compute predictions for range of ages
d.pred <- data.frame(ageOfMother=seq(0,40,0.1))
d.pred <- data.frame(ageOfFather=seq(0,40,0.1))
d.pred <- data.frame(ageOfParent=seq(0,40,0.1))



fit2 <- fitted(m2, newdata = d.pred, robust = T)
## Robust = T means use posterior median rather than mean
str(fit2)
## First col: prediction, third & fourth cols: limits CI

d.pred$y <- 40*fit2[,1]
d.pred$ymin <- 40*fit2[,3]
d.pred$ymax <- 40*fit2[,4]

ggplot(d.pred,aes(x=ageOfMother,y=y)) +
  theme_bw() +
  geom_line(col="blue",lwd=1.5) +
  geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.2)

ggplot(d.pred,aes(x=ageOfFather,y=y)) +
  theme_bw() +
  geom_line(col="blue",lwd=1.5) +
  geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.2)

ggplot(d.pred,aes(x=ageOfParent,y=y)) +
  theme_bw() +
  geom_line(col="blue",lwd=1.5) +
  geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.2)

## Somewhat less extreme than the gam fit. 
## But more nonlinear.

## Let's see if a linear model is "worse"


bf3 <- bf(y2 ~ ageOfMother) + zero_one_inflated_beta()
bf3 <- bf(y2 ~ ageOfFather) + zero_one_inflated_beta()
bf3 <- bf(y2 ~ ageOfParent) + zero_one_inflated_beta()


m3 <- brm(bf3,
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
          file = "m3.h")
m2 <- add_criterion(m2, criterion = "loo")
m3 <- add_criterion(m3, criterion = "loo")
loo_compare(m2,m3, criterion = "loo")
## m3 is slightly better but within 1 se.
## Evidence for nonlinearity weak.


##########################################################################

# lme4

d <- myLongitudinalData
d <- trackedIndividuals
males <- subset(trackedIndividuals, sexOfParent == "M")
females <- subset(trackedIndividuals, sexOfParent == "F")
d <- males
d <- females

length(unique(d$ID)) # number of parents = 981. The remainder did not have offspring 
# 478 unique fathers
# 500 unique mothers

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(ageOfParent)))
table(z$na) # indeed there are plenty with more than 1 ageOfParent!

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation (but this selects ...)
## Filter out the one with 1 observation (but this selects the individuals that die in the first time step)
# filters out the weaker individuals. probably. 


d <- d %>% filter(na>1)

d <- d %>% mutate(ID = factor(ID))

# compares 40 to the expected age at death. If expectedAgeAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath)) 
d <- d %>% mutate(y2 = y1/40)

ggplot(d,aes(x=y2)) + geom_histogram(bins = 20)

ggplot(d,aes(x=ageOfParent,y=y2)) +
  geom_point() +
  geom_smooth()
## very weak downward trend towards the right ...

library(lme4)
library(lmerTest)

m1 <- lmer(y2 ~ ageOfParent + (ageOfParent | ID), 
           data = d, 
           control=lmerControl(calc.derivs = F))
summary(m1)
# males: negative estimate with p-val = 3.92e-5
confint(m1)
# males: negative confidence intervals

plot(m1) # hmmm

## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

m1b <- lmer(y3 ~ ageOfParent + (ageOfParent | ID),data =d, control=lmerControl(calc.derivs = F))
summary(m1b)
# females: estimate is negative with a p-value of 0.017
# males: negative estimate with p-val = 1.27e-7
plot(m1b)

## Work with glmmTMB, which has beta distribution:

library(glmmTMB)
m1c <- glmmTMB(y2 ~ ageOfParent + (ageOfParent | ID), data =d, family = beta_family())
## A warning; don't know how serious
summary(m1c)
## So anyway, all models agree that there's a negative relation.
# females: estimate is negative with a p-value of 0.045 
# males: estimate is negative with p-val = 9.77e-5


## The gam model with beta family takes too long.
## I'm going to remove some data.
d2 <- d %>% filter(na>5)
## Use faster bam on logit transformed y
## bs="fs" means separate spline for each ID, same wigliness
m1i <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5),
           #family=betar(link="logit"),
           data=d,
           method = "REML")
# m1d = bam over the 500 tracked individuals!
# m1e = bam over 500 tracked individuals, based only on age-specific genes 
# m1f = bam over 500 tracked individuals, based only on binary genes. With sex implemented - but not used yet 
# m1g = females only 
# m1h = males only 

summary(m1i)
# females: edf of 2.249 with p-val = 6.59e*-6
# males: edf of 1 with p-val of 2e-16
gratia::draw(m1i, fun = logist)
## Overall curve linear (edf=1). Lots of variability

plot(m1e,all.terms = T,trans = logist)



