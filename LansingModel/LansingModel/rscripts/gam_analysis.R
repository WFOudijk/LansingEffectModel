################################################################################
# GAM ANALYSIS
################################################################################

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
d2 <- d %>% filter(na>15)

##########################
#for damage-only scenario
##########################
m1z <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) 
           + ti(mutationProbGametes, k = 5),
           #family=betar(link="logit"),
           data=d2,
           method = "REML")
# transform the predictor parameters into factors 
d2 <- d2 %>% mutate(ID = factor(ID)) 
d2 <- d2 %>% mutate(mutationProbGametes = factor(mutationProbGametes)) 
# new predicted data for determining the percentage decrease
pred_data = expand.grid(ageOfParent = c(0,40), ID = levels(d2$ID)[1], 
                        mutationProbGametes = levels(d2$mutationProbGametes)[1])

##########################
#for quality-only scenario
##########################
m1z <- bam(y3 ~ s(ageOfParent, k = 5) + s(ageOfParent, ID, bs = "fs", k = 5) 
           + ti(mutationProbAgeGenes, k = 5) 
           + ti(meanMutBias, k = 5) 
           + ti(sdMutEffectSize, k = 5),
           data=d2, 
           method="REML") 
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

#######################
## Predict, only using the first term (i.e. without the "random" effect)
z <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(ageOfParent)"))
z
str(z)
## Use logist transformed:
zl <- logist(z)

## Percentage decrease:
perc_decr <- 100*(zl[1]-zl[2])/zl[1]
perc_decr

summary(m1z)

#gratia::draw(m1z, fun = logist, main = "test")
plot(m1z, trans = logist, 
     #xlab = "Parental age",
     #ylab = "Expected age at death of offspring",
     #main = paste("% decrease = ", perc_decr)
)
