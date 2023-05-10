###############################################################################
# COMBINING OF THE MECHANISMS 
###############################################################################
library(cowplot)
library(tidyverse)
library(mgcv)
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/combiningAll/"

# 1. evolutionary consequences 
deathPop <- read.table(paste(path, "null/outputDeclineGameteQuality.txt", sep = "")) # mut rate = 0.001; mean = -0.02; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.02; sd = 0.01"
deathPop <- read.table(paste(path, "null/outputDeclineGameteQuality2.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01"
deathPop <- read.table(paste(path, "null/outputDeclineGameteQuality3.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01. Pop size = 10.000 and tEnd = 100.000
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01"
deathPop <- read.table(paste(path, "null/outputDeclineGameteQuality4.txt", sep = "")) # Null. Pop size = 10.000 and tEnd = 10.000
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01. tEnd = 10.000"
deathPop <- read.table(paste(path, "nullDamage/outputDeclineGameteQuality.txt", sep = "")) # Null + damage 
plotText <- "Null + damage"
deathPop <- read.table(paste(path, "nullDamage/outputDeclineGameteQuality2.txt", sep = "")) # Null + damage; tEnd = 10.000 
plotText <- "Null + damage, tEnd = 10.000"
deathPop <- read.table(paste(path, "nullQuality/outputDeclineGameteQuality.txt", sep = "")) # Null + quality 
plotText <- "Null + quality"
deathPop <- read.table(paste(path, "nullQuality/outputDeclineGameteQuality2.txt", sep = "")) # Null + quality; tEnd = 10.000
plotText <- "Null + quality, tEnd = 10.000"
deathPop <- read.table(paste(path, "nullResource/outputDeclineGameteQuality.txt", sep = "")) # Null + resource, tEnd = 100.000
plotText <- "Null + resource distribution"
deathPop <- read.table(paste(path, "nullResource/outputDeclineGameteQuality2.txt", sep = "")) # Null + resource; tEnd = 10.000
plotText <- "Null + resource distribution; tEnd = 10.000"
deathPop <- read.table(paste(path, "damageQuality/outputDeclineGameteQuality.txt", sep = "")) # quality + damage, tEnd = 100.000
plotText <- "Quality + damage; tEnd = 100.000"
deathPop <- read.table(paste(path, "damageResource/outputDeclineGameteQuality.txt", sep = "")) # damage + resource
plotText <- "Damage + resource; tEnd = 100.000"
deathPop <- read.table(paste(path, "qualityResource/outputDeclineGameteQuality.txt", sep = "")) # quality + resource, tEnd = 100.000
plotText <- "Quality + resource; tEnd = 100.000"
deathPop <- read.table(paste(path, "qualityResource/outputDeclineGameteQuality2.txt", sep = "")) # quality + resource, tEnd = 100.000. Again with code fix to prevent negative survival probs 
plotText <- "Quality + resource; tEnd = 100.000"
deathPop <- read.table(paste(path, "damageOnly/outputDeclineGameteQuality.txt", sep = "")) # damage only. tEnd = 100.000
plotText <- "Damage-only"
deathPop <- read.table(paste(path, "qualityOnly/outputDeclineGameteQuality.txt", sep = "")) # quality only. tEnd = 100.000
plotText <- "Quality-only"
deathPop <- read.table(paste(path, "resourceOnly/outputDeclineGameteQuality.txt", sep = "")) # resource only. tEnd = 100.000
plotText <- "Resource-only"
deathPop <- read.table(paste(path, "longrunNull/outputDeclineGameteQuality.txt", sep = "")) # Null long run, tEnd = 500.000  
plotText <- "Null, long run"

deathPop <- deathPop[, 1:6] # select only the columns that are useful 
colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb")

# subset data to get info every 100th time step 
deathPop <- deathPop[deathPop$time %in% seq(0, max(deathPop$time), 100),]

# get average age at death per time step 
averageData <- aggregate(deathPop$ageAtDeath, list(deathPop$time), mean)
colnames(averageData) <- c("time", "meanAgeAtDeath")

# get average age at death per time step 
averageDataSurvProb <- aggregate(deathPop$survivalProb, list(deathPop$time), mean)
colnames(averageDataSurvProb) <- c("time", "meanSurvProb")

# plot the average 
ggplot(averageData, aes(time, meanAgeAtDeath)) + 
  geom_line(alpha = 0.2) + 
  geom_point() + 
  theme_cowplot() + 
  labs(title = "Evolutionary consequences",
       subtitle = plotText,
       y = "average age at death") 
 # ylim(0, 15)

# plot the average survival prob over time 
ggplot(averageDataSurvProb, aes(time, meanSurvProb)) + 
  geom_line(alpha = 0.2) + 
  geom_point() + 
  theme_cowplot() + 
  labs(title = "Evolutionary consequences",
       subtitle = plotText,
       y = "average survival prob") 
# ylim(0, 15)

# use geom_smooth for the whole data 
ggplot(deathPop, aes(time, ageAtDeath)) + 
  geom_smooth() + 
  theme_cowplot() +
  labs(title = "Evolutionary consequences",
       subtitle = paste("Using ggplot gam model to smoothen. ", plotText),
       y = "age at death")  +
  coord_cartesian(ylim = c(0, 5))

# determine average survival probability to see if ageing occurs 
averageSurvProb <- aggregate(deathPop$survivalProb, list(deathPop$time, deathPop$ageAtDeath), mean)
colnames(averageSurvProb) <- c("time", "ageAtDeath", "survivalProb")
# plot the data where every panel is an age class 
ggplot(averageSurvProb, aes(time, survivalProb)) + 
  geom_line() + 
  facet_wrap(averageSurvProb$ageAtDeath) +
  labs(title = "Survival probability over time",
       subtitle = paste("Where every panel is an age class. ", plotText), 
       y = "average survival probability") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90))

# different way of plotting to see if ageing occurs 
ggplot(averageSurvProb, aes(ageAtDeath, survivalProb, group = time, colour = time)) +
  geom_line() +
  scale_color_viridis_c(option = "A") +
  labs(title = "survival probability over age classes",
       subtitle = paste("With time shown as colours. ", plotText)) +
  theme(text = element_text(size = 20))

  # 2. Expected age at death over parental age 
# make empty data frame
dataTotalLansing <- c()

myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals.txt", sep = "")) # mut rate = 0.001; mean = -0.02; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.02; sd = 0.01"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals2.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals3.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01; pop size = 10.000 and tEnd = 100.000
plotText <- "null"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals4.txt", sep = "")) # null; pop size = 10.000 and tEnd = 10.000
plotText <- "Null. tEnd = 10.000; pop size = 10.000"
myLongitudinalData <- read.table(paste(path, "nullDamage/outputLETrackedIndividuals.txt", sep = "")) # null + damage
plotText <- "Null + damage"
myLongitudinalData <- read.table(paste(path, "nullDamage/outputLETrackedIndividuals2.txt", sep = "")) # null + damage; tEnd = 10.000
plotText <- "Null + damage, tEnd = 10.000"
myLongitudinalData <- read.table(paste(path, "nullQuality/outputLETrackedIndividuals.txt", sep = "")) # null + quality
plotText <- "Null + quality"
myLongitudinalData <- read.table(paste(path, "nullQuality/outputLETrackedIndividuals2.txt", sep = "")) # null + quality; tEnd = 10.000
plotText <- "Null + quality, tEnd = 10.000"
myLongitudinalData <- read.table(paste(path, "nullResource/outputLETrackedIndividuals.txt", sep = "")) # null + resource distribution
plotText <- "Null + resource distribution"
myLongitudinalData <- read.table(paste(path, "nullResource/outputLETrackedIndividuals2.txt", sep = "")) # null + resource distribution; tEnd = 10.000
plotText <- "Null + resource distribution, tEnd = 10.000"
myLongitudinalData <- read.table(paste(path, "damageQuality/outputLETrackedIndividuals.txt", sep = "")) # quality + damage; tEnd = 100.000
plotText <- "Quality + damage"
myLongitudinalData <- read.table(paste(path, "qualityResource/outputLETrackedIndividuals.txt", sep = "")) # quality + resource; tEnd = 100.000
plotText <- "Quality + resource, tEnd = 100.000"
myLongitudinalData <- read.table(paste(path, "qualityResource/outputLETrackedIndividuals2.txt", sep = "")) # quality + resource; tEnd = 100.000; With edited code to prevent possibility of surv probs becoming negative 
plotText <- "Quality + resource"
myLongitudinalData <- read.table(paste(path, "damageOnly/outputLETrackedIndividuals.txt", sep = "")) #damage - only
plotText <- "Damage-only"
myLongitudinalData <- read.table(paste(path, "qualityOnly/outputLETrackedIndividuals.txt", sep = "")) #quality - only
plotText <- "Quality-only"
myLongitudinalData <- read.table(paste(path, "damageResource/outputLETrackedIndividuals.txt", sep = "")) # damage + resource
plotText <- "Damage + resource"
myLongitudinalData <- read.table(paste(path, "resourceOnly/outputLETrackedIndividuals.txt", sep = "")) #  resource - only
plotText <- "Resource-only, tEnd = 100.000"
myLongitudinalData <- read.table(paste(path, "longrunNull/outputLETrackedIndividuals.txt", sep = "")) #  null, tEnd = 500.000
plotText <- "Long run of null scenario, tEnd = 500.000"

myLongitudinalData <- myLongitudinalData[, 1:5] # subset only the useful columns
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

avg <- aggregate(myLongitudinalData$expectedAgeAtDeath, list(myLongitudinalData$ageOfParent), mean)
avg$scenario <- plotText
colnames(avg) <- c("ageOfParent", "averageExpectedAgeAtDeath", "scenario")
dataTotalLansing <- rbind(dataTotalLansing, avg)

ggplot(dataTotalLansing, aes(ageOfParent, averageExpectedAgeAtDeath, colour = factor(scenario), group = factor(scenario))) + 
  geom_line() +
  labs(title = "Looking at expected age of death over parental age per scenario",
       y = "average expected age at death") +
  theme_cowplot()
#+  ylim(4.8, 7)

# plot using gam from ggplot
ggplot(myLongitudinalData, aes(ageOfParent, expectedAgeAtDeath)) + 
  geom_smooth(method = "gam") + 
  theme_cowplot() + 
  labs(title = "Looking at offspring's expected age at death over parental age",
       subtitle = "Using gam from ggplot",
       x = "parental age", 
       y = "expected age at death offspring")

# perform statistical analysis 
d <- myLongitudinalData

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(ageOfParent)))
table(z$na) 

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation
d <- d %>% filter(na>1)
d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If expectedAgeAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,expectedAgeAtDeath)) 
d <- d %>% mutate(y2 = y1/40)

d <- d %>% filter(survivalProb > 0) 
## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

# remove some data otherwise it takes too long 
d2 <- d %>% filter(na>8)

Sys.time()
m1z <- bam(y3 ~ s(ageOfParent, k = 10) + s(ageOfParent, ID, bs = "fs", k = 10),
          # + s(ageOfParent, na, k = 10) + 
          #   ti(ageOfParent, na, k = c(10,10)),
           #family=betar(link="logit"),
           data=d2,
           method = "REML")
# interaction term. 
Sys.time()

gam.check(m1z)

d2 <- d2 %>% mutate(ID = factor(ID)) 

# get maximum parental age 
maxi <- max(d2$ageOfParent)
max_iterations <- maxi
index <- 0
while ((nrow(d2[d2$ageOfParent == maxi,]) / nrow(d2) * 100) <= 3 && index < max_iterations) {
  maxi <- maxi - 1 
  index = index + 1
  if (index == max_iterations) { 
    print("no age class has enough data.. Setting to the maximum age class")
    maxi = max_iterations
  }
}
pred_data = expand.grid(ageOfParent = c(0,maxi), ID = levels(d2$ID)[1])
pred_data = expand.grid(ageOfParent = c(0,max(d2$ageOfParent)), ID = levels(d2$ID)[1])

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

logist <- function(x) 40/(1+exp(-x))

plot(m1z, trans = logist, 
     xlab = "Parental age",
     ylab = "Expected age at death of offspring",
     main = paste("Expected age at death over parental age. 
                  ", plotText, "% decrease = ", perc_decr),
     # pages = 1,
     shade = TRUE, 
     shade.col = "lightgrey",
     shift = coef(m1z)[1], # adjust for intercept,
     #ylim = c(0, 5)
)
#

###############################################################################
# combining the mechanisms in one plot to be able to compare. 
###############################################################################
# create empty dataframe 
dataTotal <- c()

# first get the scenarios-only data combined 
deathPop <- read.table(paste(path, "null/outputDeclineGameteQuality3.txt", sep = "")) # Null-only (1)
txt <- "Null"
deathPop <- read.table(paste(path, "nullDamage/outputDeclineGameteQuality.txt", sep = "")) # Null + damage (2)
txt <- "Null_damage"
deathPop <- read.table(paste(path, "nullQuality/outputDeclineGameteQuality.txt", sep = "")) # Null + quality (3)
txt <- "Null_quality"
deathPop <- read.table(paste(path, "nullResource/outputDeclineGameteQuality.txt", sep = "")) # Null + resource (4)
txt <- "Null_resource"
deathPop <- read.table(paste(path, "damageQuality/outputDeclineGameteQuality.txt", sep = "")) # quality + damage (5)
txt <- "Quality_damage"
deathPop <- read.table(paste(path, "damageResource/outputDeclineGameteQuality.txt", sep = "")) # damage + resource (6)
txt <- "Damage_resource"
deathPop <- read.table(paste(path, "qualityResource/outputDeclineGameteQuality.txt", sep = "")) # quality + resource (7)
txt <- "Quality_resource"
deathPop <- read.table(paste(path, "damageOnly/outputDeclineGameteQuality.txt", sep = "")) # damage only (8)
txt <- "damageOnly"
deathPop <- read.table(paste(path, "qualityOnly/outputDeclineGameteQuality.txt", sep = "")) # quality only (9)
txt <- "qualityOnly"
deathPop <- read.table(paste(path, "resourceOnly/outputDeclineGameteQuality.txt", sep = "")) # resource only (10)
txt <- "ResourceOnly"

# subset to only get some timepoints
deathPop <- deathPop[, 1:6] # select only the columns that are useful 
colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb")

# subset data to get info every 10.000th time step 
deathPop <- deathPop[deathPop$time %in% seq(0, max(deathPop$time), 10000),]

# add column to specify which scenario 
deathPop$scenario <- txt

# add data
dataTotal <- rbind(dataTotal, deathPop)

dataTotal$scenario <- factor(dataTotal$scenario, levels = unique(dataTotal$scenario))
# plot with boxplot 
ggplot(dataTotal, aes(scenario, ageAtDeath)) + 
  geom_boxplot(lwd = 3) + geom_jitter(alpha = 0.1, size = 0.3, width = 0.3)

# only boxplot
ggplot(dataTotal, aes(scenario, ageAtDeath)) + 
  geom_boxplot() 

# boxplot plus density plot 
ggplot(dataTotal, aes(scenario, ageAtDeath)) +
  geom_boxplot() +
  geom_flat_violin(trim = T, 
                   alpha = 0.7,
                   scale = "width") +
  labs(title = "Age at death per scenario: density and boxplot portrayed")

# boxplot plus density plot with log transformed y-axis 
dataTotalAdjusted <- dataTotal 
dataTotalAdjusted$ageAtDeath <- sqrt(dataTotalAdjusted$ageAtDeath)
ggplot(dataTotalAdjusted, aes(scenario, ageAtDeath)) +
  geom_boxplot() +
  geom_flat_violin(trim = T, 
                   alpha = 0.7,
                   scale = "width") +
  labs(title = "Age at death per scenario: density and boxplot portrayed")

test <- aggregate(dataTotal$ageAtDeath, list(dataTotal$scenario), mean)
test$sd <- tapply(dataTotal$ageAtDeath, dataTotal$scenario, sd)

ggplot(test, aes(Group.1, x)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = x-sd, ymax = x+sd)) 

# get average age at death per time step 
averageDataNull <- aggregate(deathPop$ageAtDeath, list(deathPop$time), mean)
colnames(averageDataNull) <- c("time", "meanAgeAtDeath")
# add min and max 
averageDataNull$min <- tapply(deathPop$ageAtDeath, deathPop$time, min)
averageDataNull$max <- tapply(deathPop$ageAtDeath, deathPop$time, max)
averageDataNull$scenario <- "qualityOnly"

ggplot(averageDataNull, aes(time, meanAgeAtDeath)) + 
  geom_line() +
  geom_linerange(aes(ymin = min, ymax = max), 
                 linetype = 2)

#### look at alive population ####
# create empty dataframe 
dataTotalSurvPop <- c()

# first get the scenarios-only data combined 
survPop <- read.table(paste(path, "null/outputAgeAlivePop.txt", sep = "")) # Null-only (1)
txt <- "Null"
survPop <- read.table(paste(path, "nullDamage/outputAgeAlivePop.txt", sep = "")) # Null + damage (2)
txt <- "Null_damage"
survPop <- read.table(paste(path, "nullQuality/outputAgeAlivePop.txt", sep = "")) # Null + quality (3)
txt <- "Null_quality"
survPop <- read.table(paste(path, "nullResource/outputAgeAlivePop.txt", sep = "")) # Null + resource (4)
txt <- "Null_resource"
survPop <- read.table(paste(path, "damageQuality/outputAgeAlivePop.txt", sep = "")) # quality + damage (5)
txt <- "Quality_damage"
survPop <- read.table(paste(path, "damageResource/outputAgeAlivePop.txt", sep = "")) # damage + resource (6)
txt <- "Damage_resource"
survPop <- read.table(paste(path, "qualityResource/outputAgeAlivePop.txt", sep = "")) # quality + resource (7)
txt <- "Quality_resource"
survPop <- read.table(paste(path, "damageOnly/outputAgeAlivePop.txt", sep = "")) # damage only (8)
txt <- "damageOnly"
survPop <- read.table(paste(path, "qualityOnly/outputAgeAlivePop.txt", sep = "")) # quality only (9)
txt <- "qualityOnly"
survPop <- read.table(paste(path, "resourceOnly/outputAgeAlivePop.txt", sep = "")) # resource only (10)
txt <- "ResourceOnly"

colnames(survPop) <- c("maleAge", "femaleAge")

# using square root transformation 
survPop <- sqrt(survPop) # log returned infinity becasue of age 0 

average <- mean(as.matrix(survPop))
sd <- sd(as.matrix((survPop)))
tmp <- data.frame(average, sd, txt)
dataTotalSurvPop <- rbind(dataTotalSurvPop, tmp)

dataTotalSurvPop$txt <- factor(dataTotalSurvPop$txt, levels = unique(dataTotalSurvPop$txt))

ggplot(dataTotalSurvPop, aes(txt, average)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = average-sd, ymax = average+sd)) +
  labs(y = "average age",
       x = "scenarios", 
       title = "Mean and sd portrayed of the average age of the alive population",
       subtitle = "with the ages normalized by square root transformation")


## look at age-specific gene values for sampled IDs. 
age_specific_data <- read.table(paste(path, "longrunNull/outputWithAgeSpecificGenes.txt", sep = ""))
age_specific_data <- read.table(paste(path, "qualityOnly/outputWithAgeSpecificGenes.txt", sep = ""))
age_specific_data <- read.table(paste(path, "null/outputWithAgeSpecificGenes.txt", sep = "")) # pop size = 10.000, tEnd = 10.000

age_specific_data <- age_specific_data[, 1:3]
colnames(age_specific_data) <- c("ID", "age", "geneValue")
sub <- age_specific_data[age_specific_data$ID %in% sample(age_specific_data$ID, 15),]

ggplot(sub, aes(age, geneValue, group = factor(ID), colour = factor(ID))) + 
  geom_line() +
  labs(title = "looking at the age-specific gene values for 15 sampled IDs",
       y = "ag-specific gene value")

       