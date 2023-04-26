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
deathPop <- read.table(paste(path, "nullQuality/outputDeclineGameteQuality.txt", sep = "")) # Null + quality 
plotText <- "Null + quality"
deathPop <- read.table(paste(path, "nullResource/outputDeclineGameteQuality.txt", sep = "")) # Null + resource 
plotText <- "Null + resource distribution"

deathPop <- deathPop[, 1:6] # select only the columns that are useful 
colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb")

# subset data to get info every 100th time step 
deathPop <- deathPop[deathPop$time %in% seq(0, max(deathPop$time), 100),]

# get average age at death per time step 
averageData <- aggregate(deathPop$ageAtDeath, list(deathPop$time), mean)
colnames(averageData) <- c("time", "meanAgeAtDeath")

# plot the average 
ggplot(averageData, aes(time, meanAgeAtDeath)) + 
  geom_line(alpha = 0.2) + 
  geom_point() + 
  theme_cowplot() + 
  labs(title = "Evolutionary consequences",
       subtitle = plotText,
       y = "average age at death") 
 # ylim(0, 15)

# use geom_smooth for the whole data 
ggplot(deathPop, aes(time, ageAtDeath)) + 
  geom_smooth() + 
  theme_cowplot() +
  labs(title = "Evolutionary consequences",
       subtitle = paste("Using ggplot gam model to smoothen. ", plotText),
       y = "age at death") 
 # coord_cartesian(ylim = c(0, 15))

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
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals.txt", sep = "")) # mut rate = 0.001; mean = -0.02; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.02; sd = 0.01"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals2.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals3.txt", sep = "")) # mut rate = 0.001; mean = -0.01; sd = 0.01; pop size = 10.000 and tEnd = 100.000
plotText <- "Mut rate = 0.001; mean = -0.01; sd = 0.01"
myLongitudinalData <- read.table(paste(path, "null/outputLETrackedIndividuals4.txt", sep = "")) # null; pop size = 10.000 and tEnd = 10.000
plotText <- "Null. tEnd = 10.000; pop size = 10.000"
myLongitudinalData <- read.table(paste(path, "nullDamage/outputLETrackedIndividuals.txt", sep = "")) # null + damage
plotText <- "Null + damage"
myLongitudinalData <- read.table(paste(path, "nullQuality/outputLETrackedIndividuals.txt", sep = "")) # null + quality
plotText <- "Null + quality"
myLongitudinalData <- read.table(paste(path, "nullResource/outputLETrackedIndividuals.txt", sep = "")) # null + resource distribution
plotText <- "Null + resource distribution"

myLongitudinalData <- myLongitudinalData[, 1:5] # subset only the useful columns
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

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

## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))

# remove some data otherwise it takes too long 
d2 <- d %>% filter(na>6)

Sys.time()
m1z <- bam(y3 ~ s(ageOfParent, k = 10) + s(ageOfParent, ID, bs = "fs", k = 10),
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
     shift = coef(m1z)[1] # adjust for intercept,
)
#
