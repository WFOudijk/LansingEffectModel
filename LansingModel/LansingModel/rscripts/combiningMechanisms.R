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
deathPop <- read.table(paste(path, "qualityOnly/outputDeclineGameteQuality2.txt", sep = "")) # quality only. tEnd = 10.000
plotText <- "Quality-only"
deathPop <- read.table(paste(path, "resourceOnly/outputDeclineGameteQuality.txt", sep = "")) # resource only. tEnd = 100.000
plotText <- "Resource-only"
deathPop <- read.table(paste(path, "longrunNull/outputDeclineGameteQuality.txt", sep = "")) # Null long run, tEnd = 500.000  
plotText <- "Null, long run"
library(data.table)
deathPop <- fread(paste(path, "longrunNull/outputDeclineGameteQuality2.txt", sep = ""), select = 1:6) # Null long run, tEnd = 1.000.000  
plotText <- "Null, long run" 
deathPop <- fread(paste(path, "longrunNull/outputDeclineGameteQuality3.txt", sep = ""), select = 1:6) # Null long run, tEnd = 50.000 and mut rate = 0.02
plotText <- "Null, long run" 
deathPop <- fread(paste(path, "longrunNull/outputDeclineGameteQuality4.txt", sep = ""), select = 1:6) # Null long run, tEnd = 50.000 and mut rate = 0.001; bias = -0.02 and sd = 0.015
plotText <- "Null, long run" 
deathPop <- fread(paste(path, "null/outputDeclineGameteQuality6.txt", sep = ""), select = 1:6) # Null long run, tEnd = 70.000 and mut rate = 0.002; bias = -0.02 and sd = 0.015
plotText <- "Null, long run"
deathPop <- fread(paste(path, "null/outputDeclineGameteQuality7.txt", sep = ""), select = 1:6) # Null long run, tEnd = 50.000 and mut rate = 0.02; bias = -0.02 and sd = 0.015
plotText <- "Null, long run"
deathPop <- fread(paste(path, "null/outputDeclineGameteQuality8.txt", sep = ""), select = 1:6) # Null long run, tEnd = 70.000 and mut rate = 0.02; bias = -0.03 and sd = 0.015
plotText <- "Null, long run"
deathPop <- fread(paste(path, "nullDamageQuality/outputDeclineGameteQuality.txt", sep = ""), select = 1:6) # Null + damage + quality 
plotText <- "Null+damage+quality"

deathPop <- deathPop[, 1:6] # select only the columns that are useful 
colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb")

# subset data to get info every 100th time step 
deathPop <- deathPop[deathPop$time %in% seq(0, max(deathPop$time), 100),]

# get average age at death per time step 
averageData <- aggregate(deathPop$ageAtDeath, list(deathPop$time), mean)
colnames(averageData) <- c("time", "meanAgeAtDeath")

# get average survival probability per time step 
averageDataSurvProb <- aggregate(deathPop$survivalProb, list(deathPop$time), mean)
colnames(averageDataSurvProb) <- c("time", "meanSurvProb")

# plot the average age at death over time
ggplot(averageData, aes(time, meanAgeAtDeath)) + 
  geom_line(alpha = 0.2) + 
  geom_point() + 
  theme_cowplot() + 
  labs(title = "Evolutionary consequences",
       subtitle = plotText,
       y = "average age at death") 

# plot the average survival prob over time 
ggplot(averageDataSurvProb, aes(time, meanSurvProb)) + 
  geom_line(alpha = 0.2) + 
  geom_point() + 
  theme_cowplot() + 
  labs(title = "Survival prob over time",
       subtitle = plotText,
       y = "average survival prob") 

# geom_smooth of survival prob over time 
ggplot(deathPop, aes(time, survivalProb)) + 
  geom_smooth() + 
  theme_cowplot() + 
  labs(title = "Survival prob over time",
       subtitle = plotText,
       y = "average survival prob") 

# use geom_smooth for the whole data to look at age at death
ggplot(deathPop, aes(time, ageAtDeath)) + 
  geom_smooth() + 
 # theme_cowplot() +
  labs(title = "Evolutionary consequences",
       subtitle = paste("Using ggplot gam model to smoothen. ", plotText),
       y = "age at death") + 
  coord_cartesian(ylim = c(0, 12))

# 1b. Check how the ages at death are distributed 

# only select later time point 
deathPop <- subset(deathPop, deathPop$time > 9000)

# select just the sex and age at death 
deathPop <- subset(deathPop, select = c(ageAtDeath, sex)) 
deathPop <- melt(deathPop)

ggplot(deathPop, aes(x = value, fill = sex)) +
  geom_histogram(alpha=.25, position = "identity", binwidth = 1) +
  labs(title = "Histogram of age at death distribution",
       subtitle = paste("At time > 9000.", plotText), 
       x = "Age at death") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))

# 2. ageing plots. Determine average survival probability to see if ageing occurs.
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

  # 3. Offspring lifespan over parental age 
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
plotText <- "Null + resource"
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
myLongitudinalData <- read.table(paste(path, "qualityOnly/outputLETrackedIndividuals2.txt", sep = "")) #quality - only tEnd =10.000 with new implementation of lifespan simulation
plotText <- "Quality-only"
myLongitudinalData <- read.table(paste(path, "damageResource/outputLETrackedIndividuals.txt", sep = "")) # damage + resource
plotText <- "Damage + resource"
myLongitudinalData <- read.table(paste(path, "resourceOnly/outputLETrackedIndividuals.txt", sep = "")) #  resource - only
plotText <- "Resource-only"
myLongitudinalData <- read.table(paste(path, "longrunNull/outputLETrackedIndividuals.txt", sep = "")) #  null, tEnd = 500.000
plotText <- "Long run of null scenario, tEnd = 500.000"
myLongitudinalData <- read.table(paste(path, "longrunNull/outputLETrackedIndividuals3.txt", sep = "")) #  null, tEnd = 50.000; mut rate = 0.02
plotText <- "Long run of null scenario, tEnd = 50.000"
myLongitudinalData <- read.table(paste(path, "damageOnly/outputLETrackedIndividuals2.txt", sep = "")) #  damage-only again 
plotText <- "Damage only"
myLongitudinalData <- read.table(paste(path, "nullDamageQuality/outputLETrackedIndividuals.txt", sep = "")) #  null + damage + quality 
plotText <- "Null + damage + quality"
myLongitudinalData <- read.table(paste(path, "nullDamageQuality/outputLETrackedIndividuals.txt", sep = "")) #  null + damage + quality 
plotText <- "Null + damage + quality"
myLongitudinalData <- read.table(paste(path, "nullDamageResource/outputLETrackedIndividuals.txt", sep = "")) #  null + damage + resource 
plotText <- "Null + damage + resource"

myLongitudinalData <- myLongitudinalData[, 1:5] # subset only the useful columns
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

# 95th percentile of age of parent 
percentile <- quantile(myLongitudinalData$ageOfParent, probs = 0.95) 

# calculate average expected age at death from parental age of 0 until 95th percentile parental age 
avg <- aggregate(myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$expectedAgeAtDeath, 
                 list(myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$ageOfParent), mean)
colnames(avg) <- c("ageOfParent", "averageExpectedAgeAtDeath")

avg$n <- tapply(myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$expectedAgeAtDeath,
                myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$ageOfParent, length)

avg$sd <- tapply(myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$expectedAgeAtDeath,
                 myLongitudinalData[myLongitudinalData$ageOfParent <= percentile,]$ageOfParent, sd)
calcMargin = function(x, output) { # function to calculate the margin for confidence interval 
  n <- x[3] # get sample size
  sd <- x[4] # get sd 
  return(qt(0.975, df = n - 1)*sd/sqrt(n))
}
avg$margin <- apply(avg, 1, calcMargin) # calc margin per row
avg$lowerInterval <- apply(avg, 1, function(x){x[2] - x[5]}) # get lower bound per row
avg$upperInterval <- apply(avg, 1, function(x){x[2] + x[5]}) # get upper bound per row
avg$scenario <- plotText
# adjust values to between 0-1
avg$ageOfParent <- avg$ageOfParent / percentile
avg$sd <- avg$sd / avg$averageExpectedAgeAtDeath[1]
avg$lowerInterval <- avg$lowerInterval / avg$averageExpectedAgeAtDeath[1]
avg$upperInterval <- avg$upperInterval / avg$averageExpectedAgeAtDeath[1]
avg$averageExpectedAgeAtDeath <- avg$averageExpectedAgeAtDeath / avg$averageExpectedAgeAtDeath[1]
dataTotalLansing <- rbind(dataTotalLansing, avg)

# plot 
dataTotalLansing$scenario <- factor(dataTotalLansing$scenario)
ggplot(dataTotalLansing, aes(ageOfParent, offspringLifespan, colour = scenario, group = scenario, shape = scenario)) + 
  scale_shape_manual(values=1:nlevels(dataTotalLansing$scenario)) +
  geom_point() +
  geom_line() +
  #geom_ribbon(aes(ymin = lowerInterval, ymax = upperInterval), # confidence interval
  #            alpha = 0.2,
  #            size = 0) +
  labs(title = "Looking at Lansing effect per scenario",
       y = "adjusted average expected age at death",
       x = "parental age relative to the 95th percentile parental age class") +
  theme_cowplot() 

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
d2 <- d %>% filter(na>5)

Sys.time()
m1z <- bam(y3 ~ s(ageOfParent, k = 10) + s(ageOfParent, ID, bs = "fs", k = 10),
          # + s(ageOfParent, na, k = 10) + # can be used to check how na plays a part in the data
          #   ti(ageOfParent, na, k = c(10,10)),
           data=d2,
           method = "REML")
# interaction term. 
Sys.time()

gam.check(m1z)

d2 <- d2 %>% mutate(ID = factor(ID)) 

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
getwd()
#dataTotal <- read.table("desktop/dataTotal.txt", header = T)

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

# plot with boxplot and jitter
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

# create density plots with ridges
library(ggridges)
ggplot(dataTotal, aes(ageAtDeath, scenario, group = scenario)) + 
  geom_density_ridges_gradient(
    point_shape = "|",
    alpha = 0.2,
    scale = 2, 
    quantile_lines = TRUE, 
    vline_color = c("red"),
    quantile_fun = median,
    bandwidth = 0.4,
    bandwidth=1.5
    ) +
  geom_density_ridges_gradient(
    quantile_lines = TRUE, 
    vline_color = c("blue"),
    fill = NA, 
    alpha = 0.2,
    scale = 2, 
    bandwidth = 0.4,
    quantiles = c(0.025, 0.975),
    bandwidth=1.5
  ) +
  labs(title = "Distribution of age at death per scenario",
       subtitle = "Blue lines represent the lower and upper quantile and red line 
       represents the median") +
  xlim(0, 10)

# using ridges but the quantiles are coloured. 
ggplot(dataTotal, aes(ageAtDeath, scenario, group = scenario, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = T, 
    quantiles = 4, quantile_lines = T, bandwidth = 0.9) +
  xlim(0, 40)
  
# use a dot for the mean and portray error bars instead of boxplot 
df_new <- aggregate(dataTotal$ageAtDeath, list(dataTotal$scenario), mean)
df_new$sd <- tapply(dataTotal$ageAtDeath, dataTotal$scenario, sd)
colnames(df_new) <- c("ageAtDeath", "mean", "sd")

# plot the data 
ggplot(df_new, aes(ageAtDeath, mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd)) 

# Look at the data again, but first normalize the data with log 
adjustedDataTotal <- dataTotal
adjustedDataTotal$ageAtDeath <- adjustedDataTotal$ageAtDeath + 0.1 # add a constant to prevent log(0) becoming -Inf. 
adjustedDataTotal$ageAtDeath <- log(adjustedDataTotal$ageAtDeath)
# get mean 
df_new <- aggregate(adjustedDataTotal$ageAtDeath, list(adjustedDataTotal$scenario), mean)
colnames(df_new) <- c("scenario", "mean")
df_new$sd <- tapply(adjustedDataTotal$ageAtDeath, adjustedDataTotal$scenario, sd)
# plot the log-transformed data 
ggplot(df_new, aes(scenario, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd)) +
  labs(title = "Log transformed age at death data per scenario",
       subtitle = "the dots are the average and the errorbars portray the sd",
       y = "Age at death")

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
survPop <- sqrt(survPop) # log returned infinity because of age 0 

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
age_specific_data <- read.table(paste(path, "longrunNull/outputWithAgeSpecificGenes3.txt", sep = "")) # pop size = 10.000, tEnd = 50.000, mut rate = 0.02
age_specific_data <- read.table(paste(path, "null/outputWithAgeSpecificGenes7.txt", sep = "")) # Null long run, tEnd = 50.000 and mut rate = 0.02; bias = -0.02 and sd = 0.015

age_specific_data <- age_specific_data[, 1:3]
colnames(age_specific_data) <- c("ID", "age", "geneValue")
# sample a subset 
sub <- age_specific_data[age_specific_data$ID %in% sample(age_specific_data$ID, 15),]

# plot the subset. 
ggplot(sub, aes(age, geneValue, group = factor(ID), colour = factor(ID))) + 
  geom_line() +
  labs(title = "looking at the age-specific gene values for 15 sampled IDs",
       y = "ag-specific gene value") 


#######
# Offspring lifespan recalculated with simulated offspring going through 
# the mortality round an x number of times. 
# Starting from this point this is calculated by using a for loop. 
#######
dataTotalLansing <- c()
dataTotalLansing2 <- c()
dataTotalLansingPat <- c()

calcMargin = function(x, output) { # function to calculate the margin for confidence interval 
  n <- x[3] # get sample size
  sd <- x[4] # get sd 
  return(qt(0.975, df = n - 1)*sd/sqrt(n))
}

# set the path to all data files. 
parent_path <- paste(path, "combiningAll2/", sep = "")
f <- list.files(path = parent_path, pattern = "outputLifeExp", recursive = T)
dataTotalLansing <- c() # for maternal 
dataTotalLansingPat <- c() # for paternal 
totalDatasets <- c()

for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 1
  scenario <- splitted_path[[1]][nScenario]
  # read data
  local_data <- read.csv(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  # add the complete data to the totalDatasets dataframe. 
  totalDatasets <- rbind(totalDatasets, local_data) # get all data together
  
  # get data for Lansing plot 
  percentile <- quantile(local_data$maternalAge, probs = 0.95) 
  percentileP <- quantile(local_data$paternalAge, probs = 0.95) 
  
  # calculate average age at death from parental age of 0 until 95th percentile parental age 
  avg <- aggregate(local_data[local_data$maternalAge <= percentile,]$ageAtDeath, 
                   list(local_data[local_data$maternalAge <= percentile,]$maternalAge), mean)
  avgP <- aggregate(local_data[local_data$paternalAge <= percentile,]$ageAtDeath, 
                   list(local_data[local_data$paternalAge <= percentile,]$paternalAge), mean)
  colnames(avg) <- c("ageOfParent", "offspringLifespan")
  colnames(avgP) <- c("ageOfParent", "offspringLifespan")
  
  # calculate the confidence interval 
  avg$n <- tapply(local_data[local_data$maternalAge <= percentile,]$ageAtDeath,
                  local_data[local_data$maternalAge <= percentile,]$maternalAge, length)
  avgP$n <- tapply(local_data[local_data$paternalAge <= percentile,]$ageAtDeath,
                  local_data[local_data$paternalAge <= percentile,]$paternalAge, length)
  
  avg$sd <- tapply(local_data[local_data$maternalAge <= percentile,]$ageAtDeath,
                   local_data[local_data$maternalAge <= percentile,]$maternalAge, sd)
  avgP$sd <- tapply(local_data[local_data$paternalAge <= percentile,]$ageAtDeath,
                   local_data[local_data$paternalAge <= percentile,]$paternalAge, sd)

  avg$margin <- apply(avg, 1, calcMargin) # calc margin per row
  avg$lowerInterval <- apply(avg, 1, function(x){x[2] - x[5]}) # get lower bound per row
  avg$upperInterval <- apply(avg, 1, function(x){x[2] + x[5]}) # get upper bound per row
  avg$scenario <- scenario
  # paternal
  avgP$margin <- apply(avgP, 1, calcMargin) # calc margin per row
  avgP$lowerInterval <- apply(avgP, 1, function(x){x[2] - x[5]}) # get lower bound per row
  avgP$upperInterval <- apply(avgP, 1, function(x){x[2] + x[5]}) # get upper bound per row
  avgP$scenario <- scenario
  
  # adjust values to between 0-1
  avg$ageOfParent <- avg$ageOfParent / percentile
  avg$sd <- avg$sd / avg$offspringLifespan[1]
  avg$lowerInterval <- avg$lowerInterval / avg$offspringLifespan[1]
  avg$upperInterval <- avg$upperInterval / avg$offspringLifespan[1]
  avg$offspringLifespan <- avg$offspringLifespan / avg$offspringLifespan[1]
  # paternal
  avgP$ageOfParent <- avgP$ageOfParent / percentile
  avgP$sd <- avgP$sd / avgP$offspringLifespan[1]
  avgP$lowerInterval <- avgP$lowerInterval / avgP$offspringLifespan[1]
  avgP$upperInterval <- avgP$upperInterval / avgP$offspringLifespan[1]
  avgP$offspringLifespan <- avgP$offspringLifespan / avgP$offspringLifespan[1]
  # combine
  dataTotalLansing <- rbind(dataTotalLansing, avg)
  dataTotalLansingPat <- rbind(dataTotalLansingPat, avgP)
}

# plot maternal Lansing effects 
dataTotalLansing$scenario <- factor(dataTotalLansing$scenario)
ggplot(dataTotalLansing, aes(ageOfParent, offspringLifespan, colour = scenario, group = scenario, shape = scenario)) + 
  scale_shape_manual(values=1:nlevels(dataTotalLansing$scenario)) +
  geom_point() +
  geom_line() +
  #geom_ribbon(aes(ymin = lowerInterval, ymax = upperInterval), # use this for confidence interval 
  #            alpha = 0.2,
  #            size = 0) +
  labs(title = "Looking at Lansing effect per scenario",
       subtitle = "maternal",
       y = "adjusted average expected age at death",
       x = "parental age relative to the 95th percentile parental age class") +
  theme_cowplot() +
  ylim(0,1.1)

# plot paternal Lansing effects 
dataTotalLansingPat$scenario <- factor(dataTotalLansingPat$scenario)
ggplot(dataTotalLansingPat, aes(ageOfParent, offspringLifespan, colour = scenario, group = scenario, shape = scenario)) + 
  scale_shape_manual(values=1:nlevels(dataTotalLansingPat$scenario)) +
  geom_point() +
  geom_line() +
  #geom_ribbon(aes(ymin = lowerInterval, ymax = upperInterval), # use this for confidence interval. 
  #            alpha = 0.2,
  #            size = 0) +
  labs(title = "Looking at Lansing effect per scenario",
       subtitle = "paternal",
       y = "adjusted average expected age at death",
       x = "parental age relative to the 95th percentile parental age class") +
  theme_cowplot() +
  ylim(0,1.1)

###############

# perform statistical analysis 
d <- myLongitudinalData
d <- totalDatasets

length(unique(d$ID)) # number of parents = 981. The remainder did not have offspring 

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(maternalAge)))
table(z$na) # indeed there are plenty with more than 1 ageOfParent!

## Add the counter to d
d <- d %>% left_join(z,by="ID")

## Filter out the one with 1 observation (but this selects ...)

d <- d %>% filter(na>1)

d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
d <- d %>% mutate(y2 = y1/40)

## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))
d2 <- d[d$ID %in% sample(d$ID, 200),]

#d$scenario <- factor(d$scenario)


#d2 <- d %>% filter(na>8)


Sys.time()
m1z <- bam(y3 ~ s(maternalAge, k = 10, by = as.factor(scenario)), 
           data=d,
           method = "REML")

m1z <- bam(y3 ~ s(maternalAge, k = 10) + s(maternalAge, ID, bs = "fs", k = 10), 
           data=d2,
           method = "REML")


# interaction term. 
Sys.time()

# plot 
plot(m1z)

# not fully working yet. 
pred_data = expand.grid(maternalAge = c(0,max(d$maternalAge)), 
                        scenario = levels(d$scenario))
# with ID 
d$ID <- factor(d$ID)
pred_data = expand.grid(maternalAge = c(0,max(d$maternalAge)),
                        ID = levels(d$ID)[1]) # maternal


pred_data <- aggregate(d$maternalAge, list(d$scenario), min)
pred_data <- rbind(pred_data, aggregate(d$maternalAge, list(d$scenario), max))
colnames(pred_data) <- c("scenario", "maternalAge")

## Predict
z <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
z <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(paternalAge)"))

             
str(z)
z
logist(z)

# Add intercept and all terms, transform back
#pred_data$y3 <- coef(m1z)[1] + apply(z,1,sum) 
#logist <- function(x) 1/(1+exp(-x)) # transforms data back to the original data/40 (so mapped to 0-1)
#pred_data$y4 <- 40*logist(pred_data$y3)

#ggplot(pred_data, aes(maternalAge, y4, colour = as.factor(scenario))) + 
#  geom_line() +
#  ylim(0,100)

## Use logist transformed:
logist <- function(x) 40/(1+exp(-x)) # transforms data back to original age at death 
zl <- logist(z)

## Percentage decrease: Not fully working yet. 
for(i in 1:length(unique(pred_data$scenario))){
  print(z[i]-z[(i+length(unique(pred_data$scenario)))])
  perc_decr <- 100*(z[i]-z[(i+length(unique(pred_data$scenario)))])/z[i]
  print(paste(unique(pred_data$scenario[i]), "perc decre =", logist(perc_decr)))
}

perc_decr <- 100*(zl[1]-zl[2])/zl[1]
perc_decr

summary(m1z)

logist <- function(x) 40/(1+exp(-x)) # transforms data back to original age at death 

plot(m1z, trans = logist, 
     #xlab = "Parental age",
     #ylab = "Expected age at death of offspring",
     #main = paste("Expected age at death over parental age. 
    #              ", plotText, "% decrease = ", perc_decr),
     # pages = 1,
     shade = TRUE, 
     shade.col = "lightgrey",
     shift = coef(m1z)[1], # adjust for intercept,
     #ylim = c(0, 5)
)
#

###########################################
# Looking at the individuals longitudinally
##########################################

parent_path <- paste(path, "combiningAll3/", sep = "")
dataTotalLansing <- c()
f <- list.files(path = parent_path, pattern = "outputLifeExpLong", recursive = T)

# paste all data together
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 1
  scenario <- splitted_path[[1]][nScenario]
  # read data
  local_data <- read.csv(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  # subset parents 
  z <- local_data %>% group_by(ID) %>% summarize(na = length(unique(maternalAge)))
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  local_data <- local_data %>% filter(na > 6)
  # sample 200 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(local_data$ID, 200, replace = F),]
  d <-local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d2 %>% mutate(y2 = y1/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  print(d2$scenario[1])
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(local_data$maternalAge))) 
    }
  m1z <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = 4), 
             data=d2,
             method = "REML")
  # get 95th percentile 
  percentile <- quantile(d2$maternalAge, probs = 0.95)
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,percentile,1),
                          ID = levels(d2$ID)[1],
                          scenario = levels(d2$scenario)[1]) 
  
  pred_data$z <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
  # predict
  #pred_data$z <- predict(m1z, newdata = pred_data)
  colnames(pred_data)[4] <- "z"
  pred_data$z <- logist(pred_data$z) # transform back
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$z <- pred_data$z / pred_data$z[1]
  
  # add the complete data to the totalDatasets dataframe. 
  dataTotalLansing <- rbind(dataTotalLansing, pred_data)
}

ggplot(dataTotalLansing, aes(maternalAge, z, group=scenario, colour = scenario, shape = scenario)) +
  geom_line() +
  scale_shape_manual(values=1:nlevels(dataTotalLansing$scenario)) +
  geom_point() 
  

