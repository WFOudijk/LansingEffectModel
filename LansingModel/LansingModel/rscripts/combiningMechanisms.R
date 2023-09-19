###############################################################################
# COMBINING OF THE MECHANISMS 
# final Rscript
###############################################################################
library(cowplot)
library(tidyverse)
library(mgcv)
library(MetBrewer)
library(ggpubr)
library(grid)
library(gridExtra)
library(R.filesets)
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/combiningAll/"
path_to_output <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/report/figures/"

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

# do the following steps for every death pop data set 
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
  labs(title = "Survival prob over time with geom_smooth",
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

# plot the ages at death 
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

# read in the data 
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

# every time for the dataset
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
avg$scenario <- plotText # incldue the corresponding scenario to the data 
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

# plot the model. 
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

# first get the data combined 
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

###############################################################################
# look at the alive population. 
###############################################################################
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

# plot 
ggplot(dataTotalSurvPop, aes(txt, average)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = average-sd, ymax = average+sd)) +
  labs(y = "average age",
       x = "scenarios", 
       title = "Mean and sd portrayed of the average age of the alive population",
       subtitle = "with the ages normalized by square root transformation")


###############################################################################
# Sample some IDs and look at the age-specific gene values. 
###############################################################################
age_specific_data <- read.table(paste(path, "longrunNull/outputWithAgeSpecificGenes.txt", sep = ""))
age_specific_data <- read.table(paste(path, "qualityOnly/outputWithAgeSpecificGenes.txt", sep = ""))
age_specific_data <- read.table(paste(path, "null/outputWithAgeSpecificGenes.txt", sep = "")) # pop size = 10.000, tEnd = 10.000
age_specific_data <- read.table(paste(path, "longrunNull/outputWithAgeSpecificGenes3.txt", sep = "")) # pop size = 10.000, tEnd = 50.000, mut rate = 0.02
age_specific_data <- read.table(paste(path, "null/outputWithAgeSpecificGenes7.txt", sep = "")) # Null long run, tEnd = 50.000 and mut rate = 0.02; bias = -0.02 and sd = 0.015

age_specific_data <- age_specific_data[, 1:3]
colnames(age_specific_data) <- c("ID", "age", "geneValue")
# sample a subset 
sub <- age_specific_data[age_specific_data$ID %in% sample(unique(age_specific_data$ID), 15),]

# plot the subset. 
ggplot(sub, aes(age, geneValue, group = factor(ID), colour = factor(ID))) + 
  geom_line() +
  labs(title = "looking at the age-specific gene values for 15 sampled IDs",
       y = "ag-specific gene value") 


#############################################################################
# Offspring lifespan recalculated with simulated offspring going through 
# the mortality round an x number of times. 
# Starting from this point this is calculated by using a for loop. 
# CROSS-SECTIONAL
#############################################################################
dataTotalLansingLat <- c()
dataTotalLansing2 <- c()
dataTotalLansingPat <- c()

calcMargin = function(x, output) { # function to calculate the margin for confidence interval 
  n <- x[3] # get sample size
  sd <- x[4] # get sd 
  return(qt(0.975, df = n - 1)*sd/sqrt(n))
}

# set the path to all data files. 
parent_path <- paste(path, "combiningAll2/", sep = "")
f <- list.files(path = parent_path, pattern = "outputLifeExp.txt", recursive = T)
dataTotalLansingMat <- c() # for maternal 
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
  dataTotalLansingMat <- rbind(dataTotalLansingMat, avg)
  dataTotalLansingPat <- rbind(dataTotalLansingPat, avgP)
}

# plot maternal Lansing effects 
dataTotalLansingMat$scenario <- factor(dataTotalLansingMat$scenario)
ggplot(dataTotalLansingMat, aes(ageOfParent, offspringLifespan, colour = scenario, group = scenario, shape = scenario)) + 
  scale_shape_manual(values=1:nlevels(dataTotalLansingMat$scenario)) +
  geom_point() +
  scale_color_manual(values = met.brewer("Greek", 15)) +
  geom_line() +
  #geom_ribbon(aes(ymin = lowerInterval, ymax = upperInterval), # use this for confidence interval 
  #            alpha = 0.2,
  #            size = 0) +
  labs(title = "Looking at Lansing effect per scenario cross-sectionally",
       subtitle = "maternal",
       y = "adjusted average expected age at death",
       x = "Adjusted parental age") +
  theme_cowplot() +
  ylim(0,1.2)

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

length(unique(d$ID)) 

## Count how many unique ageOfParent per ID
z <- d %>% group_by(ID) %>% summarize(na = length(unique(maternalAge)))
table(z$na)

## Add the counter to d
d <- d %>% left_join(z,by="ID")

d <- d %>% filter(na>1)

d <- d %>% mutate(ID = factor(ID)) 

# compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
d <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
d <- d %>% mutate(y2 = y1/40)

## Work with logits (then (0,1) -> (-inf,+inf))
d <- d %>% mutate(y3 = car::logit(y2))
# sample 200 parents. 
d2 <- d[d$ID %in% sample(unique(d$ID), 200),]

#d$scenario <- factor(d$scenario)
#d2 <- d %>% filter(na>8)

Sys.time()
# for the combined dataset with scenario as factor. 
m1z <- bam(y3 ~ s(maternalAge, k = 10, by = as.factor(scenario)), 
           data=d,
           method = "REML")

# for one dataset. 
m1z <- bam(y3 ~ s(maternalAge, k = 10) + s(maternalAge, ID, bs = "fs", k = 10), 
           data=d2,
           method = "REML")
Sys.time()

## Use logist transformed:
logist <- function(x) 40/(1+exp(-x)) # transforms data back to original age at death 

summary(m1z)

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

############################################################################
# Looking at the individuals LONGITUDINALLY
###########################################################################
# set path 
parent_path <- paste(path, "combiningAll3/", sep = "")
dataTotalLansing <- c()
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T)

# paste all data together
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 2
  scenario <- splitted_path[[1]][nScenario]
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.csv(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  # add replicate as column 
  local_data$rep <- rep
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  # subset parents 
  z <- local_data %>% group_by(ID) %>% summarize(na = length(unique(maternalAge)))
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  local_data <- local_data %>% filter(na > 6)
  # sample 200 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 200, replace = F),]
  d <-local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d2 %>% mutate(y2 = y1/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
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
                          scenario = levels(d2$scenario)[1],
                          rep = levels(d2$rep)[1]) 
  
  tmp <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(maternalAge)"), se.fit = TRUE)
 
  pred_data$z <- tmp$fit
  pred_data$upr <- tmp$fit + (2 * tmp$se.fit)
  pred_data$lwr <- tmp$fit - (2 * tmp$se.fit)
  
  pred_data$z <- coef(m1z)[1] + pred_data$z 
  pred_data$upr <- coef(m1z)[1] + pred_data$upr 
  pred_data$lwr <- coef(m1z)[1] + pred_data$lwr 
  
  pred_data$z <- logist(pred_data$z) # transform back
  pred_data$upr <- logist(pred_data$upr) # transform upper interval back 
  pred_data$lwr <- logist(pred_data$lwr) # transform lower interval back 
  
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$upr <- pred_data$upr / pred_data$z[1]
  pred_data$lwr <- pred_data$lwr / pred_data$z[1]
  pred_data$z <- pred_data$z / pred_data$z[1]
  
  # add the complete data to the totalDatasets dataframe. 
  dataTotalLansing <- rbind(dataTotalLansing, pred_data)
}

write.table(dataTotalLansing, paste(path_to_output, "dataTotalLansingRep3.txt", sep = ""))

# the normalized Lansing data
#dataTotalLansingTMP <- read.table("/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/combiningAll/combiningAll3/dataTotalLansing.txt")
#dataTotalLansing <- read.table(paste(path_to_output, "dataTotalLansing.txt", sep = "")) # with confidence bands 
#dataTotalLansing <- read.table(paste(path_to_output, "dataTotalLansingRep2.txt", sep = "")) # longitudinal results in rep 2.  

library(MetBrewer)
library(ggrepel)

# data plotted with scenarios at the end of the lines 
dataTotalLansing %>%
  mutate(label = if_else(maternalAge == max(maternalAge), as.character(scenario), NA_character_)) %>%
  ggplot(aes(maternalAge, z, group=scenario, colour = scenario)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Cross", 15)) +
  #geom_point() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() + 
  geom_label_repel(aes(label = label),
                   nudge_x = 0.1,
                   na.rm = TRUE) +
  theme(legend.position = "none")

# data plotted with the scenarios as legend. 
ggplot(dataTotalLansing, aes(maternalAge, z, group=scenario, colour = scenario, shape = scenario)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Greek", 15)) +
  scale_shape_manual(values=1:nlevels(dataTotalLansing$scenario)) +
  geom_point() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2)

# data plotted in a grid
ggplot(dataTotalLansing, aes(maternalAge, z)) +
  geom_line() +
  #scale_color_manual(values = met.brewer("Greek", 15)) +
  #scale_shape_manual(values=1:nlevels(dataTotalLansing$scenario)) +
  geom_point() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  facet_grid(dataTotalLansing$scenario)

# plot only the singled-out scenarios 
# get data 
singleScenarios <- c("damageOnly", "qualityOnly", "resourceOnly")
singleScenarios <- c("damage", "quality", "resource")

singleSceneariosData <- dataTotalLansing[dataTotalLansing$scenario %in% singleScenarios,]
singleSceneariosData <- dataTotalLansingMat[dataTotalLansingMat$scenario %in% singleScenarios,]

# plot 
ggplot(singleSceneariosData, aes(maternalAge, z)) +
  geom_line() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_bw() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, outline.type = "full") +
  facet_wrap(singleSceneariosData$scenario) +
  theme(text = element_text(family = "Times New Roman", size = 20),
        strip.text.x = element_text(size = 25, face = "bold"))

ggplot(singleSceneariosData, aes(ageOfParent, offspringLifespan)) +
  geom_line() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_bw() +
  geom_ribbon(aes(ymin = lowerInterval, ymax = upperInterval), alpha = 0.2, outline.type = "full") +
  facet_wrap(singleSceneariosData$scenario) +
  theme(text = element_text(family = "Times New Roman", size = 20),
        strip.text.x = element_text(size = 25, face = "bold"))

ggsave(paste(path_to_output, "singleScenarios.png", sep = ""), last_plot())

############################################################################
# CROSS-SECTIONAL
############################################################################
parent_path <- paste(path, "combiningAll2/", sep = "")
f <- list.files(path = parent_path, pattern = "outputLifeExp.txt", recursive = T)
dataTotalLansingMat <- c() # for maternal 
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
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  d <-local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d2 %>% mutate(y2 = y1/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  print(d2$scenario[1])
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(local_data$maternalAge))) 
  }
  m1z <- bam(y3 ~ s(maternalAge, k = k), 
             data=d2,
             method = "REML")
  # get 95th percentile 
  percentile <- quantile(d2$maternalAge, probs = 0.95)
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,percentile,1),
                          ID = levels(d2$ID)[1],
                          scenario = levels(d2$scenario)[1]) 
  
  tmp <- predict(m1z, newdata = pred_data, type="terms",terms = c("s(maternalAge)"), se.fit = TRUE)
  
  pred_data$z <- tmp$fit
  pred_data$upr <- tmp$fit + (2 * tmp$se.fit)
  pred_data$lwr <- tmp$fit - (2 * tmp$se.fit)
  
  pred_data$z <- coef(m1z)[1] + pred_data$z 
  pred_data$upr <- coef(m1z)[1] + pred_data$upr 
  pred_data$lwr <- coef(m1z)[1] + pred_data$lwr 
  
  pred_data$z <- logist(pred_data$z) # transform back
  pred_data$upr <- logist(pred_data$upr) # transform upper interval back 
  pred_data$lwr <- logist(pred_data$lwr) # transform lower interval back 
  
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$upr <- pred_data$upr / pred_data$z[1]
  pred_data$lwr <- pred_data$lwr / pred_data$z[1]
  pred_data$z <- pred_data$z / pred_data$z[1]
  
  # add the complete data to the totalDatasets dataframe. 
  dataTotalLansingMat <- rbind(dataTotalLansingMat, pred_data)
}

# data plotted with the scenarios as legend. 
ggplot(dataTotalLansingMat, aes(maternalAge, z, group=scenario, colour = scenario, shape = scenario)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Greek", 15)) +
  scale_shape_manual(values=1:nlevels(dataTotalLansingMat$scenario)) +
  geom_point() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2)

# plot only the singled-out scenarios 
# get data 
singleScenarios <- c("damageOnly", "qualityOnly", "resourceOnly")
singleScenarios <- c("damage", "quality", "resource")

singleSceneariosData <- dataTotalLansing[dataTotalLansing$scenario %in% singleScenarios,]
singleSceneariosData <- dataTotalLansingMat[dataTotalLansingMat$scenario %in% singleScenarios,]

# plot 
ggplot(singleSceneariosData, aes(maternalAge, z)) +
  geom_line() +
  labs(title = "Offspring quality over parental age",
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_bw() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, outline.type = "full") +
  facet_wrap(singleSceneariosData$scenario) +
  theme(text = element_text(family = "Times New Roman", size = 20),
        strip.text.x = element_text(size = 25, face = "bold")) +
  ylim(0,1.5)

###############################################################################
# Using replicates as CIs for the scenarios LONGITUDINALLY.  
###############################################################################

# LONGITUDINALLY 

allData <- loadRDS("allData.rds")
models <- loadRDS("models.rds")

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = 4), data = data, method = "REML")
  return(tmp)
}

# reading all data and saving both the data as the gam models. 
parent_path <- paste(path, "combiningAll3/", sep = "")
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T)
allData <- c()
allDataComplete <- c()
models <- list()
# paste all data together
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 2
  scenario <- splitted_path[[1]][nScenario]
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.csv(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per scenario
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 

  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  print(d2$scenario[1])
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(local_data$maternalAge))) 
  }
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and scenario
  names(models)[i] <- paste("model_", rep, "_", scenario, sep = "")
}

# save the R data just in case. 
saveRDS(allData, file = "allData.rds")
saveRDS(models, file = "models.rds")
saveRDS(allDataComplete, file = "allDataComplete.rds")

allData <- loadRDS("allData.rds")
models <- loadRDS("models.rds")

# loop through the scenarios 
allData$scenario <- factor(allData$scenario)
allData$rep <- factor(allData$rep)

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()
percentiles <- data.frame(levels(allData$scenario))
colnames(percentiles) <- c("scenario")
percentiles$percentile <- 0
# go per scenario through the data 
for (x in levels(allData$scenario)){
  # make subset of scenario 
  sub <- allData[allData$scenario == x,]
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          scenario = sub$scenario[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste("model_", i, "_", x, sep = "")
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 4:13)
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)

  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data)
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$scenario == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data)
}

saveRDS(percentiles, "percentiles.RDS")

# plot the normalized data grouped by scenario
ggplot(predDataTotalNormalized, aes(maternalAge, mean, group=scenario, colour = scenario, shape = scenario)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Greek", 15)) +
  scale_shape_manual(values=1:nlevels(pred_data$scenario)) +
  geom_point() +
  labs(x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) 

###############################################################################
# Generating the plots per scenario  
###############################################################################
# get the scenarios relevant for the matrix plot
scenarios <- c()
for (x in levels(allData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) <= 2){
    scenarios <- c(scenarios, x) # get list of single scenarios and doubles. 
  }
}

# dynamically generate the plots
plots <- c()
for (i in 1:length(scenarios)){ # per scenario ..
  # .. make the ggplot 
  p <- ggplot(predDataTotalNormalized[predDataTotalNormalized$scenario == scenarios[i],], aes(maternalAge, mean)) +
    geom_line() +
    labs(x = NULL, # set labs to NULL 
         y = NULL) +
    geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) + # CIs based on reps 
    ylim(0, 1.7) + 
    theme_minimal() +
    theme(#legend.text = element_text(size=10),
      #legend.key.size = unit(0.2, "cm"),
    #legend.key.width = unit(0.1,"cm"),
      #legend.position = "bottom",
      #legend.title = element_blank(),
      axis.text = element_text(size=11,face="plain",color="black"),
    axis.title = element_text(size = 13),
    #axis.text.x = element_blank(),
    axis.line = element_line(color="black", linewidth = 0.6),
    panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.5)
  )
  plots[[i]] <- p # save the plot to a list 
  # rename the list element by scenario 
  names(plots)[i] <- paste("p_", scenarios[i], sep = "") 
}

# make matrix plot using cowplot 
plot_grid(plots$p_damage, NULL, NULL, NULL, 
          plots$p_nullDamage, plots$p_null, NULL, NULL,
          plots$p_damageQuality, plots$p_nullQuality, plots$p_quality, NULL,
          plots$p_damageResource, plots$p_nullResource, plots$p_qualityResource, plots$p_resource,
          ncol = 4)

###############################################################################
# Using replicates as CIs for the scenarios CROSS-SECTIONAL.    
############################################################################### 

allDataLat <- loadRDS("allDataLat.rds")
modelsLat <- loadRDS("modelsLat.rds")

# the gam model
run_model_lat <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k), data=data,method = "REML")
  return(tmp)
}

# reading all data and saving both the data as the gam models. 
parent_path <- paste(path, "combiningAll3/", sep = "")
f <- list.files(path = parent_path, pattern = "outputLifeExp.txt", recursive = T)
allDataLat <- c()
modelsLat <- list()
# paste all data together
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 2
  scenario <- splitted_path[[1]][nScenario]
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.csv(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per scenario
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allDataLat <- rbind(allDataLat, local_data)
  
  # GAM 
  d <- local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  print(d2$scenario[1])
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(local_data$maternalAge))) 
  }
  # run model 
  mod <- run_model_lat(d2)
  # save model in list 
  modelsLat[[i]] <- mod
  # rename model to be unique for replicate and scenario
  names(modelsLat)[i] <- paste("model_", rep, "_", scenario, sep = "")
}

# save the R data just in case. 
saveRDS(allDataLat, file = "allDataLat.rds")
saveRDS(modelsLat, file = "modelsLat.rds")

# loop through the scenarios 
allDataLat$scenario <- factor(allDataLat$scenario)
allDataLat$rep <- factor(allDataLat$rep)

# to save all data 
predDataTotalLat <- c()
predDataTotalNormalizedLat <- c()

for (x in levels(allDataLat$scenario)){
  print(x)
  # make subset of scenario 
  sub <- allDataLat[allDataLat$scenario == x,]
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          scenario = sub$scenario[1]) 
  
  sub <- sub %>% mutate(rep = factor(rep))
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste("model_", i, "_", x, sep = "")
    index <- which(names(modelsLat) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(modelsLat[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means
  repVals <- subset(pred_data, select = 4:11)
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) # transform back
  pred_data$min <- logist(pred_data$min) # transform min interval back 
  pred_data$max <- logist(pred_data$max) # transform max interval back 
  
  predDataTotalLat <- rbind(predDataTotalLat, pred_data)
  
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  predDataTotalNormalizedLat <- rbind(predDataTotalNormalizedLat, pred_data)
}

###############################################################################
# Combining longitudinal with cross-sectional in one 4x4 matrix plot.   
############################################################################### 

# add column with group to both normalized data sets 
predDataTotalNormalized$group <- "Longitudinal"
predDataTotalNormalizedLat$group <- "Cross-sectional"
# merge them together 
totalNormalizedData <- rbind(predDataTotalNormalized, predDataTotalNormalizedLat)

#totalNormalizedData <- read.table("totalNormalizedData.txt")
totalNormalizedData <- totalNormalizedData %>% mutate(scenario = factor(scenario))

# get the scenarios relevant for the matrix plot
scenarios <- c()
for (x in levels(totalNormalizedData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) <= 2){
    scenarios <- c(scenarios, x) # get list of single scenarios and doubles. 
  }
}

# dynamically generate the plots
plotsTot <- c()
for (i in 1:length(scenarios)){
  p <- ggplot(totalNormalizedData[totalNormalizedData$scenario == scenarios[i],], 
              aes(maternalAge, mean, group = group, colour = group)) +
    geom_line() +
    labs(x = NULL,
         y = NULL) +
    geom_ribbon(data = totalNormalizedData[totalNormalizedData$scenario == scenarios[i],],
                aes(ymin = min, ymax = max,  fill = group), 
                alpha = 0.2, colour = NA) +
    ylim(0, 2) + 
    theme_minimal() +
    theme(#legend.text = element_text(size=10),
      #legend.key.size = unit(0.2, "cm"),
      #legend.key.width = unit(0.1,"cm"),
      legend.position = "none",
      #legend.title = element_blank(),
      axis.text = element_text(size=13,face="plain",color="black"),
      axis.title = element_text(size = 13),
      #axis.text.x = element_blank(),
      axis.line = element_line(color="black", linewidth = 1.0),
      panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) +
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
    scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)])
  
  # save plot in list 
  plotsTot[[i]] <- p
  # adjust name to corresponding scenario run
  names(plotsTot)[i] <- paste("p_", scenarios[i], sep = "") 
}

# plot matrix using cowplot 
plot_grid(plotsTot$p_damage, NULL, NULL, NULL, 
          plotsTot$p_nullDamage, plotsTot$p_null, NULL, NULL,
          plotsTot$p_damageQuality, plotsTot$p_nullQuality, plotsTot$p_quality, NULL,
          plotsTot$p_damageResource, plotsTot$p_nullResource, plotsTot$p_qualityResource, plotsTot$p_resource,
          ncol = 4)


###############################################################################
# Adding the parental age distribution to the matrix. 
############################################################################### 

# Parental age distribution plots 
predDataTotalAgeDist <- c()

allDataComplete <- loadRDS("allDataComplete.rds")

allDataComplete <- allDataComplete %>% mutate(scenario = factor(scenario))
# loop through the scenarios. 
for (x in levels(allDataComplete$scenario)){
  print(x)
  # make subset of scenario 
  sub <- allDataComplete[allDataComplete$scenario == x,]
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          scenario = sub$scenario[1]) 
  
  sub <- sub %>% mutate(rep = factor(rep))
  # loop through the replicates
  for (i in levels(sub$rep)){
    tmp <- sub[sub$rep == i,]
    
    # use this to count the number of parents per age class 
    counts <- c()
    for (j in seq(0, max(pred_data$maternalAge), 1)) {
      counts <- c(counts, nrow(tmp[tmp$maternalAge == j,]))
    }
    new_col <- paste("rep", i, sep = "")
    # add the counts
    pred_data[[new_col]] <- counts
  }
  
  # get means
  repVals <- subset(pred_data, select = 3:12)
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)

  # combine the data per scenario
  predDataTotalAgeDist <- rbind(predDataTotalAgeDist, pred_data)
}

# add the age distribution plots to the matrix. 

#predDataTotalAgeDist <- read.table("predDataTotalAgeDist.txt")

# get the scenarios relevant for the matrix plot
scenarios <- c()
for (x in levels(totalNormalizedData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) == 2){
    scenarios <- c(scenarios, x) # get list of single scenarios and doubles. 
  }
}

# dynamically generate the plots
plotsAgeDist <- c()
for (i in 1:length(scenarios)){
  p <- ggplot(predDataTotalAgeDist[predDataTotalAgeDist$scenario == scenarios[i],], 
              aes(maternalAge, mean)) +
    geom_line() +
    labs(x = NULL,
         y = NULL) +
    geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, linewidth = 0.8) + 
    #ylim(0, 2) + 
    theme_minimal() +
    theme(#legend.text = element_text(size=10),
      #legend.key.size = unit(0.2, "cm"),
      #legend.key.width = unit(0.1,"cm"),
      legend.position = "none",
      #legend.title = element_blank(),
      axis.text = element_text(size=13,face="plain",color="black"), # size = 11
      axis.title = element_text(size = 13),
      #axis.text.x = element_blank(),
      axis.line = element_line(color="black", linewidth = 1.0),
      panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) +
    xlim(0, 40)
    #scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1))
  
  # save plot in list 
  plotsAgeDist[[i]] <- p
  # adjust name to corresponding scenario run
  names(plotsAgeDist)[i] <- paste("p_", scenarios[i], sep = "") 
}

# plot matrix using cowplot 
plot_grid(plotsTot$p_null, plotsAgeDist$p_nullDamage, plotsAgeDist$p_nullQuality, plotsAgeDist$p_nullResource, 
          plotsTot$p_nullDamage, plotsTot$p_damage, plotsAgeDist$p_damageQuality, plotsAgeDist$p_damageResource,
          plotsTot$p_nullQuality, plotsTot$p_damageQuality, plotsTot$p_quality, plotsAgeDist$p_qualityResource,
          plotsTot$p_nullResource, plotsTot$p_damageResource, plotsTot$p_qualityResource, plotsTot$p_resource,
          ncol = 4)

plotsTot$p_null <- plotsTot$p_null + theme(axis.text.x = element_blank())
plotsTot$p_nullDamage <- plotsTot$p_nullDamage + theme(axis.text.x = element_blank())
plotsTot$p_nullQuality <- plotsTot$p_nullQuality + theme(axis.text.x = element_blank())
plotsTot$p_damage <- plotsTot$p_damage + theme(axis.text.x = element_blank(),
                                               axis.text.y = element_blank())
plotsTot$p_damageQuality <- plotsTot$p_damageQuality + theme(axis.text.x = element_blank(),
                                                             axis.text.y = element_blank())
plotsTot$p_damageResource <- plotsTot$p_damageResource + theme(axis.text.y = element_blank())
plotsTot$p_quality <- plotsTot$p_quality + theme(axis.text.x = element_blank(),
                                                 axis.text.y = element_blank())
plotsTot$p_qualityResource <- plotsTot$p_qualityResource + theme(axis.text.y = element_blank())
plotsTot$p_resource <- plotsTot$p_resource + theme(axis.text.y = element_blank())
plotsAgeDist$p_nullQuality <- plotsAgeDist$p_nullQuality + theme(axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank())
plotsAgeDist$p_nullResource <- plotsAgeDist$p_nullResource + theme(axis.text.x = element_blank(),
                                                                   axis.text.y = element_blank())
plotsAgeDist$p_damageResource <- plotsAgeDist$p_damageResource + theme(axis.text.x = element_blank(),
                                                                       axis.text.y = element_blank())

plot_matrix <- plot_grid(plotsTot$p_null, plotsAgeDist$p_nullDamage, plotsAgeDist$p_nullQuality, plotsAgeDist$p_nullResource,
                         plotsTot$p_nullDamage, plotsTot$p_damage, plotsAgeDist$p_damageQuality, plotsAgeDist$p_damageResource,
                         plotsTot$p_nullQuality, plotsTot$p_damageQuality, plotsTot$p_quality, plotsAgeDist$p_qualityResource, 
                         plotsTot$p_nullResource, plotsTot$p_damageResource, plotsTot$p_qualityResource, plotsTot$p_resource,
                         ncol = 4, align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.88, label_y = 0.97)
plot_matrix

p_tmp <- plotsTot$p_resource + theme(legend.position = "bottom",
                                     legend.title = element_blank())

leg1 <- get_legend(p_tmp) 
library(ggpubr)
legend <- as_ggplot(leg1)
legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)


#null <- ggplot() + annotate("text", x = 0, y = 0, size = 4, label = "Baseline") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#damage <- ggplot() + annotate("text", x = 1, y = 0, size = 4, label = "Damage accumulation") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#quality <- ggplot() + annotate("text", x = 1, y = 0, size = 4, label = "Parental care quality") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#resource <- ggplot() + annotate("text", x = 1, y = 0, size = 4, label = "Resource allocation") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#top_labels <- plot_grid(null, damage, quality, resource, rel_widths = c(1, 1, 1, 1), nrow = 1)

#null2 <- ggplot() + annotate("text", x = 0, y = 0, size = 4, label = "Baseline", angle = "270") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#damage2 <- ggplot() + annotate("text", x = 0, y = 0, size = 4, label = "Damage accumulation", angle = "270") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#quality2 <- ggplot() + annotate("text", x = 0, y = 0, size = 4, label = "Parental care quality", angle = "270") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#resource2 <- ggplot() + annotate("text", x = 0, y = 0, size = 4, label = "Resource allocation", angle = "270") + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#right_labels <- plot_grid(null2, damage2, quality2, resource2, rel_heights = c(1, 1, 1, 1), ncol = 1)


library(grid)
library(gridExtra)
plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                         left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))

#plot_with_legend <- plot_grid(top_labels, plot_matrix2, legend_row, ncol = 1, rel_heights = c(0.02, 1, 0.01), align = "v", axis = "rl")
#plot_with_legend2 <- plot_grid(plot_with_legend, right_labels, nrow = 1, rel_widths = c(1, 0.02), align = "h", axis = "tb")

plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
plot_with_legend
ggsave("Lansing_fig.pdf", width = 12, height = 10)


###############################################################################
# plotting the three scenarios combined 
############################################################################### 

scenarios <- c()
for (x in levels(totalNormalizedData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) == 3){
    scenarios <- c(scenarios, x) # get list of single scenarios and doubles. 
  }
}

# dynamically generate the plots
plotsThreeCombs <- c()
for (i in 1:length(scenarios)){
  p <- ggplot(totalNormalizedData[totalNormalizedData$scenario == scenarios[i],], 
              aes(maternalAge, mean, group = group, colour = group)) +
    geom_line() +
    labs(x = NULL,
         y = NULL) +
    geom_ribbon(data = totalNormalizedData[totalNormalizedData$scenario == scenarios[i],],
                aes(ymin = min, ymax = max,  fill = group), 
                alpha = 0.2, colour = NA) + 
    #ylim(0, 2) + 
    theme_minimal() +
    theme(#legend.text = element_text(size=10),
      #legend.key.size = unit(0.2, "cm"),
      #legend.key.width = unit(0.1,"cm"),
      legend.position = "none",
      #legend.title = element_blank(),
      axis.text = element_text(size=13,face="plain",color="black"),
      axis.title = element_text(size = 13),
      #axis.text.x = element_blank(),
      axis.line = element_line(color="black", linewidth = 0.6),
      panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.5)) + 
  scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1)) + 
  scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
  scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
    scale_y_continuous(breaks = seq(0, 1.5, 0.2)) +
    coord_cartesian(ylim= c(0, 1.5))

  
  # save plot in list 
  plotsThreeCombs[[i]] <- p
  # adjust name to corresponding scenario run
  names(plotsThreeCombs)[i] <- paste("p_", scenarios[i], sep = "") 
}

# needs formatting 
plot_grid(plotsThreeCombs$p_nullDamageQuality,
          plotsThreeCombs$p_nullDamageResource,
          plotsThreeCombs$p_damageQualityResource,
          plotsThreeCombs$p_nullQualityResource)

plotsThreeCombs$p_nullDamageQuality <- plotsThreeCombs$p_nullDamageQuality + theme(axis.text.x = element_blank())
plotsThreeCombs$p_nullDamageResource <- plotsThreeCombs$p_nullDamageResource + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
plotsThreeCombs$p_nullQualityResource <- plotsThreeCombs$p_nullQualityResource + theme(axis.text.y = element_blank())

plot_matrix <- plot_grid(plotsThreeCombs$p_nullDamageQuality,
          plotsThreeCombs$p_nullDamageResource,
          plotsThreeCombs$p_damageQualityResource,
          plotsThreeCombs$p_nullQualityResource,
          align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.90, label_y = 0.97)

plot_matrix

p_tmp <- plotsThreeCombs$p_nullQualityResource + theme(legend.position = "bottom",
                                     legend.title = element_blank(),
                                     legend.text = element_text(size=13))

leg1 <- get_legend(p_tmp) 
library(ggpubr)
legend <- as_ggplot(leg1)
legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)


library(grid)
library(gridExtra)
plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                         left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))

plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
plot_with_legend
ggsave("Lansing_fig_suppl.pdf", width = 12, height = 10)


########## suppl figure for all 4 scenarios combined 
scenario <- "nullDamageQualityResource"
p <- ggplot(totalNormalizedData[totalNormalizedData$scenario == scenario,], 
            aes(maternalAge, mean, group = group, colour = group)) +
  geom_line() +
  labs(x = "Normalized parental age",
       y = "Normalized offspring lifespan") +
  geom_ribbon(data = totalNormalizedData[totalNormalizedData$scenario == scenario,],
              aes(ymin = min, ymax = max,  fill = group), 
              alpha = 0.2, colour = NA) + 
  theme_minimal() +
  theme(legend.text = element_text(size=13),
    #legend.key.size = unit(0.2, "cm"),
    #legend.key.width = unit(0.1,"cm"),
    legend.position = c(0.10, -0.04),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    axis.text = element_text(size=13,face="plain",color="black"),
    axis.title = element_text(size = 13),
    #axis.text.x = element_blank(),
    axis.line = element_line(color="black", linewidth = 0.6),
    panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.5)) + 
  scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1)) + 
  scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
  scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
  scale_y_continuous(breaks = seq(0, 1.3, 0.2)) +
  coord_cartesian(ylim = c(0,1.3))
p
ggsave("All_scenarios.pdf", width = 12, height = 8)

###############################################################################
# Plotting the supplementary figures with varying parameters 
############################################################################### 

# look at tracked individuals data 
parent_path <- paste0(path, "parameterSim/null/")
parent_path <- paste0(path, "parameterSim/damage/")
parent_path <- paste0(path, "parameterSim/quality/")
parent_path <- paste0(path, "parameterSim/resource/")

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
found <- c()
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nMutProb <- length(splitted_path[[1]]) - 1
  mutProb <- splitted_path[[1]][nMutProb]
  # read data
  local_data <- read.table(file_name, header = F, sep = " ")
  colnames(local_data) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  local_data$mutProb <- mutProb
  local_data$ID <- sub("^", local_data$mutProb[1], local_data$ID)
  
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  
  #allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 200, replace = F),]
  
  found <- rbind(found, local_data)
}

# average the offspring lifespan per maternal age per mutation probability 
calcMargin = function(x) { # function to calculate the margin for confidence interval 
  n <- x[4] # get sample size
  n <- as.numeric(n)
  sd <- x[5] # get sd 
  sd <- as.numeric(sd)
  return(qt(0.975, df = n - 1)*sd/sqrt(n))
}
  
found <- found %>% mutate(mutProb = factor(mutProb))
allAvgs <- c()
for (x in levels(found$mutProb)) {
  sub <- found[found$mutProb == x,]
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  sub <- sub[sub$maternalAge <= percentile,]
  # calculate average expected age at death from parental age of 0 until 95th percentile parental age 
  avg <- aggregate(sub$ageAtDeath, list(sub$maternalAge), mean)
  colnames(avg) <- c("ageOfParent", "avgOffspringLifespan")
  avg$mutProb <- x
  
  # get information for CIs. 
  avg$n <- tapply(sub$ageAtDeath, sub$maternalAge, length)
  avg$sd <- tapply(sub$ageAtDeath, sub$maternalAge, sd)
  
  avg$margin <- apply(avg, 1, calcMargin) # calc margin per row
  avg$lwr <- by(avg, seq(nrow(avg)), function(x){x[2] - x[6]}) # get lower bound per row
  avg$upr <- by(avg, seq(nrow(avg)), function(x){x[2] + x[6]}) # get upper bound per row
  
  allAvgs <- rbind(allAvgs, avg)
}

allAvgs <- allAvgs %>% mutate(lwr = as.numeric(lwr),
                              upr = as.numeric(upr))

# with 95% CIs. 
legend_title <- "Mutation probability for age-specific 
survival genes"
ggplot(allAvgs, aes(ageOfParent, avgOffspringLifespan, group = mutProb, colour = mutProb)) +
  geom_line() +
  geom_ribbon(data = allAvgs,
              aes(ymin = lwr, ymax = upr, fill = mutProb),
              alpha = 0.2, colour = NA) +
  theme_minimal() +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        #axis.text.x = element_blank(),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_colour_manual(legend_title, values = met.brewer("Egypt")[1:5]) +
  scale_fill_manual(legend_title, values = met.brewer("Egypt")[1:5]) + 
  ylim(0, 30) +
  labs(x = "Maternal age",
       y = "Average offspring lifespan") 


ggsave("SupplResource.pdf", width = 12, height = 8)


###############################################################################
# plotting the resource allocation gene values
############################################################################### 
myData <- read.table(paste0(path, "combiningAll3/resource/1/outputWithAgeSpecificGenes.txt"))
colnames(myData) <- c("ID", "age", "investmentGeneVal", "survivalGeneVal")
ggplot(myData[myData$ID %in% sample(myData$ID, 10),], aes(age, investmentGeneVal, group = as.factor(ID), colour = as.factor(ID))) +
  geom_line()



# plot to show effect of s 
funcPP <- function(x) {exp(-x[1]*x[2])}
s <- c(0, 0.01, 0.03, 0.05, 0.07, 0.1)
D <- 0:20
tmp <- expand.grid(s,D)
colnames(tmp) <- c("s", "D")
res <- apply(tmp, 1,function(x){exp(-x[1]*x[2])})
tmp$res <- res
ggplot(tmp, aes(D, res, group = as.factor(s), colour = as.factor(s))) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 30),
        strip.text.x = element_text(size = 35, face = "bold")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     limits = c(0,1)) +
  scale_color_manual("s", values = met.brewer("Egypt", 6)) +
  labs(y = "m2")
ggsave("ppEquation.pdf", width = 12, height = 8)

c <- seq(0,1,0.1)
x <- seq(0,1, 0.1)
tmp <- expand.grid(c,x)
colnames(tmp) <- c("c", "x")
res <- apply(tmp, 1, function(x){1 - x[1] * x[2]^2})
tmp$res <- res

ggplot(tmp, aes(x, res, group = as.factor(c), colour = as.factor(c))) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 30),
        strip.text.x = element_text(size = 35, face = "bold")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     limits = c(0,1)) +
  scale_color_manual("c", values = met.brewer("Egypt", 11)) +
  labs(y = "m4a")
ggsave("ppEquation2.pdf", width = 12, height = 8)

a <- 0:5
tmp <- expand.grid(a,x)
colnames(tmp) <- c("a", "x")
res <- apply(tmp, 1, function(x){1 / (1+exp(-x[1]*x[2]-1))})
tmp$res <- res

ggplot(tmp, aes(x, res, group = as.factor(a), colour = as.factor(a))) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 30),
        strip.text.x = element_text(size = 35, face = "bold")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     limits = c(0,1)) +
  scale_color_manual("a", values = met.brewer("Egypt", 11)) +
  labs(y = "m4b")
ggsave("ppEquation3.pdf", width = 12, height = 8)

b <- 0:5
tmp <- expand.grid(b,x)
colnames(tmp) <- c("b", "x")
res <- apply(tmp, 1, function(x){1 / (1+exp(-3*x[2]-x[1]))})
tmp$res <- res
funcTMP <- function(x) {1 / (1 + exp(-3*x - 1))}
x <- seq(0,1,0.1)
res <- funcTMP(x)
dfTMP <- data.frame(x, res)

ggplot(dfTMP, aes(x, res)) + 
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 30),
        strip.text.x = element_text(size = 35, face = "bold")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     limits = c(0,1)) +
  #scale_color_manual("b", values = met.brewer("Egypt", 11)) +
  labs(y = "m4b")
ggsave("ppEquation4.pdf", width = 12, height = 8)


##### plot age-specific gene values 
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/combiningAll/"

parent_path <- paste0(path, "combiningAll3/nullDamage/")
parent_path <- paste0(path, "combiningAll3/resource/")
parent_path <- paste0(path, "combiningAll3/damageResource/")
parent_path <- paste0(path, "combiningAll3/qualityResource/")
parent_path <- paste0(path, "combiningAll3/nullResource/")

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes.txt", recursive = T, all.files = T)
found <- c()
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  nScenario <- length(splitted_path[[1]]) - 2
  scenario <- splitted_path[[1]][nScenario]
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = F)
  colnames(local_data) <- c("ID", "age", "InvestmentGeneVal", "SurvivalGeneVal")
  # add the scenario as a column to the data
  local_data$scenario <- scenario
  # add replicate as column 
  local_data$rep <- rep
  local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
  local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
  
  found <- rbind(found, local_data)
  
}

found <- found %>% mutate(scenario = factor(scenario))
found <- found %>% mutate(rep = factor(rep))

# for age-specific survival genes 
tmp <- aggregate(found$SurvivalGeneVal, list(found$age), mean)
colnames(tmp) <- c("age", "mean")
tmp$min <- tapply(found$SurvivalGeneVal, found$age, min)
tmp$max <- tapply(found$SurvivalGeneVal, found$age, max)

ggplot(tmp, aes(age, mean)) + geom_line() +
  geom_ribbon(aes(ymin =min, ymax = max), alpha = 0.2) +
  theme_bw() +
  labs(x = "Age",
       y = "Averaged survival gene value") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  theme(axis.text = element_text(size=17,face="plain",color="black"), 
  axis.title = element_text(size = 17))
ggsave("ppNullDamAgeGenes.pdf", width = 10, height = 8)

# for age-specific resource genes
tmp <- aggregate(found$InvestmentGeneVal, list(found$age), mean)
colnames(tmp) <- c("age", "mean")
tmp$min <- tapply(found$InvestmentGeneVal, found$age, min)
tmp$max <- tapply(found$InvestmentGeneVal, found$age, max)

ggplot(tmp, aes(age, (1-mean))) + geom_line() +
  geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
  theme_bw() +
  labs(x = "Age",
       y = "Averaged proportion of resources 
  allocated to reproduction") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  theme(axis.text = element_text(size=17,face="plain",color="black"), 
        axis.title = element_text(size = 17))
ggsave("ppResourceNull.pdf", width = 12, height = 8)


