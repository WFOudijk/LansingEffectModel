library(ggplot2)

###############################################################################
# Decline in gamete quality analysis
###############################################################################

path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/"
# change path to where you want the plots to be outputted 
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/plots/"

# path to decline in quality with both parental effects implemented 
myData <- read.table(paste(path, "outputDeclineGameteQuality 2.csv", sep = "")) # with time 

colnames(myData) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather")
# the dataset without the age of 40 included
myData <- myData[myData$ageAtDeath != 40, ]

# averaging data 
avgDataframe <- aggregate(myData$ageAtDeath, list(myData$ageOfMother), mean)
colnames(avgDataframe) <- c("ageOfMother", "ageAtDeath")

# looking at average age of death over maternal ages
ggplot(data = avgDataframe, aes(x = ageOfMother, y = ageAtDeath)) +
  geom_point() + 
  labs(title = "Average age at death of offspring over maternal ages",
       subtitle = "paternal and maternal decline in gamete quality implemented",
       x = "Age of mother",
       y = "Average age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20))
#ggsave(paste(output_path, "plot_maternal_gamete_paternal_implemented.pdf", sep = ""), plot = last_plot()) #to save file 

### paternal effect 
# average age of death over paternal ages
avgDataframe <- aggregate(myData$ageAtDeath, list(myData$ageOfFather), mean)
colnames(avgDataframe) <- c("ageOfFather", "ageAtDeath")

ggplot(data = avgDataframe, aes(x = ageOfFather, y = ageAtDeath)) +
  geom_point() + 
  labs(title = "Average age at death of offspring over paternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells",
       x = "Age of father",
       y = "Average age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(paste(output_path, "plot_paternal_stem_cell_decline.pdf", sep = ""), plot = last_plot()) #to save file 

######### plotting final time step ######
# only look at ending time points

## data with both maternal and paternal decline in gamete quality implemented. With the stem cells identical to the genome. 
myData <- read.table(paste(path, "outputDeclineGameteQuality 3.csv", sep = "")) # with time 
myData <- read.table(paste(path, "outputDeclineGameteQuality 4.csv", sep = "")) # with a seperate small mutation prob (0.0005) for paternal decline 
myData <- read.table(paste(path, "outputDeclineGameteQuality.csv", sep = "")) # with a mutation prob for SC (0.001) for paternal decline 

colnames(myData) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC")

# look from time > 6000 
equi_timepoint <- subset(myData, time > 6000)
# get average 
avgDataframe <- aggregate(equi_timepoint$ageAtDeath, list(equi_timepoint$ageOfMother), mean)
colnames(avgDataframe) <- c("ageOfMother", "AverageAgeAtDeath")

ggplot(data = equi_timepoint, aes(x = ageOfMother, y = ageAtDeath, group = ageOfMother)) +
  geom_boxplot() + 
  labs(title = "Age at death of offspring over maternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells and maternal gametes. Time > 6000",
       x = "Age of mother",
       y = "Age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  theme_big
#ggsave(paste(output_path, "boxplot_maternal_time_interval.pdf", sep = ""), plot = last_plot()) #to save file 

# paternal
avgDataframe <- aggregate(equi_timepoint$ageAtDeath, list(equi_timepoint$ageOfFather), mean)
colnames(avgDataframe) <- c("ageOfFather", "AverageAgeAtDeath")

ggplot(data = equi_timepoint, aes(x = ageOfFather, y = ageAtDeath, group = ageOfFather)) +
  geom_boxplot() + 
  labs(title = "Age at death of offspring over paternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells and maternal gametes. Time > 6000.",
       x = "Age of father",
       y = "Age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) + 
  theme_big
#ggsave(paste(output_path, "boxplot_paternal_time_interval.pdf", sep = ""), plot = last_plot()) #to save file 

###############################################################################
# Decline in gamete quality analysis - intermediate parameter exploration
# looking at the mutation probability of the stem cells. 
###############################################################################
parent_path <- paste(path, "stemcellMutExploration/", sep = "")
f <- list.files(path = parent_path, pattern = "outputDeclineGameteQuality.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell") 

# averaging the age at deaths 
avgDataframe <- aggregate(found$ageAtDeath, list(found$time, found$mutationProbStemCell), mean)
colnames(avgDataframe) <- c("time", "mutationProbSC", "ageAtDeath")
suba <- subset(avgDataframe, mutationProbSC == 0.001) # choose 0.001 to get average ~0.001

ggplot(data = avgDataframe, aes(x = time, y = ageAtDeath, col = mutationProbSC, group = mutationProbSC)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Age of death over time with differing mutation probabilities for stem cell mutation", 
       x = "Time",
       y = "Age of death") +
  theme_classic()
#ggsave(paste(output_path, "plot_SC_mut.pdf", sep = ""), plot = last_plot()) #to save file 

###############################################################################
# Analyzing life expectancy over the surviving population 
###############################################################################
survivingPop <- read.table(paste(path, "outputLifeExpectancy.txt", sep = ""))
colnames(survivingPop) <- c("age", "expectedAgeAtDeath", "ageOfMother", "ageOfFather", "survivalProb")

# look at maternal ages 
ggplot(data = survivingPop, aes(x = ageOfMother, y = expectedAgeAtDeath, group = ageOfMother)) +
  geom_boxplot() + 
  labs(title = "Expected age at death of offspring over maternal ages after the final generation",
       subtitle = "Decline in gamete quality implemented.",
       x = "Age of mother",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)
#ggsave(paste(output_path, "life_exp_maternal.pdf", sep = ""), plot = last_plot()) #to save file 

# paternal 
ggplot(data = survivingPop, aes(x = ageOfFather, y = expectedAgeAtDeath, group = ageOfFather)) +
  geom_boxplot() + 
  #geom_point() +
  labs(title = "Expected age at death of offspring over paternal ages after the final generation",
       subtitle = "Decline in gamete quality implemented.",
       x = "Age of father",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40) 
#ggsave(paste(output_path, "life_exp_paternal.pdf", sep = ""), plot = last_plot()) #to save file 

###############################################################################
# Analyzing life expectancy over the surviving population - extrinsic mortality = 0.05
###############################################################################
myData <- read.table(paste(path, "outputDeclineGameteQuality 5.csv", sep = ""))
colnames(myData) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProbability", "mutationProbSC")

equi_timepoint <- subset(myData, time > 6000)

# get average 
avgDataframe <- aggregate(equi_timepoint$ageAtDeath, list(equi_timepoint$ageOfMother), median) 
colnames(avgDataframe) <- c("ageOfMother", "AverageAgeAtDeath")
# get minimum and maximum for every parental age class 
avgDataframe$minAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfMother, min)
avgDataframe$maxAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfMother, max)

# plotting with lines to indicate min and max 
ggplot(data = avgDataframe, aes(x = ageOfMother, y = AverageAgeAtDeath, group = ageOfMother)) +
  geom_point() +
  geom_errorbar(aes(ymin = minAge, ymax = maxAge), 
                 linetype = 2, width = 0.2) +
  labs(title = "Age at death of offspring over maternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells and maternal gametes. Time > 6000",
       x = "Age of mother",
       y = "Age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  theme_bw()

###############################################################################
# setting mutation probabilities both on 0.0045
###############################################################################
myData <- read.table(paste(path, "outputDeclineGameteQuality 6.csv", sep=""))
colnames(myData) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProbability", "mutationProbSC")

equi_timepoint <- subset(myData, time > 6000)
# get average 
avgDataframe <- aggregate(equi_timepoint$ageAtDeath, list(equi_timepoint$ageOfMother), median) # TODO: median? 
colnames(avgDataframe) <- c("ageOfMother", "AverageAgeAtDeath")
# get minimum and maximum for every parental age class 
avgDataframe$minAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfMother, min)
avgDataframe$maxAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfMother, max)

# plotting with lines to indicate min and max 
ggplot(data = avgDataframe, aes(x = ageOfMother, y = AverageAgeAtDeath, group = ageOfMother)) +
  geom_point() +
  geom_errorbar(aes(ymin = minAge, ymax = maxAge), 
                linetype = 2, width = 0.2) +
  labs(title = "Age at death of offspring over maternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells and maternal gametes. Time > 6000",
       x = "Age of mother",
       y = "Age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  theme_bw()

# paternal
# get average 
avgDataframe <- aggregate(equi_timepoint$ageAtDeath, list(equi_timepoint$ageOfFather), median) # TODO: median? 
colnames(avgDataframe) <- c("ageOfFather", "AverageAgeAtDeath")
# get minimum and maximum for every parental age class 
avgDataframe$minAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfFather, min)
avgDataframe$maxAge <- tapply(equi_timepoint$ageAtDeath, equi_timepoint$ageOfFather, max)

# plotting with lines to indicate min and max 
ggplot(data = avgDataframe, aes(x = ageOfFather, y = AverageAgeAtDeath, group = ageOfFather)) +
  geom_point() +
  geom_errorbar(aes(ymin = minAge, ymax = maxAge), 
                linetype = 2, width = 0.2) +
  labs(title = "Age at death of offspring over paternal ages",
       subtitle = "Decline in gamete quality implemented by mutating stem cells and paternal gametes. Time > 6000",
       x = "Age of father",
       y = "Age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  theme_bw()


###############################################################################
# Decline in gamete quality analysis - intermediate parameter exploration
# looking at the mutation probability of the gametes and the stem cells
###############################################################################

parent_path <- paste(path, "mutationProb/", sep = "")
f <- list.files(path = parent_path, pattern = "outputDeclineGameteQuality.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb") 

# averaging the age at deaths 
avgDataframe <- aggregate(found$ageAtDeath, list(found$time, found$mutationProb), mean)
colnames(avgDataframe) <- c("time", "mutationProb", "ageAtDeath")

ggplot(data = avgDataframe, aes(x = time, y = ageAtDeath, col = mutationProb, group = mutationProb)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Age of death over time with differing mutation probabilities for stem cell mutation", 
       x = "Time",
       y = "Age of death") +
  theme_classic()
#ggsave(paste(output_path, "plot_mut_prob.pdf", sep = ""), plot = last_plot()) #to save file 

###### looking at mutation prob for SC ###### (with mutationProb (gametes) set to 0.004)
parent_path <- paste(path, "mutationProbSC/", sep = "")
f <- list.files(path = parent_path, pattern = "outputDeclineGameteQuality.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb") 

# averaging the age at deaths 
avgDataframe <- aggregate(found$ageAtDeath, list(found$time, found$mutationProbStemCell), mean)
colnames(avgDataframe) <- c("time", "mutationProbStemCell", "ageAtDeath")
suba <- subset(avgDataframe, mutationProb == 0.004) # choose 0.001 to get average ~0.001

ggplot(data = avgDataframe, aes(x = time, y = ageAtDeath, col = mutationProbStemCell, 
                                group = mutationProbStemCell)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Age of death over time with differing mutation probabilities for stem cell mutation", 
       x = "Time",
       y = "Age of death") +
  theme_classic()
#ggsave(paste(output_path, "plot_mut_prob.pdf", sep = ""), plot = last_plot()) #to save file 

###############################################################################
# Statistical analysis
###############################################################################

###### get data ready for statistical analysis ########
output_path <- paste(path, "data_for_statistics_2/", sep = "") # set path
########## both mutation probabilities set to 0.0045 ##########
output_path_2 <- paste(output_path, "outputLifeExpectancy_same_mut_probs.txt", sep = "") # both to 0.0045 
output_path_2 <- paste(output_path, "outputLifeExpectancy_gametes_lower_SC.txt", sep = "") # gametes = 0.001; stem cells = 0.0045
output_path_2 <- paste(output_path, "outputLifeExpectancy_gametes_higher_SC.txt", sep = "") # gametes = 0.0045; stem cells = 0.001 
output_path_2 <- paste(output_path, "outputLifeExpectancy_gametes_same_SC_higherProb.txt", sep = "") # both mut probs set to 0.01
output_path_2 <- paste(output_path, "outputLifeExpectancy_gametes_higher_SC2.txt", sep = "") # gametes = 0.006; stem cells = 0.002
output_path_2 <- paste(output_path, "outputLifeExpectancy_gametes_same_SC2.txt", sep = "") # both to 0.005 
output_path_2 <- paste(output_path, "outputLifeExpectancy_same_mut_probs2.txt", sep = "") # both to 0.0045 2.0 

survivingPop <- read.table(output_path_2)
colnames(survivingPop) <- c("age", "expectedAgeAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb")

# get average
avgDataframe <- aggregate(survivingPop$expectedAgeAtDeath, list(survivingPop$ageOfMother), median) 
colnames(avgDataframe) <- c("ageOfMother", "medianAgeAtDeath")
# get minimum and maximum 
avgDataframe$minAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfMother, min)
avgDataframe$maxAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfMother, max)

# look at maternal ages 
ggplot(data = avgDataframe, aes(x = ageOfMother, y = medianAgeAtDeath, group = ageOfMother)) +
  geom_point() + 
  geom_linerange(aes(ymin = minAge, ymax = maxAge), 
                linetype = 2) +
  labs(title = "Median expected age at death of offspring over maternal ages after the final generation",
       subtitle = "Dashed lines point to min and max expected age",
       x = "Age of mother",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)

#ggsave(paste(output_path, "plots/life_exp_maternal_gametes_same_SC_higherProb.pdf", sep = ""), plot = last_plot()) #to save file 

# paternal
# get average 
avgDataframe <- aggregate(survivingPop$expectedAgeAtDeath, list(survivingPop$ageOfFather), median) # TODO: median? 
colnames(avgDataframe) <- c("ageOfFather", "medianAgeAtDeath")
# get minimum and maximum 
avgDataframe$minAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfFather, min)
avgDataframe$maxAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfFather, max)

# look at paternal ages 
ggplot(data = avgDataframe, aes(x = ageOfFather, y = medianAgeAtDeath, group = ageOfFather)) +
  geom_point() + 
  geom_linerange(aes(ymin = minAge, ymax = maxAge), 
                 linetype = 2) +
  labs(title = "Median expected age at death of offspring over paternal ages after the final generation",
       subtitle = "Dashed lines point to min and max expected age",
       x = "Age of father",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)
ggsave(paste(output_path, "life_exp_paternal_gametes_high_SC.pdf", sep = ""), plot = last_plot()) #to save file 

#make subset for statistical analysis 
sub_maternal_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfMother))
sub_paternal_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfFather))
sub_both_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfMother, ageOfFather))

###############################################################################
# To gain more results - running the simulation 10 times. Both mutation probabilities = 0.006
###############################################################################

parent_path <- paste(output_path, "lifeExpectancies10simulations/", sep = "")
parent_path <- paste(output_path, "lifeExp10Sim/", sep = "")

f <- list.files(path = parent_path, pattern = "outputLifeExpectancy.txt", recursive = T)
survivingPop <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.table(file_name)
  survivingPop <- rbind(survivingPop, local_data)
}
colnames(survivingPop) <- c("age", "expectedAgeAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb")

# maternal effects 
avgDataframe <- aggregate(survivingPop$expectedAgeAtDeath, list(survivingPop$ageOfMother), median) 
colnames(avgDataframe) <- c("ageOfMother", "medianAgeAtDeath")
avgDataframe$minAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfMother, min)
avgDataframe$maxAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfMother, max)
avgDataframe <- avgDataframe %>% mutate(y1 = pmin(40,maxAge))

# look at maternal ages 
ggplot(data = avgDataframe, aes(x = ageOfMother, y = medianAgeAtDeath, group = ageOfMother)) +
  geom_point() + 
  geom_linerange(aes(ymin = minAge, ymax = y1), 
                 linetype = 2) +
  labs(title = "Median expected age at death of offspring over maternal ages after the final generation",
       subtitle = "Dashed lines point to min and max expected age",
       x = "Age of mother",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)
#output_path_2 <- paste(output_path, "data_for_statistics_2/", sep = "")
#ggsave(paste(output_path, "life_exp_maternal_gametes_same_SC_higherProb.pdf", sep = ""), plot = last_plot()) #to save file 

# paternal
avgDataframe <- aggregate(survivingPop$expectedAgeAtDeath, list(survivingPop$ageOfFather), median)  
colnames(avgDataframe) <- c("ageOfFather", "medianAgeAtDeath")
avgDataframe$minAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfFather, min)
avgDataframe$maxAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfFather, max)
avgDataframe <- avgDataframe %>% mutate(y1 = pmin(40,maxAge))

# look at paternal ages 
ggplot(data = avgDataframe, aes(x = ageOfFather, y = medianAgeAtDeath, group = ageOfFather)) +
  geom_point() + 
  geom_linerange(aes(ymin = minAge, ymax = y1), 
                 linetype = 2) +
  labs(title = "Median expected age at death of offspring over paternal ages after the final generation",
       subtitle = "Dashed lines point to min and max expected age",
       x = "Age of father",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)
#ggsave(paste(output_path, "life_exp_paternal_gametes_high_SC.pdf", sep = ""), plot = last_plot()) #to save file 

# make subset for statistical analysis 
sub_maternal_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfMother))
sub_paternal_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfFather))

###############################################################################
# Looking at sex specific effects 
###############################################################################
output_path_2 <- paste(output_path, "outputLifeExpectancy_sexSpecific.txt", sep = "")
survivingPop <- read.table(output_path_2)
colnames(survivingPop) <- c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

sub_sex_spec <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfParent, sexOfParent))

###############################################################################
# Looking at the individuals longitudinally 
###############################################################################
##### looking at longitudinal effects #####
path_to_read <- paste(output_path, "longitudinal/", sep = "")
myLongitudinalData <- read.table(paste(path_to_read, "outputLETrackedIndividuals.txt", sep = ""))
myLongitudinalData <- read.table(paste(path_to_read, "outputLETrackedIndividuals_500_tracked.txt", sep = ""))

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "survivalProb", "expectedAgeAtDeath")
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

myData <- read.table(paste(path_to_read, "outputLifeExpectancy.txt", sep = ""))
colnames(myData) <- c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
subSexSpec <- subset(myData, select = c(expectedAgeAtDeath, ageOfParent, sexOfParent))

###############################################################################
# The final data - sex of parent in the datasets. 
# both mutation probabilities are set to 0.0045. 
# 500 individuals are tracked after the first simulation 
###############################################################################
allDeadIndividuals <- read.table(paste(path_to_read, "outputDeclineGameteQualityFINAL.txt", sep = ""))
colnames(allDeadIndividuals) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb") 

survivingPop <- read.table(paste(path_to_read, "outputLifeExpectancyFINAL.txt", sep = ""))
colnames(survivingPop) <- c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

trackedIndividuals <- read.table(paste(path_to_read, "outputLETrackedIndividualsFINAL.txt", sep = ""))
colnames(trackedIndividuals) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")






