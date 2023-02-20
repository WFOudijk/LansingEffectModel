###############################################################################
# AGEING ADDED - ageing is implemented by adding an age-specific gene array
# represented by survival probabilities. The individuals overall survival 
# probability is now calculated by multiplying the binary-based survival 
# probability with the age-specific survival probability 
###############################################################################
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/AGEING_ADDED/"
survivingPop <- read.table(paste(path, "outputLifeExpectancyCorrect.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

avgDataframe <- aggregate(survivingPop$expectedAgeAtDeath, list(survivingPop$ageOfParent), median) 
colnames(avgDataframe) <- c("ageOfParent", "medianAgeAtDeath")
avgDataframe$minAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfParent, min)
avgDataframe$maxAge <- tapply(survivingPop$expectedAgeAtDeath, survivingPop$ageOfParent, max)

# look at data
ggplot(data = avgDataframe, aes(x = ageOfParent, y = medianAgeAtDeath, group = ageOfParent)) +
  geom_point() + 
  geom_linerange(aes(ymin = minAge, ymax = maxAge), 
                 linetype = 2) +
  labs(title = "Median expected age at death of offspring over parental ages after the final generation",
       subtitle = "Dashed lines point to min and max expected age",
       x = "Age of parent",
       y = "Expected age of death") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) +
  ylim(0, 40)
#output_path_2 <- paste(output_path, "data_for_statistics_2/", sep = "")
#ggsave(paste(output_path, "life_exp_parental_gametes_same_SC_higherProb.pdf", sep = ""), plot = last_plot()) #to save file 

#make subset for statistical analysis 
sub_parental_LE <- subset(survivingPop, select = c(expectedAgeAtDeath, ageOfParent))

# look at tracked individuals 
myLongitudinalData <- read.table(paste(path, "outputLETrackedIndividualsCorrect.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

# look at age of death for every time step 
myData <- read.table(paste(path, "outputDeclineGameteQuality.txt", sep = ""))
colnames(myData) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb") 

equi_timepoint <- subset(myData, time > 6000)

# maternal results 
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

# paternal effects 
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
# Parameter exploration of ratio between mutation probability of stem cells and gametes
# and mutation probability of the age-specific genes 
###############################################################################

# {mutationProbAgeSpecificGenes; mutationProbStemCell; mutationProb}
# first try = {0.003; 0.004; 0.004} > no significant decrease 
# second try = {0.002; 0.0035; 0.0035} > no significant decrease and more randomeness. Too low? 
# third try: {0.0025; 0.0035; 0.0035} > first lmer model not significant, second and third are significant 
# fourth try = {0.0025; 0.0037; 0.0037} > no singularity issue! 
survivingPop <- read.table(paste(path, "outputLifeExpectancySmallMuts3.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

myLongitudinalData <- read.table(paste(path, "outputLETrackedIndividualsSmallMuts4.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

deathIndividuals <- read.table(paste(path, "outputDeclineGameteQualitySmallMuts3.txt", sep = ""))
colnames(deathIndividuals) <- c("time", "ageAtDeath", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbStemCell", "mutationProb") 

# fourth try parameter values but with 5000 individuals and 3000 tracked. 

