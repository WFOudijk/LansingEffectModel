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
survivingPop <- read.table(paste(path, "outputLifeExpectancy5000Indv.txt", sep = ""))
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

myLongitudinalData <- read.table(paste(path, "outputLETrackedIndividuals5000Indv.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

# test to see if ageing is implemented correct by turning the age-specific mutation rate to 0 
# both mutation probabilities for gametes and stem cells are set to 0.0045
survivingPop <- read.table(paste(path, "outputLifeExpectancyNoAgeing.txt", sep = ""))
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

myLongitudinalData <- read.table(paste(path, "outputLETrackedIndividualsNoAgeing.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

######### old model ########
survivingPopOld <- read.table(paste(path, "outputLifeExpectancyOLD.txt", sep = ""))
colnames(survivingPopOld) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

avgDataframe <- aggregate(survivingPopOld$expectedAgeAtDeath, list(survivingPopOld$ageOfParent), median) 
colnames(avgDataframe) <- c("ageOfParent", "medianAgeAtDeath")
avgDataframe$minAge <- tapply(survivingPopOld$expectedAgeAtDeath, survivingPopOld$ageOfParent, min)
avgDataframe$maxAge <- tapply(survivingPopOld$expectedAgeAtDeath, survivingPopOld$ageOfParent, max)

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

###### new model ##########
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancyFirst.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancySecond.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancyThird.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancyFourth.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancyBigPop.txt", sep = ""))
colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")

myLongitudinalData <- read.table(paste(path, "correctModel/outputLETrackedIndividualsBigPop.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 



#### strength of selection = 0 ####
survivingPop <- read.table(paste(path, "correctModel/outputLifeExpectancyNoStrengthSelec.txt", sep = ""))
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

myLongitudinalData <- read.table(paste(path, "correctModel/outputLETrackedIndividualsNoStrengthSelec.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

### {0.002; 0.002; 0.0002}
survivingPop <- read.table(paste(path, "correctModel/smallProbs/outputLifeExpectancy002_0002.txt", sep = ""))
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

myLongitudinalData <- read.table(paste(path, "correctModel/smallProbs/outputLETrackedIndividuals002_0002.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

####### ONLY AGE-SPECIFIC GENES #########
# 1 = {0.0004; 0.002; 0.002}
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy1.txt", sep = ""))
# 2 = {0.0005; 0.002; 0.002}
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy2.txt", sep = ""))
# 3 = {0.0005; 0.002; 0.002} with mutbias = -0.03
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy3.txt", sep = ""))
# 4 = {0.0007; 0.002; 0.002} with mutbias = -0.01
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy4.txt", sep = ""))
# 5 = {0.0009; 0.002; 0.002} with mutbias = -0.01
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy5.txt", sep = ""))
# 6 = {0.0009; 0.002; 0.002} with mutbias = -0.02
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy6.txt", sep = ""))
# 7 = {0.001; 0.002; 0.002} with mutbias = -0.02
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy7.txt", sep = ""))
# 8 = {0.0015; 0.002; 0.002} with mutbias = -0.02
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy8.txt", sep = ""))
# 9 = {0.002; 0.002; 0.002} with mutbias = -0.02
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy9.txt", sep = ""))
# 10 = big pop {0.001; 0.002; 0.002} mutbias = -0.02
survivingPop <- read.table(paste(path, "correctModel/age_spec_genes/outputLifeExpectancy10.txt", sep = ""))

colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
sum(survivingPop$survivalProb > 0.95)

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals1.txt", sep = "")) # 1
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals2.txt", sep = "")) # 2
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals3.txt", sep = "")) # 3
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals4.txt", sep = "")) # 4
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals5.txt", sep = "")) # 5
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals6.txt", sep = "")) # 6
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals7.txt", sep = "")) # 7
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals8.txt", sep = "")) # 8
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals9.txt", sep = "")) # 9
myLongitudinalData <- read.table(paste(path, "correctModel/age_spec_genes/outputLETrackedIndividuals10.txt", sep = "")) # 10

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

######### ONLY BINARY GENES #########
# 1 = {0.002; 0.002}
survivingPop <- read.table(paste(path, "correctModel/binary_genes/outputLifeExpectancy1.txt", sep = ""))
# 2 = {0.004; 0.004}
survivingPop <- read.table(paste(path, "correctModel/binary_genes/outputLifeExpectancy2.txt", sep = ""))

colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
sum(survivingPop$survivalProb > 0.95)

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

myLongitudinalData <- read.table(paste(path, "correctModel/binary_genes/outputLETrackedIndividuals1.txt", sep = "")) # 1
myLongitudinalData <- read.table(paste(path, "correctModel/binary_genes/outputLETrackedIndividuals2.txt", sep = "")) # 2

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

######### QUALITY + AGE-SPECIFIC EFFECTS #########
# 1 = {0.0004}
survivingPop <- read.table(paste(path, "correctModel/quality_age_spec/outputLifeExpectancy.txt", sep = ""))

colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
sum(survivingPop$survivalProb > 0.95)

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

myLongitudinalData <- read.table(paste(path, "correctModel/quality_age_spec/outputLETrackedIndividuals.txt", sep = "")) # 1
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

  
######### BIG POPULATION #########
# 0.002 all mut probs 
survivingPop <- read.table(paste(path, "correctModel/bigpop/outputLifeExpectancy.txt", sep = ""))

colnames(survivingPop) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
sum(survivingPop$survivalProb > 0.95)

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

myLongitudinalData <- read.table(paste(path, "correctModel/bigpop/outputLETrackedIndividuals.txt", sep = "")) # 1
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")

######## QUALITY ##########
# only quality is on. To look at effect and find parameter. 
survivalData <- read.table(paste(path, "correctModel/survival_data/outputWithSurvivalProbs.txt", sep = "")) 
colnames(survivalData) <- c("ID", "Age", "SurvivalProb")
sub <- subset(survivalData, survivalData$ID > 90)
sub <- survivalData[0:4000,]
sub <- subset(sub, sub$ID < 10)

ggplot(data = survivalData, aes(x = Age, y = SurvivalProb, group = ID, color = factor(ID))) +
  geom_line() +
  labs(title = "",
       subtitle = "",
       x = "Age of parent",
       y = "Survival probability") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39)  +
  ylim(0, 1)


######## TEST: TO REMOVE ##########

test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Release/outputLifeExpectancy.txt")
test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Debug/outputLifeExpectancy.txt")


colnames(test) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProbStemCell", "mutationProb")
survivingPop <- test

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
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Release/outputLETrackedIndividuals.txt")
test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Debug/outputLETrackedIndividuals.txt")

colnames(test) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath")
myLongitudinalData <- test

