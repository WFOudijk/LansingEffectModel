###############################################################################
# Resource budget added. Investment in repair/ reproduction 
# only mechanism 'on' in model 
###############################################################################
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/resourceBudget/"

##### AGE DISTRIBUTION ######
ageDist <- read.table(paste(path, "outputAgeAlivePop.txt", sep = ""))
ageDist <- read.table(paste(path, "0.2InitInvestment/outputAgeAlivePop.txt", sep = "")) # with init investment in repair = 0.2
ageDist <- read.table(paste(path, "0.8InitInvestment/outputAgeAlivePop.txt", sep = "")) # with init investment in repair = 0.8
ageDist <- read.table(paste(path, "resource-only0.004and0.02/outputAgeAlivePop.txt", sep = "")) # init = 0.5; sd = 0.02 mut prob = 0.004

# run time = 200.000
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop.txt", sep = ""))
colnames(ageDist) <- c("maleAge", "femaleAge")

#library(reshape2)
ageDist <- melt(ageDist)
ggplot(ageDist, aes(x = value, fill = variable)) + 
  geom_density(alpha=.25) + 
  labs(title = "Density plot of the age distribution at the final time point",
       subtitle = "Quality-only scenario.",
       x = "age") + 
  theme_big + 
  scale_fill_manual(values = c("darkblue", "pink")) + 
  xlim(0, 40)
# right-skewed 

##### AGE AT DEATH PLOT ######
deathPop <- read.table(paste(path, "outputDeclineGameteQuality.txt", sep = "")) # mut prob of investment genes = 0.001
deathPop <- read.table(paste(path, "0.2InitInvestment/outputDeclineGameteQuality.txt", sep = "")) # init genes = 0.2
deathPop <- read.table(paste(path, "0.8InitInvestment/outputDeclineGameteQuality.txt", sep = "")) # init genes = 0.8
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality.txt", sep = "")) # mut prob of investment genes = 0.001
deathPop <- read.table(paste(path, "resource-only0.004and0.02/outputDeclineGameteQuality.txt", sep = ""))

colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC", "mutationProbGamete")

# only select later time point 
deathPop <- subset(deathPop, deathPop$time > 9000)

# select just the sex and age at death 
deathPop <- subset(deathPop, select = c(ageAtDeath, sex)) 
deathPop <- melt(deathPop)

ggplot(deathPop, aes(x = value, fill = sex)) +
  geom_histogram(alpha=.25, position = "identity") +
  labs(title = "Histogram of age at death distribution",
       subtitle = "At time > 9000. Quality-only scenario.", 
       x = "Age at death") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
# right-skewed 

##### TRACKED INDIVIDUALS #####
myLongitudinalData <- read.table(paste(path, "outputLETrackedIndividuals.txt", sep = ""))
myLongitudinalData <- read.table(paste(path, "0.2InitInvestment/outputLETrackedIndividuals.txt", sep = ""))
myLongitudinalData <- read.table(paste(path, "0.8InitInvestment/outputLETrackedIndividuals.txt", sep = ""))
myLongitudinalData <- read.table(paste(path, "resource-only0.004and0.02/outputLETrackedIndividuals.txt", sep = ""))
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")
# use this for gam analysis in gam_analysis Rscript. 

##### RESOURCE DISTRIBUTION #####
parentalInvestment <- read.table(paste(path, "outputWithInvestmentDistribution.txt", sep = "")) 
parentalInvestment <- read.table(paste(path, "0.2InitInvestment/outputWithInvestment.txt", sep = "")) # with init investment in repair = 0.2 
parentalInvestment <- read.table(paste(path, "0.8InitInvestment/outputWithInvestment.txt", sep = "")) # with init investment in repair = 0.8 
parentalInvestment <- read.table(paste(path, "resource-only0.004and0.02/outputWithInvestment.txt", sep = "")) # 

# a runtime of 200.000
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestmentDistribution.txt", sep = "")) 
colnames(parentalInvestment) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")

# get 10 IDs to examine
#parentalInvestmentSub <- subset(parentalInvestment, as.numeric(sapply(strsplit(parentalInvestment$ID, "_"), getElement, 1)) <= 10)
parentalInvestmentSub <- subset(parentalInvestment, parentalInvestment$ID %in% sample(parentalInvestment$ID, 10)) 

# plot ten sampled IDs
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class with init = 0.8",
       #subtitle = "With differing mean mutation bias",
       x = "Age",
       y = "Age-specific investment in repair") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 7)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39)  +
  ylim(0, 1)

# population average 
parentalInvestmentSubPopAvg <- aggregate(parentalInvestment$investmentInRepair, list(parentalInvestment$age), mean)
parentalInvestmentSubPopAvg$min <- tapply(parentalInvestment$investmentInRepair, list(parentalInvestment$age), min)
parentalInvestmentSubPopAvg$max <- tapply(parentalInvestment$investmentInRepair, list(parentalInvestment$age), max)
colnames(parentalInvestmentSubPopAvg) <- c("age", "investmentInRepair", "min", "max")

# plot population average 
ggplot(data = parentalInvestmentSubPopAvg, aes( x = age, y = investmentInRepair)) +
  geom_line() +
  ylim(0,1) + 
  geom_linerange(aes(ymin = min, ymax = max),
                 linetype = "dotted") +
  labs(title = "Population average",
       subtitle = "with dotted lines being minimum and maximum",
       x = "Age",
       y = "Investment in repair") +
  theme(text = element_text(size = 20))

###############################################################################
# parameter simulations
###############################################################################
# mutation probability 
parent_path <- paste(path, "mutationProbVary/", sep = "")

# look at tracked individuals data 
f <- list.files(path = parent_path, pattern = "outputLETrackedIndividuals.txt", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestmentGenes", "sdInvestmentGenes")

# plot this data using geom_smooth function 
ggplot(data = found, aes(x = ageOfParent, y = expectedAgeAtDeath, color = as.factor(mutationProbInvestmentGenes))) +
  geom_smooth() + 
  labs(title = "Expected age at death of offspring over parental ages",
       x = "Age of parent",
       y = "Expected age of death offspring",
       color = "Mutation probability") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

# sample from the data
sampled_data_investment <- c()
for (x in unique(found$mutationProbInvestmentGenes)){
  # get subset of this mutation probability
  sub <- subset(found, found$mutationProbInvestmentGenes == x)
  # get 100 randomly sampled and unique IDs 
  tmp <- unique(sub$ID)[sample(length(unique(sub$ID)), 25)]
  # subset the data of these 100 IDs
  tmp2 <- subset(sub, sub$ID %in% tmp)
  # add them to the dataframe 
  sampled_data_investment <- rbind(sampled_data_investment, tmp2)
}

# standard deviation 
parent_path <- paste(path, "sdVary/", sep = "")

# look at tracked individuals data 
f <- list.files(path = parent_path, pattern = "outputLETrackedIndividuals.txt", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestmentGenes", "sdInvestmentGenes")

# plot this data using geom_smooth function 
ggplot(data = found, aes(x = ageOfParent, y = expectedAgeAtDeath, color = as.factor(sdInvestmentGenes))) +
  geom_smooth() + 
  labs(title = "Expected age at death of offspring over parental ages",
       x = "Age of parent",
       y = "Expected age of death offspring",
       color = "Mutational effect size") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

# sample from the data
sampled_data_sd_investment <- c()
for (x in unique(found$sdInvestmentGenes)){
  # get subset of this mutation probability
  sub <- subset(found, found$sdInvestmentGenes == x)
  # get 100 randomly sampled and unique IDs 
  tmp <- unique(sub$ID)[sample(length(unique(sub$ID)), 25)]
  # subset the data of these 100 IDs
  tmp2 <- subset(sub, sub$ID %in% tmp)
  # add them to the dataframe 
  sampled_data_sd_investment <- rbind(sampled_data_sd_investment, tmp2)
}

sampled_data_investment_tot <- rbind(sampled_data_investment, sampled_data_sd_investment)


###############################################################################
# Sampling ten individuals 
###############################################################################
parent_path <- paste(path, "mutationProbVary/", sep = "")

f <- list.files(path = parent_path, pattern = "outputWithInvestmentDistribution.txt", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}

colnames(found) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestmentGenes")

# resource-only scenario with initial investment in repair set to 0.2. 
parentalInvestment <- read.table(paste(path, "0.2InitInvestment/outputWithInvestmentDistribution.txt", sep = ""))

# get 10 IDs to examine
#parentalInvestmentSub <- subset(found, as.numeric(sapply(strsplit(found$ID, "_"), getElement, 1)) <= 10)
parentalInvestmentSub <- subset(found, found$ID %in% sample(found$ID, 10)) 

# plot ten sampled IDs
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class",
       #subtitle = "With differing mean mutation bias",
       x = "Age",
       y = "Age-specific investment in repair") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 7)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39)  +
  ylim(0, 1) + 
  facet_wrap(as.factor(parentalInvestmentSub$mutationProbInvestment))

# population average 
parentalInvestmentSubPopAvg <- aggregate(found$investmentInRepair, list(found$age, found$mutationProbInvestment), mean)
#parentalInvestmentSubPopAvg$min <- tapply(found$investmentInRepair, list(found$age, found$sdInvestmentGenes), min)
#parentalInvestmentSubPopAvg$max <- tapply(found$investmentInRepair, list(found$age, found$sdInvestmentGenes), max)
#colnames(parentalInvestmentSubPopAvg) <- c("age", "investmentInRepair", "min", "max")
colnames(parentalInvestmentSubPopAvg) <- c("age", "mutationProb", "investmentInRepair")

# plot population average 
ggplot(data = parentalInvestmentSubPopAvg, aes( x = age, y = investmentInRepair)) +
  geom_line() +
  ylim(0,1) + 
  #geom_linerange(aes(ymin = min, ymax = max),
  #               linetype = "dotted") +
  labs(title = "Population average",
       #subtitle = "with dotted lines being minimum and maximum",
       x = "Age",
       y = "Investment in repair") +
  theme(text = element_text(size = 20)) +
  facet_wrap(parentalInvestmentSubPopAvg$mutationProb)

# varying sd and looking at plots 
parent_path <- paste(path, "sdVary/", sep = "")

# look at tracked individuals data 
f <- list.files(path = parent_path, pattern = "outputWithInvestmentDistribution.txt", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}

colnames(found) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestmentGenes")

# get 10 IDs to examine
#parentalInvestmentSub <- subset(found, as.numeric(sapply(strsplit(found$ID, "_"), getElement, 1)) <= 10)
parentalInvestmentSub <- subset(found, found$ID %in% sample(found$ID, 10)) 

# plot ten sampled IDs
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class",
       #subtitle = "With differing mean mutation bias",
       x = "Age",
       y = "Age-specific investment in repair") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 7)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39)  +
  ylim(0, 1) + 
  facet_wrap(as.factor(parentalInvestmentSub$sdInvestmentGenes))

# population average 
parentalInvestmentSubPopAvg <- aggregate(found$investmentInRepair, list(found$age, found$sdInvestmentGenes), mean)
colnames(parentalInvestmentSubPopAvg) <- c("age", "sd", "investmentInRepair")

# plot population average 
ggplot(data = parentalInvestmentSubPopAvg, aes( x = age, y = investmentInRepair)) +
  geom_line() +
  ylim(0,1) + 
  #geom_linerange(aes(ymin = min, ymax = max),
  #               linetype = "dotted") +
  labs(title = "Population average",
       x = "Age",
       y = "Investment in repair") +
  theme(text = element_text(size = 20)) +
  facet_wrap(parentalInvestmentSubPopAvg$sd)


###############################################################################
# COMBINING RESOURCE DISTRIBUTION WITH THE OTHER MECHANISMS
###############################################################################
### AGE DISTRIBUTION ###
# resource and damage
alivePop <- read.table(paste(path, "resource_damage/outputAgeAlivePop.txt", sep = ""))
mechanisms <- "resource + damage"

# resource and quality
alivePop <- read.table(paste(path, "resource-quality/outputAgeAlivePop.txt", sep = ""))
mechanisms <- "resource + quality"

# resource + quality + damage
alivePop <- read.table(paste(path, "resource-damage-quality/outputAgeAlivePop.txt", sep = ""))
mechanisms <- "resource + quality + damage"

colnames(alivePop) <- c("maleAge", "femaleAge")

#library(reshape2)
alivePop <- melt(alivePop)
ggplot(alivePop, aes(x = value, fill = variable)) + 
  geom_density(alpha=.25) + 
  labs(title = "Density plot of the age distribution at the final time point",
       subtitle = mechanisms, 
       x = "age") + 
  theme_big + 
  scale_fill_manual(values = c("darkblue", "pink")) + 
  xlim(0, 40)

### age at death ### 
deathPop <- read.table(paste(path, "resource_damage/outputDeclineGameteQuality.txt", sep = "")) 
mechanisms <- "resource + damage"

deathPop <- read.table(paste(path, "resource-quality/outputDeclineGameteQuality.txt", sep = "")) 
mechanisms <- "resource + quality"

deathPop <- read.table(paste(path, "resource-damage-quality/outputDeclineGameteQuality.txt", sep = "")) 
mechanisms <- "resource + quality + damage"

colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC", "mutationProbGamete")

# only select later time point 
deathPop <- subset(deathPop, deathPop$time > 9000)

# select just the sex and age at death 
deathPop <- subset(deathPop, select = c(ageAtDeath, sex)) 
deathPop <- melt(deathPop)

ggplot(deathPop, aes(x = value, fill = sex)) +
  geom_histogram(alpha=.25, position = "identity") +
  labs(title = "Histogram of age at death distribution",
       subtitle = paste("At time > 9000. ", mechanisms), 
       x = "Age at death") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
# right-skewed

### tracked individuals ###
myLongitudinalData <- read.table(paste(path, "resource_damage/outputLETrackedIndividuals.txt", sep = "")) # resource + damage
myLongitudinalData <- read.table(paste(path, "resource-quality/outputLETrackedIndividuals.txt", sep = "")) # resource + quality 
myLongitudinalData <- read.table(paste(path, "resource-damage-quality/outputLETrackedIndividuals.txt", sep = "")) # resource + quality + damage

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")
# use this for gam analysis in gam_analysis Rscript. 

### Resource distribution ###
parentalInvestment <- read.table(paste(path, "resource_damage/outputWithInvestment.txt", sep = "")) # resource + damage
mechanisms <- "resource + damage" 

parentalInvestment <- read.table(paste(path, "resource-quality/outputWithInvestment.txt", sep = "")) # resource + quality
mechanisms <- "resource + quality"

colnames(parentalInvestment) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")

# get 10 IDs to examine
parentalInvestmentSub <- subset(parentalInvestment, parentalInvestment$ID %in% sample(parentalInvestment$ID, 10)) 

# plot ten sampled IDs
# mut prob = 0.004 and sd = 0.02
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class",
       subtitle = mechanisms,
       x = "Age",
       y = "Age-specific investment in repair") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 7)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39)  +
  ylim(0, 1)

# population average 
parentalInvestmentSubPopAvg <- aggregate(parentalInvestment$investmentInRepair, list(parentalInvestment$age), mean)
parentalInvestmentSubPopAvg$min <- tapply(parentalInvestment$investmentInRepair, list(parentalInvestment$age), min)
parentalInvestmentSubPopAvg$max <- tapply(parentalInvestment$investmentInRepair, list(parentalInvestment$age), max)
colnames(parentalInvestmentSubPopAvg) <- c("age", "investmentInRepair", "min", "max")

# plot population average 
ggplot(data = parentalInvestmentSubPopAvg, aes( x = age, y = investmentInRepair)) +
  geom_line() +
  ylim(0,1) + 
  geom_linerange(aes(ymin = min, ymax = max),
                 linetype = "dotted") +
  labs(title = "Population average",
       subtitle = paste("with dotted lines being minimum and maximum. ", mechanisms),
       x = "Age",
       y = "Investment in repair") +
  theme(text = element_text(size = 20))



######## TEST: TO REMOVE ##########

test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Release/outputLifeExpectancy.txt")
test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Debug/outputLifeExpectancy.txt")
test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Release/outputWithParentalQuality.txt")
test <- read.table("/Users/willemijnoudijk/Library/Developer/Xcode/DerivedData/LansingModel-bfhrejexgadtgjexzmzzbuoxivtr/Build/Products/Release/outputLETrackedIndividuals.txt")


colnames(test) <-  c("age", "expectedAgeAtDeath", "ageOfParent", "sexOfParent", "survivalProb", "mutationProb")
colnames(test) <- c("ID", "age", "survivalProb", "mutationProb") 
colnames(test) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProb")
colnames(test) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes")



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

# geom_smooth for smooth. Use group and color as differing parameter values. Method = gam. 
plot(funct(seq(0,1,0.01)), x = seq(0,1,0.01), ylab = "number of offspring", xlab = "investment in reproduction", main = expression(paste("Number of offspring calculated by y = 5", italic(x)^(0.5+italic(x)))))

a <- seq(0,1,0.1) # a = investment in repair 
c3 <- 0.3
1 - c3 * (1-a)^2
