###############################################################################
# Resource budget added. Investment in repair/ reproduction 
# only mechanism 'on' in model 

# INVESTMENT IN REPAIR 2.0 = investment affects offspring quality instead of 
# number of offspring. 
###############################################################################
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/resourceBudget/"

##### AGE DISTRIBUTION ######
ageDist <- read.table(paste(path, "outputAgeAlivePop.txt", sep = "")) # mut prob = 0.001
ageDist <- read.table(paste(path, "0.2InitInvestment/outputAgeAlivePop.txt", sep = "")) # with init investment in repair = 0.2
ageDist <- read.table(paste(path, "0.8InitInvestment/outputAgeAlivePop.txt", sep = "")) # with init investment in repair = 0.8
ageDist <- read.table(paste(path, "resource-only0.004and0.02/outputAgeAlivePop.txt", sep = "")) # init = 0.5; sd = 0.02 mut prob = 0.004
ageDist <- read.table(paste(path, "resource-only0.004and0.03/outputAgeAlivePop.txt", sep = "")) # init = 0.5; sd = 0.03 mut prob = 0.004
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop.txt", sep = "")) # run time = 200.000
ageDist <- read.table(paste(path, "resource-only/outputAgeAlivePop.txt", sep = "")) # sd = 0.03 and init = 0.3 
ageDist <- read.table(paste(path, "resource-only/outputAgeAlivePopWeight0.9.txt", sep = "")) # sd = 0.03 and init = 0.3 and weight = 0.9
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop2.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.4. tEnd = 100.000
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop3.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.7. tEnd = 200.000
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop4.txt", sep = "")) #  tEnd = 200.000; weight = 0.3; init = 0.5 
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop5.txt", sep = "")) #  tEnd = 500.000; weight = 0.3; init = 0.5 
ageDist <- read.table(paste(path, "longrun/outputAgeAlivePop6.txt", sep = "")) #  tEnd = 500.000; weight = 0.3; init = 0.5. Again same thing

#### INVESTMENT 2.0 #####
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop.txt", sep = "")) 
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop2.txt", sep = "")) # init = 0.2; sd = 0.02; mut rate = 0.004
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop3.txt", sep = "")) # init = 0.8; sd = 0.02; mut rate = 0.004
ageDist <- read.table(paste(path, "investment2/longrun/outputAgeAlivePop2.txt", sep = "")) # init = 0.5; sd = 0.04; mut rate = 0.004 and tEnd = 200.000
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop4.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 10.000
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop5.txt", sep = "")) # sd = 0.02; mut rate = 0.002 and tEnd = 10.000
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop6.txt", sep = "")) # sd = 0.01; mut rate = 0.004 and tEnd = 10.000
ageDist <- read.table(paste(path, "investment2/outputAgeAlivePop7.txt", sep = "")) # sd = 0.02; mut rate = 0.004 and tEnd = 10.000
ageDist <- read.table(paste(path, "investment2/longrun/outputAgeAlivePop3.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 30.000; pop size = 10.000
ageDist <- read.table(paste(path, "investment2/longrun/outputAgeAlivePop4.txt", sep = "")) # sd = 0.02; mut rate = 0.004 and tEnd = 30.000; pop.size = 10.000; num of offspring per female = 4; num of individuals to follow = 1000; num of SC = 70. 

colnames(ageDist) <- c("maleAge", "femaleAge")

#library(reshape2)
ageDist <- melt(ageDist)
ggplot(ageDist, aes(x = value, fill = variable)) + 
  geom_density(alpha=.25) + 
  labs(title = "Density plot of the age distribution at the final time point",
       #subtitle = "Quality-only scenario.",
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
deathPop <- read.table(paste(path, "resource-only0.004and0.02/outputDeclineGameteQuality.txt", sep = "")) # sd = 0.02 
deathPop <- read.table(paste(path, "resource-only0.004and0.03/outputDeclineGameteQuality.txt", sep = "")) # sd = 0.03
deathPop <- read.table(paste(path, "resource-only/outputDeclineGameteQuality.txt", sep = "")) # sd = 0.03 and init = 0.3
deathPop <- read.table(paste(path, "resource-only/outputDeclineGameteQualityWeight0.9.txt", sep = "")) # sd = 0.03, init = 0.3 and weight = 0.9
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality2.txt", sep = "")) # sd = 0.03, init = 0.5 and weight = 0.4, tEnd= 100.000
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality3.txt", sep = "")) # sd = 0.03, init = 0.5 and weight = 0.7, tEnd= 200.000
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality4.txt", sep = "")) # sd = 0.03, init = 0.5 and weight = 0.3, tEnd= 200.000
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality5.txt", sep = "")) # sd = 0.03, init = 0.5 and weight = 0.3, tEnd= 500.000
deathPop <- read.table(paste(path, "longrun/outputDeclineGameteQuality6.txt", sep = "")) # sd = 0.03, init = 0.5 and weight = 0.3, tEnd= 500.000

### INVESTMENT 2.0 ###
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality.txt", sep = "")) # init = 0.5; sd = 0.02; rate = 0.004
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality2.txt", sep = "")) # init = 0.2; sd = 0.02; rate = 0.004
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality3.txt", sep = "")) # init = 0.8; sd = 0.02; rate = 0.004
deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality2.txt", sep = "")) # init = 0.5; sd = 0.04; rate = 0.004. tEnd = 200.000
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality4.txt", sep = "")) #sd = 0.01; rate = 0.002. tEnd = 10.000
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality5.txt", sep = "")) #sd = 0.02; rate = 0.002. tEnd = 10.000
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality6.txt", sep = "")) #sd = 0.01; rate = 0.004. tEnd = 10.000
deathPop <- read.table(paste(path, "investment2/outputDeclineGameteQuality7.txt", sep = "")) #sd = 0.02; rate = 0.004. tEnd = 10.000
deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality3.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 30.000; pop size = 10.000
deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality4.txt", sep = "")) # sd = 0.02; mut rate = 0.004 and tEnd = 30.000; pop.size = 10.000; num of offspring per female = 4; num of individuals to follow = 1000; num of SC = 70. 

colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC", "mutationProbGamete")

# only select later time point 
deathPop <- subset(deathPop, deathPop$time > 9000)

# select just the sex and age at death 
deathPop <- subset(deathPop, select = c(ageAtDeath, sex)) 
deathPop <- melt(deathPop)

ggplot(deathPop, aes(x = value, fill = sex)) +
  geom_histogram(alpha=.25, position = "identity", binwidth = 1) +
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
myLongitudinalData <- read.table(paste(path, "resource-only0.004and0.03/outputLETrackedIndividuals.txt", sep = ""))
myLongitudinalData <- read.table(paste(path, "resource-only/outputLETrackedIndividuals.txt", sep = "")) # with sd = 0.03 and init = 0.3 
myLongitudinalData <- read.table(paste(path, "resource-only/outputLETrackedIndividualsWeight0.9.txt", sep = "")) # with sd = 0.03 and init = 0.3 and weight = 0.9
myLongitudinalData <- read.table(paste(path, "longrun/outputLETrackedIndividuals2.txt", sep = "")) # with sd = 0.03 and init = 0.5 and weight = 0.4; tEnd= 100.000
myLongitudinalData <- read.table(paste(path, "longrun/outputLETrackedIndividuals3.txt", sep = "")) # with sd = 0.03 and init = 0.5 and weight = 0.7; tEnd= 200.000
myLongitudinalData <- read.table(paste(path, "longrun/outputLETrackedIndividuals4.txt", sep = "")) # with sd = 0.03 and init = 0.5 and weight = 0.3; tEnd= 200.000
myLongitudinalData <- read.table(paste(path, "longrun/outputLETrackedIndividuals5.txt", sep = "")) # with sd = 0.03 and init = 0.5 and weight = 0.3; tEnd= 500.000
myLongitudinalData <- read.table(paste(path, "longrun/outputLETrackedIndividuals6.txt", sep = "")) # with sd = 0.03 and init = 0.5 and weight = 0.3; tEnd= 500.000

### INVESTMENT 2.0 ###
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals.txt", sep = "")) 
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals2.txt", sep = "")) # init = 0.2; sd = 0.02; rate = 0.004 
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals3.txt", sep = "")) # init = 0.8; sd = 0.02; rate = 0.004 
myLongitudinalData <- read.table(paste(path, "investment2/longrun/outputLETrackedIndividuals2.txt", sep = "")) # init = 0.5; sd = 0.04; rate = 0.004; tEnd = 200.000
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals4.txt", sep = "")) # sd = 0.01; rate = 0.002; tEnd = 10.000
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals5.txt", sep = "")) # sd = 0.02; rate = 0.002; tEnd = 10.000
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals6.txt", sep = "")) # sd = 0.01; rate = 0.004; tEnd = 10.000
myLongitudinalData <- read.table(paste(path, "investment2/outputLETrackedIndividuals7.txt", sep = "")) # sd = 0.02; rate = 0.004; tEnd = 10.000
myLongitudinalData <- read.table(paste(path, "investment2/longrun/outputLETrackedIndividuals3.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 30.000; pop size = 10.000
myLongitudinalData <- read.table(paste(path, "investment2/longrun/outputLETrackedIndividuals4.txt", sep = ""))  # sd = 0.02; mut rate = 0.004 and tEnd = 30.000; pop.size = 10.000; num of offspring per female = 4; num of individuals to follow = 1000; num of SC = 70. 

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")
#colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment", "weightInvestment")

# use this for gam analysis in gam_analysis Rscript. 

##### RESOURCE DISTRIBUTION #####
parentalInvestment <- read.table(paste(path, "outputWithInvestmentDistribution.txt", sep = "")) 
parentalInvestment <- read.table(paste(path, "0.2InitInvestment/outputWithInvestment.txt", sep = "")) # with init investment in repair = 0.2 
parentalInvestment <- read.table(paste(path, "0.8InitInvestment/outputWithInvestment.txt", sep = "")) # with init investment in repair = 0.8 
parentalInvestment <- read.table(paste(path, "resource-only0.004and0.02/outputWithInvestment.txt", sep = "")) # sd = 0.02
parentalInvestment <- read.table(paste(path, "resource-only0.004and0.03/outputWithInvestment.txt", sep = "")) # sd = 0.03 
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestmentDistribution.txt", sep = "")) # a runtime of 200.000
parentalInvestment <- read.table(paste(path, "resource-only/outputWithInvestment.txt", sep = "")) # sd = 0.03 and init = 0.3
parentalInvestment <- read.table(paste(path, "resource-only/outputWithInvestmentWeight0.9.txt", sep = "")) # sd = 0.03 and init = 0.3 and weight = 0.9
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestment2.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.4; tEnd = 100.000
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestment3.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.7; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestment5.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.7; tEnd = 500.000
parentalInvestment <- read.table(paste(path, "longrun/outputWithInvestment6.txt", sep = "")) # sd = 0.03 and init = 0.5 and weight = 0.7; tEnd = 500.000

### INVESTMENT 2.0 ###
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment.txt", sep = "")) 
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment2.txt", sep = "")) # init = 0.2; sd = 0.02; rate = 0.004
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment3.txt", sep = "")) # init = 0.8; sd = 0.02; rate = 0.004
parentalInvestment <- read.table(paste(path, "investment2/longrun/outputWithInvestment2.txt", sep = "")) # init = 0.5; sd = 0.04; rate = 0.004; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment4.txt", sep = "")) # sd = 0.01; rate = 0.002; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment5.txt", sep = "")) # sd = 0.02; rate = 0.002; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment6.txt", sep = "")) # sd = 0.01; rate = 0.004; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/outputWithInvestment7.txt", sep = "")) # sd = 0.02; rate = 0.004; tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/longrun/outputWithInvestment3.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 30.000; pop size = 10.000
parentalInvestment <- read.table(paste(path, "investment2/longrun/outputWithInvestment4.txt", sep = "")) # sd = 0.02; mut rate = 0.004 and tEnd = 30.000; pop.size = 10.000; num of offspring per female = 4; num of individuals to follow = 1000; num of SC = 70.

colnames(parentalInvestment) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")

# get 10 IDs to examine
#parentalInvestmentSub <- subset(parentalInvestment, as.numeric(sapply(strsplit(parentalInvestment$ID, "_"), getElement, 1)) <= 10)
parentalInvestmentSub <- subset(parentalInvestment, parentalInvestment$ID %in% sample(parentalInvestment$ID, 15)) 

# plot ten sampled IDs
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class",
       subtitle = paste("sd = ", parentalInvestmentSub$sdInvestment[1], "; mut rate = ", parentalInvestmentSub$mutationProbInvestment[1], "; tEnd = 30.000"),
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
##### varying mutation probability with the resource 2.0 #####
parent_path <- paste(path, "investment2/mutProb/", sep = "")
##### varying mutation probability with the resource 2.0 for a longer time. tEnd = 30.000 #####
parent_path <- paste(path, "investment2/mutProb30/", sep = "")

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
  tmp <- unique(sub$ID)[sample(length(unique(sub$ID)), 40)] 
  # subset the data of these 100 IDs
  tmp2 <- subset(sub, sub$ID %in% tmp)
  # add them to the dataframe 
  sampled_data_investment <- rbind(sampled_data_investment, tmp2)
}

# standard deviation for resource 
parent_path <- paste(path, "sdVary/", sep = "")
# standard deviation for resource 2.0
parent_path <- paste(path, "investment2/sd/", sep = "")

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
  tmp <- unique(sub$ID)[sample(length(unique(sub$ID)), 50)]
  # subset the data of these 100 IDs
  tmp2 <- subset(sub, sub$ID %in% tmp)
  # add them to the dataframe 
  sampled_data_sd_investment <- rbind(sampled_data_sd_investment, tmp2)
}

sampled_data_investment_tot <- rbind(sampled_data_investment, sampled_data_sd_investment)


## varying the weight of the investment genes 
parent_path <- paste(path, "resource-only/varyWeight/", sep = "")

# look at tracked individuals data 
f <- list.files(path = parent_path, pattern = "outputLETrackedIndividuals.txt", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestmentGenes", "sdInvestmentGenes", "weightInvestment")

# plot this data using geom_smooth function 
ggplot(data = found, aes(x = ageOfParent, y = expectedAgeAtDeath, color = as.factor(weightInvestment))) +
  geom_smooth() + 
  labs(title = "Expected age at death of offspring over parental ages",
       x = "Age of parent",
       y = "Expected age of death offspring",
       color = "weight investment genes") +
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15)) +
  scale_x_continuous(labels = as.character(0:39), breaks = 0:39) 

# sample from data 
sampled_data_weight_investment <- c()
for (x in unique(found$weightInvestment)){
  # get subset of this mutation probability
  sub <- subset(found, found$weightInvestment == x)
  # get 100 randomly sampled and unique IDs 
  tmp <- unique(sub$ID)[sample(length(unique(sub$ID)), 50)]
  # subset the data of these 100 IDs
  tmp2 <- subset(sub, sub$ID %in% tmp)
  # add them to the dataframe 
  sampled_data_weight_investment <- rbind(sampled_data_weight_investment, tmp2)
}


###############################################################################
# Sampling ten individuals 
###############################################################################
parent_path <- paste(path, "mutationProbVary/", sep = "")
parent_path <- paste(path, "investment2/mutProb/", sep = "") # resource-only 2.0 
parent_path <- paste(path, "investment2/mutProb30/", sep = "") # run for 30.000 time steps 

f <- list.files(path = parent_path, pattern = "outputWithInvestmentDistribution.txt", recursive = T)
f <- list.files(path = parent_path, pattern = "outputWithInvestment.txt", recursive = T)

found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}

colnames(found) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestmentGenes")

# get 10 IDs to examine
parentalInvestmentSub <- subset(found, found$ID %in% sample(found$ID, 10)) 

# plot ten sampled IDs
ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
  #geom_smooth(se = F) +
  geom_line() +
  labs(title = "Looking at parental investment in repair per age class",
       subtitle = "Where every panel is a muation probability and tEnd = 30.000",
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
# resource-only
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputAgeAlivePop.txt", sep = ""))
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputAgeAlivePop2.txt", sep = "")) # pop.size = 10.000; tEnd = 50.000
mechanisms <- "resource-only"

# resource and damage
alivePop <- read.table(paste(path, "resource_damage/outputAgeAlivePop.txt", sep = ""))
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputAgeAlivePop.txt", sep = "")) # resource 2.0
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputAgeAlivePop2.txt", sep = "")) # resource 2.0 and tEnd = 200.000
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputAgeAlivePop3.txt", sep = "")) # resource 2.0 with higher mutation rate for damage accumulation (0.004)
mechanisms <- "resource + damage"

# resource and quality
alivePop <- read.table(paste(path, "resource-quality/outputAgeAlivePop.txt", sep = ""))
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputAgeAlivePop.txt", sep = "")) # resource 2.0
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputAgeAlivePop2.txt", sep = "")) # resource 2.0 with tEnd = 200.000
mechanisms <- "resource + quality"

# resource + quality + damage
alivePop <- read.table(paste(path, "resource-damage-quality/outputAgeAlivePop.txt", sep = ""))
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputAgeAlivePop.txt", sep = "")) # resource 2.0 
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputAgeAlivePop2.txt", sep = "")) # resource 2.0 with tEnd = 200.000
alivePop <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputAgeAlivePop3.txt", sep = "")) # resource 2.0
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
# resource-only
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputDeclineGameteQuality.txt", sep = "")) 
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputDeclineGameteQuality2.txt", sep = "")) # pop.size = 10.000; tEnd = 50.000
mechanisms <- "resource"

deathPop <- read.table(paste(path, "resource_damage/outputDeclineGameteQuality.txt", sep = "")) 
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputDeclineGameteQuality.txt", sep = "")) # resource 2.0
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputDeclineGameteQuality2.txt", sep = "")) # resource 2.0 with tEnd = 200.000
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputDeclineGameteQuality3.txt", sep = "")) # resource 2.0; mut rate gametes/SCs = 0.004
mechanisms <- "resource + damage"

deathPop <- read.table(paste(path, "resource-quality/outputDeclineGameteQuality.txt", sep = "")) 
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputDeclineGameteQuality.txt", sep = "")) # resource 2.0
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputDeclineGameteQuality2.txt", sep = "")) # resource 2.0 with tEnd = 200.000
mechanisms <- "resource + quality"

deathPop <- read.table(paste(path, "resource-damage-quality/outputDeclineGameteQuality.txt", sep = "")) 
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputDeclineGameteQuality.txt", sep = "")) # resource 2.0 
deathPop <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputDeclineGameteQuality2.txt", sep = "")) # resource 2.0 with tEnd = 200.000
mechanisms <- "resource + quality + damage"

colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC", "mutationProbGamete", "mutationProbInvestment", "sdInvestment", "investmentInRepair")

# only select later time point 
deathPop <- subset(deathPop, deathPop$time > 9000)

# select just the sex and age at death 
deathPop <- subset(deathPop, select = c(ageAtDeath, sex)) 
deathPop <- melt(deathPop)

ggplot(deathPop, aes(x = value, fill = sex)) +
  geom_histogram(alpha=.25, position = "identity", binwidth = 1) +
  labs(title = "Histogram of age at death distribution",
       subtitle = paste("At time > 9000. ", mechanisms), 
       x = "Age at death") + 
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))
# right-skewed

### tracked individuals ###
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputLETrackedIndividuals.txt", sep = "")) # resource
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputLETrackedIndividuals2.txt", sep = "")) # resource; pop.size = 10.000; tEnd = 50.000

myLongitudinalData <- read.table(paste(path, "resource_damage/outputLETrackedIndividuals.txt", sep = "")) # resource + damage
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputLETrackedIndividuals.txt", sep = "")) # resource 2.0 + damage
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputLETrackedIndividuals2.txt", sep = "")) # resource 2.0 + damage with tEnd = 200.000
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputLETrackedIndividuals3.txt", sep = "")) # resource 2.0 + damage; mut rate gametes/SCs = 0.004

myLongitudinalData <- read.table(paste(path, "resource-quality/outputLETrackedIndividuals.txt", sep = "")) # resource + quality 
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputLETrackedIndividuals.txt", sep = "")) # resource 2.0 + quality 
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputLETrackedIndividuals2.txt", sep = "")) # resource 2.0 + quality with tEnd = 200.000

myLongitudinalData <- read.table(paste(path, "resource-damage-quality/outputLETrackedIndividuals.txt", sep = "")) # resource + quality + damage
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputLETrackedIndividuals.txt", sep = "")) # resource + quality + damage
myLongitudinalData <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputLETrackedIndividuals2.txt", sep = "")) # resource + quality + damage with tEnd = 200.000

colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes", "mutationProbInvestment", "sdInvestment")
# use this for gam analysis in gam_analysis Rscript. 

ggplot(myLongitudinalData, aes(ageOfParent, expectedAgeAtDeath)) + geom_smooth(method = "gam")

### Resource distribution ###
parentalInvestment <- read.table(paste(path, "resource-only-baseline/outputWithInvestment.txt", sep = "")) # resource 
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-only-baseline/outputWithInvestment2.txt", sep = "")) # resource; pop.size = 10.000; tEnd = 50.000
mechanisms <- "resource"

parentalInvestment <- read.table(paste(path, "resource_damage/outputWithInvestment.txt", sep = "")) # resource + damage
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputWithInvestment.txt", sep = "")) # resource 2.0 + damage
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputWithInvestment2.txt", sep = "")) # resource 2.0 + damage with tEnd = 200.000
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-damage/outputWithInvestment3.txt", sep = "")) # resource 2.0 + damage; mut rate gametes/SCs = 0.004
mechanisms <- "resource + damage" 

parentalInvestment <- read.table(paste(path, "resource-quality/outputWithInvestment.txt", sep = "")) # resource + quality
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputWithInvestment.txt", sep = "")) # resource 2.0 + quality
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-quality/outputWithInvestment2.txt", sep = "")) # resource 2.0 + quality with tEnd = 200.000
mechanisms <- "resource + quality"

parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputWithInvestment.txt", sep = "")) # resource 2.0 + quality + damage
parentalInvestment <- read.table(paste(path, "investment2/CombiningResource/resource-quality-damage/outputWithInvestment2.txt", sep = "")) # resource 2.0 + quality + damage with tEnd = 200.000
mechanisms <- "resource + quality + damage"

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


###############################################################################
# INVESTMENT IN REPAIR OVER TIME 
###############################################################################
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/resourceBudget/"

deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality2.txt", sep = "")) # sd = 0.04, init = 0.5, rate= 0.004, tEnd= 200.000
deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality3.txt", sep = "")) # sd = 0.01; mut rate = 0.002 and tEnd = 30.000; pop size = 10.000
deathPop <- read.table(paste(path, "investment2/longrun/outputDeclineGameteQuality4.txt", sep = "")) # sd = 0.02; mut rate = 0.004 and tEnd = 30.000; pop.size = 10.000; num of offspring per female = 4; num of individuals to follow = 1000; num of SC = 70.

colnames(deathPop) <- c("time", "ageAtDeath", "sex", "ageOfMother", "ageOfFather", "survivalProb", "mutationProbSC", "mutationProbGamete", "mutationProbInvestmentGenes", "sdInvestment", "investmentInRepair")

deathPopSub <- deathPop[deathPop$ageAtDeath %in% c(0, 10, 20, 30, 32, 35, 37, 39),]
deathPopSub <- deathPop[deathPop$ageAtDeath %in% c(25, 28, 30, 32, 35, 37),]

deathPopSub <- deathPop[deathPop$ageAtDeath %in% c(0, 10, 20, 30, 39),]
deathPopSub <- deathPopSub[deathPopSub$time %in% seq(0, 200000, 100),]

# use geom_smooth
ggplot(data = deathPopSub, aes(time, investmentInRepair, colour = factor(ageAtDeath), group = factor(ageAtDeath))) + 
  #geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) + 
  labs(title = "Investment in repair over time using smooth from ggplot") + 
  theme(text = element_text(size = 20))

# fit gam model and plot this
mod <- bam(investmentInRepair ~ s(time, k = 10, by = as.factor(ageAtDeath)),
           data = deathPopSub)
z <- predict(mod)
deathPopSub$z <- z

#library(cowplot)
ggplot(deathPopSub,aes(x=time,y=z,color=as.factor(ageAtDeath))) +
  labs(title = "resource-only scenario") +
  theme_cowplot() +
  geom_line() +
  background_grid(major = "xy")

# calculate average investment per age class and plot
averageInvestmentPerTimeStep <- aggregate(deathPopSub$investmentInRepair, list(deathPopSub$ageAtDeath, deathPopSub$time), mean)
colnames(averageInvestmentPerTimeStep) <- c("ageAtDeath", "time", "averageInvestmentInRepair")
ggplot(data = averageInvestmentPerTimeStep, aes(time, averageInvestmentInRepair, colour = factor(ageAtDeath), group = factor(ageAtDeath))) + 
  geom_point() + 
  geom_smooth() + 
  labs(title = "Average investment in repair over time") + 
  theme(text = element_text(size = 20))


###############################################################################
  
###############################################################################
