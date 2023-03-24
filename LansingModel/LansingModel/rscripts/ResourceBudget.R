###############################################################################
# Resource budget added. Investment in repair/ reproduction 
# only mechanism 'on' in model 
###############################################################################
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/resourceBudget/"

##### AGE DISTRIBUTION ######
ageDist <- read.table(paste(path, "outputAgeAlivePop.txt", sep = ""))
colnames(ageDist) <- c("maleAge", "femaleAge")

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
colnames(myLongitudinalData) <- c("ID", "ageOfParent", "sexOfParent", "survivalProb", "expectedAgeAtDeath", "mutationProbGametes", "mutationProbSC", "meanMutBias", "sdMutEffectSize", "mutationProbAgeGenes")
# use this for gam analysis in gam_analysis Rscript. 

##### RESOURCE DISTRIBUTION #####
parentalInvestment <- read.table(paste(path, "outputWithInvestmentDistribution.txt", sep = "")) 
colnames(parentalInvestment) <- c("ID", "age", "investmentInRepair", "mutationProbGametes", "mutationProbStemCell", "meanMutationBias", "sdMutationalEffectSize")

# get 10 IDs to examine
parentalInvestmentSub <- subset(parentalInvestment, parentalInvestment$ID >= 990)

ggplot(data = parentalInvestmentSub, aes(age, investmentInRepair, group = ID, color = factor(ID))) +
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
  ylim(0, 1)




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
