###############################################################################


# MANUSCRIPT R CODE 


###############################################################################

###############################################################################
# Supplementary materials: looking at parameter simulations longitudinal. 
# with ten replicates for the confidence intervals. 
###############################################################################
# standard format for images. 
# plot the normalized data grouped by mutation prob
ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                    group=mutationProbAgeSpecificGenes, 
                                    colour = mutationProbAgeSpecificGenes, 
                                    shape = mutationProbAgeSpecificGenes)) +
  geom_line() +
  labs( x = "Normalized parental age",
        y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) +
  scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
  scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.84, 0.86),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_y_continuous(breaks = seq(0,1.5,0.2),
                     limits = c(0,1.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1))

# save param sim (ps) image. 
ggsave("baseline_ps.pdf", width = 12, height = 8)


path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/manuscript/"

# BASELINE 

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
  return(tmp)
}

# reading all data and saving both the data as the gam models. 
parent_path <- paste0(path, "baseline/short_time/") # pop.size = 1000 and tEnd = 10.000
# pop.size = 10.000 and tEnd = 100.000: 
name = "setting_zero_short_time"
parent_path <- paste0(path, "baseline/setting_one/") # mean = -0.02; sd = 0.02
name = "setting_one"
parent_path <- paste0(path, "baseline/setting_two/") # mean = 0; sd = 0.02
name = "setting_two"
parent_path <- paste0(path, "baseline/setting_three/") # mean = 0; sd = 0.08
name = "setting_three"
parent_path <- paste0(path, "baseline/figure_s1/") # final parameter sim for in paper
name = "Figure_s1"

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allDataBaseline <- c()
allDataCompleteBaseline <- c()
modelsBaseline <- list()
counter = 0
# paste all data together
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate and scenario as column 
  local_data$rep <- rep
  #local_data$scenario <- scenario
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$mutationProbAgeSpecificGenes[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(mutationProbAgeSpecificGenes = factor(mutationProbAgeSpecificGenes))
  
  allDataCompleteBaseline <- rbind(allDataCompleteBaseline, local_data) # for the maternal age distribution plot 
  
  # perform gam
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  to_sample <- 100 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
    #to_sample <- length(unique(local_data$ID))
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allDataBaseline <- rbind(allDataBaseline, local_data)
    
  # GAM 
  d <- local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  #d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  modelsBaseline[[i]] <- mod
  # rename model to be unique for replicate and scenario
  names(modelsBaseline)[i] <- paste0("model_", rep)
}
print(paste0("number of skipped replicates: ", counter))

# save the R data just in case. 
saveRDS(allDataBaseline, file = paste0("allDataBaseline_", name, ".rds"))
saveRDS(modelsBaseline, file = paste0("modelsBaseline_", name, ".rds"))
saveRDS(allDataCompleteBaseline, file = paste0("allDataCompleteBaseline_", name, ".rds"))

allDataBaseline$rep <- factor(allDataBaseline$rep)
logist <- function(x) 40/(1+exp(-x))

allDataBaseline <- loadRDS(paste0("allDataBaseline_", name, ".rds"))
modelsBaseline <- loadRDS(paste0("modelsBaseline_", name, ".rds"))
allDataCompleteBaseline <- loadRDS(paste0("allDataCompleteBaseline_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()
# go per scenario through the data 
for (x in levels(allDataBaseline$mutationProbAgeSpecificGenes)){
  # make subset of scenario 
  sub <- allDataBaseline[allDataBaseline$mutationProbAgeSpecificGenes == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          mutationProbAgeSpecificGenes = sub$mutationProbAgeSpecificGenes[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i)
    index <- which(names(modelsBaseline) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(modelsBaseline[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 4:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0("predDataTotalNormalized_", name, ".rds"))

predDataTotalNormalized_final <- predDataTotalNormalized
predDataTotal_final <- predDataTotal

# load from saved extra in environment 
predDataTotalNormalized <- predDataTotalNormalized_one
predDataTotal <- predDataTotal_one
name = "setting_one"
predDataTotalNormalized <- predDataTotalNormalized_two
predDataTotal <- predDataTotal_two
name = "setting_two"
predDataTotalNormalized <- predDataTotalNormalized_three
predDataTotal <- predDataTotal_three
name = "setting_three"
# TO RECREATE FIGURE 1 
predDataTotalNormalized <- predDataTotalNormalized_final
predDataTotal <- predDataTotal_final
name = "figure_s1"

# load from RDS
#predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".rds"))

# plot the normalized data grouped by mutation prob
ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                    group=mutationProbAgeSpecificGenes, 
                                    colour = mutationProbAgeSpecificGenes, 
                                    shape = mutationProbAgeSpecificGenes)) +
  geom_line() +
  labs( x = "Normalized parental age",
        y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) +
  scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
  scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.84, 0.86),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_y_continuous(breaks = seq(0,1.5,0.2),
                     limits = c(0,1.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1))

# save param sim (ps) image. 
ggsave("baseline_ps.pdf", width = 12, height = 8)


# plot the normalized data grouped by mutation prob
ggplot(predDataTotal, aes(maternalAge, mean, 
                                  group=mutationProbAgeSpecificGenes, 
                                  colour = mutationProbAgeSpecificGenes, 
                                  shape = mutationProbAgeSpecificGenes)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(predDataTotal$mutationProbAgeSpecificGenes)) +
  geom_point() +
  labs(title = name, 
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() 
  #geom_ribbon(data = predDataTotal,
  #          aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
  #          alpha = 0.2, colour = NA) 

###############################################################################

# GAMETE DAMAGE ACCUMULATION  

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
  return(tmp)
}

# reading all data and saving both the data as the gam models. 
parent_path <- paste0(path, "gamete_damage/figure_s2/") # final parameter sim for in paper
name = "Figure_s2"
parent_path <- paste0(path, "gamete_damage/gamete_small/") # gamete damage accumulation with smaller steps between the parameters
name <- "Gamete_small_steps"

# GAMETE + BASELINE
parent_path <- paste0(path, "gamete_baseline/defaults/") # gamete damage accumulation + baseline
name <- "Gamete_baseline"
parent_path <- paste0(path, "gamete_baseline/no_bias/") # no mutation bias
name = "no_bias"


f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allDataBaseline <- c()
allDataCompleteBaseline <- c()
modelsBaseline <- list()
counter = 0
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract scenario from file name 
  splitted_path <- strsplit(file_name, "/")
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate and scenario as column 
  local_data$rep <- rep
  #local_data$scenario <- scenario
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$mutationProb[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(mutationProb = factor(mutationProb))
  
  allDataCompleteBaseline <- rbind(allDataCompleteBaseline, local_data) # for the maternal age distribution plot 

  # perform gam
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  to_sample <- 50 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
    #to_sample <- length(unique(local_data$ID))
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allDataBaseline <- rbind(allDataBaseline, local_data)
  
  # GAM 
  d <- local_data
  # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
  #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  #d2 <- d2 %>% mutate(scenario = factor(scenario)) 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  modelsBaseline[[i]] <- mod
  # rename model to be unique for replicate and scenario
  names(modelsBaseline)[i] <- paste0("model_", rep, "_", d2$mutationProb)
}
Sys.time()
print(paste0("number of skipped replicates: ", counter))

# save the R data just in case. 
saveRDS(allDataBaseline, file = paste0("allDataBaseline_", name, ".rds"))
saveRDS(modelsBaseline, file = paste0("modelsBaseline_", name, ".rds"))
saveRDS(allDataCompleteBaseline, file = paste0("allDataCompleteBaseline_", name, ".rds"))

allDataBaseline$rep <- factor(allDataBaseline$rep)
logist <- function(x) 40/(1+exp(-x))

allDataBaseline <- loadRDS(paste0("allDataBaseline_", name, ".rds"))
modelsBaseline <- loadRDS(paste0("modelsBaseline_", name, ".rds"))
allDataCompleteBaseline <- loadRDS(paste0("allDataCompleteBaseline_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()
# go per scenario through the data 
for (x in levels(allDataBaseline$mutationProb)){
  # make subset of scenario 
  sub <- allDataBaseline[allDataBaseline$mutationProb == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          mutationProb = sub$mutationProb[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$mutationProb[1])
    index <- which(names(modelsBaseline) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(modelsBaseline[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 4:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0("predDataTotalNormalized_", name, ".rds"))

predDataTotalNormalized_gamete_baseline_no_bias <- predDataTotalNormalized
predDataTotal_gamete_baseline_no_bias <- predDataTotal

# load from saved extra in environment 

# gamete damage only 
predDataTotalNormalized <- predDataTotalNormalized_figure_s2
predDataTotal <- predDataTotal_figure_s2
name = "figure_s2"

# gamete + baseline broad
predDataTotalNormalized <- predDataTotalNormalized_gamete_baseline
predDataTotal <- predDataTotal_gamete_baseline
name = "gamete_baseline"
predDataTotalNormalized <- predDataTotalNormalized_gamete_baseline_no_bias
predDataTotal <- predDataTotal_gamete_baseline_no_bias
name = "gamete_baseline_no_bias"

# gamete + baseline small steps 
predDataTotalNormalized <- predDataTotalNormalized_gamete_small
predDataTotal <- predDataTotal_gamete_small
name = "gamete_small_steps"

# load from RDS
#predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".rds"))

predDataTotalNormalized_tmp <- predDataTotalNormalized[order(as.numeric(predDataTotalNormalized$mutationProb)),]

# plot the normalized data grouped by mutation prob
ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                    group=mutationProb, 
                                    colour = mutationProb, 
                                    shape = mutationProb)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 14)) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$mutationProb)) +
  geom_point() +
  labs(title = name,
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() 
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProb), 
              alpha = 0.2, colour = NA) 

# plot the normalized data grouped by mutation prob
ggplot(predDataTotal, aes(maternalAge, mean, 
                          group=mutationProb, 
                          colour = mutationProb, 
                          shape = mutationProb)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(predDataTotal$mutationProb)) +
  geom_point() +
  labs(title = name, 
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() 
  geom_ribbon(data = predDataTotal,
          aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
          alpha = 0.2, colour = NA) 

  ###############################################################################
  
# PARENTAL CARE QUALITY 
  
  # the gam model
  run_model <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
    return(tmp)
  }
  
  # reading all data and saving both the data as the gam models. 
  parent_path <- paste0(path, "parental_care/defaults/") # with mutation bias
  name = "parental_care_default"
  
  # reading all data and saving both the data as the gam models. 
  parent_path <- paste0(path, "parental_care/no_bias/") # with mutation bias
  name = "parental_care_no_bias"
  
  f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  # paste all data together
  Sys.time()
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path, x)
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    # extract replicate from file name 
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate and scenario as column 
    local_data$rep <- rep
    #local_data$scenario <- scenario
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbAgeSpecificGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbAgeSpecificGenes = factor(mutationProbAgeSpecificGenes))
    
    allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
    
    # perform gam
    local_data <- local_data %>% filter(na > 6)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    to_sample <- 50 
    if (length(unique(local_data$ID)) < to_sample){
      counter = counter + 1
      next
      #to_sample <- length(unique(local_data$ID))
    }
    # sample 100 parents by their IDs
    local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
    # save the data
    allData <- rbind(allData, local_data)
    
    # GAM 
    d <- local_data
    # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
    #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
    #d2 <- d2 %>% mutate(scenario = factor(scenario)) 
    k = 10; 
    if (length(levels(as.factor(d2$maternalAge))) < k) { 
      k = length(levels(as.factor(d2$maternalAge))) 
    }
    z = 4
    if (length(levels(as.factor(d2$ID))) < z) { 
      z = length(levels(as.factor(d2$ID))) 
    }
    print(d2$rep[1])
    # run model 
    mod <- run_model(d2)
    # save model in list 
    models[[i]] <- mod
    # rename model to be unique for replicate and scenario
    names(models)[i] <- paste0("model_", rep, "_", d2$mutationProbAgeSpecificGenes)
  }
  Sys.time()
  print(paste0("number of skipped replicates: ", counter))
  
  # save the R data just in case. 
  saveRDS(allData, file = paste0("allData_", name, ".rds"))
  saveRDS(models, file = paste0("models_", name, ".rds"))
  saveRDS(allDataComplete, file = paste0("allDataComplete_", name, ".rds"))
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))
  
  allData <- loadRDS(paste0("allData_", name, ".rds"))
  models <- loadRDS(paste0("models_", name, ".rds"))
  allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per scenario through the data 
  for (x in levels(allData$mutationProbAgeSpecificGenes)){
    # make subset of scenario 
    sub <- allData[allData$mutationProbAgeSpecificGenes == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProbAgeSpecificGenes = sub$mutationProbAgeSpecificGenes[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      mod <- paste0("model_", i, "_", sub$mutationProbAgeSpecificGenes[1])
      index <- which(names(models) == mod)
      new_col <- paste("rep", i, sep = "")
      tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
      # add to predicted data with adjusting for the intercept. 
      pred_data[[new_col]] <- tmp + attr(tmp, "constant")
    }
    
    # get means, minimum and maximum values 
    repVals <- subset(pred_data, select = 4:ncol(pred_data))
    pred_data$mean <- rowMeans(repVals)
    pred_data$min <- apply(repVals, 1, min)
    pred_data$max <- apply(repVals, 1, max)
    
    # transform back 
    pred_data$mean <- logist(pred_data$mean) 
    pred_data$min <- logist(pred_data$min) 
    pred_data$max <- logist(pred_data$max) 
    
    # save by binding to one big dataframe 
    predDataTotal <- rbind(predDataTotal, pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
    #predDataTotal$nRep <- length(levels(sub$rep)) 
    
    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  saveRDS(predDataTotalNormalized, file = paste0("predDataTotalNormalized_", name, ".rds"))
  
  predDataTotalNormalized_parental_care_no_bias <- predDataTotalNormalized
  predDataTotal_parental_care_no_bias <- predDataTotal
  
  # load from saved extra in environment 
  predDataTotalNormalized <- predDataTotalNormalized_parental_care_defaults
  predDataTotal <- predDataTotal_parental_care_defaults
  name = "parental_care_defaults"
  
  # load from saved extra in environment 
  predDataTotalNormalized <- predDataTotalNormalized_parental_care_no_bias
  predDataTotal <- predDataTotal_parental_care_no_bias
  name = "parental_care_no_bias"
  
  # load from RDS
  #predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProbAgeSpecificGenes, 
                                      colour = mutationProbAgeSpecificGenes, 
                                      shape = mutationProbAgeSpecificGenes)) +
    geom_line() +
    scale_color_manual(values = met.brewer("Egypt", 14)) +
    scale_shape_manual(values=1:nlevels(predDataTotalNormalized$mutationProbAgeSpecificGenes)) +
    geom_point() +
    labs(title = name,
         x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_cowplot() 
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) 
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProbAgeSpecificGenes, 
                            colour = mutationProbAgeSpecificGenes, 
                            shape = mutationProbAgeSpecificGenes)) +
    geom_line() +
    scale_color_manual(values = met.brewer("Egypt", 11)) +
    scale_shape_manual(values=1:nlevels(predDataTotal$mutationProbAgeSpecificGenes)) +
    geom_point() +
    labs(title = name, 
         x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_cowplot() 
  geom_ribbon(data = predDataTotal,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) 
  
  ###############################################################################
  
# RESOURCE DISTRIBUTION 
  
  # the gam model
  run_model <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
    return(tmp)
  }
  
  # reading all data and saving both the data as the gam models. 
  parent_path <- paste0(path, "resource_distribution/defaults/") # with mutation bias
  name = "resource_distribution_defaults"
  
  f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  # paste all data together
  Sys.time()
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path, x)
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    # extract replicate from file name 
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate and scenario as column 
    local_data$rep <- rep
    #local_data$scenario <- scenario
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbInvestmentGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbInvestmentGenes = factor(mutationProbInvestmentGenes))
    
    allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
    
    # perform gam
    local_data <- local_data %>% filter(na > 6)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    to_sample <- 50 
    if (length(unique(local_data$ID)) < to_sample){
      counter = counter + 1
      next
      #to_sample <- length(unique(local_data$ID))
    }
    # sample 100 parents by their IDs
    local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
    # save the data
    allData <- rbind(allData, local_data)
    
    # GAM 
    d <- local_data
    # compares 40 to the expected age at death. If ageAtDeath is lower > that will be y1; else it becomes 40 
    #d2 <- d %>% mutate(y1 = pmin(40,ageAtDeath)) 
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
    #d2 <- d2 %>% mutate(scenario = factor(scenario)) 
    k = 10; 
    if (length(levels(as.factor(d2$maternalAge))) < k) { 
      k = length(levels(as.factor(d2$maternalAge))) 
    }
    z = 4
    if (length(levels(as.factor(d2$ID))) < z) { 
      z = length(levels(as.factor(d2$ID))) 
    }
    print(d2$rep[1])
    # run model 
    mod <- run_model(d2)
    # save model in list 
    models[[i]] <- mod
    # rename model to be unique for replicate and scenario
    names(models)[i] <- paste0("model_", rep, "_", d2$mutationProbInvestmentGenes)
  }
  Sys.time()
  print(paste0("number of skipped replicates: ", counter))
  
  # save the R data just in case. 
  saveRDS(allData, file = paste0("allData_", name, ".rds"))
  saveRDS(models, file = paste0("models_", name, ".rds"))
  saveRDS(allDataComplete, file = paste0("allDataComplete_", name, ".rds"))
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))
  
  allData <- loadRDS(paste0("allData_", name, ".rds"))
  models <- loadRDS(paste0("models_", name, ".rds"))
  allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per scenario through the data 
  for (x in levels(allData$mutationProbInvestmentGenes)){
    # make subset of scenario 
    sub <- allData[allData$mutationProbInvestmentGenes == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProbInvestmentGenes = sub$mutationProbInvestmentGenes[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      mod <- paste0("model_", i, "_", sub$mutationProbInvestmentGenes[1])
      index <- which(names(models) == mod)
      new_col <- paste("rep", i, sep = "")
      tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
      # add to predicted data with adjusting for the intercept. 
      pred_data[[new_col]] <- tmp + attr(tmp, "constant")
    }
    
    # get means, minimum and maximum values 
    repVals <- subset(pred_data, select = 4:ncol(pred_data))
    pred_data$mean <- rowMeans(repVals)
    pred_data$min <- apply(repVals, 1, min)
    pred_data$max <- apply(repVals, 1, max)
    
    # transform back 
    pred_data$mean <- logist(pred_data$mean) 
    pred_data$min <- logist(pred_data$min) 
    pred_data$max <- logist(pred_data$max) 
    
    # save by binding to one big dataframe 
    predDataTotal <- rbind(predDataTotal, pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
    #predDataTotal$nRep <- length(levels(sub$rep)) 
    
    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  saveRDS(predDataTotalNormalized, file = paste0("predDataTotalNormalized_", name, ".rds"))
  
  predDataTotalNormalized_resource_dist <- predDataTotalNormalized
  predDataTotal_resource_dist <- predDataTotal
  
  # load from saved extra in environment 
  predDataTotalNormalized <- predDataTotalNormalized_resource_dist
  predDataTotal <- predDataTotal_resource_dist
  name = "resource_distribution_defaults"
  
  # load from RDS
  #predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProbInvestmentGenes, 
                                      colour = mutationProbInvestmentGenes, 
                                      shape = mutationProbInvestmentGenes)) +
    geom_line() +
    scale_color_manual(values = met.brewer("Egypt", 14)) +
    scale_shape_manual(values=1:nlevels(predDataTotalNormalized$mutationProbInvestmentGenes)) +
    geom_point() +
    labs(title = name,
         x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_cowplot() 
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbInvestmentGenes), 
              alpha = 0.2, colour = NA) 
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProbInvestmentGenes, 
                            colour = mutationProbInvestmentGenes, 
                            shape = mutationProbInvestmentGenes)) +
    geom_line() +
    scale_color_manual(values = met.brewer("Egypt", 11)) +
    scale_shape_manual(values=1:nlevels(predDataTotal$mutationProbInvestmentGenes)) +
    geom_point() +
    labs(title = name, 
         x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_cowplot() 
  geom_ribbon(data = predDataTotal,
              aes(ymin = min, ymax = max,  fill = mutationProbInvestmentGenes), 
              alpha = 0.2, colour = NA) 
  

###############################################################################
# Supplementary materials: Normalizing figures S5 - S7
###############################################################################
path_s <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/combiningAll/"

parent_path <- paste0(path_s, "combiningAll3/resource/")
parent_path <- paste0(path_s, "combiningAll3/nullResource/")
parent_path <- paste0(path_s, "combiningAll3/qualityResource/")

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes.txt", recursive = T, all.files = T)
found <- c()
found_2 <- c()
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
  
  # read longitudinal data to determine the 95th percentile. 
  file_name_2 <- paste0(parent_path, rep, "/outputLifeExpLong.txt")
  local_data_2 <- read.table(file_name_2, header = F)
  colnames(local_data_2) <- c("ID", "ageAtDeath", "maternalAge", "paternalAge")
  found_2 <- rbind(found_2, local_data_2)
}

# determine the mean per replicate. 
age <- c(0:max(found$age))
data_per_rep <- data.frame(age)
for (i in levels(as.factor(found$rep))){
  sub <- found[found$rep == i,]
  new_col <- paste0("rep", i)
  data_per_rep[[new_col]] <- tapply(sub$InvestmentGeneVal, sub$age, mean)
}

# get means
repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
data_per_rep$mean <- rowMeans(repVals)
data_per_rep$min <- apply(repVals, 1, min)
data_per_rep$max <- apply(repVals, 1, max)
cut_off <- quantile(found_2$maternalAge, probs = 0.95) # ageAtDeath > maternalAge for nullResource. maternalAge > ageAtDeath for qualityResource. 
data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] # cut off at 95th percentile

ggplot(data_per_rep, aes(age, (1-mean))) + geom_line() +
  geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
  theme_bw() +
  labs(x = "Age",
       y = "Averaged proportion of resources 
  allocated to reproduction") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  theme(axis.text = element_text(size=17,face="plain",color="black"), 
        axis.title = element_text(size = 17)) 
  scale_x_continuous(breaks = seq(0, 9, 3),
                     limits = c(0, 9))
ggsave("ppResourceNull.pdf", width = 12, height = 8)

  
###############################################################################
# Adjusting upper right corner of matrix by cutting off at 95th percentile. 
###############################################################################

# Parental age distribution plots 
predDataTotalAgeDist <- c()
  
allDataComplete <- loadRDS("allDataComplete.rds")
  
allDataComplete <- allDataComplete %>% mutate(scenario = factor(scenario))
total_age_dist_data <- c()
# loop through the scenarios. 
for (x in levels(allDataComplete$scenario)){
    print(x)
    # make subset of scenario 
    sub <- allDataComplete[allDataComplete$scenario == x,]
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,quantile(sub$maternalAge,probs = 0.95),1),
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
    total_age_dist_data <- rbind(total_age_dist_data, pred_data)
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
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) 
      xlim(0, 40)
    #scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1))
    
    # save plot in list 
    plotsAgeDist[[i]] <- p
    # adjust name to corresponding scenario run
    names(plotsAgeDist)[i] <- paste("p_", scenarios[i], sep = "") 
  }
