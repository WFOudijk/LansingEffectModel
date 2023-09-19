###############################################################################


# MANUSCRIPT R CODE 
# Created at 13-09-2023


###############################################################################

path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/data/"
output_path_data <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/output_data/"
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/Figures/"

###############################################################################
# Resource distribution matrix. 
###############################################################################

###############################################################################
# EQUATION 1
###############################################################################
theme_plots <- theme(axis.text = element_text(size=15,face="plain",color="black"),
                     axis.title = element_text(size = 16),
                     axis.line = element_line(color="black", linewidth = 1.0),
                     panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
                     text = element_text(size = 13),
                     legend.position = c(0.84, 0.86),
                     legend.text = element_text(size = 15),
                     legend.title = element_text(size = 16))


# define equation 1
eq1 <- function(x, c) {
  1 - c * x^2
}
c <- c(0.1, 0.3, 0.7)
dummy_data <- data.frame(x = seq(0, 1, 0.1))
dummy_data$c <- c[1]
dummy_data_2 <- data.frame(x = seq(0, 1, 0.1))
dummy_data_2$c <- c[2]
dummy_data_3 <- data.frame(x = seq(0, 1, 0.1))
dummy_data_3$c <- c[3]

dummy_data <- rbind(dummy_data, dummy_data_2, dummy_data_3)

for (i in 1:nrow(dummy_data)){
  dummy_data[i,"y"] <- eq1(dummy_data[i,]$x, dummy_data[i,]$c)
}

p1 <- ggplot(dummy_data, aes(x, y, group = as.factor(c), colour = as.factor(c))) + 
  geom_line() +
  scale_colour_manual("c", values = met.brewer("Archambault")[c(1, 4, 7)]) +
  theme_minimal() +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.1, 0.8),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1))

p1

# define equation 2
eq2 <- function(x, a, d) {
  1 / (1 + exp(-a * x - d))
}

a <- c(4.5, 3, 10)
d <- c(0, 2, -4)
dummy_data_eq2 <- data.frame(x = seq(0, 1, 0.1))
dummy_data_eq2$a <- a[1]
dummy_data_eq2$d <- d[1]
dummy_data__eq22 <- data.frame(x = seq(0, 1, 0.1))
dummy_data__eq22$a <- a[2]
dummy_data__eq22$d <- d[2]
dummy_data__eq23 <- data.frame(x = seq(0, 1, 0.1))
dummy_data__eq23$a <- a[3]
dummy_data__eq23$d <- d[3]

dummy_data_eq2 <- rbind(dummy_data_eq2, dummy_data__eq22, dummy_data__eq23)

for (i in 1:nrow(dummy_data_eq2)){
  dummy_data_eq2[i,"y"] <- eq2(dummy_data_eq2[i,]$x, dummy_data_eq2[i,]$a, dummy_data_eq2[i,]$d)
}

p2 <- ggplot(dummy_data_eq2, aes(x, y, group = as.factor(a), colour = as.factor(a))) + 
  geom_line() +
  scale_colour_manual("a", values = met.brewer("Archambault")[c(1, 4, 7)]) +
  theme_minimal() +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.1, 0.8),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1))

p2

# set path to data for the resource matrix. 
path_rema <- paste0(path, "resource_matrix/")

# EQUATION 1

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
  return(tmp)
}

# set path 
path_rema_eq1 <- paste0(path_rema, "eq1/")

####### copy from here #################

# resource-only 
parent_path <- paste0(path_rema_eq1, "resource_eq1/") 
name = "resource_eq1"

# resource + parental care quality
parent_path <- paste0(path_rema, "resource_parentalQuality_eq1/") 
name = "resource_quality_eq1"

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
counter = 0
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$weightInvestment[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  # filter parents that had offspring at at least 6 different ages 
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # check if enough IDs are present to sample 
  to_sample <- 100 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  # normalize between 0 and 1 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # for gam > determine number of k
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  # to keep track of where you are in the run 
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and varying parameter
  names(models)[i] <- paste0("model_", rep, "_", d2$weightInvestment)
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  # TODO: remove this row if data is complete. 
  local_data_age$weightInvestment <- local_data$weightInvestment[1]
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}
Sys.time()

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))
  
  #allData <- loadRDS(paste0("allData_", name, ".rds"))
  #models <- loadRDS(paste0("models_", name, ".rds"))
  #allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
  percentiles <- data.frame(levels(allData$weightInvestment))
  colnames(percentiles) <- c("weightInvestment")
  percentiles$percentile <- 0
  
  # go per parameter value through the data 
  for (x in levels(allData$weightInvestment)){
    # make subset of parameter
    sub <- allData[allData$weightInvestment == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            weightInvestment = sub$weightInvestment[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      # get the name for the model 
      mod <- paste0("model_", i, "_", sub$weightInvestment[1])
      # get index of model 
      index <- which(names(models) == mod)
      new_col <- paste("rep", i, sep = "")
      # predict using the corresponding model 
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

    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save the percentile for this subset 
    percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the data 
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }

Sys.time() 

saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p3 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                         group = weightInvestment, 
                                         colour = weightInvestment, 
                                         shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$weightInvestment)) +
  geom_point() +
  labs(title = name ) +
       #x = "Normalized parental age",
       #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) +
  theme_plots 
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))

# make age-specific figure 

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$weightInvestment))){
  sub <- found[found$weightInvestment == i,]
  #tmp <- expand.grid(age = c(0:max(found$age)),
                     # ID = levels(sub$ID)[1],
  #                   weightInvestment = sub$weightInvestment[1])
  percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes between 0 and 1 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))

saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p4 <- ggplot(combine_data, aes(age, mean, 
                         group= weightInvestment, 
                         colour = weightInvestment, 
                         shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$weightInvestment)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 

# RESOURCE + BASELINE 

# resource-only 
parent_path <- paste0(path_rema_eq1, "resource_baseline_eq1/") 
name = "resource_baseline_eq1"

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
counter = 0
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$weightInvestment[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  # filter parents that had offspring at at least 6 different ages 
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # check if enough IDs are present to sample 
  to_sample <- 100 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  # normalize between 0 and 1 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # for gam > determine number of k
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  # to keep track of where you are in the run 
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and varying parameter
  names(models)[i] <- paste0("model_", rep, "_", d2$weightInvestment)
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  # TODO: remove this row if data is complete. 
  local_data_age$weightInvestment <- local_data$weightInvestment[1]
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}
Sys.time()

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))


allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

#allData <- loadRDS(paste0("allData_", name, ".rds"))
#models <- loadRDS(paste0("models_", name, ".rds"))
#allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$weightInvestment))
colnames(percentiles) <- c("weightInvestment")
percentiles$percentile <- 0

# go per parameter value through the data 
for (x in levels(allData$weightInvestment)){
  # make subset of parameter
  sub <- allData[allData$weightInvestment == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          weightInvestment = sub$weightInvestment[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    # get the name for the model 
    mod <- paste0("model_", i, "_", sub$weightInvestment[1])
    # get index of model 
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    # predict using the corresponding model 
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
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # save the percentile for this subset 
  percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  # save the data 
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

Sys.time() 

saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p7 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                          group = weightInvestment, 
                                          colour = weightInvestment, 
                                          shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$weightInvestment)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))

# make age-specific figure 

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$weightInvestment))){
  sub <- found[found$weightInvestment == i,]
  #tmp <- expand.grid(age = c(0:max(found$age)),
  # ID = levels(sub$ID)[1],
  #                   weightInvestment = sub$weightInvestment[1])
  percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p8 <- ggplot(combine_data, aes(age, mean, 
                               group= weightInvestment, 
                               colour = weightInvestment, 
                               shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$weightInvestment)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 

# RESOURCE + GAMETE 

# resource-only 
parent_path <- paste0(path_rema_eq1, "resource_gamete_eq1/") 
name = "resource_gamete_eq1"

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
counter = 0
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$weightInvestment[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  # filter parents that had offspring at at least 6 different ages 
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # check if enough IDs are present to sample 
  to_sample <- 100 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  # normalize between 0 and 1 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # for gam > determine number of k
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  # to keep track of where you are in the run 
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and varying parameter
  names(models)[i] <- paste0("model_", rep, "_", d2$weightInvestment)
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  # TODO: remove this row if data is complete. 
  local_data_age$weightInvestment <- local_data$weightInvestment[1]
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}
Sys.time()

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))


allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

#allData <- loadRDS(paste0("allData_", name, ".rds"))
#models <- loadRDS(paste0("models_", name, ".rds"))
#allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$weightInvestment))
colnames(percentiles) <- c("weightInvestment")
percentiles$percentile <- 0

# go per parameter value through the data 
for (x in levels(allData$weightInvestment)){
  # make subset of parameter
  sub <- allData[allData$weightInvestment == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          weightInvestment = sub$weightInvestment[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    # get the name for the model 
    mod <- paste0("model_", i, "_", sub$weightInvestment[1])
    # get index of model 
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    # predict using the corresponding model 
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
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # save the percentile for this subset 
  percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  # save the data 
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p11 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                          group = weightInvestment, 
                                          colour = weightInvestment, 
                                          shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$weightInvestment)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))

# make age-specific figure 

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$weightInvestment))){
  sub <- found[found$weightInvestment == i,]
  #tmp <- expand.grid(age = c(0:max(found$age)),
  # ID = levels(sub$ID)[1],
  #                   weightInvestment = sub$weightInvestment[1])
  percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p12 <- ggplot(combine_data, aes(age, mean, 
                               group= weightInvestment, 
                               colour = weightInvestment, 
                               shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$weightInvestment)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 
# RESOURCE + GAMETE 

# resource-only 
parent_path <- paste0(path_rema_eq1, "resource_parentalQuality_eq1/") 
name = "resource_quality_eq1"

f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
counter = 0
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$weightInvestment[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  # filter parents that had offspring at at least 6 different ages 
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # check if enough IDs are present to sample 
  to_sample <- 100 
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
  }
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), to_sample, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  # normalize between 0 and 1 
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # for gam > determine number of k
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  # to keep track of where you are in the run 
  print(d2$rep[1])
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and varying parameter
  names(models)[i] <- paste0("model_", rep, "_", d2$weightInvestment)
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  # TODO: remove this row if data is complete. 
  local_data_age$weightInvestment <- local_data$weightInvestment[1]
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}
Sys.time()

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

#allData <- loadRDS(paste0("allData_", name, ".rds"))
#models <- loadRDS(paste0("models_", name, ".rds"))
#allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$weightInvestment))
colnames(percentiles) <- c("weightInvestment")
percentiles$percentile <- 0

# go per parameter value through the data 
for (x in levels(allData$weightInvestment)){
  # make subset of parameter
  sub <- allData[allData$weightInvestment == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          weightInvestment = sub$weightInvestment[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    # get the name for the model 
    mod <- paste0("model_", i, "_", sub$weightInvestment[1])
    # get index of model 
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    # predict using the corresponding model 
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
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # save the percentile for this subset 
  percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  # save the data 
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p15 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                           group = weightInvestment, 
                                           colour = weightInvestment, 
                                           shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$weightInvestment)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))

# make age-specific figure 

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$weightInvestment))){
  sub <- found[found$weightInvestment == i,]
  #tmp <- expand.grid(age = c(0:max(found$age)),
  # ID = levels(sub$ID)[1],
  #                   weightInvestment = sub$weightInvestment[1])
  percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p16 <- ggplot(combine_data, aes(age, mean, 
                                group= weightInvestment, 
                                colour = weightInvestment, 
                                shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$weightInvestment)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 

########## EQUATION 2 ##########

# set path 
path_rema_eq2 <- paste0(path_rema, "eq2/")

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
  return(tmp)
}

####### copy from here #################

# resource-only 
parent_path <- paste0(path_rema_eq2, "resource_eq2/") 
name = "resource_eq2"

# list the files 
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
# to save the data
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
# for every listed file 
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  nVal <- length(splitted_path[[1]]) - 2
  val <- splitted_path[[1]][nVal]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(scalingStrengthOfAllocationToReproduce = factor(scalingStrengthOfAllocationToReproduce))
  local_data <- local_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # set k for the gam model 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  
  # to keep track 
  print(d2$rep[1])
  
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and parameter setting 
  names(models)[i] <- paste0("model_", rep, "_", d2$steepnessAllocationToReproduce[1])
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, val, "/", rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
colnames(percentiles) <- c("steepnessAllocationToReproduce")
percentiles$percentile <- 0

# go per scenario through the data 
for (x in levels(allData$steepnessAllocationToReproduce)){
  # make subset of scenario 
  sub <- allData[allData$steepnessAllocationToReproduce == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1],
                          scalingStrengthOfAllocationToReproduce = sub$scalingStrengthOfAllocationToReproduce[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$steepnessAllocationToReproduce[1])
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 5:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p5 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                           group = steepnessAllocationToReproduce, 
                                           colour = steepnessAllocationToReproduce, 
                                           shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))


# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
  sub <- found[found$steepnessAllocationToReproduce == i,]
  percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p6 <- ggplot(combine_data, aes(age, mean, 
                                group= steepnessAllocationToReproduce, 
                                colour = steepnessAllocationToReproduce, 
                                shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 

# RESOURCE + BASELINE

# resource-only 
parent_path <- paste0(path_rema_eq2, "resource_baseline_eq2/") 
name = "resource_baseline_eq2"

# list the files 
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
# to save the data
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
# for every listed file 
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  nVal <- length(splitted_path[[1]]) - 2
  val <- splitted_path[[1]][nVal]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(scalingStrengthOfAllocationToReproduce = factor(scalingStrengthOfAllocationToReproduce))
  local_data <- local_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # set k for the gam model 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  
  # to keep track 
  print(d2$rep[1])
  
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and parameter setting 
  names(models)[i] <- paste0("model_", rep, "_", d2$steepnessAllocationToReproduce[1])
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, val, "/", rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
colnames(percentiles) <- c("steepnessAllocationToReproduce")
percentiles$percentile <- 0

# go per scenario through the data 
for (x in levels(allData$steepnessAllocationToReproduce)){
  # make subset of scenario 
  sub <- allData[allData$steepnessAllocationToReproduce == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1],
                          scalingStrengthOfAllocationToReproduce = sub$scalingStrengthOfAllocationToReproduce[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$steepnessAllocationToReproduce[1])
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 5:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p9 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                          group = steepnessAllocationToReproduce, 
                                          colour = steepnessAllocationToReproduce, 
                                          shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))


# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
  sub <- found[found$steepnessAllocationToReproduce == i,]
  percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p10 <- ggplot(combine_data, aes(age, mean, 
                               group= steepnessAllocationToReproduce, 
                               colour = steepnessAllocationToReproduce, 
                               shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 

# RESOURCE + GAMETE

# resource-only 
parent_path <- paste0(path_rema_eq2, "resource_gamete_eq2/") 
name = "resource_gamete_eq2"

# list the files 
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
# to save the data
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
# for every listed file 
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  nVal <- length(splitted_path[[1]]) - 2
  val <- splitted_path[[1]][nVal]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(scalingStrengthOfAllocationToReproduce = factor(scalingStrengthOfAllocationToReproduce))
  local_data <- local_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # set k for the gam model 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  
  # to keep track 
  print(d2$rep[1])
  
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and parameter setting 
  names(models)[i] <- paste0("model_", rep, "_", d2$steepnessAllocationToReproduce[1])
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, val, "/", rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
colnames(percentiles) <- c("steepnessAllocationToReproduce")
percentiles$percentile <- 0

# go per scenario through the data 
for (x in levels(allData$steepnessAllocationToReproduce)){
  # make subset of scenario 
  sub <- allData[allData$steepnessAllocationToReproduce == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1],
                          scalingStrengthOfAllocationToReproduce = sub$scalingStrengthOfAllocationToReproduce[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$steepnessAllocationToReproduce[1])
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 5:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p13 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                          group = steepnessAllocationToReproduce, 
                                          colour = steepnessAllocationToReproduce, 
                                          shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))


# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
  sub <- found[found$steepnessAllocationToReproduce == i,]
  percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p14 <- ggplot(combine_data, aes(age, mean, 
                                group= steepnessAllocationToReproduce, 
                                colour = steepnessAllocationToReproduce, 
                                shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 

# RESOURCE + PARENTAL CARE QUALITY

# resource-only 
parent_path <- paste0(path_rema_eq2, "resource_and_parentalQuality_eq2/") 
name = "resource_and_parentalQuality_eq2"

# list the files 
f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
# to save the data
allData <- c()
allDataComplete <- c()
models <- list()
found <- c()
counter <- 0
# for every listed file 
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  nVal <- length(splitted_path[[1]]) - 2
  val <- splitted_path[[1]][nVal]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  # rename ID so it is unique per replicate per parameter setting
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(scalingStrengthOfAllocationToReproduce = factor(scalingStrengthOfAllocationToReproduce))
  local_data <- local_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
  
  # perform gam
  na_filter = 6
  if (max(local_data$na) <= na_filter){
    na_filter = 3
  }
  local_data <- local_data %>% filter(na > 6)
  local_data <- local_data %>% mutate(ID = factor(ID)) 
  
  to_sample <- 100
  if (length(unique(local_data$ID)) < to_sample){
    counter = counter + 1
    next
  }
  
  # sample 100 parents by their IDs
  local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
  # save the data
  allData <- rbind(allData, local_data)
  
  # GAM 
  d <- local_data
  d2 <- d %>% mutate(y2 = ageAtDeath/40)
  ## Work with logits (then (0,1) -> (-inf,+inf))
  d2 <- d2 %>% mutate(y3 = car::logit(y2))
  d2 <- d2 %>% mutate(ID = factor(ID)) 
  d2 <- d2 %>% mutate(rep = factor(rep)) 
  
  # set k for the gam model 
  k = 10; 
  if (length(levels(as.factor(d2$maternalAge))) < k) { 
    k = length(levels(as.factor(d2$maternalAge))) 
  }
  z = 4
  if (length(levels(as.factor(d2$ID))) < z) { 
    z = length(levels(as.factor(d2$ID))) 
  }
  if (k <= z){
    z = (k - 1)
  }
  # to keep track 
  print(d2$rep[1])
  
  # run model 
  mod <- run_model(d2)
  # save model in list 
  models[[i]] <- mod
  # rename model to be unique for replicate and parameter setting 
  names(models)[i] <- paste0("model_", rep, "_", d2$steepnessAllocationToReproduce[1])
  
  # COLLECT DATA FOR AGE-SPECIFIC FIGURE
  # get path 
  file_name <- paste0(parent_path, val, "/", rep, "/outputWithAgeSpecificGenes.txt")
  # read data
  local_data_age <- read.table(file_name, header = T)
  # add replicate as column 
  local_data_age$rep <- rep
  local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
  local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
  found <- rbind(found, local_data_age)
}

# save the R data just in case. 
saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
colnames(percentiles) <- c("steepnessAllocationToReproduce")
percentiles$percentile <- 0

# go per scenario through the data 
for (x in levels(allData$steepnessAllocationToReproduce)){
  # make subset of scenario 
  sub <- allData[allData$steepnessAllocationToReproduce == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1],
                          scalingStrengthOfAllocationToReproduce = sub$scalingStrengthOfAllocationToReproduce[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$steepnessAllocationToReproduce[1])
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 5:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))

# plot the normalized data grouped by weight investment
p17 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                           group = steepnessAllocationToReproduce, 
                                           colour = steepnessAllocationToReproduce, 
                                           shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name ) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_minimal() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) +
  theme_plots 
scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0,1.6, 0.2), limits = c(0, 1.6))


# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
  sub <- found[found$steepnessAllocationToReproduce == i,]
  percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
  tmp <- expand.grid(age = c(0:percentile),
                     steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
  sub <- sub[sub$age <= percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  # recalculate to make it investment in reproduction instead of repair. 
  tmp$mean <- (1 - tmp$mean)
  tmp$min <- (1 - tmp$min)
  tmp$max <- (1 - tmp$max)
  
  # normalize the axes 
  tmp$age <- tmp$age / percentile
  tmp$min <- tmp$min / tmp$mean[1]
  tmp$max <- tmp$max / tmp$mean[1]
  tmp$mean <- tmp$mean / tmp$mean[1]
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))

# plot 
p18 <- ggplot(combine_data, aes(age, mean, 
                                group= steepnessAllocationToReproduce, 
                                colour = steepnessAllocationToReproduce, 
                                shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 

############ make grid ##################

resource_matrix <- plot_grid(p1, p2,
                             p3, p4, p5, p6,
                             p7, p8, p9, p10,
                             p11, p12, p13, p14,
                             p15, p16, p17, p18,
                         ncol = 4, align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.88, label_y = 0.97)
resource_matrix

############ until here #################

allData$rep <- factor(allData$rep)
logist <- function(x) 40/(1+exp(-x))

allData <- loadRDS(paste0("allData_", name, ".rds"))
models <- loadRDS(paste0("models_", name, ".rds"))
allDataComplete <- loadRDS(paste0("allDataComplete_", name, ".rds"))

# to save all data 
predDataTotal <- c()
predDataTotalNormalized <- c()

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$weightInvestment))
colnames(percentiles) <- c("weightInvestment")
percentiles$percentile <- 0

# go per parameter value through the data 
for (x in levels(allData$weightInvestment)){
  # make subset of parameter
  sub <- allData[allData$weightInvestment == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          weightInvestment = sub$weightInvestment[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    # get the name for the model 
    mod <- paste0("model_", i, "_", sub$weightInvestment[1])
    # get index of model 
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    # predict using the corresponding model 
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
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  # save the percentile for this subset 
  percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  # save the data 
  pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

# save in environment
predDataTotalNormalized_resource_quality_eq1 <- predDataTotalNormalized
predDataTotal_resource_quality_eq1 <- predDataTotal
percentiles_resource_quality_eq1 <- percentiles

# load from saved extra in environment 
# RESOURCE
predDataTotalNormalized <- predDataTotalNormalized_resource_eq1
predDataTotal <- predDataTotal_resource_eq1
name = "resource_eq1"

# load from saved extra in environment 
# RESOURCE - old 
predDataTotalNormalized <- predDataTotalNormalized_resource_vary_eq1
predDataTotal <- predDataTotal_resource_vary_eq1
name = "resource_eq1"

# load from saved extra in environment 
# RESOURCE + QUALITY
predDataTotalNormalized <- predDataTotalNormalized_resource_quality_eq1
predDataTotal <- predDataTotal_resource_quality_eq1
percentiles <- percentiles_resource_quality_eq1
name = "resource_quality_eq1"

# load from saved extra in environment 
# RESOURCE + QUALITY - old 
predDataTotalNormalized <- predDataTotalNormalized_resource_quality_varying_weight_bias_0.01
predDataTotal <- predDataTotal_resource_quality_varying_weight_bias_0.01
percentiles <- percentiles_resource_quality_vary_eq2
name = "resource_quality_vary_eq2"

# plot the normalized data grouped by weight investment
ggplot(predDataTotalNormalized_resource_eq1, aes(maternalAge, mean, 
                                    group = weightInvestment, 
                                    colour = weightInvestment, 
                                    shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 14)) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$weightInvestment)) +
  geom_point() +
  labs(title = name,
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 

# plot the not-normalized data grouped by weight investment
ggplot(predDataTotal, aes(maternalAge, mean, 
                          group=weightInvestment, 
                          colour = weightInvestment, 
                          shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(predDataTotal$weightInvestment)) +
  geom_point() +
  labs(title = name, 
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() 
geom_ribbon(data = predDataTotal,
            aes(ymin = min, ymax = max,  fill = weightInvestment), 
            alpha = 0.2, colour = NA) 


# making figure with age-specific genes 
# AGE-SPECIFIC FIGURE 
parent_path <- paste0(path_rema, "resource_eq1/") 
name = "resource_eq1"

parent_path <- paste0(path, "resource_parentalQuality/vary_eq1/") 
name = "resource_parentalQuality_vary_eq1"

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes.txt", recursive = T, all.files = T)
found <- c()
# paste all data together
Sys.time()
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  splitted_path <- strsplit(file_name, "/")
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
  local_data$ID <- sub("^", local_data$weightInvestment[1], local_data$ID)
  
  found <- rbind(found, local_data)
}
# TODO: place within one for loop. 

saveRDS(found, paste0(name, ".RDS"))

found <- loadRDS(paste0(name, ".RDS"))

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$weightInvestment))){
  sub <- found[found$weightInvestment == i,]
  tmp <- expand.grid(age = c(0:max(found$age)),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  tmp <- expand.grid(age = c(0:percentiles[percentiles$weightInvestment == i,]$percentile),
                     # ID = levels(sub$ID)[1],
                     weightInvestment = sub$weightInvestment[1])
  sub <- sub[sub$age <= percentiles[percentiles$weightInvestment == i,]$percentile,]
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  
  # get means
  repVals <- subset(tmp, select = 3:ncol(tmp))
  tmp$mean <- rowMeans(repVals)
  tmp$min <- apply(repVals, 1, min)
  tmp$max <- apply(repVals, 1, max)
  
  
  combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
}

combine_data_resource_eq1 <- combine_data

# RESOURCE
combine_data <- combine_data_resource_eq1
name <- "resource_eq1"

# RESOURCE - old
combine_data <- combine_data_resource_vary_eq1
name <- "resource_eq1"

# RESOURCE + QUALITY
combine_data <- combine_data_resource_parentalQuality_vary_eq1
name <- "age-specfic_resource_parentalQuality_vary_eq1"

# RESOURCE + QUALITY - old
combine_data <- combine_data_resource_parentalQuality_vary_eq1
name <- "age-specfic_resource_parentalQuality_vary_eq1"

combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))

# plot 
ggplot(combine_data, aes(age, (1-mean), 
                         group= weightInvestment, 
                         colour = weightInvestment, 
                         shape = weightInvestment)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$weightInvestment)) +
  geom_point() +
  labs(title = name) +
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = (1-min), ymax = (1-max),  fill = weightInvestment), 
              alpha = 0.2, colour = NA) 

# RESOURCE DISTRIBUTION 3 > varying equation 2 

# the gam model
run_model <- function(data) {
  tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
  return(tmp)
}

# reading all data and saving both the data as the gam models. 
# RESOURCE
parent_path <- paste0(path, "resource_distribution/vary_eq2/") 
name = "resource_vary_eq2"

# reading all data and saving both the data as the gam models. 
# RESOURCE + QUALITY
parent_path <- paste0(path, "resource_parentalQuality/vary_eq2/") 
name = "resource_parentalQuality_vary_eq2"

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
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], paste0("_", local_data$ID))
  local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
  # subset parents 
  z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
  ## Add the counter to data
  local_data <- local_data %>% left_join(z,by="ID")
  local_data <- local_data %>% mutate(scalingStrengthOfAllocationToReproduce = factor(scalingStrengthOfAllocationToReproduce))
  local_data <- local_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
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
  names(models)[i] <- paste0("model_", rep, "_", d2$steepnessAllocationToReproduce)
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

# to save the percentiles where the data is cut-off per scenario
percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
colnames(percentiles) <- c("steepnessAllocationToReproduce")
percentiles$percentile <- 0

# go per scenario through the data 
for (x in levels(allData$steepnessAllocationToReproduce)){
  # make subset of scenario 
  sub <- allData[allData$steepnessAllocationToReproduce == x,]
  # refactor the replicates of the subset 
  sub <- sub %>% mutate(rep = factor(rep))
  # make data expansion 
  pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                          ID = levels(sub$ID)[1],
                          steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1],
                          scalingStrengthOfAllocationToReproduce = sub$scalingStrengthOfAllocationToReproduce[1]) 
  
  # loop through the replicates
  for (i in levels(sub$rep)){
    mod <- paste0("model_", i, "_", sub$steepnessAllocationToReproduce[1])
    index <- which(names(models) == mod)
    new_col <- paste("rep", i, sep = "")
    tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
    # add to predicted data with adjusting for the intercept. 
    pred_data[[new_col]] <- tmp + attr(tmp, "constant")
  }
  
  # get means, minimum and maximum values 
  repVals <- subset(pred_data, select = 5:ncol(pred_data))
  pred_data$mean <- rowMeans(repVals)
  pred_data$min <- apply(repVals, 1, min)
  pred_data$max <- apply(repVals, 1, max)
  
  # transform back 
  pred_data$mean <- logist(pred_data$mean) 
  pred_data$min <- logist(pred_data$min) 
  pred_data$max <- logist(pred_data$max) 
  
  # save by binding to one big dataframe 
  predDataTotal <- rbind(predDataTotal, pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
  #predDataTotal$nRep <- length(levels(sub$rep)) 
  
  # get 95th percentile 
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
  # cut-off the data at that percentile 
  pred_data <- pred_data[pred_data$maternalAge <= percentile,]
  # normalize the axes between 0 and 1 
  pred_data$maternalAge <- pred_data$maternalAge / percentile
  pred_data$min <- pred_data$min / pred_data$mean[1]
  pred_data$max <- pred_data$max / pred_data$mean[1]
  pred_data$mean <- pred_data$mean / pred_data$mean[1]
  
  pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
  pred_data_tmp$nRep <- length(levels(sub$rep))
  
  # save the data 
  predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
}

saveRDS(predDataTotalNormalized, file = paste0("predDataTotalNormalized_", name, ".rds"))

predDataTotalNormalized_resource_parentalQuality_vary_eq2 <- predDataTotalNormalized
predDataTotal_resource_parentalQuality_vary_eq2 <- predDataTotal
percentiles_resource_parentalQuality_vary_eq2 <- percentiles 

# load from saved extra in environment 
# RESOURCE
predDataTotalNormalized <- predDataTotalNormalized_resource_vary_eq2
predDataTotal <- predDataTotal_resource_vary_eq2
percentiles <- percentiles_resource_vary_eq2
name = "resource_vary_eq2"

# load from saved extra in environment 
# RESOURCE - old 
predDataTotalNormalized <- predDataTotalNormalized_resource_vary_eq2
predDataTotal <- predDataTotal_resource_vary_eq2
percentiles <- percentiles_resource_vary_eq2
name = "resource_vary_eq2"

# load from saved extra in environment 
# RESOURCE + QUALITY
predDataTotalNormalized <- predDataTotalNormalized_resource_parentalQuality_vary_eq2
predDataTotal <- predDataTotal_resource_parentalQuality_vary_eq2
percentiles <- percentiles_resource_parentalQuality_vary_eq2
name = "resource_parentalQuality_vary_eq2"

# load from saved extra in environment 
# RESOURCE + QUALITY - old 
predDataTotalNormalized <- predDataTotalNormalized_resource_parentalQuality_vary_eq2
predDataTotal <- predDataTotal_resource_parentalQuality_vary_eq2
percentiles <- percentiles_resource_parentalQuality_vary_eq2
name = "resource_parentalQuality_vary_eq2"

# load from RDS
#predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".rds"))

# plot the normalized data grouped by mutation prob
ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                    group=steepnessAllocationToReproduce, 
                                    colour = steepnessAllocationToReproduce, 
                                    shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 14)) +
  scale_shape_manual(values=1:nlevels(predDataTotalNormalized$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name,
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 

# plot the normalized data grouped by mutation prob
ggplot(predDataTotal, aes(maternalAge, mean, 
                          group=steepnessAllocationToReproduce, 
                          colour = steepnessAllocationToReproduce, 
                          shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(predDataTotal$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name, 
       x = "Normalized parental age",
       y = "Normalized offspring age at death") +
  theme_cowplot() 
geom_ribbon(data = predDataTotal,
            aes(ymin = min, ymax = max,  fill = weightInvestment), 
            alpha = 0.2, colour = NA) 

# AGE-SPECIFIC FIGURE 
# RESOURCE
parent_path <- paste0(path, "resource_distribution/vary_eq2/") 
name = "resource_vary_eq2"

# RESOURCE + QUALITY
parent_path <- paste0(path, "resource_parentalQuality/vary_eq2/") 
name = "resource_parentalQuality_vary_eq2"

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes.txt", recursive = T, all.files = T)
found <- c()
counter = 0
# paste all data together
Sys.time()
for (x in f) {
  # get path 
  file_name <- paste0(parent_path, x)
  splitted_path <- strsplit(file_name, "/")
  # extract replicate from file name 
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  # add replicate as column 
  local_data$rep <- rep
  local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
  local_data$ID <- sub("^", local_data$steepnessAllocationToReproduce[1], local_data$ID)
  
  found <- rbind(found, local_data)
}

# determine the mean per replicate. 
combine_data <- c()
for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
  sub <- found[found$steepnessAllocationToReproduce == i,]
  tmp <- expand.grid(age = c(0:max(found$age)),
                     # ID = levels(sub$ID)[1],
                     steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
  
  for (x in levels(as.factor(sub$rep))){
    sub_2 <- sub[sub$rep == x,]
    new_col <- paste0("rep", x)
    tmp[[new_col]] <- tapply(sub_2$investmentGeneVal, sub_2$age, mean)
  }
  combine_data <- rbind(combine_data, tmp)
}

# get means
repVals <- subset(combine_data, select = 3:ncol(combine_data))
combine_data$mean <- rowMeans(repVals)
combine_data$min <- apply(repVals, 1, min)
combine_data$max <- apply(repVals, 1, max)

combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))

# save data
combine_data_resource_parentalQuality_vary_eq2 <- combine_data

# load data
# RESOURCE
combine_data <- combine_data_resource_vary_eq2
name = "resource_vary_eq2"

# RESOURCE + QUALITY
combine_data <- combine_data_resource_parentalQuality_vary_eq2
name = "resource_parentalQuality_vary_eq2"

# plot 
ggplot(combine_data, aes(age, (1-mean), 
                         group= steepnessAllocationToReproduce, 
                         colour = steepnessAllocationToReproduce, 
                         shape = steepnessAllocationToReproduce)) +
  geom_line() +
  scale_color_manual(values = met.brewer("Egypt", 11)) +
  scale_shape_manual(values=1:nlevels(combine_data$steepnessAllocationToReproduce)) +
  geom_point() +
  labs(title = name) + 
  #x = "Normalized parental age",
  #y = "Normalized offspring age at death") +
  theme_cowplot() +
  geom_ribbon(data = combine_data,
              aes(ymin = (1-min), ymax = (1-max),  fill = steepnessAllocationToReproduce), 
              alpha = 0.2, colour = NA) 


#### remove: #######

# resource + parental care quality
parent_path <- paste0(path_rema_eq1, "resource_gamete_eq1/") 
name = "resource_quality_eq1"

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes.txt", recursive = T, all.files = T)
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  # read data
  local_data <- read.table(file_name, header = T)
  parameters <- read.table(paste0(parent_path, rep, "/parameters.txt"))
  nVal <- which(strsplit(as.character(parameters), " ") == "weightInvestment")[1] + 1
  local_data$weightInvestment <- as.numeric(strsplit(as.character(parameters), " ")[nVal])
  write.table(local_data, paste0(parent_path, x))
}

# resource + parental care quality
parent_path <- paste0(path_rema_eq2, "resource_and_parentalQuality_eq2/") 

f <- list.files(path = parent_path, pattern = "outputWithAgeSpecificGenes", recursive = T, all.files = T)
# paste all data together
Sys.time()
for (i in 1:length(f)) {
  x <- f[i]
  # get path 
  file_name <- paste0(parent_path, x)
  # extract replicate from file name 
  splitted_path <- strsplit(file_name, "/")
  nRep <- length(splitted_path[[1]]) - 1
  rep <- splitted_path[[1]][nRep]
  nType <- length(splitted_path[[1]]) - 2
  type <- splitted_path[[1]][nType]
  
  # read data
  local_data <- read.table(file_name, header = T)
  parameters <- read.table(paste0(parent_path, type, "/", rep, "/parameters.txt"))
  nVal <- which(strsplit(as.character(parameters), " ") == "scalingStrengthOfAllocationToReproduce")[1] + 1
  local_data$scalingStrengthOfAllocationToReproduce <- as.numeric(strsplit(as.character(parameters), " ")[nVal])
  nVal2 <- which(strsplit(as.character(parameters), " ") == "steepnessAllocationToReproduce")[1] + 1
  local_data$steepnessAllocationToReproduce <- as.numeric(strsplit(as.character(parameters), " ")[nVal2])
  
  write.table(local_data, paste0(parent_path, x))
}
