###############################################################################


# MANUSCRIPT R CODE 
# Created at 13-09-2023


###############################################################################
# SET PATHS 

# path to the data 
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/data/"

# set save_output to TRUE or FALSE. If TRUE; set output_path_data to where the data should be saved
save_output <- TRUE # set to false if you do not want to save all the data
output_path_data <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/output_data/new/"

# output path to the figures. 
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/Figures/"

# libraries 
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(mgcv)
library(MetBrewer)
library(R.filesets) # remove 

###############################################################################
# Figure 3 - the matrix
###############################################################################
  # user-defined function to transform the ages at death back 
  logist <- function(x) 40/(1+exp(-x))

# 1. MAKING GAM MODELS FOR EVERY REPLICATE PER SCENARIO 

  # the gam model
  run_model <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = 4), data = data, method = "REML")
    return(tmp)
  }

  parent_path <- paste(path, "matrix_data/", sep = "")
  f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  # paste all data together
  for (i in 1:length(f)) {
    x <- f[i]
    # get path to file 
    file_name <- paste0(parent_path, x)
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    nScenario <- length(splitted_path[[1]]) - 2
    scenario <- splitted_path[[1]][nScenario]
    # extract replicate from file name 
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    
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
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # filter data using na; number of ages a parent gets offspring
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # sample 100 parents by their IDs
    local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
    # save the data
    allData <- rbind(allData, local_data)
    
    # Perform GAM on the data
    d <- local_data
    # normalize age at death 
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(scenario = factor(scenario)) 
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
    print(d2$scenario[1]) # to keep track 
    # set k to 10 
    k = 10;
    # if there are less unique maternal ages than k, set to max 
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_matrix.rds"))
    saveRDS(models, file = paste0(output_path_data, "models_matrix.rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_matrix.rds"))
  }
  
  #allData <- loadRDS(paste0(output_path_data, "allData_matrix.rds"))
  #models <- loadRDS(paste0(output_path_data, "models_matrix.rds"))
  
# 2. PREDICTING AND NORMALIZING THE DATA 

  # make rep and scenario a factor 
  allData$scenario <- factor(allData$scenario)
  allData$rep <- factor(allData$rep)
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
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
      # find corresponding model 
      mod <- paste("model_", i, "_", x, sep = "")
      index <- which(names(models) == mod)
      new_col <- paste("rep", i, sep = "")
      # predict using the corresponding model 
      tmp <- predict(models[[index]], newdata = pred_data, type="terms",terms = c("s(maternalAge)"))
      # add to predicted data with adjusting for the intercept. 
      pred_data[[new_col]] <- tmp + attr(tmp, "constant")
    }
    
    # get means, minimum and maximum values per age 
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
    
    # get 95th percentile and add to the saved percentiles df
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    percentiles[percentiles$scenario == x,]$percentile <- percentile
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data)
  }
  
  if (save_output){
    # save the percentiles 
    saveRDS(percentiles, paste0(output_path_data, "percentiles_matrix.RDS"))
    # save the normalized dataset 
    saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_matrix.RDS"))
  }
  
# 3. PERFORMING CROSS-SECTIONAL GAM ANALYSIS PER REPLICATE PER SCENARIO 
  
  # the gam model
  run_model_lat <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k), data=data,method = "REML")
    return(tmp)
  }
  
  parent_path <- paste(path, "matrix_data/", sep = "")
  f <- list.files(path = parent_path, pattern = "outputLifeExp.txt", recursive = T)
  allDataLat <- c()
  modelsLat <- list()
  for (i in 1:length(f)) {
    x <- f[i]
    # get path to file 
    file_name <- paste0(parent_path, x)
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    nScenario <- length(splitted_path[[1]]) - 2
    scenario <- splitted_path[[1]][nScenario]
    # extract replicate from file name 
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
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
    # normalize the data 
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(scenario = factor(scenario)) 
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
    print(d2$scenario[1])
    
    # set k to 10 
    k = 10; 
    # if there are less unique maternal ages than k, set to max
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allDataLat, file = paste0(output_path_data, "allDataLat_matrix.rds"))
    saveRDS(modelsLat, file = paste0(output_path_data, "modelsLat_matrix.rds"))
  }
  
  # make factor 
  allDataLat$scenario <- factor(allDataLat$scenario)
  allDataLat$rep <- factor(allDataLat$rep)
  
# 4. PREDICTING AND NORMALIZING DATA   
  
  # to save all data 
  predDataTotalLat <- c()
  predDataTotalNormalizedLat <- c()
  
  # loop through the data per scenario 
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
      # find corresponding model 
      mod <- paste("model_", i, "_", x, sep = "")
      index <- which(names(modelsLat) == mod)
      new_col <- paste("rep", i, sep = "")
      # predict using the corresponding model 
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
    pred_data$mean <- logist(pred_data$mean)
    pred_data$min <- logist(pred_data$min) 
    pred_data$max <- logist(pred_data$max) 
    
    # save the data 
    predDataTotalLat <- rbind(predDataTotalLat, pred_data)
    
    # determine percentile and cut the data off at this point 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    
    # normalize the axes 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the normalized data 
    predDataTotalNormalizedLat <- rbind(predDataTotalNormalizedLat, pred_data)
  }
  
# 5. DYNAMICALLY GENERATING PLOTS 
  
  # add column with group to both normalized data sets 
  predDataTotalNormalized$group <- "Longitudinal"
  predDataTotalNormalizedLat$group <- "Cross-sectional"
  # merge them together 
  totalNormalizedData <- rbind(predDataTotalNormalized, predDataTotalNormalizedLat)
  
  totalNormalizedData <- totalNormalizedData %>% mutate(scenario = factor(scenario))
  totalNormalizedData <- totalNormalizedData %>% mutate(group = factor(group))
  
  # dynamically generate the plots per scenario 
  plotsTot <- c()
  for (i in 1:length(levels(totalNormalizedData$scenario))){
    p <- ggplot(totalNormalizedData[totalNormalizedData$scenario == scenarios[i],], 
                aes(maternalAge, mean, group = group, colour = group)) +
      geom_line() +
      labs(x = NULL,
           y = NULL) +
      geom_ribbon(data = totalNormalizedData[totalNormalizedData$scenario == scenarios[i],],
                  aes(ymin = min, ymax = max,  fill = group), 
                  alpha = 0.2, colour = NA) +
      theme_minimal() +
      theme(legend.position = "none",
        axis.text = element_text(size=13,face="plain",color="black"),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) +
      scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
      scale_y_continuous(breaks = seq(0, 1.5, 0.5), limits = c(0, 1.5), 
                         oob = scales::oob_keep) +
      scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
      scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)])
    
    # save plot in list 
    plotsTot[[i]] <- p
    # adjust name to corresponding scenario run
    names(plotsTot)[i] <- paste("p_", scenarios[i], sep = "") 
  }
  
# 7. MAKE MATRIX 
  
  # reset the axis numbering that should not be portrayed 
  plotsTot$p_baseline <- plotsTot$p_baseline + theme(axis.text.x = element_blank())
  plotsTot$p_gamete_and_baseline <- plotsTot$p_gamete_and_baseline + theme(axis.text.x = element_blank())
  plotsTot$p_baseline_quality <- plotsTot$p_baseline_quality + theme(axis.text.x = element_blank())
  plotsTot$p_gamete_damage <- plotsTot$p_gamete_damage + theme(axis.text.x = element_blank(),
                                                 axis.text.y = element_blank())
  plotsTot$p_gamete_quality <- plotsTot$p_gamete_quality + theme(axis.text.x = element_blank(),
                                                               axis.text.y = element_blank())
  plotsTot$p_gamete_resource <- plotsTot$p_gamete_resource + theme(axis.text.y = element_blank())
  plotsTot$p_parental_care_quality <- plotsTot$p_parental_care_quality + theme(axis.text.x = element_blank(),
                                                   axis.text.y = element_blank())
  plotsTot$p_resource_and_parentalQuality <- plotsTot$p_resource_and_parentalQuality + theme(axis.text.y = element_blank())
  plotsTot$p_resource_distribution <- plotsTot$p_resource_distribution + theme(axis.text.y = element_blank())
  
  # name the columns and rows 
  plotsTot$p_baseline <- plotsTot$p_baseline + 
    labs(tag = "Mutation accumulation",
          title = "Mutation accumulation") + 
    theme(plot.tag = element_text(angle = -90, size = 16), 
          plot.tag.position = c(1.08, 0.5), 
          legend.box.margin = margin(l = 20), 
          plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 16))
  
  plotsTot$p_gamete_damage <- plotsTot$p_gamete_damage + 
    labs(tag = "Gamete quality",
        title = "Gamete quality") + 
    theme(plot.tag = element_text(angle = -90, size = 16), 
          plot.tag.position = c(1.08, 0.5), 
          legend.box.margin = margin(l = 20), 
          plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 16))
  
  plotsTot$p_parental_care_quality <- plotsTot$p_parental_care_quality + 
    labs(tag = "Parental care",
         title = "Parental care") + 
    theme(plot.tag = element_text(angle = -90, size = 16), 
          plot.tag.position = c(1.08, 0.5), 
          legend.box.margin = margin(l = 20), 
          plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 16))
  
  plotsTot$p_resource_distribution <- plotsTot$p_resource_distribution + 
    labs(tag = "Resource allocation",
        title = "Resource allocation") + 
    theme(plot.tag = element_text(angle = -90, size = 16), 
          plot.tag.position = c(1.08, 0.5), 
          legend.box.margin = margin(l = 20), 
          plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 16))
  
  # set the label numbering for the matrix 
  labels <- c("A", NA, NA, NA,
              "B", "C", NA, NA,
              "D", "E", "F", NA,
              "G", "H", "I", "J")
  
  # make the matrix 
  plot_matrix <- plot_grid(plotsTot$p_baseline, NULL, NULL, NULL,
                           plotsTot$p_gamete_and_baseline, plotsTot$p_gamete_damage, NULL, NULL,
                           plotsTot$p_baseline_quality, plotsTot$p_gamete_quality, plotsTot$p_parental_care_quality, NULL, 
                           plotsTot$p_baseline_resource, plotsTot$p_gamete_resource, plotsTot$p_resource_and_parentalQuality, plotsTot$p_resource_distribution,
                           ncol = 4, align = "hv", axis = "brlt", labels = labels, label_x = 0.82, label_y = 0.9)
  
  #plot_matrix
  
  # make legend 
  p_tmp <- plotsTot$p_resource_distribution + theme(legend.position = "bottom",
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 17))
  leg1 <- get_legend(p_tmp) 
  legend <- as_ggplot(leg1)
  legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
  blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
  legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)
  
  # make x and y labeling 
  plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                           left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))
  
  # make figure with legend 
  plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
  plot_with_legend
  # save figure 
  ggsave(paste0(output_path, "Lansing_fig.pdf"), width = 15, height = 15)
  
###############################################################################
# Supplementary materials: looking at parameter simulations longitudinal. 
# with ten replicates for the confidence intervals. 
# for figures S1 - S4
###############################################################################
  theme_plots_s <- theme(axis.text = element_text(size=15,face="plain",color="black"),
                         axis.title = element_text(size = 16),
                         axis.line = element_line(color="black", linewidth = 1.0),
                         panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
                         text = element_text(size = 13),
                         legend.position = c(0.84, 0.86),
                         legend.text = element_text(size = 15),
                         legend.title = element_text(size = 16)) 
  
  parent_path_sup <- paste(path, "sup_figures/", sep = "")
  
  
# BASELINE FIGURE S1
  parent_path_s1 <- paste0(parent_path_sup, "s1/")
  
  f <- list.files(path = parent_path_s1, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  # paste all data together
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path_s1, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate as column 
    local_data$rep <- rep
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbAgeSpecificGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbAgeSpecificGenes = factor(mutationProbAgeSpecificGenes))
    
    # gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough parents are present to sample 
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
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
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
    # rename model to be unique for replicate
    names(models)[i] <- paste0("model_", rep)
  }
  
  if (save_output){
    saveRDS(allData, paste0(output_path_data, "allData_s1.RDS"))
    saveRDS(allDataComplete, paste0(output_path_data, "allDataComplete_s1.RDS"))
    saveRDS(models, paste0(output_path_data, "models_s1.RDS"))
    
  }
  
  # NORMALIZE AND PREDICT 
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per mutation rate through the data 
  for (x in levels(allData$mutationProbAgeSpecificGenes)){
    # make subset of mutation rate  
    sub <- allData[allData$mutationProbAgeSpecificGenes == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProbAgeSpecificGenes = sub$mutationProbAgeSpecificGenes[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      # find corresponding model 
      mod <- paste0("model_", i)
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- round(quantile(sub$maternalAge, probs = 0.95))
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }  
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s1.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotal_s1.RDS"))
  }  
  
  #predDataTotalNormalized <- paste0(output_path_data, "predDataTotalNormalized_s1.RDS")
  #predDataTotal <- paste0(output_path_data, "predDataTotal_s1.RDS")
  
  # figure s1
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProbAgeSpecificGenes, 
                                      colour = mutationProbAgeSpecificGenes)) +
    geom_line() +
    labs( x = "Normalized parental age",
          y = "Normalized offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
                alpha = 0.2, colour = NA) +
    scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    theme_plots_s + 
    scale_y_continuous(breaks = seq(0,1.5,0.2),
                       limits = c(0,1.5)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0,1))
  
  # save image. 
  ggsave(paste0(output_path, "figure_s1.pdf"), width = 12, height = 8)
  
  # figure s1 - without normalization 
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProbAgeSpecificGenes, 
                            colour = mutationProbAgeSpecificGenes)) +
    geom_line() +
    labs( x = "Parental age",
          y = "Offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
                alpha = 0.2, colour = NA) +
    scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    theme_plots_s

# GAMETE DAMAGE ACCUMULATION - S2
  
  parent_path_s2 <- paste0(parent_path_sup, "s2/")
  
  f <- list.files(path = parent_path_s2, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path_s2, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate as column 
    local_data$rep <- rep
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProb[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProb = factor(mutationProb))
    
    # gather all data
    allDataComplete <- rbind(allDataComplete, local_data)  
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough parents are present to sample
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
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
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
    # rename model to be unique for replicate
    names(models)[i] <- paste0("model_", rep)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_s2.RDS"))
    saveRDS(models, file = paste0(output_path_data, "models_s2.RDS"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_s2.RDS"))
  }
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per mutation rate through the data 
  for (x in levels(allData$mutationProb)){
    # make subset of scenario 
    sub <- allData[allData$mutationProb == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProb = sub$mutationProb[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      mod <- paste0("model_", i) 
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

    # get 95th percentile 
    percentile <- round(quantile(sub$maternalAge, probs = 0.95))
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output) {
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s2.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotal_s2.RDS"))
  }
  
  #predDataTotalNormalized <- paste0(output_path_data, "predDataTotalNormalized_s2.RDS")
  #predDataTotal <- paste0(output_path_data, "predDataTotal_s2.RDS")
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProb, 
                                      colour = mutationProb)) +
    geom_line() +
    labs(x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProb), 
              alpha = 0.2, colour = NA) +
    scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    theme_plots_s +
    scale_y_continuous(breaks = seq(0,1.5,0.2),
                       limits = c(0,1.5)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0,1))
  
  # save image. 
  ggsave(paste0(output_path, "figure_s2.pdf"), width = 12, height = 8)
  
  # plot without normalization
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProb, 
                            colour = mutationProb)) +
    geom_line() +
    labs(x = "Parental age",
         y = "Offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotal,
              aes(ymin = min, ymax = max,  fill = mutationProb), 
              alpha = 0.2, colour = NA) +
    scale_colour_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    theme_plots_s
  
# PARENTAL CARE QUALITY - S3
  
  parent_path_s3 <- paste0(parent_path_sup, "s3/")
  
  f <- list.files(path = parent_path_s3, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path_s3, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate as column 
    local_data$rep <- rep
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbAgeSpecificGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbAgeSpecificGenes = factor(mutationProbAgeSpecificGenes))
    
    # tp gather all data
    allDataComplete <- rbind(allDataComplete, local_data)  
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough parents are present to sample from 
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
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
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
    # rename model to be unique for replicate
    names(models)[i] <- paste0("model_", rep)
  }

  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_s3.RDS"))
    saveRDS(models, file = paste0(output_path_data, "models_s3.RDS"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_s3.RDS"))
  }
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))

  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per parameter setting through the data 
  for (x in levels(allData$mutationProbAgeSpecificGenes)){
    # make subset of this parameter setting 
    sub <- allData[allData$mutationProbAgeSpecificGenes == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProbAgeSpecificGenes = sub$mutationProbAgeSpecificGenes[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      mod <- paste0("model_", i)
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- round(quantile(sub$maternalAge, probs = 0.95))
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes between 0 and 1 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s3.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotalNormalized_s3.RDS"))
  }
  
  #predDataTotalNormalized <- paste0(output_path_data, "predDataTotalNormalized_s3.RDS")
  #predDataTotal <- paste0(output_path_data, "predDataTotalNormalized_s3.RDS")
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProbAgeSpecificGenes, 
                                      colour = mutationProbAgeSpecificGenes)) +
    geom_line() +
    scale_color_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:5]) +
    labs(x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) +
    theme_plots_s +
    scale_y_continuous(breaks = seq(0,1.5,0.2),
                       limits = c(0,1.5)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0,1))
  
  ggsave(paste0(output_path, "figure_s3.pdf"), width = 12, height = 8)
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProbAgeSpecificGenes, 
                            colour = mutationProbAgeSpecificGenes)) +
    geom_line() +
    scale_color_manual("Mutation rate", values = met.brewer("Egypt")[1:4]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:4]) +
    labs(x = "Parental age",
         y = "Offspring age at death") +
    theme_minimal() +
    geom_ribbon(data = predDataTotal,
              aes(ymin = min, ymax = max,  fill = mutationProbAgeSpecificGenes), 
              alpha = 0.2, colour = NA) +
    theme_plots_s
  
# RESOURCE DISTRIBUTION - S4
  
  parent_path_s4 <- paste0(parent_path_sup, "s4/")
  
  f <- list.files(path = parent_path_s4, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  counter = 0
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path_s4, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data$rep <- rep
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbInvestmentGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbInvestmentGenes = factor(mutationProbInvestmentGenes))
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough parents are present to sample from 
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
    d2 <- d %>% mutate(y2 = ageAtDeath/40)
    ## Work with logits (then (0,1) -> (-inf,+inf))
    d2 <- d2 %>% mutate(y3 = car::logit(y2))
    d2 <- d2 %>% mutate(ID = factor(ID)) 
    d2 <- d2 %>% mutate(rep = factor(rep)) 
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
    # rename model to be unique for replicate
    names(models)[i] <- paste0("model_", rep)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_s4.RDS"))
    saveRDS(models, file = paste0(output_path_data, "models_s4.RDS"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_s4.RDS"))
  }
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))

  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  # go per mutation rate through the data 
  for (x in levels(allData$mutationProbInvestmentGenes)){
    # make subset of mutation rate  
    sub <- allData[allData$mutationProbInvestmentGenes == x,]
    # refactor the replicates of the subset 
    sub <- sub %>% mutate(rep = factor(rep))
    # make data expansion 
    pred_data = expand.grid(maternalAge =seq(0,max(sub$maternalAge),1),
                            ID = levels(sub$ID)[1],
                            mutationProbInvestmentGenes = sub$mutationProbInvestmentGenes[1]) 
    
    # loop through the replicates
    for (i in levels(sub$rep)){
      mod <- paste0("model_", i)
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- round(quantile(sub$maternalAge, probs = 0.95))
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s4.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotal_s4.RDS"))
  }
  
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_s4.RDS"))
  #predDataTotal <- loadRDS(paste0(output_path_data, "predDataTotal_s4.RDS"))
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                      group=mutationProbInvestmentGenes, 
                                      colour = mutationProbInvestmentGenes)) +
    geom_line() +
    scale_color_manual("Mutation rate", values = met.brewer("Egypt")[1:4]) +
    scale_fill_manual("Mutation rate", values = met.brewer("Egypt")[1:4]) +
    labs(x = "Normalized parental age",
         y = "Normalized offspring age at death") +
    theme_minimal() +  
    geom_ribbon(data = predDataTotalNormalized,
              aes(ymin = min, ymax = max,  fill = mutationProbInvestmentGenes), 
              alpha = 0.2, colour = NA) +
    theme_plots_s +
    scale_y_continuous(breaks = seq(0,1.9,0.2),
                       limits = c(0,1.9)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0,1))
  
  ggsave(paste0(output_path, "figure_s4.pdf"), width = 12, height = 8)
  
  # plot the normalized data grouped by mutation prob
  ggplot(predDataTotal, aes(maternalAge, mean, 
                            group=mutationProbInvestmentGenes, 
                            colour = mutationProbInvestmentGenes)) +
    geom_line() +
    scale_color_manual(values = met.brewer("Egypt", 11)) +
    labs(x = "Parental age",
         y = "Offspring age at death") +
    theme_minimal() +  
    geom_ribbon(data = predDataTotal,
              aes(ymin = min, ymax = max,  fill = mutationProbInvestmentGenes), 
              alpha = 0.2, colour = NA) +
    theme_plots_s

  
###############################################################################
# Supplementary materials.
# for figures S5 - S7
###############################################################################
  
# RESOURCE-ONLY - S5
  
  parent_path <- paste0(path, "matrix_data/resource_distribution/")

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
    local_data <- read.table(file_name, header = TRUE)
    # add the scenario as a column to the data
    local_data$scenario <- scenario
    # add replicate as column 
    local_data$rep <- rep
    # make ID unique
    local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
    local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
    
    # paste the data together 
    found <- rbind(found, local_data)
  }
  
  # determine the mean per replicate for every age.  
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$investmentGeneVal, sub$age, mean)
  }
  
  # get means over all replicates 
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  
  # cut off at 95th percentile
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] 
  
  # make plot 
  ggplot(data_per_rep, aes(age, (1-mean))) + geom_line() +
    geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
    theme_minimal() +
    labs(x = "Age",
         y = "Averaged proportion of resources 
  allocated to reproduction") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    theme(axis.text = element_text(size=17,face="plain",color="black"), 
          axis.title = element_text(size = 17)) +
    theme(axis.text = element_text(size=15,face="plain",color="black"),
          axis.title = element_text(size = 16),
          axis.line = element_line(color="black", linewidth = 1.0),
          panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
          text = element_text(size = 13),
          legend.position = c(0.84, 0.86),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 16)) 
  
  scale_x_continuous(breaks = seq(0, 9, 3),
                     limits = c(0, 9))
  
  ggsave(paste0(output_path, "figure_s5.pdf"), width = 12, height = 8)

# RESOURCE + BASELINE - S6
  
  parent_path <- paste0(path, "matrix_data/baseline_resource/")
  
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
    local_data <- read.table(file_name, header = TRUE)
    # add the scenario as a column to the data
    local_data$scenario <- scenario
    # add replicate as column 
    local_data$rep <- rep
    # make ID unique
    local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
    local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
    
    # paste the data together 
    found <- rbind(found, local_data)
  }
  
  # determine the mean per replicate per age 
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$investmentGeneVal, sub$age, mean)
  }
  
  # get means over all replicates 
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  
  # cut off at 95th percentile
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] 
  
  # make plot 
  ggplot(data_per_rep, aes(age, (1-mean))) + geom_line() +
    geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
    theme_minimal() +
    labs(x = "Age",
         y = "Averaged proportion of resources 
  allocated to reproduction") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    theme(axis.text = element_text(size=17,face="plain",color="black"), 
          axis.title = element_text(size = 17)) +
    theme(axis.text = element_text(size=15,face="plain",color="black"),
          axis.title = element_text(size = 16),
          axis.line = element_line(color="black", linewidth = 1.0),
          panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
          text = element_text(size = 13),
          legend.position = c(0.84, 0.86),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 16)) 
  
  ggsave(paste0(output_path, "figure_s6.pdf"), width = 12, height = 8)
  
# RESOURCE + QUALITY
  
  parent_path <- paste0(path, "matrix_data/resource_and_parentalQuality/")
  
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
    local_data <- read.table(file_name, header = TRUE)
    # add the scenario as a column to the data
    local_data$scenario <- scenario
    # add replicate as column 
    local_data$rep <- rep
    # make ID unique
    local_data$ID <- sub("^", local_data$scenario[1], local_data$ID)
    local_data$ID <- sub("^", local_data$rep[1], local_data$ID)
    
    # paste the data together 
    found <- rbind(found, local_data)
  }
  
  # determine the mean per replicate per age 
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$survivalGeneVal, sub$age, mean)
  }
  
  # get means over all replicates 
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  
  # cut off at 95th percentile
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] 
  
  # make plot 
  ggplot(data_per_rep, aes(age, (1-mean))) + geom_line() +
    geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
    theme_minimal() +
    labs(x = "Age") +
         #y = "Averaged proportion of resources 
  #allocated to reproduction") +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    theme(axis.text = element_text(size=17,face="plain",color="black"), 
          axis.title = element_text(size = 17)) +
    theme(axis.text = element_text(size=15,face="plain",color="black"),
          axis.title = element_text(size = 16),
          axis.line = element_line(color="black", linewidth = 1.0),
          panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
          text = element_text(size = 13),
          legend.position = c(0.84, 0.86),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 16)) 
  
  ggsave(paste0(output_path, "figure_s8.pdf"), width = 12, height = 8)
  
###############################################################################
# Supplementary materials.
# for figure S8
###############################################################################
  
  # set axes 
  plotsTot$p_baseline_gamete_quality <- plotsTot$p_baseline_gamete_quality + theme(axis.text.x = element_blank())
  plotsTot$p_baseline_gamete_resource <- plotsTot$p_baseline_gamete_resource + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  plotsTot$p_baseline_quality_resource <- plotsTot$p_baseline_quality_resource + theme(axis.text.y = element_blank())
  
  # make matrix
  plot_matrix <- plot_grid(plotsTot$p_baseline_gamete_quality,
                           plotsTot$p_baseline_gamete_resource,
                           plotsTot$p_gamete_quality_resource,
                           plotsTot$p_baseline_quality_resource,
                           align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.90, label_y = 0.97)
  
  # make legend 
  p_tmp <- plotsTot$p_baseline_quality_resource + theme(legend.position = "bottom",
                                                         legend.title = element_blank(),
                                                         legend.text = element_text(size=13))
  leg1 <- get_legend(p_tmp) 
  legend <- as_ggplot(leg1)
  legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
  blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
  legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)
  
  # make matrix with axes labeling
  plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                           left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))
  
  # make matrix with axes and legend 
  plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
  plot_with_legend
  
  # save figure
  ggsave(paste0(output_path, "Figure_S8.pdf"), width = 12, height = 10)
  
###############################################################################
# Supplementary materials.
# for figure S9
###############################################################################
  
  plotsTot$p_baseline_gamete_quality_resource
  ggsave(paste0(output_path, "Figure_S9.pdf"), width = 12, height = 10)
  
###############################################################################
# Resource distribution matrix 
###############################################################################
  # define equation 1
  eq1 <- function(x, c) {
    1 - c * x^2
  }
  
  # make dummy data
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
  geom_line(linewidth = 1) +
  scale_colour_manual(name = "", values = met.brewer("Archambault")[c(1, 4, 7)],
                      labels = c("c = 0.1", "c = 0.3", "c = 0.7")) +
  theme_minimal() +
  labs(x = expression(paste("Allocation to reproduction (", italic("x"),")")),
       y = expression(paste("Effect on parental survival (", m["4a"], ")"))) +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.key.width = unit(1, "cm")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1)) 
  
  # define equation 2
  eq2 <- function(x, a, d) {
    1 / (1 + exp(-a * x - d))
  }
  
  # make dummy data 
  a <- c(1, 3, 10)
  d <- c(2.5, 1, -4)
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
  geom_line(linewidth = 1) +
  scale_colour_manual("", values = met.brewer("Archambault")[c(1, 4, 7)],
                      labels = c("a = 1, d = 2.5", "a = 3, d = 1", "a = 10, d = -4")) +
  theme_minimal() +
  labs(x = expression(paste("Allocation to reproduction (", italic("x"),")")),
       y = expression(paste("Effect on offspring's survival (", m["4b"], ")"))) +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.key.width = unit(1, "cm")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1)) 
  
################################################################################
# EQUATION 1
################################################################################
  # general theme 
  theme_plots <- theme(axis.text = element_text(size=15,face="plain",color="black"),
                       axis.text.x = element_blank(),
                       axis.title = element_text(size = 16),
                       axis.line = element_line(color="black", linewidth = 1.0),
                       panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
                       text = element_text(size = 13),
                       legend.position = "none",
                       plot.title = element_text(hjust = 0.5))

  # set path to data for the resource matrix. 
  path_rema <- paste0(path, "resource_matrix/")
  
  # the gam model
  run_model <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
    return(tmp)
  }
  
  # general ggplot format 
  color_eq1 <- scale_color_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) 
  fill_eq1 <- scale_fill_manual(values = met.brewer("Archambault")[c(1, 4, 7)]) 
  scale_y_eq1_lansing <- scale_y_continuous(breaks = seq(0, 2, 1), limits = c(0, 2))
  scale_y_eq1_gene_vals <- scale_y_continuous(breaks = seq(0, 1), limits = c(0, 1))
  scale_x <- scale_x_continuous(breaks = seq(0,1, 0.2), limits = c(0,1))
  
  # set path to eq1 
  path_rema_eq1 <- paste0(path_rema, "eq1/")
  
# RESOURCE-ONLY
  parent_path <- paste0(path_rema_eq1, "resource_eq1/") 
  name = "resource_eq1"
  
  f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  found <- c()
  counter = 0
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
    
    # to gather all data
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    # filter parents that had offspring at at least 5 different ages 
    local_data <- local_data %>% filter(na >= 5)
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
    # rename model to be unique for replicate
    names(models)[i] <- paste0("model_", rep)
    
    # COLLECT DATA FOR AGE-SPECIFIC FIGURE
    # get path 
    file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
    # read data
    local_data_age <- read.table(file_name, header = T)
    # add replicate as column 
    local_data_age$rep <- rep
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }

  if (save_output) {
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
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
      mod <- paste0("model_", i)
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
    
    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save the percentile for this subset 
    percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the normalized data 
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  name = "resource_eq1"
  predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
  
  # plot the normalized data grouped by weight investment
  p3 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                            group = weightInvestment, 
                                            colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL,
         title = "Normalized offspring lifespan") 
  
  # make age-specific figure 
  
  # loop through every parameter value  
  combine_data <- c()
  for (i in levels(as.factor(found$weightInvestment))){
    sub <- found[found$weightInvestment == i,]
    # get the 95th percentile. 
    percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
    # expand data from age 0 to percentile age 
    tmp <- expand.grid(age = c(0:percentile),
                       weightInvestment = sub$weightInvestment[1])
    sub <- sub[sub$age <= percentile,]
    
    # loop through every replicate 
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
    
    # normalize the x-axis
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p4 <- ggplot(combine_data, aes(age, mean, 
                                 group= weightInvestment, 
                                 colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL,
         title = "Allocation to reproduction") 
  
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
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    # filter parents that had offspring at at least 5 different ages 
    local_data <- local_data %>% filter(na >= 5)
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
    names(models)[i] <- paste0("model_", rep)
    
    # COLLECT DATA FOR AGE-SPECIFIC FIGURE
    # get path 
    file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
    # read data
    local_data_age <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data_age$rep <- rep
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }

  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
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
      mod <- paste0("model_", i)
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:3,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 
    
    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save the percentile for this subset 
    percentiles[percentiles$weightInvestment == x,]$percentile <- percentile
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the normalized data 
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_baseline_eq1"
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
  
  # plot the normalized data grouped by weight investment
  p7 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                            group = weightInvestment, 
                                            colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL) 
  
  # make age-specific figure 
  
  # determine the mean per replicate per parameter value. 
  combine_data <- c()
  for (i in levels(as.factor(found$weightInvestment))){
    sub <- found[found$weightInvestment == i,]
    # get percentile for cut-off
    percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       weightInvestment = sub$weightInvestment[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through every replicate and determine mean gene value per age
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p8 <- ggplot(combine_data, aes(age, mean, 
                                 group= weightInvestment, 
                                 colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL) 
  
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
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data
    allDataComplete <- rbind(allDataComplete, local_data)  
    
    # perform gam
    # filter parents that had offspring at at least 5 different ages 
    local_data <- local_data %>% filter(na >= 5)
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
    names(models)[i] <- paste0("model_", rep)
    
    # COLLECT DATA FOR AGE-SPECIFIC FIGURE
    # get path 
    file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
    # read data
    local_data_age <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data_age$rep <- rep
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
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
      mod <- paste0("model_", i)
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
    # normalize the axes
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    # save the data 
    pred_data_tmp <- pred_data[c(1:3, ((ncol(pred_data)-2):ncol(pred_data)))] # select only relevant info 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_gamete_eq1"
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
  
  # plot the normalized data grouped by weight investment
  p11 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                             group = weightInvestment, 
                                             colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL)
  
  # make age-specific figure 
  
  # determine the mean per replicate. 
  combine_data <- c()
  for (i in levels(as.factor(found$weightInvestment))){
    sub <- found[found$weightInvestment == i,]
    # get percentile for cut-off
    percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       weightInvestment = sub$weightInvestment[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through every replicate and determine mean per age
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p12 <- ggplot(combine_data, aes(age, mean, 
                                  group= weightInvestment, 
                                  colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL)
  
# RESOURCE + QUALITY
  
  # resource-only 
  parent_path <- paste0(path_rema_eq1, "resource_parentalQuality_eq1/") 
  name = "resource_quality_eq1"
  
  f <- list.files(path = parent_path, pattern = "outputLifeExpLong.txt", recursive = T, all.files = T)
  allData <- c()
  allDataComplete <- c()
  models <- list()
  found <- c()
  counter = 0
  for (i in 1:length(f)) {
    x <- f[i]
    # get path 
    file_name <- paste0(parent_path, x)
    # extract replicate from file name 
    splitted_path <- strsplit(file_name, "/")
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    # filter parents that had offspring at at least 5 different ages 
    local_data <- local_data %>% filter(na >= 5)
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
    names(models)[i] <- paste0("model_", rep)
    
    # COLLECT DATA FOR AGE-SPECIFIC FIGURE
    # get path 
    file_name <- paste0(parent_path, rep, "/outputWithAgeSpecificGenes.txt")
    # read data
    local_data_age <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data_age$rep <- rep
    # TODO: remove this row if data is complete. 
    local_data_age$weightInvestment <- local_data$weightInvestment[1]
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$weightInvestment[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }

  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
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
      mod <- paste0("model_", i)
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
  
  if (save_output){
    saveRDS(predDataTotalNormalized, paste0(output_path_data, "predDataTotalNormalized_", name, ".RDS"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_quality_eq1"
  #predDataTotalNormalized <- loadRDS(paste0("predDataTotalNormalized_", name, ".RDS"))
  
  # plot the normalized data grouped by weight investment
  p15 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                             group = weightInvestment, 
                                             colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_lansing +
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    labs(x = NULL,
         y = NULL) +
    theme(axis.text.x = element_text(colour = "black", size =16))
  
  
  # make age-specific figure 
  
  # determine the mean per replicate. 
  combine_data <- c()
  for (i in levels(as.factor(found$weightInvestment))){
    sub <- found[found$weightInvestment == i,]
    # get percentile as cut-off
    percentile <- percentiles[percentiles$weightInvestment == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       weightInvestment = sub$weightInvestment[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through every replicate 
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(weightInvestment = factor(weightInvestment))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p16 <- ggplot(combine_data, aes(age, mean, 
                                  group= weightInvestment, 
                                  colour = weightInvestment)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = weightInvestment), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq1 +
    fill_eq1 +
    scale_y_eq1_gene_vals +
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    labs(x = NULL,
         y = NULL) +
    theme(axis.text.x = element_text(colour = "black", size =16))
  
  
################################################################################
# EQUATION 2 
################################################################################
  
  # set path 
  path_rema_eq2 <- paste0(path_rema, "eq2/")
  
  # the gam model
  run_model <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k) + s(maternalAge, ID, bs = "fs", k = z), data = data, method = "REML")
    return(tmp)
  }
  
  # for the lansing figures (column 3)
  color_eq2_lansing <- scale_color_manual(values = met.brewer("Archambault")[c(7, 4, 1)]) 
  fill_eq2_lansing <- scale_fill_manual(values = met.brewer("Archambault")[c(7, 4, 1)]) 
  # for the gene values figures (column 4)
  color_eq2_gene_vals <- scale_color_manual(values = met.brewer("Archambault")[c(4, 1, 7)]) 
  fill_eq2_gene_vals <- scale_fill_manual(values = met.brewer("Archambault")[c(4, 1, 7)])
  # scale y for Lansing figures (column 3)
  scale_y_eq2_lansing <- scale_y_continuous(breaks = seq(0, 2, 1), limits = c(0, 2))
  # scale y for the gene value figures (column 4)
  scale_y_eq2_gene_vals <- scale_y_continuous(breaks = seq(0, 1), limits = c(0, 1))
  
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
    # extract which setting folder the file is in 
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
    
    # to gather all data
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough IDs are present to sample 
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
  percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
  colnames(percentiles) <- c("steepnessAllocationToReproduce")
  percentiles$percentile <- 0
  
  # go per parameter setting through the data 
  for (x in levels(allData$steepnessAllocationToReproduce)){
    # make subset of parameter setting 
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save the percentile value 
    percentiles[percentiles$steepnessAllocationToReproduce == x,]$percentile <- percentile
    # cut-off the data at that percentile 
    pred_data <- pred_data[pred_data$maternalAge <= percentile,]
    # normalize the axes 
    pred_data$maternalAge <- pred_data$maternalAge / percentile
    pred_data$min <- pred_data$min / pred_data$mean[1]
    pred_data$max <- pred_data$max / pred_data$mean[1]
    pred_data$mean <- pred_data$mean / pred_data$mean[1]
    
    pred_data_tmp <- pred_data[c(1:4, ((ncol(pred_data)-2):ncol(pred_data)))]
    pred_data_tmp$nRep <- length(levels(sub$rep))
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_eq2"
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by weight investment
  p5 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean,  
                                            group = steepnessAllocationToReproduce, 
                                            colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max, fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_lansing +
    fill_eq2_lansing +
    scale_y_eq2_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL,
         title = "Lansing effect") 
  
  # determine the mean per replicate per parameter setting. 
  combine_data <- c()
  for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
    sub <- found[found$steepnessAllocationToReproduce == i,]
    # get the corresponding percentile to use as cut-off
    percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through every replicate and determine mean per age 
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p6 <- ggplot(combine_data, aes(age, mean, 
                                 group= steepnessAllocationToReproduce, 
                                 colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_gene_vals +
    fill_eq2_gene_vals +
    scale_y_eq2_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL,
         title = "Allocation to reproduction") 
  
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
    # get the parameter folder setting 
    nVal <- length(splitted_path[[1]]) - 2
    val <- splitted_path[[1]][nVal]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data
    allDataComplete <- rbind(allDataComplete, local_data)  
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough IDs are present to sample 
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
  percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
  colnames(percentiles) <- c("steepnessAllocationToReproduce")
  percentiles$percentile <- 0
  
  # go per parameter setting through the data 
  for (x in levels(allData$steepnessAllocationToReproduce)){
    # make subset of data with this parameter setting  
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save the percentile value 
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
    
    # save the normalized data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_baseline_eq2"
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by weight investment
  p9 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                            group = steepnessAllocationToReproduce, 
                                            colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_lansing +
    fill_eq2_lansing +
    scale_y_eq2_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL) 
  
  # determine the mean per replicate per parameter setting. 
  combine_data <- c()
  for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
    sub <- found[found$steepnessAllocationToReproduce == i,]
    # get the percentile to use as cut-off 
    percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through the data per replicate 
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p10 <- ggplot(combine_data, aes(age, mean, 
                                  group= steepnessAllocationToReproduce, 
                                  colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_gene_vals +
    fill_eq2_gene_vals +
    scale_y_eq2_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL) 
  
  
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
    # extract the parameter setting folder 
    nVal <- length(splitted_path[[1]]) - 2
    val <- splitted_path[[1]][nVal]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough IDs are present to sample 
    to_sample <- 100 
    if (length(unique(local_data$ID)) < to_sample){
      counter = counter + 1
      next
    }
    
    # sample 100 parents by their IDs
    local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = FALSE),]
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
    local_data_age <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data_age$rep <- rep
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
  percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
  colnames(percentiles) <- c("steepnessAllocationToReproduce")
  percentiles$percentile <- 0
  
  # go per parameter setting through the data 
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save percentile 
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
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  #name = "resource_gamete_eq2"
  #predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by weight investment
  p13 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                             group = steepnessAllocationToReproduce, 
                                             colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_lansing +
    fill_eq2_lansing +
    scale_y_eq2_lansing +
    scale_x +
    labs(x = NULL,
         y = NULL)
  
  # determine the mean per replicate per parameter setting. 
  combine_data <- c()
  for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
    sub <- found[found$steepnessAllocationToReproduce == i,]
    # get corresponding percentile to use as cut-off
    percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through the data per replicate 
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
    
    # normalize the x-axis 
    tmp$age <- tmp$age / percentile
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p14 <- ggplot(combine_data, aes(age, mean, 
                                  group= steepnessAllocationToReproduce, 
                                  colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_gene_vals +
    fill_eq2_gene_vals +
    scale_y_eq2_gene_vals +
    scale_x +
    labs(x = NULL,
         y = NULL) 
  
# RESOURCE + PARENTAL CARE QUALITY
  
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
    # extract the parameter setting folder 
    nVal <- length(splitted_path[[1]]) - 2
    val <- splitted_path[[1]][nVal]
    # read data
    local_data <- read.table(file_name, header = TRUE)
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
    
    # to gather all data 
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # check if enough individuals are present to sample 
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
    local_data_age <- read.table(file_name, header = TRUE)
    # add replicate as column 
    local_data_age$rep <- rep
    local_data_age$ID <- sub("^", local_data_age$rep[1], local_data_age$ID)
    local_data_age$ID <- sub("^", local_data_age$steepnessAllocationToReproduce[1], local_data_age$ID)
    found <- rbind(found, local_data_age)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0(output_path_data, "allData_", name, ".rds"))
    saveRDS(models, file = paste0(output_path_data, "models_", name, ".rds"))
    saveRDS(allDataComplete, file = paste0(output_path_data, "allDataComplete_", name, ".rds"))
    saveRDS(found, file = paste0(output_path_data, "found_", name, ".rds"))
  }
  
  allData$rep <- factor(allData$rep)
  
  # to save all data 
  predDataTotal <- c()
  predDataTotalNormalized <- c()
  
  # to save the percentiles where the data is cut-off per scenario
  percentiles <- data.frame(levels(allData$steepnessAllocationToReproduce))
  colnames(percentiles) <- c("steepnessAllocationToReproduce")
  percentiles$percentile <- 0
  
  # go per parameter setting through the data 
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
    predDataTotal <- rbind(predDataTotal, 
                           pred_data[c(1:4,(ncol(pred_data)-2):ncol(pred_data))]) # select only columns with the relevant information 

    # get 95th percentile 
    percentile <- quantile(sub$maternalAge, probs = 0.95)
    # save percentile 
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
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
    saveRDS(predDataTotal, paste0(output_path_data, "predDataTotal_", name, ".RDS"))
    saveRDS(percentiles, paste0(output_path_data, "percentiles_", name, ".RDS"))
  }
  
  name = "resource_and_parentalQuality_eq2"
  predDataTotalNormalized <- loadRDS(paste0(output_path_data, "predDataTotalNormalized_", name, ".rds"))
  
  # plot the normalized data grouped by weight investment
  p17 <- ggplot(predDataTotalNormalized, aes(maternalAge, mean, 
                                             group = steepnessAllocationToReproduce, 
                                             colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = predDataTotalNormalized,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_lansing +
    fill_eq2_lansing +
    scale_y_eq2_lansing +
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    labs(x = NULL,
         y = NULL) +
    theme(axis.text.x = element_text(colour = "black", size =16))
  
  # determine the mean per replicate. 
  combine_data <- c()
  for (i in levels(as.factor(found$steepnessAllocationToReproduce))){
    sub <- found[found$steepnessAllocationToReproduce == i,]
    # get the percentile to use as cut-off 
    percentile <- percentiles[percentiles$steepnessAllocationToReproduce == i,]$percentile
    tmp <- expand.grid(age = c(0:percentile),
                       steepnessAllocationToReproduce = sub$steepnessAllocationToReproduce[1])
    sub <- sub[sub$age <= percentile,]
    
    # go through the data per replicate 
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
    
    combine_data <- rbind(combine_data, tmp[c(1:2,(ncol(tmp)-2):ncol(tmp))])
  }
  
  combine_data <- combine_data %>% mutate(steepnessAllocationToReproduce = factor(steepnessAllocationToReproduce))
  
  if (save_output){
    saveRDS(combine_data, paste0(output_path_data, "combine_data_", name, ".RDS"))
  }
  
  #combine_data <- loadRDS(paste0(output_path_data, "combine_data_", name, ".RDS"))
  
  # plot 
  p18 <- ggplot(combine_data, aes(age, mean, 
                                  group= steepnessAllocationToReproduce, 
                                  colour = steepnessAllocationToReproduce)) +
    geom_line() +
    theme_minimal() +
    geom_ribbon(data = combine_data,
                aes(ymin = min, ymax = max,  fill = steepnessAllocationToReproduce), 
                alpha = 0.2, colour = NA) +
    theme_plots +
    color_eq2_gene_vals +
    fill_eq2_gene_vals +
    scale_y_eq2_gene_vals +
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    labs(x = NULL,
         y = NULL) +
    theme(axis.text.x = element_text(colour = "black", size =16))

  
  ############ make grid ##################
  
  top_row <- plot_grid(p1, p2, ncol = 2, labels = "AUTO")
  labels <- c("C", "D", "E", "F",
              "G", "H", "I", "J",
              "K", "L", "M", "N",
              "O", "P", "Q", "R")
  bottom_part <- plot_grid(p3, p4, p5, p6,
                           p7, p8, p9, p10,
                           p11, p12, p13, p14,
                           p15, p16, p17, p18,
                           ncol = 4, labels = labels, 
                           align = "hv") #,
                          # label_x = 0.1, label_y = 0.9, 
                          # align = "hv")
  
  resource_matrix <- plot_grid(top_row, bottom_part, ncol = 1, rel_heights = c(1, 2.5))
  resource_matrix
  
  plot_matrix2 <- grid.arrange(arrangeGrob(resource_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)))) 

  plot_matrix2
  ggsave(paste0(output_path, "resource_matrix.pdf"), width = 15, height = 20)
  