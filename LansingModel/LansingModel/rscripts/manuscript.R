###############################################################################


# MANUSCRIPT R CODE 
# Created at 13-09-2023


###############################################################################

path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/data/"
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/Manuscript/Figures/"
save_output = TRUE

###############################################################################
# Figure S2 - the matrix
###############################################################################

# 1. MAKING GAM MODELS FOR EVERY REPLICATE PER SCENARIO 

  # load the elements - lines TO_FILL to TO_FILL not necessary 
  allData <- loadRDS("allData.rds")
  models <- loadRDS("models.rds")
  
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
    
    allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
    
    local_data <- local_data %>% filter(na >= 5)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
    # sample 100 parents by their IDs
    local_data <- local_data[local_data$ID %in% sample(unique(local_data$ID), 100, replace = F),]
    # save the data
    allData <- rbind(allData, local_data)
    
    # Perform GAM on the data
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = "allData_matrix.rds")
    saveRDS(models, file = "models_matrix.rds")
    saveRDS(allDataComplete, file = "allDataComplete_matrix.rds")
  }
  
  #allData <- loadRDS("allData_matrix.rds")
  #models <- loadRDS("models_matrix.rds")
  
# 2. PREDICTING AND NORMALIZING THE DATA 

  # loop through the scenarios 
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
    
    # save the data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data)
  }
  
  if (save_output){
    # save the percentiles 
    saveRDS(percentiles, "percentiles_matrix.RDS")
    # save the normalized dataset 
    saveRDS(predDataTotalNormalized, "predDataTotalNormalized_matrix.RDS")
  }
  
# 3. PERFORMING CROSS-SECTIONAL GAM ANALYSIS 
  
  allDataLat <- loadRDS("allDataLat.rds")
  modelsLat <- loadRDS("modelsLat.rds")
  
  # the gam model
  run_model_lat <- function(data) {
    tmp <- bam(y3 ~ s(maternalAge, k = k), data=data,method = "REML")
    return(tmp)
  }
  
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
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allDataLat, file = "allDataLat_matrix.rds")
    saveRDS(modelsLat, file = "modelsLat_matrix.rds")
  }
  
  # loop through the scenarios 
  allDataLat$scenario <- factor(allDataLat$scenario)
  allDataLat$rep <- factor(allDataLat$rep)
  
# 4. PREDICTING AND NORMALIZING DATA   
  
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
  
# 5. DYNAMICALLY GENERATING PLOTS 
  
  # add column with group to both normalized data sets 
  predDataTotalNormalized$group <- "Longitudinal"
  predDataTotalNormalizedLat$group <- "Cross-sectional"
  # merge them together 
  totalNormalizedData <- rbind(predDataTotalNormalized, predDataTotalNormalizedLat)
  
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
      ylim(0, 1.5) + 
      theme_minimal() +
      theme(#legend.text = element_text(size=10),
        #legend.key.size = unit(0.2, "cm"),
        #legend.key.width = unit(0.1,"cm"),
        legend.position = "none",
        #legend.title = element_blank(),
        axis.text = element_text(size=13,face="plain",color="black"),
        axis.title = element_text(size = 16),
        #axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
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
  
# 6. MAKING THE PARENTAL AGE DISTRIBUTION PLOTS
 #  
 #  # load the saved complete data 
 #  allDataComplete <- loadRDS("allDataComplete_matrix.rds")
 #  # load the saved percentiles 
 #  percentiles <- loadRDS("percentiles_matrix.RDS")
 #  
 #  allDataComplete <- allDataComplete %>% mutate(scenario = factor(scenario))
 #  
 #  predDataTotalAgeDist <- c()
 #  # loop through the scenarios. 
 #  for (x in levels(allDataComplete$scenario)){
 #    print(x)
 #    # make subset of scenario 
 #    sub <- allDataComplete[allDataComplete$scenario == x,]
 #    # get percentile 
 #    cut_off <- percentiles[percentiles$scenario == x,]$percentile
 #    
 #    # make data expansion 
 #    pred_data = expand.grid(maternalAge =seq(0,cut_off,1),
 #                            scenario = sub$scenario[1]) 
 #    
 #    sub <- sub %>% mutate(rep = factor(rep))
 #    # loop through the replicates
 #    for (i in levels(sub$rep)){
 #      tmp <- sub[sub$rep == i,]
 #      
 #      # use this to count the number of parents per age class 
 #      counts <- c()
 #      for (j in seq(0, max(pred_data$maternalAge), 1)) {
 #        counts <- c(counts, nrow(tmp[tmp$maternalAge == j,]))
 #      }
 #      new_col <- paste("rep", i, sep = "")
 #      # add the counts
 #      pred_data[[new_col]] <- counts
 #    }
 #    
 #    # get means
 #    repVals <- subset(pred_data, select = 3:12)
 #    pred_data$mean <- rowMeans(repVals)
 #    pred_data$min <- apply(repVals, 1, min)
 #    pred_data$max <- apply(repVals, 1, max)
 #    
 #    # combine the data per scenario
 #    predDataTotalAgeDist <- rbind(predDataTotalAgeDist, pred_data)
 #  }
 #  
 #  # plots 
 #  # get the scenarios relevant for the matrix plot
 #  scenarios <- c("baseline_quality", "baseline_resource", "gamete_and_baseline", "gamete_quality",
 #               "gamete_resource", "resource_and_parentalQuality")
 #  
 #  # dynamically generate the plots
 #  plotsAgeDist <- c()
 #  for (i in 1:length(scenarios)){
 #    p <- ggplot(predDataTotalAgeDist[predDataTotalAgeDist$scenario == scenarios[i],], 
 #                aes(maternalAge, mean)) +
 #      geom_line() +
 #      labs(x = NULL,
 #           y = NULL) +
 #      geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, linewidth = 0.8) + 
 #      #ylim(0, 2) + 
 #      theme_minimal() +
 #      theme(#legend.text = element_text(size=10),
 #        #legend.key.size = unit(0.2, "cm"),
 #        #legend.key.width = unit(0.1,"cm"),
 #        legend.position = "none",
 #        #legend.title = element_blank(),
 #        axis.text = element_text(size=13,face="plain",color="black"), # size = 11
 #        axis.title = element_text(size = 13),
 #        #axis.text.x = element_blank(),
 #        axis.line = element_line(color="black", linewidth = 1.0),
 #        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) +
 #      scale_x_continuous(breaks = c(0, max(predDataTotalAgeDist[predDataTotalAgeDist$scenario == scenarios[i],]$maternalAge)), 
 #                         limits = c(0, max(predDataTotalAgeDist[predDataTotalAgeDist$scenario == scenarios[i],]$maternalAge)))
 #    
 # #   scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1))
 #    
 #    # save plot in list 
 #    plotsAgeDist[[i]] <- p
 #    # adjust name to corresponding scenario run
 #    names(plotsAgeDist)[i] <- paste("p_", scenarios[i], sep = "") 
 #  }
 #  
# 7. MAKE MATRIX 
  
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
  
  labels <- c("A", NA, NA, NA,
              "B", "C", NA, NA,
              "D", "E", "F", NA,
              "G", "H", "I", "J")
  
  plot_matrix <- plot_grid(plotsTot$p_baseline, NULL, NULL, NULL,
                           plotsTot$p_gamete_and_baseline, plotsTot$p_gamete_damage, NULL, NULL,
                           plotsTot$p_baseline_quality, plotsTot$p_gamete_quality, plotsTot$p_parental_care_quality, NULL, 
                           plotsTot$p_baseline_resource, plotsTot$p_gamete_resource, plotsTot$p_resource_and_parentalQuality, plotsTot$p_resource_distribution,
                           ncol = 4, align = "hv", axis = "brlt", labels = labels, label_x = 0.82, label_y = 0.9)
  
  plot_matrix
  
  p_tmp <- plotsTot$p_resource_distribution + theme(legend.position = "bottom",
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
  ggsave(paste0(output_path, "Lansing_fig_with_labels.pdf"), width = 15, height = 15)

  
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
    
    to_sample <- 100 
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
    
    # save the data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }  
  
  if (save_output){
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s1.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotal_s1.RDS"))
  }  
  
  predDataTotalNormalized <- paste0(output_path_data, "predDataTotalNormalized_s1.RDS")
  predDataTotal <- paste0(output_path_data, "predDataTotal_s1.RDS")
  
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
  
  # save image. 
  ggsave(paste0(output_path, "figure_s1.pdf"), width = 12, height = 8)

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
    
    allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
    
    # perform gam
    local_data <- local_data %>% filter(na > 6)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
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
  # go per scenario through the data 
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
    
    # save the data 
    predDataTotalNormalized <- rbind(predDataTotalNormalized, pred_data_tmp)
  }
  
  if (save_output) {
    saveRDS(predDataTotalNormalized, file = paste0(output_path_data, "predDataTotalNormalized_s2.RDS"))
    saveRDS(predDataTotal, file = paste0(output_path_data, "predDataTotal_s2.RDS"))
  }
  
  predDataTotalNormalized <- paste0(output_path_data, "predDataTotalNormalized_s2.RDS")
  predDataTotal <- paste0(output_path_data, "predDataTotal_s2.RDS")
  
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
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    # extract replicate from file name 
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
    
    allDataComplete <- rbind(allDataComplete, local_data) # for the maternal age distribution plot 
    
    # perform gam
    local_data <- local_data %>% filter(na > 6)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
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
    
    # save the data 
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
    # extract scenario from file name 
    splitted_path <- strsplit(file_name, "/")
    # extract replicate from file name 
    nRep <- length(splitted_path[[1]]) - 1
    rep <- splitted_path[[1]][nRep]
    # read data
    local_data <- read.table(file_name, header = T)
    # add replicate and scenario as column 
    local_data$rep <- rep
    # rename ID so it is unique per replicate per parameter setting
    local_data$ID <- sub("^", local_data$mutationProbInvestmentGenes[1], paste0("_", local_data$ID))
    local_data$ID <- sub("^", local_data$rep[1], paste("_", local_data$ID, sep = ""))
    # subset parents 
    z <- local_data %>% dplyr::group_by(ID) %>% dplyr::summarize(na = length(unique(maternalAge))) 
    ## Add the counter to data
    local_data <- local_data %>% left_join(z,by="ID")
    local_data <- local_data %>% mutate(mutationProbInvestmentGenes = factor(mutationProbInvestmentGenes))
    
    allDataComplete <- rbind(allDataComplete, local_data) 
    
    # perform gam
    local_data <- local_data %>% filter(na > 6)
    local_data <- local_data %>% mutate(ID = factor(ID)) 
    
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
    # rename model to be unique for replicate and scenario
    names(models)[i] <- paste0("model_", rep)
  }
  
  if (save_output){
    # save the R data just in case. 
    saveRDS(allData, file = paste0("allData_s4.RDS"))
    saveRDS(models, file = paste0("models_s4.RDS"))
    saveRDS(allDataComplete, file = paste0("allDataComplete_s4.RDS"))
  }
  
  allData$rep <- factor(allData$rep)
  logist <- function(x) 40/(1+exp(-x))

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
    
    # save the data 
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
  
  # determine the mean per replicate. 
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$investmentGeneVal, sub$age, mean)
  }
  
  # get means
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] # cut off at 95th percentile
  
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
  
  # determine the mean per replicate. 
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$investmentGeneVal, sub$age, mean)
  }
  
  # get means
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] # cut off at 95th percentile
  
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
  
  # determine the mean per replicate. 
  age <- c(0:max(found$age))
  data_per_rep <- data.frame(age)
  for (i in levels(as.factor(found$rep))){
    sub <- found[found$rep == i,]
    new_col <- paste0("rep", i)
    data_per_rep[[new_col]] <- tapply(sub$survivalGeneVal, sub$age, mean)
  }
  
  # get means
  repVals <- subset(data_per_rep, select = 2:ncol(data_per_rep))
  data_per_rep$mean <- rowMeans(repVals)
  data_per_rep$min <- apply(repVals, 1, min)
  data_per_rep$max <- apply(repVals, 1, max)
  cut_off <- percentiles[percentiles$scenario == found$scenario[1],]$percentile
  data_per_rep <- data_per_rep[data_per_rep$age <= cut_off,] # cut off at 95th percentile
  
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
  
  # to remove: to check the survival genes
  ggplot(data_per_rep, aes(age, mean)) + geom_line() +
    geom_ribbon(aes(ymin =min, ymax = max), alpha = 0.2) +
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
  