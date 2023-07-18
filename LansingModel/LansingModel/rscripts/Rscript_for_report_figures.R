###############################################################################

# Willemijn Oudijk - s4995805
# Masters project
# 18/07/23
# Parental age and offspring lifespan: the Lansing effect and its underlying mechanisms
# R script to generate the figures in the report 

###############################################################################
# install if necessary, then load
library(cowplot)
library(tidyverse)
library(mgcv)
library(MetBrewer)
library(ggpubr)
library(grid) # check if necessary 
library(gridExtra)
#library(R.filesets) # remove 

# set path with downloaded data
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/dataForReport/"
# set path to where the figures should be saved
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/figuresForReport/"

###############################################################################
# FIGURE 1
###############################################################################

funct <- function(x) {1/(1 + exp(a * (x - b)))}
a <- 10 # FIGURE 1A
b <- 0.5 # FIGURE 1A
a <- 2.7 # FIGURE 1B
b <- 1 # FIGURE 1B
df <- as.data.frame(funct(seq(0,1,0.026)))
df$age <- 1:39
colnames(df)[1] <- "survivalProb"
ggplot(df, aes(age, survivalProb)) + 
  geom_point() +
  labs(#y = "Resource allocation to repair", # FIGURE 1B
    y = "Survival probability, parental care quality", # FIGURE 1A
    x = "Age") +
  theme_bw() +
  theme(text = element_text(size = 13),
        strip.text.x = element_text(size = 13, face = "bold")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     limits = c(0,1))

# save figure 
ggsave(paste0(output_path, "Figure1.pdf"), width = 12, height = 8)

###############################################################################
# FIGURE 3
###############################################################################

# contains the predicted values based on the gam models for the longitudinal data
predDataTotalNormalized <- read.table(paste0(path, "predDataTotalNormalized.txt"))

# contains the predicted values based on the gam models for the latitudinal data
predDataTotalNormalizedLat <- read.table(paste0(path, "predDataTotalNormalizedLat.txt"))

# add column with group to both normalized data sets 
predDataTotalNormalized$group <- "Longitudinal"
predDataTotalNormalizedLat$group <- "Cross-sectional"

# merge them together 
totalNormalizedData <- rbind(predDataTotalNormalized, predDataTotalNormalizedLat)

# GENERATE PLOTS FOR LOWER TRIANGULAR
scenarios <- c() # get the scenarios that should be included in the matrix
for (x in levels(totalNormalizedData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) <= 2){
    # get list of single scenarios and the set combinations
    scenarios <- c(scenarios, x) 
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
    ylim(0, 2) + 
    theme_minimal() +
    theme(legend.position = "none",
      axis.text = element_text(size=13,face="plain",color="black"),
      axis.title = element_text(size = 13),
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

# GENERATE PLOTS OF PARENTAL AGE DISTRIBUTION (UPPER TRIANGULAR)
predDataTotalAgeDist <- read.table(paste0(path, "predDataTotalAgeDist.txt"))

# get the scenarios relevant for the matrix plot
scenarios <- c()
for (x in levels(predDataTotalAgeDist$scenario)) {
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
    theme_minimal() +
    theme(legend.position = "none",
      axis.text = element_text(size=13,face="plain",color="black"), # size = 11
      axis.title = element_text(size = 13),
      axis.line = element_line(color="black", linewidth = 1.0),
      panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7)) +
  scale_x_continuous(breaks = seq(0,40,10), limits = c(0,40))
  
  # save plot in list 
  plotsAgeDist[[i]] <- p
  # adjust name to corresponding scenario run
  names(plotsAgeDist)[i] <- paste("p_", scenarios[i], sep = "") 
}

# remove some of the axes numbering to make the matrix better readable
plotsTot$p_null <- plotsTot$p_null + theme(axis.text.x = element_blank())
plotsTot$p_nullDamage <- plotsTot$p_nullDamage + theme(axis.text.x = element_blank())
plotsTot$p_nullQuality <- plotsTot$p_nullQuality + theme(axis.text.x = element_blank())
plotsTot$p_damage <- plotsTot$p_damage + theme(axis.text.x = element_blank(),
                                               axis.text.y = element_blank())
plotsTot$p_damageQuality <- plotsTot$p_damageQuality + theme(axis.text.x = element_blank(),
                                                             axis.text.y = element_blank())
plotsTot$p_damageResource <- plotsTot$p_damageResource + theme(axis.text.y = element_blank())
plotsTot$p_quality <- plotsTot$p_quality + theme(axis.text.x = element_blank(),
                                                 axis.text.y = element_blank())
plotsTot$p_qualityResource <- plotsTot$p_qualityResource + theme(axis.text.y = element_blank())
plotsTot$p_resource <- plotsTot$p_resource + theme(axis.text.y = element_blank())
plotsAgeDist$p_nullQuality <- plotsAgeDist$p_nullQuality + theme(axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank())
plotsAgeDist$p_nullResource <- plotsAgeDist$p_nullResource + theme(axis.text.x = element_blank(),
                                                                   axis.text.y = element_blank())
plotsAgeDist$p_damageResource <- plotsAgeDist$p_damageResource + theme(axis.text.x = element_blank(),
                                                                       axis.text.y = element_blank())

# plot the matrix
plot_matrix <- plot_grid(plotsTot$p_null, plotsAgeDist$p_nullDamage, plotsAgeDist$p_nullQuality, plotsAgeDist$p_nullResource,
                         plotsTot$p_nullDamage, plotsTot$p_damage, plotsAgeDist$p_damageQuality, plotsAgeDist$p_damageResource,
                         plotsTot$p_nullQuality, plotsTot$p_damageQuality, plotsTot$p_quality, plotsAgeDist$p_qualityResource, 
                         plotsTot$p_nullResource, plotsTot$p_damageResource, plotsTot$p_qualityResource, plotsTot$p_resource,
                         ncol = 4, align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.88, label_y = 0.97)
plot_matrix

# add legend
p_tmp <- plotsTot$p_resource + theme(legend.position = "bottom",
                                     legend.title = element_blank(),
                                     legend.text = element_text(size = 13))
leg1 <- get_legend(p_tmp) 
legend <- as_ggplot(leg1)
legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)

# arrange grid again, with legend and axis titles
plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                         left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))

plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
# visualize the result
plot_with_legend

# save figure 
ggsave(paste0(output_path, "Figure3.pdf"), width = 12, height = 10)

###############################################################################
# FIGURE S1 - S4
###############################################################################

found <- read.table(paste0(path, "baseline_ParameterSim.txt")) # FIGURE S1
found <- read.table(paste0(path, "mechanism1_ParameterSim.txt")) # FIGURE S2
found <- read.table(paste0(path, "mechanism2_ParameterSim.txt")) # FIGURE S3
found <- read.table(paste0(path, "mechanism3_ParameterSim.txt")) # FIGURE S4

# average the offspring lifespan per maternal age per mutation probability 
calcMargin = function(x) { # function to calculate the margin for confidence interval 
  n <- x[4] # get sample size
  n <- as.numeric(n)
  sd <- x[5] # get sd 
  sd <- as.numeric(sd)
  return(qt(0.975, df = n - 1)*sd/sqrt(n))
}

found <- found %>% mutate(mutProb = factor(mutProb))
allAvgs <- c()
for (x in levels(found$mutProb)) {
  sub <- found[found$mutProb == x,]
  percentile <- quantile(sub$maternalAge, probs = 0.95)
  sub <- sub[sub$maternalAge <= percentile,]
  # calculate average expected age at death from parental age of 0 until 95th percentile parental age 
  avg <- aggregate(sub$ageAtDeath, list(sub$maternalAge), mean)
  colnames(avg) <- c("ageOfParent", "avgOffspringLifespan")
  avg$mutProb <- x
  
  # get information for CIs. 
  avg$n <- tapply(sub$ageAtDeath, sub$maternalAge, length)
  avg$sd <- tapply(sub$ageAtDeath, sub$maternalAge, sd)
  
  avg$margin <- apply(avg, 1, calcMargin) # calc margin per row
  avg$lwr <- by(avg, seq(nrow(avg)), function(x){x[2] - x[6]}) # get lower bound per row
  avg$upr <- by(avg, seq(nrow(avg)), function(x){x[2] + x[6]}) # get upper bound per row
  
  allAvgs <- rbind(allAvgs, avg)
}

allAvgs <- allAvgs %>% mutate(lwr = as.numeric(lwr),
                              upr = as.numeric(upr))

# with 95% CIs. 
# FIGURE S1 & FIGURE S3
legend_title <- "Mutation probability for 
age-specific survival genes" 

# FIGURE S2
legend_title <- "Rate of gamete decline" 

# FIGURE S4
legend_title <- "Mutation probability for age-specific 
resource allocation genes"

ggplot(allAvgs, aes(ageOfParent, avgOffspringLifespan, group = mutProb, colour = mutProb)) +
  geom_line() +
  geom_ribbon(data = allAvgs,
              aes(ymin = lwr, ymax = upr, fill = mutProb),
              alpha = 0.2, colour = NA) +
  theme_minimal() +
  theme(axis.text = element_text(size=15,face="plain",color="black"),
        axis.title = element_text(size = 16),
        axis.line = element_line(color="black", linewidth = 1.0),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.7),
        text = element_text(size = 13),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16)) +
  scale_colour_manual(legend_title, values = met.brewer("Egypt")[1:5]) +
  scale_fill_manual(legend_title, values = met.brewer("Egypt")[1:5]) + 
  ylim(0, 30) +
  labs(x = "Maternal age",
       y = "Average offspring lifespan") 

# save figure
ggsave(paste0(output_path, "figure_S1.pdf"), width = 12, height = 8)
ggsave(paste0(output_path, "figure_S2.pdf"), width = 12, height = 8)
ggsave(paste0(output_path, "figure_S3.pdf"), width = 12, height = 8)
ggsave(paste0(output_path, "figure_S4.pdf"), width = 12, height = 8)

###############################################################################
# FIGURE S5 & S6
###############################################################################

found <- read.table(paste0(path, "resourceGeneVals.txt")) # FIGURE S5
found <- read.table(paste0(path, "resourceGeneVals.txt")) # FIGURE S6

# get the mean gene value per age class
tmp <- aggregate(found$InvestmentGeneVal, list(found$age), mean)
colnames(tmp) <- c("age", "mean")
tmp$min <- tapply(found$InvestmentGeneVal, found$age, min)
tmp$max <- tapply(found$InvestmentGeneVal, found$age, max)

ggplot(tmp, aes(age, (1-mean))) + geom_line() +
  geom_ribbon(aes(ymin =(1-min), ymax = (1-max)), alpha = 0.2) +
  theme_bw() +
  labs(x = "Age",
       y = "Averaged proportion of resources 
  allocated to reproduction") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  theme(axis.text = element_text(size=17,face="plain",color="black"), 
        axis.title = element_text(size = 17))

# save figure
ggsave(paste0(output_path, "Figure_S5.pdf"), width = 12, height = 8)
ggsave(paste0(output_path, "Figure_S6.pdf"), width = 12, height = 8)

###############################################################################
# FIGURE S7
###############################################################################

scenarios <- c()
for (x in levels(totalNormalizedData$scenario)) {
  if (length(strsplit(x, "(?<=.)(?=[A-Z])", perl = TRUE)[[1]]) == 3){
    scenarios <- c(scenarios, x) # get list of three scenarios combined
  }
}

# dynamically generate the plots
plotsThreeCombs <- c()
for (i in 1:length(scenarios)){
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
      axis.title = element_text(size = 13),
      axis.line = element_line(color="black", linewidth = 0.6),
      panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.5)) + 
    scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1)) + 
    scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
    scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
    scale_y_continuous(breaks = seq(0, 1.5, 0.2)) +
    coord_cartesian(ylim= c(0, 1.5))
  
  # save plot in list 
  plotsThreeCombs[[i]] <- p
  # adjust name to corresponding scenario run
  names(plotsThreeCombs)[i] <- paste("p_", scenarios[i], sep = "") 
}

plotsThreeCombs$p_nullDamageQuality <- plotsThreeCombs$p_nullDamageQuality + theme(axis.text.x = element_blank())
plotsThreeCombs$p_nullDamageResource <- plotsThreeCombs$p_nullDamageResource + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
plotsThreeCombs$p_nullQualityResource <- plotsThreeCombs$p_nullQualityResource + theme(axis.text.y = element_blank())

plot_matrix <- plot_grid(plotsThreeCombs$p_nullDamageQuality,
                         plotsThreeCombs$p_nullDamageResource,
                         plotsThreeCombs$p_damageQualityResource,
                         plotsThreeCombs$p_nullQualityResource,
                         align = "hv", axis = "brlt", labels = "AUTO", label_x = 0.90, label_y = 0.97)

p_tmp <- plotsThreeCombs$p_nullQualityResource + theme(legend.position = "bottom",
                                                       legend.title = element_blank(),
                                                       legend.text = element_text(size=13))

leg1 <- get_legend(p_tmp) 
legend <- as_ggplot(leg1)
legend <- legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
blank <- ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 1, 0), "cm"))
legend_row <- plot_grid(legend, blank, blank, blank, nrows = 1)

plot_matrix2 <- grid.arrange(arrangeGrob(plot_matrix, bottom = textGrob("Normalized parental age", gp=gpar(fontsize=15)), 
                                         left = textGrob("Normalized offspring lifespan", gp=gpar(fontsize=15), rot = 90)))

plot_with_legend <- plot_grid(plot_matrix2, legend_row, ncol = 1, rel_heights = c(1, 0.01))
plot_with_legend

# save figure
ggsave(paste0(output_path, "Figure_S7.pdf"), width = 12, height = 10)

###############################################################################
# FIGURE S8
###############################################################################

scenario <- "nullDamageQualityResource"
ggplot(totalNormalizedData[totalNormalizedData$scenario == scenario,], 
            aes(maternalAge, mean, group = group, colour = group)) +
  geom_line() +
  labs(x = "Normalized parental age",
       y = "Normalized offspring lifespan") +
  geom_ribbon(data = totalNormalizedData[totalNormalizedData$scenario == scenario,],
              aes(ymin = min, ymax = max,  fill = group), 
              alpha = 0.2, colour = NA) + 
  theme_minimal() +
  theme(legend.text = element_text(size=13),
        legend.position = c(0.10, -0.04),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text = element_text(size=13,face="plain",color="black"),
        axis.title = element_text(size = 13),
        axis.line = element_line(color="black", linewidth = 0.6),
        panel.border = element_rect(colour = "darkgray", fill=NA, linewidth=0.5)) + 
  scale_x_continuous(breaks = seq(0,40,0.2), limits = c(0,1)) + 
  scale_colour_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
  scale_fill_manual(values = met.brewer("Egypt", 4)[c(2,4)]) +
  scale_y_continuous(breaks = seq(0, 1.3, 0.2)) +
  coord_cartesian(ylim = c(0,1.3))

# save figure
ggsave(paste0(output_path, "Figure_S8.pdf"), width = 12, height = 8)