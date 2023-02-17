###############################################################################
# Parameter exploration 
###############################################################################

# change path to your data folder 
path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/"
# change path to where you want the plots to be outputted 
output_path <- "/Users/willemijnoudijk/Documents/STUDY/Master Biology/ResearchProject1/data/plots/"

theme_big <- theme(axis.text.x = element_text(angle = 90),
                   axis.title = element_text(size = 30),
                   axis.text = element_text(size = 20),
                   title = element_text(size = 25),
                   legend.key.size = unit(2, 'cm'),
                   legend.text = element_text(size = 30))

########### varying mutation probability #########
parent_path <- paste(path, "mutationProbVaried/", sep = "")
####### smaller variation: 0.0015 - 0.015 
parent_path <- paste(path, "mutationProbVariedSmaller/", sep = "")

f <- list.files(path = parent_path, pattern = "outputDeathAge.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "mutationProb", "extrinsicMort", "deathAge") 

ggplot(data = found, aes(x = time, y = deathAge, col = mutationProb, group = mutationProb)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Average age of death over time with differing mutation probabilities", 
       x = "Time",
       y = "Average age of death") 
ggsave(paste(output_path, "plot_mut_prob_small.pdf", sep = ""), plot = last_plot()) #to save file 

############## varying extrinsic mortality ##############
parent_path <- paste(path, "extrinsicMortVaried/", sep = "")

f <- list.files(path = parent_path, pattern = "outputDeathAge.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "mutationProb", "extrinsicMort", "strengthOfSelection", "populationSize", "deathAge") 

ggplot(data = found, aes(x = time, y = deathAge, col = extrinsicMort, group = extrinsicMort)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Average age of death over time with differing extrinsic mortality probabilities", 
       x = "Time",
       y = "Average age of death") +
  xlim(0, 1000)
ggsave(paste(output_path, "plot_vary_extr_mort.pdf", sep = ""), plot = last_plot()) #to save file 

########## varying strength of selection #############
parent_path <- paste(path, "strength_of_selection/", sep = "")

f <- list.files(path = parent_path, pattern = "outputDeathAge.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "mutationProb", "extrinsicMort", "strengthOfSelection", "populationSize", "deathAge") 

ggplot(data = found, aes(x = time, y = deathAge, col = strengthOfSelection, group = strengthOfSelection)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Average age of death over time with differing strength of selection values", 
       x = "Time",
       y = "Average age of death") +
  theme_big
ggsave(paste(output_path, "plot_stregth_of_selection.pdf", sep = ""), plot = last_plot()) #to save file 

######### varying population size ##########
parent_path <- paste(path, "popSize/", sep = "")

f <- list.files(path = parent_path, pattern = "outputDeathAge.csv", recursive = T)
found <- c()
for (x in f) {
  file_name <- paste0(parent_path, x)
  local_data <- read.csv(file_name, header = F, sep = " ")
  found <- rbind(found, local_data)
}
colnames(found) <- c("time", "mutationProb", "extrinsicMort", "strengthOfSelection", "populationSize", "deathAge") 

ggplot(data = found, aes(x = time, y = deathAge, col = populationSize, group = populationSize)) + 
  geom_line() + scale_color_viridis_c(option = "A") + 
  labs(title = "Average age of death over time with differing population sizes", 
       x = "Time",
       y = "Average age of death") +
  theme_big
ggsave(paste(output_path, "plot_popsize.pdf", sep = ""), plot = last_plot()) #to save file 
