# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(diceplot)

# Define data
cell_types <- c("TypeA", "TypeB", "TypeC")
pathways <- c("Pathway1", "Pathway2", "Pathway3")
pathway_groups <- data.frame(
  Pathway = pathways,
  Group = c("Group1", "Group2", "Group3"),
  stringsAsFactors = FALSE
)

set.seed(123)
replicates <- 10
data <- expand.grid(
  CellType = cell_types,
  Pathway = pathways,
  Replicate = 1:replicates,
  stringsAsFactors = FALSE
)

data <- data %>%
  mutate(
    treat_time = sample(c("0", "2", "6", "18"), size = n(), replace = TRUE)
  )

missing_treat_times <- setdiff(c("0", "2", "6", "18"), unique(data$treat_time))
if(length(missing_treat_times) > 0){
  data$treat_time[1:length(missing_treat_times)] <- missing_treat_times
}

data <- data %>%
  mutate(
    treat_time = factor(treat_time, levels =c( "0","2","6","18"))
  ) %>%
  left_join(pathway_groups, by = "Pathway") %>%
  select(-Replicate)

# Create and plot the dice plot
p <- dice_plot(
  data = data, 
  cat_a = "treat_time", 
  cat_b = "Pathway", 
  cat_c = "CellType", 
  group = "Group",
  group_alpha = 0.6,
  title = "Dice Plot with 6 Pathology Variables and Treatment Time",
  custom_theme = theme_minimal(),
  min_dot_size = 2,
  max_dot_size = 4,
  cluster_by_column = F,
  cluster_by_row = FALSE,
)

print(p)
