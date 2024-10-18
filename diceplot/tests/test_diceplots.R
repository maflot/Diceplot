# Load necessary libraries
library(diceplot)
library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
library(cowplot)

plot_path = "./"

# Define the variables and their colors for 3 variables
pathology_variables <- c("Stroke", "Cancer", "Flu")
cat_c_colors <- c(
  "Stroke" = "#d5cccd",
  "Cancer" = "#cb9992",
  "Flu" = "#ad310f"
)

# Define cell types (cat_a)
cell_types <- c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial")

# Define pathways (cat_b) and groups
pathways <- c(
  "Apoptosis", "Inflammation", "Metabolism", "Signal Transduction", "Synaptic Transmission",
  "Cell Cycle", "DNA Repair", "Protein Synthesis", "Lipid Metabolism", "Neurotransmitter Release"
)

# Assign groups to pathways
pathway_groups <- data.frame(
  Pathway = pathways,
  Group = c(
    "Linked", "UnLinked", "Other", "Linked", "UnLinked",
    "UnLinked", "Other", "Other", "Other", "Linked"
  ),
  stringsAsFactors = FALSE
)

# Define group colors
group_colors <- c(
  "Linked" = "#333333",
  "UnLinked" = "#888888",
  "Other" = "#DDDDDD"
)

# Create dummy data
set.seed(123)
data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)

# Assign random pathology variables to each combination
data <- data %>%
  rowwise() %>%
  mutate(
    PathologyVariable = list(sample(pathology_variables, size = sample(1:3, 1)))
  ) %>%
  unnest(cols = c(PathologyVariable))

# Merge the group assignments into the data
data <- data %>%
  left_join(pathway_groups, by = c("Pathway" = "Pathway"))

# Use the dice_plot function
dice_plot(data = data, 
          cat_a = "CellType", 
          cat_b = "Pathway", 
          cat_c = "PathologyVariable", 
          group = "Group",
          plot_path = plot_path, 
          output_str = "dice_plot_3_example", 
          group_alpha = 0.6,
          title = "Dice Plot with 3 Pathology Variables",
          cat_c_colors = cat_c_colors, 
          group_colors = group_colors, 
          format = ".png",
          custom_theme = theme_minimal())



# Define the pathology variables and their colors
pathology_variables <- c("Stroke", "Cancer", "Flu", "ADHD")
cat_c_colors <- c(
  "Stroke" = "#d5cccd",
  "Cancer" = "#cb9992",
  "Flu" = "#ad310f",
  "ADHD" = "#7e2a20"
)

# Define cell types (cat_a)
cell_types <- c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial")

# Define pathways (cat_b) and add 10 more to make a total of 15
pathways <- c(
  "Apoptosis", "Inflammation", "Metabolism", "Signal Transduction", "Synaptic Transmission",
  "Cell Cycle", "DNA Repair", "Protein Synthesis", "Lipid Metabolism", "Neurotransmitter Release",
  "Oxidative Stress", "Energy Production", "Calcium Signaling", "Synaptic Plasticity", "Immune Response"
)

# Assign groups to pathways (ensuring each pathway has only one group)
pathway_groups <- data.frame(
  Pathway = pathways,
  Group = c(
    "Linked", "UnLinked", "Other", "Linked", "UnLinked",
    "UnLinked", "Other", "Other", "Other", "Linked",
    "Other", "Other", "Linked", "UnLinked", "Other"
  ),
  stringsAsFactors = FALSE
)

# Update group colors to shades of greys
group_colors <- c(
  "Linked" = "#333333",
  "UnLinked" = "#888888",
  "Other" = "#DDDDDD"
)

# Create dummy data
set.seed(123)
data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)

# Assign random pathology variables to each combination
data <- data %>%
  rowwise() %>%
  mutate(
    PathologyVariable = list(sample(pathology_variables, size = sample(1:4, 1)))
  ) %>%
  unnest(cols = c(PathologyVariable))

# Merge the group assignments into the data
data <- data %>%
  left_join(pathway_groups, by = "Pathway")

group_colors <- c(
  "Linked" = "#333333",
  "UnLinked" = "#888888",
  "Other" = "#DDDDDD"
)
# Use the modified dice_plot function
dice_plot(data = data, 
          cat_a = "CellType", 
          cat_b = "Pathway", 
          cat_c = "PathologyVariable", 
          group = "Group",
          plot_path = plot_path, 
          output_str = "dice_plot_4_example", 
          group_alpha = 0.6,
          title = "Dummy Dice Plot with Pathology Variables",
          cat_c_colors = cat_c_colors, 
          group_colors = group_colors, 
          format = ".png",
          custom_theme = theme_minimal())


# Define the variables and their colors for 5 variables
pathology_variables <- c("Stroke", "Cancer", "Flu", "ADHD", "Lymphom")
cat_c_colors <- c(
  "Stroke" = "#d5cccd",
  "Cancer" = "#cb9992",
  "Flu" = "#ad310f",
  "ADHD" = "#7e2a20",
  "Lymphom" = "#FFD700"  # Gold color for Lymphom
)

# Create dummy data
set.seed(123)
data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)

# Assign random pathology variables to each combination
data <- data %>%
  rowwise() %>%
  mutate(
    PathologyVariable = list(sample(pathology_variables, size = sample(1:5, 1)))
  ) %>%
  unnest(cols = c(PathologyVariable))

# Merge the group assignments into the data
data <- data %>%
  left_join(pathway_groups, by = c("Pathway" = "Pathway"))

# Use the dice_plot function
dice_plot(data = data, 
          cat_a = "CellType", 
          cat_b = "Pathway", 
          cat_c = "PathologyVariable", 
          group = "Group",
          plot_path = plot_path, 
          output_str = "dice_plot_5_example", 
          group_alpha = 0.6,
          title = "Dice Plot with 5 Pathology Variables",
          cat_c_colors = cat_c_colors, 
          group_colors = group_colors, 
          format = ".png",
          custom_theme = theme_minimal())


# Define the variables and their colors for 6 variables
pathology_variables <- c("Alzheimer's disease", "Cancer", "Flu", "ADHD", "Age", "Weight")
cat_c_colors <- c(
  "Stroke" = "#d5cccd",
  "Cancer" = "#cb9992",
  "Flu" = "#ad310f",
  "ADHD" = "#7e2a20",
  "Age" = "#FFD700",  # Gold color for Lymphom
  "Weight" = "#FF6622"  # Cyan color for Var6
)

# Create dummy data
set.seed(123)
data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)

# Assign random pathology variables to each combination
data <- data %>%
  rowwise() %>%
  mutate(
    PathologyVariable = list(sample(pathology_variables, size = sample(1:6, 1)))
  ) %>%
  unnest(cols = c(PathologyVariable))

# Merge the group assignments into the data
data <- data %>%
  left_join(pathway_groups, by = c("Pathway" = "Pathway"))

# Use the dice_plot function
dice_plot(data = data, 
          cat_a = "CellType", 
          cat_b = "Pathway", 
          cat_c = "PathologyVariable", 
          group = "Group",
          plot_path = plot_path, 
          output_str = "dice_plot_6_example", 
          group_alpha = 0.6,
          title = "Dice Plot with 6 Pathology Variables",
          cat_c_colors = cat_c_colors, 
          group_colors = group_colors, 
          format = ".png",
          custom_theme = theme_minimal())




