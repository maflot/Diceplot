# Load necessary libraries
library(diceplot)
library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
library(cowplot)
library(RColorBrewer)

# Define common variables
cell_types <- c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial")
pathways <- c(
  "Apoptosis", "Inflammation", "Metabolism", "Signal Transduction", "Synaptic Transmission",
  "Cell Cycle", "DNA Repair", "Protein Synthesis", "Lipid Metabolism", "Neurotransmitter Release",
  "Oxidative Stress", "Energy Production", "Calcium Signaling", "Synaptic Plasticity", "Immune Response"
)

# Assign groups to pathways
pathway_groups <- data.frame(
  Pathway = pathways,
  Group = c(
    "Linked", "UnLinked", "Other", "Linked", "UnLinked",
    "UnLinked", "Other", "Other", "Other", "Linked",
    "Other", "Other", "Linked", "UnLinked", "Other"
  ),
  stringsAsFactors = FALSE
)

pathology_variables <- c("Alzheimer's disease", "Cancer", "Flu", "ADHD", "Age", "Weight")

n_colors <- length(pathology_variables)
colors <- brewer.pal(n = n_colors, name = "Set1")
cat_c_colors <- setNames(colors, pathology_variables)


# Function to create and plot dice plots
create_and_plot_dice <- function(pathology_variables, cat_c_colors, title, cell_types, pathways, pathway_groups, min_dot_size=3, max_dot_size=6) {
  # Create dummy data
  set.seed(123)
  data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)

  data <- data %>%
    rowwise() %>%
    mutate(
      PathologyVariable = list(sample(pathology_variables, size = sample(1:length(pathology_variables), 1)))
    ) %>%
    unnest(cols = c(PathologyVariable))
  
  # Merge the group assignments into the data
  data <- data %>%
    left_join(pathway_groups, by = "Pathway")
  
  # Use the dice_plot function
  dice_plot(data = data, 
            cat_a = "CellType", 
            cat_b = "Pathway", 
            cat_c = "PathologyVariable", 
            group = "Group",
            group_alpha = 0.6,
            title = title,
            cat_c_colors = cat_c_colors, 
            custom_theme = theme_minimal(),
            min_dot_size = min_dot_size,
            max_dot_size = max_dot_size
            )
}


# First plot with 3 pathology variables
pathology_variables_3 <- c("Stroke", "Cancer", "Flu")
create_and_plot_dice(
  pathology_variables = pathology_variables_3,
  cat_c_colors = cat_c_colors[pathology_variables_3],
  title = "Dice Plot with 3 Pathology Variables",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Second plot with 4 pathology variables
pathology_variables_4 <- c("Stroke", "Cancer", "Flu", "ADHD")
create_and_plot_dice(
  pathology_variables = pathology_variables_4,
  cat_c_colors = cat_c_colors[pathology_variables_4],
  title = "Dice Plot with 4 Pathology Variables",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Third plot with 5 pathology variables
pathology_variables_5 <- c("Stroke", "Cancer", "Flu", "ADHD", "Lymphom")
create_and_plot_dice(
  pathology_variables = pathology_variables_5,
  cat_c_colors = cat_c_colors[pathology_variables_5],
  title = "Dice Plot with 5 Pathology Variables",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Fourth plot with 6 pathology variables
pathology_variables_6 <- c("Alzheimer's disease", "Cancer", "Flu", "ADHD", "Age", "Weight")
create_and_plot_dice(
  pathology_variables = pathology_variables_6,
  cat_c_colors = cat_c_colors[pathology_variables_6],
  title = "Dice Plot with 6 Pathology Variables",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Example of a large dice plot with adjusted dot sizes
# Define more cell types and pathways for the large plot
cell_types_large <- c(
  "Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "Endothelial",
  "Ependymal", "Pericyte", "Satellite Cell", "Schwann Cell", "Macrophage"
)

pathways_large <- c(
  "Apoptosis", "Inflammation", "Metabolism", "Signal Transduction", "Synaptic Transmission",
  "Cell Cycle", "DNA Repair", "Protein Synthesis", "Lipid Metabolism", "Neurotransmitter Release",
  "Oxidative Stress", "Energy Production", "Calcium Signaling", "Synaptic Plasticity", "Immune Response",
  "Gene Expression", "Membrane Transport", "Cell Migration", "Cell Adhesion", "Cell Differentiation",
  "Angiogenesis", "Neurogenesis", "Protein Folding", "Autophagy", "Endocytosis"
)

# Assign groups to the larger set of pathways
set.seed(456)
pathway_groups_large <- data.frame(
  Pathway = pathways_large,
  Group = sample(c("Linked", "UnLinked", "Other"), size = length(pathways_large), replace = TRUE),
  stringsAsFactors = FALSE
)

# Large dice plot with adjusted dot sizes
create_and_plot_dice(
  pathology_variables = pathology_variables_6,
  cat_c_colors = cat_c_colors[pathology_variables_6],
  title = "Large Dice Plot with Adjusted Dot Sizes",
  cell_types = cell_types_large,
  pathways = pathways_large,
  pathway_groups = pathway_groups_large,
  min_dot_size = 1,
  max_dot_size = 3
)

