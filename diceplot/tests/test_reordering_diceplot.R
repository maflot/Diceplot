# Load necessary libraries
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

# Define all pathology variables used
all_pathology_variables <- c("Alzheimer's disease", "Cancer", "Flu", "ADHD", "Age", "Weight", "Stroke", "Lymphoma")

# Create color mapping for all pathology variables
n_colors <- length(all_pathology_variables)
colors <- brewer.pal(n = n_colors, name = "Set1")
cat_c_colors <- setNames(colors, all_pathology_variables)

# Original create_and_plot_dice function (unchanged)
create_and_plot_dice <- function(pathology_variables, cat_c_colors, title, cell_types, pathways, pathway_groups, min_dot_size=3, max_dot_size=6) {
  # Create dummy data
  set.seed(123)  # For reproducibility
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
  
  # Filter data to include only specified pathology variables
  data <- data %>% filter(PathologyVariable %in% pathology_variables)
  
  # Use the dice_plot function and store the plot
  p <- dice_plot(data = data, 
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
  
  return(p)
}

# Modified create_and_plot_dice function with varying counts
create_and_plot_dice_varying_counts <- function(pathology_variables, cat_c_colors, title, cell_types, 
                                                pathways, pathway_groups, min_dot_size=3, max_dot_size=6, 
                                                show_legend = TRUE) {
  # Create base data frame with all combinations
  data <- expand.grid(CellType = cell_types, Pathway = pathways, stringsAsFactors = FALSE)
  
  # Merge the group assignments into the data
  data <- data %>%
    left_join(pathway_groups, by = "Pathway")
  
  # Assign differing binomial parameters based on the group
  data <- data %>%
    mutate(
      # Assign probability based on the group
      prob = case_when(
        Group == "Linked" ~ 0.8,
        Group == "UnLinked" ~ 0.4,
        Group == "Other" ~ 0.2
      ),
      # Assign size parameter (number of trials) for binomial distribution
      size = 10
    )
  
  # Generate counts using binomial distribution
  set.seed(123)  # For reproducibility
  data <- data %>%
    rowwise() %>%
    mutate(
      Count = rbinom(1, size = size, prob = prob)
    ) %>%
    ungroup()
  
  # Remove rows with zero count to avoid empty data
  data <- data %>% filter(Count > 0)
  
  # Assign pathology variables randomly to each row
  data <- data %>%
    mutate(
      PathologyVariable = sample(pathology_variables, size = nrow(data), replace = TRUE)
    )
  
  # Expand the data according to the counts
  data_expanded <- data %>%
    uncount(weights = Count)
  
  # Use the dice_plot function and store the plot
  p <- dice_plot(data = data_expanded, 
            cat_a = "CellType", 
            cat_b = "Pathway", 
            cat_c = "PathologyVariable", 
            group = "Group",
            group_alpha = 0.6,
            title = title,
            cat_c_colors = cat_c_colors, 
            custom_theme = theme_minimal(),
            min_dot_size = min_dot_size,
            max_dot_size = max_dot_size,
            show_legend = show_legend
            )
  
  return(p)
}

# Define the pathology variables and colors
pathology_variables_4 <- c("Alzheimer's disease", "Cancer", "Flu", "ADHD")
cat_c_colors_4 <- cat_c_colors[pathology_variables_4]

# Generate and store the original plot
plot_original <- create_and_plot_dice(
  pathology_variables = pathology_variables_4,
  cat_c_colors = cat_c_colors_4,
  title = "Original Dice Plot",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Generate and store the plot with varying counts
plot_varying_counts <- create_and_plot_dice_varying_counts(
  pathology_variables = pathology_variables_4,
  cat_c_colors = cat_c_colors_4,
  title = "Dice Plot with Varying Counts",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups
)

# Display both plots side by side
plot_grid(plot_original, plot_varying_counts, ncol = 2)


no_legend <- create_and_plot_dice_varying_counts(
  pathology_variables = pathology_variables_4,
  cat_c_colors = cat_c_colors_4,
  title = "Dice Plot with Varying Counts",
  cell_types = cell_types,
  pathways = pathways,
  pathway_groups = pathway_groups,
  show_legend = F
)
no_legend
