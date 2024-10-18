# Load necessary libraries
library(diceplot)
library(dplyr)
library(ggplot2)
library(tidyr)

# Define genes
gene_list <- c("GeneA", "GeneB", "GeneC")

# Define cell types
cell_types <- c("Neuron", "Astrocyte", "Microglia")

# Define Contrasts
contrasts <- c("Type1", "Type2")  # Changed for demonstration

# Define vars for each Contrast
vars_type1 <- c("MCI-NCI", "AD-MCI", "AD-NCI")
vars_type2 <- c("Amyloid", "Plaq N", "Tangles", "NFT")

# Create a data frame with all combinations
data <- expand.grid(
  gene = gene_list,
  Cell_Type = cell_types,  # Renamed column
  Group = contrasts,       # Renamed column
  stringsAsFactors = FALSE
)

# Add the appropriate vars to each Contrast
set.seed(123) # Ensure reproducibility
data_type1 <- data %>% 
  filter(Group == "Type1") %>% 
  mutate(var = sample(vars_type1, n(), replace = TRUE))

data_type2 <- data %>% 
  filter(Group == "Type2") %>% 
  mutate(var = sample(vars_type2, n(), replace = TRUE))

# Combine the data
data <- bind_rows(data_type1, data_type2)

# Assign random values for logFC and adjusted p-values
data <- data %>%
  mutate(
    logFC = runif(n(), min = -2, max = 2),  # Renamed column
    adj_p_value = runif(n(), min = 0.0001, max = 0.05)
  )

# Use the modified function with custom parameters
p <- domino_plot(
  data = data,
  gene_list = gene_list,
  feature_col = "gene",
  celltype_col = "Cell_Type",
  contrast_col = "Group",
  contrast_levels = c("Type1", "Type2"),
  contrast_labels = c("Type 1", "Type 2"),
  logfc_col = "logFC",
  pval_col = "adj_p_value",
  switch_axis = FALSE,
  min_dot_size = 1,
  max_dot_size = 5,
  output_file = "domino_plot_example.png"
)

# Display the plot
print(p)