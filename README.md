[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/diceplot)](https://CRAN.R-project.org/package=diceplot)

> [!Note]
> This repository is in active development
> Fudging around with the paramters for a proper legend might be necessary

# DicePlot

The **DicePlot** package allows you to create visualizations (dice plots) for datasets with more than two categorical variables and additional continuous variables. This tool is particularly useful for exploring complex categorical data and their relationships with continuous variables.

## Installation

To install the **DicePlot** package, follow these steps:

### 1. Install R

Ensure that you have R installed on your system. You can download it from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

### 2. Install Required Packages

The `DicePlot` package depends on several other R packages. Install them by running:

```r
install.packages(c(
    "devtools",
    "dplyr",
    "ggplot2",
    "tidyr",
    "data.table",
    "ggdendro"
))
```

### 3.1 Install DicePlot from GitHub

You can install the `DicePlot` package directly from GitHub using the `devtools` package

```r
# Install devtools if you haven't already
install.packages("devtools")
# Install DicePlot from GitHub
devtools::install_github("maflot/DicePlot/diceplot")
```
### 3.2 Install DicePlot from files
Download the repository and run following code to install the package
```r
install.packages("$path on your local machine$/DicePlot/diceplot",repos = NULL, type="source")
```

### 4. Load the Package

After installation, load the `DicePlot` package into your R session:

```r
library(DicePlot)
```

## Example Usage

Here is a simple example of how to use the `DicePlot` package:

```r
# Load necessary libraries
library(diceplot)
library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
library(cowplot)
library(ggplotify)

plot_path = "sample_plots"

# Define the variables and their colors for 3 variables
pathology_variables <- c("Amyloid", "NFT", "Tangles")
cat_c_colors <- c(
  "Amyloid" = "#d5cccd",
  "NFT" = "#cb9992",
  "Tangles" = "#ad310f"
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
    "BBB-linked", "Cell-proliferation", "Other", "BBB-linked", "Cell-proliferation",
    "Cell-proliferation", "Other", "Other", "Other", "BBB-linked"
  ),
  stringsAsFactors = FALSE
)

# Define group colors
group_colors <- c(
  "BBB-linked" = "#333333",
  "Cell-proliferation" = "#888888",
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
          switch_axis = FALSE,
          group_alpha = 0.6,
          title = "Dice Plot with 3 Pathology Variables",
          cat_c_colors = cat_c_colors, 
          group_colors = group_colors, 
          format = ".png",
          custom_theme = theme_minimal())
```

This code will generate a dice plot visualizing the relationships between the categorical variables `Category1`, `Category2`, `Category3`, and the continuous variable `Value`.

## Use the dice plots in a real programming language
for using dice plots in python please refer to [pyDicePlot](https://github.com/maflot/pyDicePlot/tree/main)

## Sample Output

![Sample Dice with 3 categories Plot](./sample_plots/dice_plot_3_example_dice_plot.png)
![Sample Dice with 6 categories Plot](./sample_plots/dice_plot_5_example_dice_plot.png)

*Figure: A sample dice plot generated using the `DicePlot` package.*

## Documentation

For full documentation and additional examples, please refer to the [documentation](https://dice-and-domino-plot.readthedocs.io/en/latest/index.html#)

- **Visualize Complex Data:** Easily create plots for datasets with multiple categorical variables.
- **Customization:** Customize plots with titles, labels, and themes.
- **Integration with ggplot2:** Leverages the power of `ggplot2` for advanced plotting capabilities.

## Contributing

We welcome contributions from the community! If you'd like to contribute:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with a detailed description of your changes.

## Contact

If you have any questions, suggestions, or issues, please open an issue on GitHub.
