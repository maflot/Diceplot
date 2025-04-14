# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)
library(RColorBrewer)
library(UpSetR)
library(ggplot2)
library(diceplot)

# Set your file path
file_path <- "pnas.2023216118.sd05.xlsx"

# Function to create the properly formatted CSV
process_excel_to_csv <- function(file_path) {
  # Read Excel file with detailed options to ensure proper data reading
  raw_data <- read_excel(file_path, col_names = FALSE, na = "", trim_ws = TRUE)
  
  # Extract cell types from row 2
  cell_types_row <- raw_data[2,]
  
  # Extract demographic info from row 3
  demo_row <- raw_data[3,]
  
  # Create a list to store all transformed data
  all_data <- list()
  
  # Define cell type mapping
  cell_type_map <- c(
    "NK" = "Natural Killer (NK) cell",
    "TC" = "T cell (TC)",
    "BC" = "B cell (BC)",
    "DC" = "Dendritic cell (DC)",
    "MC" = "Monocyte (MC)"
  )
  
  # Special handling for the staggered format
  # First find the cell type columns
  cell_type_columns <- c()
  for (i in 1:ncol(raw_data)) {
    if (!is.na(cell_types_row[[i]]) && cell_types_row[[i]] != "") {
      cell_type_columns <- c(cell_type_columns, i)
    }
  }
  
  # Print debug information
  print(paste("Found cell type columns:", paste(cell_type_columns, collapse = ", ")))
  
  # Process each cell type column and its associated demographic columns
  for (col_idx in cell_type_columns) {
    cell_type <- cell_types_row[[col_idx]]
    cell_type_full <- cell_type_map[cell_type]
    
    # Look at the next 4 columns (OM, OF, YM, YF)
    for (offset in 0:3) {
      demo_col <- col_idx + offset
      
      # Check if this column exists and has a valid demographic
      if (demo_col <= ncol(raw_data) && !is.na(demo_row[[demo_col]]) && demo_row[[demo_col]] != "") {
        demo_info <- demo_row[[demo_col]]
        
        # Print debug info
        print(paste("Processing column", demo_col, "- Cell type:", cell_type, "- Demo info:", demo_info))
        
        # Extract demographic information
        age <- case_when(
          substr(demo_info, 4, 4) == "O" ~ "old",
          substr(demo_info, 4, 4) == "Y" ~ "young",
          TRUE ~ NA_character_
        )
        
        sex <- case_when(
          substr(demo_info, 5, 5) == "M" ~ "male",
          substr(demo_info, 5, 5) == "F" ~ "female",
          TRUE ~ NA_character_
        )
        
        # Print the parsed demographic info for debugging
        print(paste("  Demographic parsed as:", age, sex))
        
        # Process each gene in this column
        gene_count <- 0
        for (row_idx in 4:nrow(raw_data)) {
          gene <- raw_data[row_idx, demo_col][[1]]
          
          # Skip empty genes
          if (is.na(gene) || gene == "") {
            next
          }
          
          gene_count <- gene_count + 1
          
          # Create a row for this gene
          gene_row <- data.frame(
            id = paste0(cell_type, "_", demo_info, "_", gene),
            gene = gene,
            cell_type_code = cell_type,
            cell_type = cell_type_full,
            age_code = substr(demo_info, 4, 4),
            age = age,
            sex_code = substr(demo_info, 5, 5),
            sex = sex,
            demo_code = demo_info
          )
          
          # Add to our list
          all_data[[length(all_data) + 1]] <- gene_row
        }
        
        print(paste("  Processed", gene_count, "genes in this column"))
      }
    }
  }
  
  # Combine all data frames
  if (length(all_data) == 0) {
    stop("No data was processed. Check the Excel file structure.")
  }
  
  final_data <- bind_rows(all_data)
  
  # Return the final data frame
  return(final_data)
}

# Process the data
processed_data <- process_excel_to_csv(file_path)

# Save as CSV
write.csv(processed_data, "processed_gene_data.csv", row.names = FALSE, quote = TRUE)

# Display first few rows
head(processed_data)

# Summary of the data
summary_stats <- processed_data %>%
  group_by(cell_type, age, sex) %>%
  summarise(gene_count = n(), .groups = 'drop')

print("Summary statistics by cell type, age, and sex:")
print(summary_stats)

# Verify all demographic combinations are present
demo_validation <- processed_data %>%
  group_by(age, sex) %>%
  summarise(tmp_count = n(), .groups = 'drop')

print("Demographic validation:")
print(demo_validation)

# Create a nicely formatted report of gene counts by cell type and demographic
gene_counts_report <- processed_data %>%
  group_by(cell_type, age, sex) %>%
  summarise(
    gene_count = n(),
    top_genes = paste(head(gene, 5), collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(cell_type, age, sex)

# Print the report
print("Gene count report by cell type and demographics:")
print(gene_counts_report)

# Now prepare the data for diceplot visualization
# Create a demographic combination column that combines age and sex
processed_data <- processed_data %>%
  mutate(demo_combination = case_when(
    age == "old" & sex == "male" ~ "Old Male",
    age == "old" & sex == "female" ~ "Old Female",
    age == "young" & sex == "male" ~ "Young Male",
    age == "young" & sex == "female" ~ "Young Female",
    TRUE ~ paste(age, sex) # Fallback for any unexpected combinations
  ))

# Print unique combinations to verify
print("Unique demographic combinations:")
print(unique(processed_data$demo_combination))

# Order the demographic combinations factor according to the requested order
processed_data$demo_combination <- factor(
  processed_data$demo_combination,
  levels = c("Old Male", "Old Female", "Young Male", "Young Female")
)

# Make sure cell types are also properly ordered
processed_data$cell_type <- factor(
  processed_data$cell_type,
  levels = c(
    "Natural Killer (NK) cell",
    "T cell (TC)",
    "B cell (BC)",
    "Dendritic cell (DC)",
    "Monocyte (MC)"
  )
)

# For visualization, create a summary table with gene counts by cell type and demographic combination
gene_counts <- processed_data %>%
  group_by(gene, cell_type, demo_combination) %>%
  summarize(tmp_count = n(), .groups = "drop")

# Check the structure of gene_counts 
print("Structure of gene_counts data frame:")
str(gene_counts)

# Check if we have genes in all demographic combinations
print("Genes by cell type and demographic combination:")
print(table(gene_counts$cell_type, gene_counts$demo_combination))

# Define colors for the demographic combinations
demo_colors <- c(
  "Old Male" = "#E41A1C",     # Red
  "Old Female" = "#377EB8",   # Blue
  "Young Male" = "#4DAF4A",   # Green
  "Young Female" = "#984EA3"  # Purple
)

# Get the top 25 most frequent genes instead of all genes with >=10 occurrences
top_25_genes <- processed_data %>%
  count(gene) %>%
  arrange(desc(n)) %>%
  head(25) %>%
  pull(gene)

# Print the top 25 genes
print("Top 25 most frequent genes:")
print(top_25_genes)

# Filter the gene_counts data frame to include only top 25 genes
filtered_gene_counts <- gene_counts %>%
  filter(gene %in% top_25_genes)

# Print the number of rows before and after filtering
print(paste("Original gene_counts rows:", nrow(gene_counts)))
print(paste("Filtered gene_counts rows (top 25 genes):", nrow(filtered_gene_counts)))

# Get the actual combinations present in the filtered data
actual_combinations <- unique(as.character(filtered_gene_counts$demo_combination))
actual_colors <- demo_colors[actual_combinations]

# Add the default group column
filtered_gene_counts$default = ""

# Create the diceplot using the proper parameters with filtered data (top 25 genes)
p_dice_filtered <- dice_plot(
  data = filtered_gene_counts,
  x = "gene",                    # x-axis: genes
  y = "cell_type",               # y-axis: cell types
  z = "demo_combination",        # z parameter: demographic combinations
  cluster_by_column = T,
  cluster_by_row = F,
  title = "Gene Expression across Cell Types and Demographics\n(Top 25 Genes)",
  z_colors = actual_colors,      # Use the proper color palette
  max_dot_size = 6,
  min_dot_size = 3,
  legend_width = 0.2,
  legend_height = 0.25,
  show_legend = T
)

# Display the filtered diceplot
print(p_dice_filtered)

# Save the filtered diceplot
pdf("top_25_genes_diceplot.pdf", width = 12, height = 10)
print(p_dice_filtered)
dev.off()

# Save the filtered diceplot
pdf("top_25_genes_diceplot.pdf", width = 12, height = 10)
print(p_dice_filtered)
dev.off()

# We can also create a table to show which genes appear most frequently
top_genes_table <- processed_data %>%
  count(gene, sort = TRUE) %>%
  head(25) %>%
  mutate(
    cell_types = sapply(gene, function(g) {
      paste(sort(unique(processed_data$cell_type_code[processed_data$gene == g])), collapse = ", ")
    }),
    demographics = sapply(gene, function(g) {
      paste(sort(unique(processed_data$demo_code[processed_data$gene == g])), collapse = ", ")
    })
  )

# Print the table of top genes
print("Top 25 genes:")
print(top_genes_table)

# You can also create a CSV file with the top genes information
write.csv(top_genes_table, "top_genes_table.csv", row.names = FALSE)

# Create UpSet plots for different aspects of the data

# 1. UpSet plot for cell types - Which genes are expressed in which cell types
# First, create a presence/absence matrix - FIXED APPROACH
cell_type_lists <- list()
for (ct in c("NK", "TC", "BC", "DC", "MC")) {
  cell_type_lists[[ct]] <- unique(
    processed_data$gene[processed_data$cell_type_code == ct & processed_data$gene %in% top_25_genes]
  )
}

# Create the cell type upset plot for top 25 genes using fromList
cell_type_upset_plot <- upset(
  fromList(cell_type_lists),
  order.by = "freq",
  keep.order = TRUE,
  sets = c("NK", "TC", "BC", "DC", "MC"),
  sets.bar.color = brewer.pal(5, "Set2"),
  main.bar.color = "black",
  matrix.color = "darkblue",
  sets.x.label = "Number of Genes",
  mainbar.y.label = "Intersection Size",
  text.scale = 1.2,
  point.size = 3,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

# Display the cell type upset plot
print(cell_type_upset_plot)

# Save the cell type upset plot
pdf("top_25_genes_cell_type_upset_plot.pdf", width = 3, height = 3)
print(cell_type_upset_plot)
dev.off()

# 2. UpSet plot for demographic combinations for top 25 genes - FIXED APPROACH
demo_lists <- list()
for (demo in c("in OM", "in OF", "in YM", "in YF")) {
  demo_lists[[demo]] <- unique(
    processed_data$gene[processed_data$demo_code == demo & processed_data$gene %in% top_25_genes]
  )
}

# Create the demographic upset plot for top 25 genes
demo_upset_plot <- upset(
  fromList(demo_lists),
  order.by = "freq",
  keep.order = TRUE,
  sets = c("in OM", "in OF", "in YM", "in YF"),
  sets.bar.color = brewer.pal(4, "Set1"),
  main.bar.color = "black",
  matrix.color = "darkblue",
  sets.x.label = "Number of Genes",
  mainbar.y.label = "Intersection Size",
  text.scale = 1.2,
  point.size = 3,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

# Display the demographic upset plot
print(demo_upset_plot)

# Save the demographic upset plot
pdf("top_25_genes_demographic_upset_plot.pdf", width = 3, height = 3)
print(demo_upset_plot)
dev.off()

# Summary of the top 25 genes data
summary_by_cell_demo <- processed_data %>%
  filter(gene %in% top_25_genes) %>%
  group_by(cell_type, demo_combination) %>%
  summarise(gene_count = n(), .groups = 'drop') %>%
  pivot_wider(
    names_from = demo_combination,
    values_from = gene_count
  )

print("Summary of top 25 genes by cell type and demographic combination:")
print(summary_by_cell_demo)
