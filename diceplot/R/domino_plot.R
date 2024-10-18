utils::globalVariables(c("label_x", "label_y", "adj_logfc"))

#' Domino Plot Visualization
#'
#' This function generates a plot to visualize gene expression levels for a given list of genes. The size of the dots can be customized, and the plot can be saved to an output file if specified.
#'
#' @param data A data frame containing gene expression data.
#' @param gene_list A character vector of gene names to include in the plot.
#' @param switch_axis A logical value indicating whether to switch the x and y axes. Default is `FALSE`.
#' @param min_dot_size A numeric value indicating the minimum dot size in the plot. Default is `1`.
#' @param max_dot_size A numeric value indicating the maximum dot size in the plot. Default is `5`.
#' @param spacing_factor A numeric value indicating the spacing between gene pairs. Default is `3`.
#' @param var_id A string representing the column name in `data` for the variable identifier. Default is `"var"`.
#' @param feature_col A string representing the column name in `data` for the feature variable (e.g., genes). Default is `"gene"`.
#' @param celltype_col A string representing the column name in `data` for the cell type variable. Default is `"Celltype"`.
#' @param contrast_col A string representing the column name in `data` for the contrast variable. Default is `"Contrast"`.
#' @param contrast_levels A character vector specifying the levels of the contrast variable. Default is `c("Clinical", "Pathological")`.
#' @param contrast_labels A character vector specifying the labels for the contrasts in the plot. Default is `c("Clinical", "Pathological")`.
#' @param logfc_col A string representing the column name in `data` for the log fold change values. Default is `"avg_log2FC"`.
#' @param pval_col A string representing the column name in `data` for the adjusted p-values. Default is `"p_val_adj"`.
#' @param logfc_limits A numeric vector of length 2 specifying the limits for the log fold change color scale. Default is `c(-1.5, 1.5)`.
#' @param logfc_colors A named vector specifying the colors for the low, mid, and high values in the color scale. Default is `c(low = "blue", mid = "white", high = "red")`.
#' @param color_scale_name A string specifying the name of the color scale in the legend. Default is `"Log2 Fold Change"`.
#' @param size_scale_name A string specifying the name of the size scale in the legend. Default is `"-log10(adj. p-value)"`.
#' @param axis_text_size A numeric value specifying the size of the axis text. Default is `8`.
#' @param aspect_ratio A numeric value specifying the aspect ratio of the plot. If `NULL`, it's calculated automatically. Default is `NULL`.
#' @param base_width A numeric value specifying the base width for saving the plot. Default is `5`.
#' @param base_height A numeric value specifying the base height for saving the plot. Default is `4`.
#' @param output_file An optional string specifying the path to save the plot. If `NULL`, the plot is not saved. Default is `NULL`.
#'
#' @return A ggplot object representing the domino plot.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_gradient2 scale_size_continuous scale_x_continuous expansion scale_y_continuous labs theme_minimal theme element_text element_line element_blank element_rect annotate coord_flip coord_cartesian ggsave
#' @importFrom dplyr select mutate filter left_join arrange desc bind_rows
#' @importFrom utils globalVariables
#' @export
domino_plot <- function(data, 
                        gene_list, 
                        switch_axis = FALSE, 
                        min_dot_size = 1, 
                        max_dot_size = 5, 
                        spacing_factor = 3,
                        var_id = "var",
                        feature_col = "gene",
                        celltype_col = "Celltype",
                        contrast_col = "Contrast",
                        contrast_levels = c("Clinical", "Pathological"),
                        contrast_labels = c("Clinical", "Pathological"),
                        logfc_col = "avg_log2FC",
                        pval_col = "p_val_adj",
                        logfc_limits = c(-1.5, 1.5),
                        logfc_colors = c(low = "blue", mid = "white", high = "red"),
                        color_scale_name = "Log2 Fold Change",
                        size_scale_name = "-log10(adj. p-value)",
                        axis_text_size = 8,
                        aspect_ratio = NULL,
                        base_width = 5,
                        base_height = 4,
                        output_file = NULL) {
  # Ensure Contrast is a factor with correct levels
  data[[contrast_col]] <- factor(data[[contrast_col]], levels = contrast_levels)
  
  # Add feature_col if it doesn't exist
  if(!feature_col %in% colnames(data)) {
    data[[feature_col]] <- rownames(data)
  }
  
  # Filter data for specified genes
  data <- data %>% dplyr::filter(!!sym(feature_col) %in% gene_list)
  
  # Create a complete dataset
  all_celltypes <- unique(data[[celltype_col]])
  all_contrasts <- contrast_levels
  all_vars <- unique(data[[var_id]])
  
  # Create a complete data frame
  complete_data <- expand.grid(
    temp_gene = gene_list,
    temp_celltype = all_celltypes,
    temp_contrast = all_contrasts,
    temp_var = all_vars,
    stringsAsFactors = FALSE
  )
  names(complete_data) <- c(feature_col, celltype_col, contrast_col, var_id)
  
  # Join the complete dataset with the original data
  data <- dplyr::left_join(complete_data, data, by = c(feature_col, celltype_col, contrast_col, var_id))
  
  # Revised var_positions creation
  var_info <- data %>% select(!!sym(var_id), !!sym(contrast_col)) %>% distinct()
  
  var_positions <- var_info %>%
    mutate(
      x_offset = ifelse(!!sym(contrast_col) == contrast_levels[1], -0.2, 0.2),
      y_offset = ifelse(!!sym(contrast_col) == contrast_levels[1], 0.2, -0.2)
    )
  
  # Prepare the data
  data[[celltype_col]] <- factor(data[[celltype_col]], levels = rev(all_celltypes))
  
  plot_data <- data %>%
    dplyr::left_join(var_positions, by = c(var_id, contrast_col)) %>%
    dplyr::mutate(
      gene_index = match(!!sym(feature_col), gene_list),
      x_pos = dplyr::case_when(
        !!sym(contrast_col) == contrast_levels[1] ~ (gene_index - 1) * spacing_factor + 1 + x_offset,
        !!sym(contrast_col) == contrast_levels[2] ~ (gene_index - 1) * spacing_factor + 2 + x_offset
      ),
      y_pos = as.numeric(!!sym(celltype_col)) + y_offset
    )
  
  # Calculate the aspect ratio to make rectangles square
  n_celltypes <- length(all_celltypes)
  n_genes <- length(gene_list)
  
  if (is.null(aspect_ratio)) {
    aspect_ratio <- (n_celltypes) / (n_genes * spacing_factor)
  }
  
  # Prepare plot data
  plot_data <- plot_data %>%
    mutate(
      adj_logfc = pmax(pmin(!!sym(logfc_col), logfc_limits[2]), logfc_limits[1]),
      log_p_val_adj = -log10(!!sym(pval_col))
    ) 
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x_pos, y = y_pos)) +
    geom_rect(
      aes(
        xmin = dplyr::case_when(
          !!sym(contrast_col) == contrast_levels[1] ~ (gene_index - 1) * spacing_factor + 1 - 0.4,
          !!sym(contrast_col) == contrast_levels[2] ~ (gene_index - 1) * spacing_factor + 2 - 0.4
        ),
        xmax = dplyr::case_when(
          !!sym(contrast_col) == contrast_levels[1] ~ (gene_index - 1) * spacing_factor + 1 + 0.4,
          !!sym(contrast_col) == contrast_levels[2] ~ (gene_index - 1) * spacing_factor + 2 + 0.4
        ),
        ymin = as.numeric(!!sym(celltype_col)) - 0.4,
        ymax = as.numeric(!!sym(celltype_col)) + 0.4
      ),
      fill = "white", color = "grey", alpha = 0.5, linewidth = 0.5
    ) +
    geom_point(aes(color = adj_logfc, size = log_p_val_adj)) +
    geom_point(aes(size = log_p_val_adj), color = "black", shape = 1) + 
    scale_color_gradient2(
      name = color_scale_name,
      low = logfc_colors["low"], mid = logfc_colors["mid"], high = logfc_colors["high"],
      limits = logfc_limits
    ) +
    scale_size_continuous(name = size_scale_name, range = c(min_dot_size, max_dot_size)) +
    scale_x_continuous(
      breaks = seq(1.5, by = spacing_factor, length.out = length(gene_list)),
      labels = gene_list,
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      breaks = seq(1, length(unique(plot_data[[celltype_col]]))),
      labels = levels(plot_data[[celltype_col]]),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.text.y = element_text(size = axis_text_size),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(face = "bold"),
      aspect.ratio = aspect_ratio  # Set aspect ratio to make rectangles square
    ) +
    annotate("text", 
             x = seq(1, by = spacing_factor, length.out = length(gene_list)), 
             y = n_celltypes + 1, 
             label = contrast_labels[1], 
             angle = 90, hjust = 0, size = 3) +
    annotate("text", 
             x = seq(2, by = spacing_factor, length.out = length(gene_list)), 
             y = n_celltypes + 1, 
             label = contrast_labels[2], 
             angle = 90, hjust = 0, size = 3)
  
  if (switch_axis) {
    p <- p + coord_flip()
  } else {
    p <- p + coord_cartesian(clip = "off")  # Allow annotations outside plot area
  }
  
  if (!is.null(output_file)) {
    # Calculate dimensions based on number of genes and cell types
    n_genes <- length(gene_list)
    n_celltypes <- length(unique(plot_data[[celltype_col]]))
    
    # Base width and height
    width <- base_width + (n_genes * 0.5)  # Add 0.5 inch per gene
    height <- base_height + (n_celltypes * 0.2)  # Add 0.2 inch per cell type
    
    # Ensure minimum dimensions
    width <- max(width, base_width)  
    height <- max(height, base_height)  
    
    # Save the plot
    ggsave(output_file, plot = p, width = width, height = height, dpi = 300)
  }
  
  return(p)
}