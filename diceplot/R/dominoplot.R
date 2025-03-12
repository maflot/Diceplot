utils::globalVariables(c("label_x", "label_y", "adj_logfc", "feature_id", "gene_index", "celltype_numeric", "x_offset", "y_offset", "aes", "case_when"))

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
#' @param x A string representing the column name in `data` for the feature variable (e.g., genes). Default is `"gene"`.
#' @param y A string representing the column name in `data` for the cell type variable. Default is `"Celltype"`.
#' @param contrast A string representing the column name in `data` for the contrast variable. Default is `"Contrast"`.
#' @param contrast_levels A character vector specifying the levels of the contrast variable. Default is `c("Clinical", "Pathological")`.
#' @param contrast_labels A character vector specifying the labels for the contrasts in the plot. Default is `c("Clinical", "Pathological")`.
#' @param log_fc A string representing the column name in `data` for the log fold change values. Default is `"avg_log2FC"`.
#' @param p_val A string representing the column name in `data` for the adjusted p-values. Default is `"p_val_adj"`.
#' @param logfc_limits A numeric vector of length 2 specifying the limits for the log fold change color scale. Default is `c(-1.5, 1.5)`.
#' @param logfc_colors A named vector specifying the colors for the low, mid, and high values in the color scale. Default is `c(low = "blue", mid = "white", high = "red")`.
#' @param color_scale_name A string specifying the name of the color scale in the legend. Default is `"Log2 Fold Change"`.
#' @param size_scale_name A string specifying the name of the size scale in the legend. Default is `"-log10(adj. p-value)"`.
#' @param axis_text_size A numeric value specifying the size of the axis text. Default is `8`.
#' @param aspect_ratio A numeric value specifying the aspect ratio of the plot. If `NULL`, it's calculated automatically. Default is `NULL`.
#' @param base_width A numeric value specifying the base width for saving the plot. Default is `5`.
#' @param base_height A numeric value specifying the base height for saving the plot. Default is `4`.
#' @param output_file An optional string specifying the path to save the plot. If `NULL`, the plot is not saved. Default is `NULL`.
#' @param cluster_method The clustering method to use. Default is `"complete"`.
#' @param feature_col Deprecated. Use `x` instead.
#' @param celltype_col Deprecated. Use `y` instead.
#' @param contrast_col Deprecated. Use `contrast` instead.
#' @param logfc_col Deprecated. Use `log_fc` instead.
#' @param pval_col Deprecated. Use `p_val` instead.
#'
#' @return A ggplot object representing the domino plot.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_gradient2 scale_size_continuous scale_x_continuous expansion scale_y_continuous labs theme_minimal theme element_text element_line element_blank element_rect annotate coord_flip coord_cartesian ggsave
#' @importFrom dplyr select mutate filter left_join arrange desc bind_rows case_when
#' @importFrom utils globalVariables
#' @export
domino_plot <- function(data, 
                        gene_list, 
                        switch_axis = FALSE, 
                        min_dot_size = 1, 
                        max_dot_size = 5, 
                        spacing_factor = 3,
                        var_id = "var",
                        x = "gene",
                        y = "Celltype",
                        contrast = "Contrast",
                        contrast_levels = c("Clinical", "Pathological"),
                        contrast_labels = c("Clinical", "Pathological"),
                        log_fc = "avg_log2FC",
                        p_val = "p_val_adj",
                        logfc_limits = c(-1.5, 1.5),
                        logfc_colors = c(low = "blue", mid = "white", high = "red"),
                        color_scale_name = "Log2 Fold Change",
                        size_scale_name = "-log10(adj. p-value)",
                        axis_text_size = 8,
                        aspect_ratio = NULL,
                        base_width = 5,
                        base_height = 4,
                        output_file = NULL,
                        cluster_method = "complete",
                        feature_col = NULL,
                        celltype_col = NULL,
                        contrast_col = NULL,
                        logfc_col = NULL,
                        pval_col = NULL) {
  
  # Handle deprecated parameters with warnings
  if (!is.null(feature_col)) {
    warning("The argument 'feature_col' is deprecated and will be removed in a future version >v1.5. Please use 'x' instead.", 
            call. = FALSE, immediate. = TRUE)
    x <- feature_col
  }
  
  if (!is.null(celltype_col)) {
    warning("The argument 'celltype_col' is deprecated and will be removed in a future version >v1.5. Please use 'y' instead.", 
            call. = FALSE, immediate. = TRUE)
    y <- celltype_col
  }
  
  if (!is.null(contrast_col)) {
    warning("The argument 'contrast_col' is deprecated and will be removed in a future version >v1.5. Please use 'contrast' instead.", 
            call. = FALSE, immediate. = TRUE)
    contrast <- contrast_col
  }
  
  if (!is.null(logfc_col)) {
    warning("The argument 'logfc_col' is deprecated and will be removed in a future version >v1.5. Please use 'log_fc' instead.", 
            call. = FALSE, immediate. = TRUE)
    log_fc <- logfc_col
  }
  
  if (!is.null(pval_col)) {
    warning("The argument 'pval_col' is deprecated and will be removed in a future version >v1.5. Please use 'p_val' instead.", 
            call. = FALSE, immediate. = TRUE)
    p_val <- pval_col
  }
  
  if (length(unique(data[[contrast]])) != 2){
    warning("Only two contrast levels are supported right now", 
            call. = FALSE, immediate. = TRUE)
    return()
  }
  
  if (!is.null(contrast_levels)) {
    warning("The argument 'contrast_levels' is deprecated and will be removed in a future version >v1.5", 
            call. = FALSE, immediate. = TRUE)
  } 
  if (!is.null(contrast_labels)) {
    warning("The argument 'contrast_labels' is deprecated and will be removed in a future version >v1.5", 
            call. = FALSE, immediate. = TRUE)
  }
  
  if(!x %in% colnames(data)) {
    data$gene <- rownames(data)
  }
  
  data[[contrast]] <- factor(data[[contrast]], levels = unique(data[[contrast]]))
  
  if(!x %in% colnames(data)) {
    warning("feature column not present!", 
            call. = FALSE, immediate. = TRUE)
    return()
  }
  
  cluster_data <- data %>%
    # Select the columns we need
    select(!!sym(var_id), !!sym(y), !!sym(contrast), !!sym(x), !!sym(log_fc)) %>%
    # Handle NAs in logFC by replacing with 0 (or another appropriate value)
    mutate(!!sym(log_fc) := ifelse(is.na(!!sym(log_fc)), 0, !!sym(log_fc)))
  
  # For each var_id, create a feature vector combining all genes, celltypes, and contrasts
  # This gives us a "profile" for each variable that we can cluster
  cluster_matrix <- cluster_data %>%
    # Create a unique identifier for each gene-celltype-contrast combination
    mutate(feature_id = paste(!!sym(x), !!sym(y), !!sym(contrast), sep = "_")) %>%
    # Pivot to wide format with var_id as rows and feature combinations as columns
    tidyr::pivot_wider(
      id_cols = !!sym(var_id),
      names_from = feature_id,
      values_from = !!sym(log_fc),
      values_fill = 0
    )
  
  # Extract the matrix for clustering (everything except the var_id column)
  if(ncol(cluster_matrix) <= 1) {
    warning("Not enough data for clustering", call. = FALSE, immediate. = TRUE)
    return()
  }
  
  cluster_var_ids <- cluster_matrix[[var_id]]
  cluster_matrix_values <- as.matrix(cluster_matrix[, -1])
  rownames(cluster_matrix_values) <- cluster_var_ids
  
  # Compute distance matrix
  dist_matrix <- dist(cluster_matrix_values)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = cluster_method)
  
  # Get the ordering from the dendrogram
  var_order <- cluster_var_ids[hc$order]
  
  # Update the var_id factor in the original data to reflect the clustering order
  data[[var_id]] <- factor(data[[var_id]], levels = var_order)
  
  
  data <- data %>% filter(!!sym(x) %in% gene_list)
  
  all_celltypes <- unique(data[[y]])
  # Ensure celltype is a factor for proper ordering
  data[[y]] <- factor(data[[y]], levels = all_celltypes)
  
  all_contrasts <- contrast_levels
  all_vars <- unique(data[[var_id]])
  
  complete_data <- expand.grid(
    temp_gene = gene_list,
    temp_celltype = all_celltypes,
    temp_contrast = all_contrasts,
    temp_var = all_vars,
    stringsAsFactors = FALSE
  )
  names(complete_data) <- c(x, y, contrast, var_id)
  
  data <- left_join(complete_data, data, by = c(x, y, contrast, var_id))
  
  # Create one row per (var, contrast) and define custom offsets:
  var_info <- data %>%
    distinct(!!sym(var_id), !!sym(contrast)) 
  
  
  # Create named lists for var_positions
  var_list_left <- setNames(
    as.list(unique(data[[var_id]][data[[contrast]] == unique(data[[contrast]])[1]])),
    unique(data[[var_id]][data[[contrast]] == unique(data[[contrast]])[1]])
  )
  
  var_list_right <- setNames(
    as.list(unique(data[[var_id]][data[[contrast]] == unique(data[[contrast]])[2]])),
    unique(data[[var_id]][data[[contrast]] == unique(data[[contrast]])[2]])
  )
  
  # Assuming create_var_positions is defined elsewhere in your code
  # Generate variable positions using named lists
  var_positions_left <- create_var_positions(
    var_list_left, 
    length(var_list_left)
  ) %>%
    mutate(
      !!sym(var_id) := names(var_list_left),
      !!sym(contrast) := contrast_levels[1],
      x_offset = 1 + x_offset
    )
  
  var_positions_right <- create_var_positions(
    var_list_right, 
    length(var_list_right)
  ) %>%
    mutate(
      !!sym(var_id) := names(var_list_right),
      !!sym(contrast) := contrast_levels[2],
      x_offset = 2 + x_offset 
    )
  
  spacing_factor <- 3
  var_positions <- dplyr::bind_rows(var_positions_left, var_positions_right)
  
  mini_plot <- ggplot(var_positions, aes(x = "x_offset", y = "y_offset", 
                                                label = var_id, color = contrast)) +
    geom_point() +
    geom_text(nudge_x = 0.05, nudge_y = 0.05, size = 3) +
    theme_minimal() +
    labs(title = "Dice Arrangement in var Positions", x = "x_offset", y = "y_offset")
  print(mini_plot)
  
  # Here we join the position information
  plot_data <- data %>%
    dplyr::left_join(var_positions, by = c(var_id, contrast)) %>%
    dplyr::mutate(
      gene_index = match(!!sym(x), gene_list),
      x_pos = (gene_index - 1) * spacing_factor + x_offset,
      # Make sure celltype_numeric is properly assigned based on the factor levels
      celltype_numeric = as.numeric(factor(!!sym(y), levels = all_celltypes)),
      y_pos = celltype_numeric + y_offset
    )
  
  n_celltypes <- length(all_celltypes)
  n_genes <- length(gene_list)
  if (is.null(aspect_ratio)) {
    aspect_ratio <- n_celltypes / (n_genes * spacing_factor)
  }
  
  plot_data <- plot_data %>%
    mutate(
      adj_logfc = pmax(pmin(!!sym(log_fc), logfc_limits[2]), logfc_limits[1]),
      log_p_val_adj = -log10(!!sym(p_val))
    )
  
  # Fix the rectangle drawing to properly use celltype_numeric
  p <- ggplot(plot_data, aes(x = x_pos, y = y_pos)) +
    geom_rect(
      aes(
        xmin = (gene_index - 1) * spacing_factor + case_when(
          !!sym(contrast) == contrast_levels[1] ~ 1 - 0.4,
          !!sym(contrast) == contrast_levels[2] ~ 2 - 0.4,
          TRUE ~ 0
        ),
        xmax = (gene_index - 1) * spacing_factor + case_when(
          !!sym(contrast) == contrast_levels[1] ~ 1 + 0.4,
          !!sym(contrast) == contrast_levels[2] ~ 2 + 0.4,
          TRUE ~ 0
        ),
        ymin = celltype_numeric - 0.4,  # Use the numeric value directly
        ymax = celltype_numeric + 0.4   # Use the numeric value directly
      ),
      fill = "white", color = "black", alpha = 0.5, linewidth = 0.5
    ) +
    geom_point(aes(color = adj_logfc, size = log_p_val_adj)) +
    geom_point(aes(size = log_p_val_adj), color = "black", shape = 1) + 
    scale_color_gradient2(
      name = color_scale_name,
      low = logfc_colors["low"],
      mid = logfc_colors["mid"],
      high = logfc_colors["high"],
      limits = logfc_limits
    ) +
    scale_size_continuous(name = size_scale_name, range = c(min_dot_size, max_dot_size)) +
    scale_x_continuous(
      breaks = seq(1.5, by = spacing_factor, length.out = length(gene_list)),
      labels = gene_list,
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    # Fix the y-axis scale to properly display celltype labels
    scale_y_continuous(
      breaks = 1:length(all_celltypes),
      labels = all_celltypes,
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
      aspect.ratio = aspect_ratio
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
  
  print(p)
  if (switch_axis) {
    p <- p + coord_flip()
  } else {
    p <- p + coord_cartesian(clip = "off")
  }
  
  if (!is.null(output_file)) {
    n_genes <- length(gene_list)
    n_celltypes <- length(unique(plot_data[[y]]))
    width <- base_width + (n_genes * 0.5)
    height <- base_height + (n_celltypes * 0.2)
    width <- max(width, base_width)  
    height <- max(height, base_height)  
    ggsave(output_file, plot = p, width = width, height = height, dpi = 300)
  }
  
  p
}