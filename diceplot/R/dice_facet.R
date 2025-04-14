utils::globalVariables(c("label_x", "label_y", "adj_logfc", "feature_id", "gene_index", "celltype_numeric","mean_logfc", "x_offset", "y_offset", "aes", "case_when"))

#' Domino Plot Visualization with Categorical Colors
#'
#' This function generates a plot to visualize categorical data in a domino plot format. The size of the dots is fixed, and the plot can be saved to an output file if specified.
#' This version supports categorical colors and allows setting colors for left and right rectangle plots.
#'
#' @param data A data frame containing the categorical data.
#' @param gene_list A character vector of gene names to include in the plot.
#' @param x A string representing the column name in `data` for the feature variable (e.g., genes). Default is `"gene"`.
#' @param y A string representing the column name in `data` for the cell type variable. Default is `"Celltype"`.
#' @param contrast A string representing the column name in `data` for the contrast variable. Default is `"Contrast"`.
#' @param var_id A string representing the column name in `data` for the variable identifier. Default is `"var"`.
#' @param spacing_factor A numeric value indicating the spacing between gene pairs. Default is `3`.
#' @param categorical_colors A named vector of colors to use for categorical values in the data. Default is NULL.
#' @param color_scale_name A string specifying the name of the color scale in the legend. Default is `"Category"`.
#' @param left_rect_color A string specifying the color for the left rectangles. Default is `"lightblue"`.
#' @param right_rect_color A string specifying the color for the right rectangles. Default is `"lightpink"`.
#' @param rect_alpha A numeric value between 0 and 1 indicating the transparency of the rectangles. Default is `0.5`.
#' @param axis_text_size A numeric value specifying the size of the axis text. Default is `8`.
#' @param x_axis_text_size A numeric value specifying the size of the x-axis text. If NULL, uses `axis_text_size`. Default is `NULL`.
#' @param y_axis_text_size A numeric value specifying the size of the y-axis text. If NULL, uses `axis_text_size`. Default is `NULL`.
#' @param legend_text_size A numeric value specifying the size of the legend text. Default is `8`.
#' @param cluster_method The clustering method to use. Default is `"complete"`.
#' @param cluster_y_axis A logical value indicating whether to cluster the y-axis (cell types). Default is `TRUE`.
#' @param cluster_var_id A logical value indicating whether to cluster the var_id. Default is `TRUE`.
#' @param base_width A numeric value specifying the base width for saving the plot. Default is `5`.
#' @param base_height A numeric value specifying the base height for saving the plot. Default is `4`.
#' @param show_legend A logical value indicating whether to show the legend. Default is `TRUE`.
#' @param legend_width A numeric value specifying the relative width of the legend. Default is `0.25`.
#' @param legend_height A numeric value specifying the relative height of the legend. Default is `0.5`.
#' @param custom_legend A logical value indicating whether to use a custom legend. Default is `TRUE`.
#' @param aspect_ratio A numeric value specifying the aspect ratio of the plot. If `NULL`, it's calculated automatically. Default is `NULL`.
#' @param switch_axis A logical value indicating whether to switch the x and y axes. Default is `FALSE`.
#' @param reverse_y_ordering A logical value indicating whether to reverse the y-axis ordering after clustering. Default is `FALSE`.
#' @param show_var_positions A logical value indicating whether to show the intermediate variable positions plot. Default is `FALSE`.
#'   When `output_file` is specified with a PDF extension, both plots will be saved to a multi-page PDF if this is `TRUE`.
#'   A warning will be shown if `show_var_positions` is `TRUE` but the output file is not a PDF.
#' @param output_file An optional string specifying the path to save the plot. If `NULL`, the plot is not saved. Default is `NULL`.
#' @param feature_col Deprecated. Use `x` instead.
#' @param celltype_col Deprecated. Use `y` instead.
#' @param contrast_col Deprecated. Use `contrast` instead.
#'
#' @return A list containing the domino plot and optionally the variable positions plot.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid ggdraw draw_plot
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_manual scale_x_continuous expansion scale_y_continuous labs theme_minimal theme element_text element_line element_blank element_rect annotate coord_flip coord_cartesian ggsave
#' @importFrom dplyr select mutate filter left_join arrange desc bind_rows case_when
#' @importFrom utils globalVariables
#' @importFrom grDevices pdf dev.off
#' @export
dice_facet_plot <- function(data, 
                            gene_list, 
                            x = "gene",
                            y = "Celltype",
                            contrast = "Contrast",
                            var_id = "var",
                            spacing_factor = 3,
                            categorical_colors = NULL,
                            color_scale_name = "Category",
                            left_rect_color = "lightblue",
                            right_rect_color = "lightpink",
                            rect_alpha = 0.5,
                            axis_text_size = 8,
                            x_axis_text_size = NULL,
                            y_axis_text_size = NULL,
                            legend_text_size = 8,
                            cluster_method = "complete",
                            cluster_y_axis = TRUE,
                            cluster_var_id = TRUE,
                            base_width = 5,
                            base_height = 4,
                            show_legend = TRUE,
                            legend_width = 0.25,
                            legend_height = 0.5,
                            custom_legend = TRUE,
                            aspect_ratio = NULL,
                            switch_axis = FALSE,
                            reverse_y_ordering = FALSE,
                            show_var_positions = FALSE,
                            output_file = NULL,
                            feature_col = NULL,
                            celltype_col = NULL,
                            contrast_col = NULL) {
  
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
  
  # Set axis text sizes if NULL
  if (is.null(x_axis_text_size)) {
    x_axis_text_size <- axis_text_size
  }
  
  if (is.null(y_axis_text_size)) {
    y_axis_text_size <- axis_text_size
  }
  
  # Get unique contrast levels from the data
  contrast_levels <- unique(data[[contrast]])
  
  if (length(contrast_levels) != 2){
    warning("Only two contrast levels are supported right now", 
            call. = FALSE, immediate. = TRUE)
    return()
  }
  
  # Use the actual contrast levels from the data for labels
  contrast_labels <- contrast_levels
  
  if(!x %in% colnames(data)) {
    data$gene <- rownames(data)
  }
  
  data[[contrast]] <- factor(data[[contrast]], levels = contrast_levels)
  
  if(!x %in% colnames(data)) {
    warning("feature column not present!", 
            call. = FALSE, immediate. = TRUE)
    return()
  }
  
  # Filter data to only include genes in gene_list before clustering
  data <- data %>% filter(!!sym(x) %in% gene_list)
  
  # Implement cell type (y-axis) clustering if requested
  if (cluster_y_axis) {
    # Create a profile matrix for each cell type with gene-contrast combinations
    celltype_profile_matrix <- data %>%
      # Create a unique identifier for each gene-contrast combination
      mutate(feature_id = paste(!!sym(x), !!sym(contrast), sep = "_")) %>%
      # Group by y value and feature_id
      group_by(!!sym(y), feature_id) %>%
      # Calculate the count for each y-feature combination
      summarise(count = n(), .groups = "drop") %>%
      # Pivot to wide format with y values as rows and feature combinations as columns
      tidyr::pivot_wider(
        id_cols = !!sym(y),
        names_from = feature_id,
        values_from = count,
        values_fill = 0L  # Use integer 0 instead of double
      )
    
    # Extract row names for y categories
    celltype_categories <- celltype_profile_matrix[[y]]
    
    # Remove the y column to get just the numeric matrix
    celltype_profile_matrix <- as.matrix(celltype_profile_matrix[, -1])
    rownames(celltype_profile_matrix) <- celltype_categories
    
    # Compute distance matrix
    dist_matrix <- dist(celltype_profile_matrix)
    
    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = cluster_method)
    
    # Get the ordering from the dendrogram
    y_order <- celltype_categories[hc$order]
    
    # Reverse the order if requested
    if (reverse_y_ordering) {
      y_order <- rev(y_order)
    }
    
    # Set the ordered levels for the y factor
    all_celltypes <- y_order
    
    # Ensure celltype is a factor with the clustered order
    data[[y]] <- factor(data[[y]], levels = all_celltypes)
  } else {
    # Check if y is already a factor
    if (is.factor(data[[y]])) {
      # If it's already a factor, keep its levels
      all_celltypes <- levels(data[[y]])
    } else {
      # If not a factor, convert it to factor with default ordering
      all_celltypes <- unique(data[[y]])
      data[[y]] <- factor(data[[y]], levels = all_celltypes)
    }
  }
  
  # Ensure celltype is a factor for proper ordering
  data[[y]] <- factor(data[[y]], levels = all_celltypes)
  
  # Handle var_id clustering or factorization
  if (cluster_var_id) {
    # Cluster var_id based on their profiles
    cluster_data <- data %>%
      # Select the columns we need
      select(!!sym(var_id), !!sym(y), !!sym(contrast), !!sym(x)) %>%
      # Create a unique identifier for each gene-celltype-contrast combination
      mutate(feature_id = paste(!!sym(x), !!sym(y), !!sym(contrast), sep = "_")) %>%
      # Group by var_id and feature_id to get counts
      group_by(!!sym(var_id), feature_id) %>%
      summarise(count = n(), .groups = "drop") %>%
      # Pivot to wide format with var_id as rows and feature combinations as columns
      tidyr::pivot_wider(
        id_cols = !!sym(var_id),
        names_from = feature_id,
        values_from = count,
        values_fill = 0L  # Use integer 0 instead of double
      )
    
    # Extract the matrix for clustering (everything except the var_id column)
    if(ncol(cluster_data) <= 1) {
      warning("Not enough data for clustering", call. = FALSE, immediate. = TRUE)
      return()
    }
    
    cluster_var_ids <- cluster_data[[var_id]]
    cluster_matrix_values <- as.matrix(cluster_data[, -1])
    rownames(cluster_matrix_values) <- cluster_var_ids
    
    # Compute distance matrix
    dist_matrix <- dist(cluster_matrix_values)
    
    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = cluster_method)
    
    # Get the ordering from the dendrogram
    var_order <- cluster_var_ids[hc$order]
    
    # Update the var_id factor in the original data to reflect the clustering order
    data[[var_id]] <- factor(data[[var_id]], levels = var_order)
  } else {
    # Check if var_id is already a factor
    if (is.factor(data[[var_id]])) {
      # If it's already a factor, keep its levels
      # No need to do anything
    } else {
      # If not a factor, convert it to factor with default ordering
      data[[var_id]] <- factor(data[[var_id]], levels = unique(data[[var_id]]))
    }
  }
  
  all_contrasts <- contrast_levels
  all_vars <- levels(data[[var_id]])
  
  # Only work with actual data points
  data <- data %>%
    filter(!is.na(!!sym(var_id))) %>%  # Remove NA values
    mutate(
      gene_index = match(!!sym(x), gene_list),
      celltype_numeric = as.numeric(factor(!!sym(y), levels = all_celltypes))
    )
  
  # Create one row per (var, contrast) and define custom offsets:
  var_info <- data %>%
    distinct(!!sym(var_id), !!sym(contrast)) 
  
  # Create named lists for var_positions
  var_list_left <- setNames(
    as.list(unique(data[[var_id]][data[[contrast]] == contrast_levels[1]])),
    unique(data[[var_id]][data[[contrast]] == contrast_levels[1]])
  )
  
  var_list_right <- setNames(
    as.list(unique(data[[var_id]][data[[contrast]] == contrast_levels[2]])),
    unique(data[[var_id]][data[[contrast]] == contrast_levels[2]])
  )
  
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
  
  var_positions <- dplyr::bind_rows(var_positions_left, var_positions_right)
  
  # Create the variable positions plot
  var_pos_plot <- ggplot(var_positions, aes(x = x_offset, y = y_offset)) +
    geom_point(size = 5) +
    geom_text(aes(label = !!sym(var_id)), hjust = -0.3, vjust = 0) +
    scale_color_manual(values = c(var_list_left, var_list_right)) +
    theme_minimal() +
    labs(
      title = "Variable Positions Visualization",
      subtitle = "Showing positions for two groups with variables",
      x = "X Position",
      y = "Y Position"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = legend_text_size)
    ) +
    coord_fixed()
  
  # Join the position information with the actual data
  plot_data <- data %>%
    dplyr::left_join(var_positions, by = c(var_id, contrast)) %>%
    dplyr::mutate(
      x_pos = (gene_index - 1) * spacing_factor + x_offset,
      y_pos = celltype_numeric + y_offset
    )
  
  n_celltypes <- length(all_celltypes)
  n_genes <- length(gene_list)
  if (is.null(aspect_ratio)) {
    aspect_ratio <- n_celltypes / (n_genes * spacing_factor)
  }
  
  # If categorical_colors is NULL, we need to create default colors
  if (is.null(categorical_colors)) {
    # Get unique categories from the data
    categories <- unique(na.omit(data[[var_id]]))
    if (length(categories) == 0) {
      # If there are no categories (all NAs), create a default category
      categories <- c("No data")
      data[[var_id]] <- "No data"
      # Default color for "No data" category
      categorical_colors <- c("No data" = "grey")
    } else {
      # Create default colors using RColorBrewer or other color palette
      default_colors <- c("red", "blue", "green", "purple", "orange", "yellow", "brown", "pink", "cyan", "magenta")
      categorical_colors <- setNames(
        default_colors[1:min(length(categories), length(default_colors))],
        as.character(sort(categories))
      )
    }
  }
  
  # Create a complete grid for rectangles only
  rect_data <- expand.grid(
    gene_index = 1:length(gene_list),
    celltype_numeric = 1:length(all_celltypes),
    contrast = contrast_levels
  ) %>%
    mutate(
      x_pos = (gene_index - 1) * spacing_factor + ifelse(contrast == contrast_levels[1], 1, 2),
      y_pos = celltype_numeric
    )
  
  # Create the plot with rectangles and points
  p <- ggplot() +
    # Add rectangles for all possible positions
    geom_rect(
      data = rect_data,
      aes(
        xmin = x_pos - 0.4,
        xmax = x_pos + 0.4,
        ymin = y_pos - 0.4,
        ymax = y_pos + 0.4,
        fill = contrast
      ),
      color = "darkgrey",
      alpha = rect_alpha,
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = setNames(
        c(left_rect_color, right_rect_color),
        contrast_levels
      )
    ) +
    # Add points only for actual data
    geom_point(
      data = plot_data,
      aes(x = x_pos, y = y_pos, color = !!sym(var_id)),
      size = 3
    ) +
    geom_point(
      data = plot_data,
      aes(x = x_pos, y = y_pos),
      size = 3.5,
      shape = 1,
      color = "black"
    ) +
    scale_color_manual(
      name = color_scale_name,
      values = categorical_colors,
      na.value = "grey"
    ) +
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
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = x_axis_text_size),
      axis.text.y = element_text(size = y_axis_text_size),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_text_size + 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(face = "bold"),
      aspect.ratio = aspect_ratio
    ) +
    annotate("text", 
             x = seq(1, by = spacing_factor, length.out = length(gene_list)), 
             y = n_celltypes + 1, 
             label = contrast_labels[1], 
             angle = 90, hjust = 0, size = x_axis_text_size/3) +
    annotate("text", 
             x = seq(2, by = spacing_factor, length.out = length(gene_list)), 
             y = n_celltypes + 1, 
             label = contrast_labels[2], 
             angle = 90, hjust = 0, size = x_axis_text_size/3)
  
  if (switch_axis) {
    p <- p + coord_flip()
  } else {
    p <- p + coord_cartesian(clip = "off")
  }
  
  legend_var_positions <- create_var_positions(var_list_left, length(var_list_left)) %>%
    mutate(!!sym(var_id) := names(var_list_left),
           !!sym(contrast) := contrast_levels[1])
  
  # Handle legend display
  if (show_legend && custom_legend) {
    # Create custom legends for categorical data
    custom_legend_plot <- create_custom_domino_legends_categorical(
      contrast_levels = contrast_labels,
      var_positions = legend_var_positions,
      var_id = var_id,
      contrast = contrast,
      categorical_colors = categorical_colors,
      color_scale_name = color_scale_name,
      legend_text_size = legend_text_size,
      left_rect_color = left_rect_color,
      right_rect_color = right_rect_color
    )
    
    # Hide the default legends from the main plot
    p <- p + theme(legend.position = "none")
    
    # Combine the main plot and custom legends
    final_plot <- cowplot::ggdraw() +
      cowplot::draw_plot(
        p, 
        x = 0, 
        y = 0, 
        width = 1 - legend_width, 
        height = 1
      ) +
      cowplot::draw_plot(
        custom_legend_plot, 
        x = 1 - legend_width, 
        y = (1 - legend_height) / 2,  # Center vertically
        width = legend_width, 
        height = legend_height
      )
  } else if (show_legend && !custom_legend) {
    # Use default legends
    final_plot <- p
  } else {
    # No legends
    final_plot <- p + theme(legend.position = "none")
  }
  
  # Calculate dimensions for saving
  n_genes <- length(gene_list)
  n_celltypes <- length(unique(plot_data[[y]]))
  width <- base_width + (n_genes)
  height <- base_height + (n_celltypes * 0.5)
  width <- max(width, base_width)  
  height <- max(height, base_height)  
  
  # If showing legend, adjust width to accommodate it
  if(show_legend ) {
    width <- width / (1 - legend_width)
  }
  
  # Save plots as requested
  if (!is.null(output_file)) {
    # Check if output_file has a PDF extension for multi-page support
    is_pdf <- tolower(tools::file_ext(output_file)) == "pdf"
    
    if (show_var_positions && !is_pdf) {
      warning("The show_var_positions=TRUE option is only fully supported with PDF output files. ",
              "Only the main plot will be saved to '", output_file, "'. ",
              "Use a .pdf extension to save both plots in a multi-page document.",
              call. = FALSE, immediate. = TRUE)
      
      # Save just the main plot to the non-PDF file
      ggsave(output_file, plot = final_plot, width = width, height = height, dpi = 300)
    } else if (is_pdf && show_var_positions) {
      # Let the user know we're saving a multi-page PDF
      message("Saving both plots to a multi-page PDF: ", output_file)
      
      # Open PDF device
      grDevices::pdf(output_file, width = width, height = height)
      print(final_plot)
      print(var_pos_plot)
      grDevices::dev.off()
    } else {
      # Regular single plot output
      ggsave(output_file, plot = final_plot, width = width, height = height, dpi = 300)
    }
  }
  
  # Return the plots in a list
  result <- list(domino_plot = final_plot)
  if (show_var_positions) {
    result$var_positions_plot <- var_pos_plot
  }
  return(result)
}
