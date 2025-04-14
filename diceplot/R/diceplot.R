utils::globalVariables(c(
  "unique_groups", "x_min", "x_max", "y_min", "y_max","var",
  "x_pos", "y_pos", "x_offset", "y_offset", "x", "y",
  "present", "combined", "value", "count", "gene", "Contrast",
  "Celltype", "avg_log2FC", "p_val_adj", "log_p_val_adj"
))

#' Dice Plot Visualization
#'
#' This function generates a custom plot based on three categorical variables and a group variable. It adapts to the number of unique categories in `z` and allows customization of various plot aesthetics.
#'
#' @param data A data frame containing the categorical and group variables for plotting.
#' @param x A string representing the column name in `data` for the first categorical variable.
#' @param y A string representing the column name in `data` for the second categorical variable.
#' @param z A string representing the column name in `data` for the third categorical variable.
#' @param group A string representing the column name in `data` for the grouping variable.
#' @param group_alpha A numeric value for the transparency level of the group rectangles. Default is `0.5`.
#' @param title An optional string for the plot title. Defaults to `NULL`.
#' @param z_colors A named vector of colors for `z` categories or a string to chose a colorbrewer palette. Defaults to `NULL` using the first suitable colorbrewer palette to use.
#' @param group_colors A named vector of colors for the group variableor a string to chose a colorbrewer palette. Defaults to `NULL` using the first suitable colorbrewer palette to use.
#' @param custom_theme A ggplot2 theme for customizing the plot's appearance. Defaults to `theme_minimal()`.
#' @param max_dot_size Maximal dot size for the plot to scale the dot sizes.
#' @param min_dot_size Minimal dot size for the plot to scale the dot sizes.
#' @param legend_width Relative width of your legend. Default is 0.25.
#' @param legend_height Relative width of your legend. Default is 0.5.
#' @param base_width_per_x Used for dynamically scaling the width. Default is 0.5.
#' @param base_height_per_y Used for dynamically scaling the height. Default is 0.3.
#' @param reverse_ordering Should the cluster ordering be reversed?. Default is FALSE.
#' @param cluster_by_row Cluster rows, defaults to TRUE
#' @param cluster_by_column Cluster columns, defaults to TRUE
#' @param show_legend Do you want to show the legend? Default is TRUE
#' @param cat_a Deprecated. Use `x` instead.
#' @param cat_b Deprecated. Use `y` instead.
#' @param cat_c Deprecated. Use `z` instead.
#' @param cat_c_colors Deprecated. Use `z_colors` instead.
#' @param cat_b_order Deprecated. Use `cluster_by_row` instead. Will be removed in a future version.
#' @param base_width_per_cat_a Deprecated. Use `base_width_per_x` instead.
#' @param base_height_per_cat_b Deprecated. Use `base_height_per_y` instead.
#'
#' @return A ggplot object representing the dice plot.
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_manual scale_fill_manual scale_x_discrete scale_y_discrete theme element_text element_blank unit labs coord_fixed ggtitle guides ggsave theme_minimal
#' @importFrom cowplot ggdraw draw_plot
#' @importFrom ggplot2 theme_minimal ggsave
#' @importFrom grid unit
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @importFrom RColorBrewer brewer.pal
#' @export
dice_plot <- function(data, 
                      x = NULL, 
                      y = NULL, 
                      z = NULL, 
                      group = NULL, 
                      group_alpha = 0.5,
                      title = NULL,
                      z_colors = NULL, 
                      group_colors = NULL, 
                      custom_theme = theme_minimal(),
                      max_dot_size = 5,
                      min_dot_size = 2,
                      legend_width = 0.25,
                      legend_height = 0.5,
                      base_width_per_x = 0.5,  
                      base_height_per_y = 0.3,
                      reverse_ordering = FALSE,
                      cluster_by_row = TRUE,
                      cluster_by_column = TRUE,
                      show_legend = TRUE,
                      cat_a = NULL,
                      cat_b = NULL,
                      cat_c = NULL,
                      cat_c_colors = NULL,
                      cat_b_order = NULL,
                      base_width_per_cat_a = NULL,
                      base_height_per_cat_b = NULL
) {
  
  # Handle deprecated parameters with warnings
  if (!is.null(cat_a)) {
    warning("The argument 'cat_a' is deprecated and will be removed in a future version >v1.5. Please use 'x' instead.", 
            call. = FALSE, immediate. = TRUE)
    x <- cat_a
  }
  
  if (!is.null(cat_b)) {
    warning("The argument 'cat_b' is deprecated and will be removed in a future version >v1.5. Please use 'y' instead.", 
            call. = FALSE, immediate. = TRUE)
    y <- cat_b
  }
  
  if (!is.null(cat_c)) {
    warning("The argument 'cat_c' is deprecated and will be removed in a future version >v1.5. Please use 'z' instead.", 
            call. = FALSE, immediate. = TRUE)
    z <- cat_c
  }
  
  if (!is.null(cat_c_colors)) {
    warning("The argument 'cat_c_colors' is deprecated and will be removed in a future version >v1.5. Please use 'z_colors' instead.", 
            call. = FALSE, immediate. = TRUE)
    z_colors <- cat_c_colors
  }
  
  if (!is.null(base_width_per_cat_a)) {
    warning("The argument 'base_width_per_cat_a' is deprecated and will be removed in a future version >v1.5. Please use 'base_width_per_x' instead.", 
            call. = FALSE, immediate. = TRUE)
    base_width_per_x <- base_width_per_cat_a
  }
  
  if (!is.null(base_height_per_cat_b)) {
    warning("The argument 'base_height_per_cat_b' is deprecated and will be removed in a future version >v1.5. Please use 'base_height_per_y' instead.", 
            call. = FALSE, immediate. = TRUE)
    base_height_per_y <- base_height_per_cat_b
  }
  
  if (!is.null(cat_b_order)) {
    warning("The argument 'cat_b_order' is deprecated and will be removed in a future version >v1.5. Please use 'cluster_by_row' instead.", 
            call. = FALSE, immediate. = TRUE)
  }
  
  num_vars <- length(unique(data[[z]]))
  z_levels = unique(data[[z]])
  if (is.null(z_colors)) {
    available_colors <- RColorBrewer::brewer.pal.info
    suitable_palettes <- rownames(available_colors[available_colors$category == "qual" & available_colors$maxcolors >= num_vars, ])
    if (length(suitable_palettes) == 0) {
      stop("No suitable ColorBrewer palette found for the number of categories in z.")
    }
    palette_name <- suitable_palettes[1]  
    default_z_colors <- RColorBrewer::brewer.pal(n = num_vars, name = palette_name)
    names(default_z_colors) <- unique(data[[z]])
    z_colors <- default_z_colors
  } else if (is.character(z_colors) && length(z_colors) == 1) {
    if (!z_colors %in% rownames(RColorBrewer::brewer.pal.info)) {
      stop(paste("The specified palette '", z_colors, "' is not a valid ColorBrewer palette.", sep = ""))
    }
    max_colors_in_palette <- RColorBrewer::brewer.pal.info[z_colors, "maxcolors"]
    if (num_vars > max_colors_in_palette) {
      stop(paste("The specified palette '", z_colors, "' does not have enough colors (needs ", num_vars, ", but only has ", max_colors_in_palette, ").", sep = ""))
    }
    palette_colors <- RColorBrewer::brewer.pal(n = num_vars, name = z_colors)
    names(palette_colors) <- z_levels
    z_colors <- palette_colors
  } else {
    # z_colors is assumed to be an explicit color palette
    if (length(z_colors) != num_vars) {
      stop("The length of z_colors does not match the number of categories in z.")
    }
    if (is.null(names(z_colors))) {
      names(z_colors) <- z_levels
    }
    z_colors <- z_colors[z_levels]
    
  }
  
  if (!is.null(group)) {
    if (is.null(group_colors)) {
      unique_groups <- unique(data[[group]])
      num_groups <- length(unique_groups)
      if (num_groups > 9) {
        stop("The number of groups exceeds the maximum colors available in ColorBrewer palettes.")
      }
      available_colors <- RColorBrewer::brewer.pal.info
      suitable_palettes <- rownames(available_colors[
        available_colors$category == "qual" & available_colors$maxcolors >= num_groups, 
      ])
      if (length(suitable_palettes) == 0) {
        stop("No suitable ColorBrewer palette found for the number of groups.")
      }
      palette_name <- suitable_palettes[2]  # Select a palette (e.g., the second one)
      group_colors_palette <- RColorBrewer::brewer.pal(n = num_groups, name = palette_name)
      names(group_colors_palette) <- unique_groups
      group_colors <- group_colors_palette
    } else if (is.character(group_colors) && length(group_colors) == 1) {
      if (!group_colors %in% rownames(RColorBrewer::brewer.pal.info)) {
        stop(paste("The specified palette '", group_colors, "' is not a valid ColorBrewer palette.", sep = ""))
      }
      max_colors_in_palette <- RColorBrewer::brewer.pal.info[group_colors, "maxcolors"]
      if (num_groups > max_colors_in_palette) {
        stop(paste("The specified palette '", group_colors, "' does not have enough colors (needs ", num_groups, ", but only has ", max_colors_in_palette, ").", sep = ""))
      }
      group_colors_palette <- RColorBrewer::brewer.pal(n = num_groups, name = group_colors)
      unique_groups <- unique(data[[group]])
      names(group_colors_palette) <- unique_groups
      group_colors <- group_colors_palette
    } else {
      unique_groups <- unique(data[[group]])
      num_groups <- length(unique_groups)
      if (length(group_colors) != num_groups) {
        stop("The length of group_colors does not match the number of groups.")
      }
      if (is.null(names(group_colors))) {
        names(group_colors) <- unique_groups
      } else {
        if (!all(unique_groups %in% names(group_colors))) {
          stop("The names of group_colors do not match the groups in the data.")
        }
        group_colors <- group_colors[unique_groups]
      }
    }
    
    # Ensure group is a factor with levels matching group_colors
    data[[group]] <- factor(data[[group]], levels = names(group_colors))
  }
  
  # Ensure consistent ordering of factors
  if(!is.factor(data[[x]])) {
    data[[x]] <- factor(data[[x]], levels = unique(data[[x]]))
  }
  
  if(!is.factor(data[[y]])) {
    data[[y]] <- factor(data[[y]], levels = unique(data[[y]]))
  }
  
  if(!is.factor(data[[z]])) {
    data[[z]] <- factor(data[[z]], levels = names(z_colors))
  }
  if (!is.null(group)) {
    # Check for unique group per y
    group_check <- data %>%
      group_by(!!sym(y)) %>%
      summarise(unique_groups = n_distinct(!!sym(group)), .groups = "drop") %>%
      filter(unique_groups > 1)
    
    if (nrow(group_check) > 0) {
      warning("Warning: The following y categories have multiple groups assigned:\n",
              paste(group_check[[y]], collapse = ", "))
    }
  }
  
  # Define variable positions dynamically
  var_positions <- create_var_positions(z_colors, num_vars)
  
  if(cluster_by_row){
    print("cluster by row")
    print(data)
    y_order <- order_cat_b(data, group, y, group_colors, reverse_ordering)
  } else {
    y_order <- levels(data[[y]])
  }
  
  if(cluster_by_column){
    x_order <- perform_clustering(data, x, y, z)
  } else {
    x_order <- levels(data[[x]])
  }
  
  plot_data <- prepare_plot_data(data, x, y, z, group, var_positions, x_order, y_order)
  
  # Always create box_data, but in different ways depending on group
  if (!is.null(group)) {
    box_data <- prepare_box_data(data, x, y, group, x_order, y_order)
  } else {
    box_data <- prepare_simple_box_data(data, x, y, x_order, y_order)
  }
  
  dot_size <- calculate_dot_size(num_vars, max_dot_size, min_dot_size)
  
  p <- ggplot()
  
  # Add rectangles with or without group fill
  if (!is.null(group)) {
    p <- p +
      geom_rect(data = box_data, 
                aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = !!sym(group)),
                color = "grey", alpha = group_alpha, linewidth = 0.5)
  } else {
    # When group is NULL, use white boxes
    p <- p +
      geom_rect(data = box_data, 
                aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
                fill = "white", color = "grey", alpha = group_alpha, linewidth = 0.5)
  }
  
  # Add points for z
  p <- p +
    geom_point(data = plot_data, 
               aes(x = x_pos, y = y_pos, color = !!sym(z)), 
               size = dot_size, show.legend = FALSE) +
    geom_point(data = plot_data, 
               aes(x = x_pos, y = y_pos), 
               size = dot_size + 0.5, shape = 1, color = "black", show.legend = FALSE) +
    scale_color_manual(values = z_colors, name = z, breaks = names(z_colors), guide = "none")
  
  # Add fill scale for groups if group is provided
  if (!is.null(group)) {
    p <- p +
      scale_fill_manual(values = group_colors, name = group)
  }
  
  p <- p +
    scale_x_discrete(limits = levels(plot_data[[x]])) +
    scale_y_discrete(limits = levels(plot_data[[y]])) +
    custom_theme +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    ) +
    labs(x = "", y = "") +
    coord_fixed(ratio = 1) +
    ggtitle(title)
  
  # Add guides based on whether group is provided
  if (!is.null(group) && show_legend) {
    p <- p + guides(fill = "none")
  }
  
  # Remove legends if show_legend is FALSE
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  # Create custom legends only if group is provided and show legend is true
  if (show_legend) {
    combined_legend_plot <- create_custom_legends(
      data, z, group, z_colors, group_colors, var_positions, num_vars, dot_size
    )
    
    # Combine the main plot and legends without 'preserve = "aspect"'
    combined_plot <- ggdraw() +
      draw_plot(
        p, 
        x = 0, 
        y = 0, 
        width = 1 - legend_width, 
        height = 1
      ) +
      draw_plot(
        combined_legend_plot, 
        x = 1 - legend_width, 
        y = (1 - legend_height) / 2,  # Center vertically
        width = legend_width, 
        height = legend_height
      )
  } else {
    combined_plot <- p
  }
  
  return(combined_plot)
}