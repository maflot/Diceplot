utils::globalVariables(c(
  "unique_groups", "x_min", "x_max", "y_min", "y_max","var",
  "x_pos", "y_pos", "x_offset", "y_offset", "x", "y",
  "present", "combined", "value", "count", "gene", "Contrast",
  "Celltype", "avg_log2FC", "p_val_adj", "log_p_val_adj"
))

#' Dice Plot Visualization
#'
#' This function generates a custom plot based on three categorical variables and a group variable. It adapts to the number of unique categories in `cat_c` and allows customization of various plot aesthetics.
#'
#' @param data A data frame containing the categorical and group variables for plotting.
#' @param cat_a A string representing the column name in `data` for the first categorical variable.
#' @param cat_b A string representing the column name in `data` for the second categorical variable.
#' @param cat_c A string representing the column name in `data` for the third categorical variable.
#' @param group A string representing the column name in `data` for the grouping variable.
#' @param group_alpha A numeric value for the transparency level of the group rectangles. Default is `0.5`.
#' @param title An optional string for the plot title. Defaults to `NULL`.
#' @param cat_c_colors A named vector of colors for `cat_c` categories or a string to chose a colorbrewer palette. Defaults to `NULL` using the first suitable colorbrewer palette to use.
#' @param group_colors A named vector of colors for the group variableor a string to chose a colorbrewer palette. Defaults to `NULL` using the first suitable colorbrewer palette to use.
#' @param custom_theme A ggplot2 theme for customizing the plot's appearance. Defaults to `theme_minimal()`.
#' @param max_dot_size Maximal dot size for the plot to scale the dot sizes.
#' @param min_dot_size Minimal dot size for the plot to scale the dot sizes.
#' @param legend_width Relative width of your legend. Default is 0.25.
#' @param legend_height Relative width of your legend. Default is 0.5.
#' @param base_width_per_cat_a Used for dynamically scaling the width. Default is 0.5.
#' @param base_height_per_cat_b Used for dynamically scaling the height. Default is 0.3.
#' @param reverse_ordering Should the cluster ordering be reversed?. Default is FALSE.
#' @param cat_b_order Do you want to pass an explicit order?. Default is NULL.
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
                      cat_a, 
                      cat_b, 
                      cat_c, 
                      group = NULL, 
                      group_alpha = 0.5,
                      title = NULL,
                      cat_c_colors = NULL, 
                      group_colors = NULL, 
                      custom_theme = theme_minimal(),
                      max_dot_size = 5,
                      min_dot_size = 2,
                      legend_width = 0.25,
                      legend_height = 0.5,
                      base_width_per_cat_a = 0.5,  
                      base_height_per_cat_b = 0.3,
                      reverse_ordering = FALSE,
                      cat_b_order = NULL
                      ) {
  
  num_vars <- length(unique(data[[cat_c]]))
  cat_c_levels = unique(data[[cat_c]])
  if (is.null(cat_c_colors)) {
    available_colors <- RColorBrewer::brewer.pal.info
    suitable_palettes <- rownames(available_colors[available_colors$category == "qual" & available_colors$maxcolors >= num_vars, ])
    if (length(suitable_palettes) == 0) {
      stop("No suitable ColorBrewer palette found for the number of categories in cat_c.")
    }
    palette_name <- suitable_palettes[1]  
    default_cat_c_colors <- RColorBrewer::brewer.pal(n = num_vars, name = palette_name)
    names(default_cat_c_colors) <- unique(data[[cat_c]])
    cat_c_colors <- default_cat_c_colors
  } else if (is.character(cat_c_colors) && length(cat_c_colors) == 1) {
    if (!cat_c_colors %in% rownames(RColorBrewer::brewer.pal.info)) {
      stop(paste("The specified palette '", cat_c_colors, "' is not a valid ColorBrewer palette.", sep = ""))
    }
    max_colors_in_palette <- RColorBrewer::brewer.pal.info[cat_c_colors, "maxcolors"]
    if (num_vars > max_colors_in_palette) {
      stop(paste("The specified palette '", cat_c_colors, "' does not have enough colors (needs ", num_vars, ", but only has ", max_colors_in_palette, ").", sep = ""))
    }
    palette_colors <- RColorBrewer::brewer.pal(n = num_vars, name = cat_c_colors)
    names(palette_colors) <- cat_c_levels
    cat_c_colors <- palette_colors
  } else {
    # cat_c_colors is assumed to be an explicit color palette
    if (length(cat_c_colors) != num_vars) {
      stop("The length of cat_c_colors does not match the number of categories in cat_c.")
    }
    if (is.null(names(cat_c_colors))) {
      names(cat_c_colors) <- cat_c_levels
    }
    cat_c_colors <- cat_c_colors[cat_c_levels]
  
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
  data[[cat_a]] <- factor(data[[cat_a]], levels = unique(data[[cat_a]]))
  data[[cat_b]] <- factor(data[[cat_b]], levels = unique(data[[cat_b]]))
  data[[cat_c]] <- factor(data[[cat_c]], levels = names(cat_c_colors))
  
  if (!is.null(group)) {
    # Check for unique group per cat_b
    group_check <- data %>%
      group_by(!!sym(cat_b)) %>%
      summarise(unique_groups = n_distinct(!!sym(group)), .groups = "drop") %>%
      filter(unique_groups > 1)
    
    if (nrow(group_check) > 0) {
      warning("Warning: The following cat_b categories have multiple groups assigned:\n",
              paste(group_check[[cat_b]], collapse = ", "))
    }
  }
  
  # Define variable positions dynamically
  var_positions <- create_var_positions(cat_c_colors, num_vars)
  cat_a_order <- perform_clustering(data, cat_a, cat_b, cat_c)
  if (!is.null(group) & is.null(cat_b_order)) {
    cat_b_order <- order_cat_b(data, group, cat_b, group_colors, reverse_ordering)
  } else if (!is.null(cat_b_order)){
    cat_b_order = cat_b_order
  } else {
    cat_b_order <- levels(data[[cat_b]])
  }
  
  plot_data <- prepare_plot_data(data, cat_a, cat_b, cat_c, group, var_positions, cat_a_order, cat_b_order)
  if (!is.null(group)) {
    box_data <- prepare_box_data(data, cat_a, cat_b, group, cat_a_order, cat_b_order)
  }
  
  dot_size <- calculate_dot_size(num_vars,max_dot_size,min_dot_size)
  
  p <- ggplot()
  
  if (!is.null(group)) {
    p <- p +
      geom_rect(data = box_data, 
                aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = !!sym(group)),
                color = "grey", alpha = group_alpha, linewidth = 0.5)
  }
  
  # Add points for cat_c
  p <- p +
    geom_point(data = plot_data, 
               aes(x = x_pos, y = y_pos, color = !!sym(cat_c)), 
               size = dot_size, show.legend = FALSE) +
    geom_point(data = plot_data, 
               aes(x = x_pos, y = y_pos), 
               size = dot_size + 0.5, shape = 1, color = "black", show.legend = FALSE) +
    scale_color_manual(values = cat_c_colors, name = cat_c, breaks = names(cat_c_colors), guide = "none")
  
  # Add fill scale for groups if group is provided
  if (!is.null(group)) {
    p <- p +
      scale_fill_manual(values = group_colors, name = group)
  }
  
  p <- p +
    scale_x_discrete(limits = levels(plot_data[[cat_a]])) +
    scale_y_discrete(limits = levels(plot_data[[cat_b]])) +
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
  if (!is.null(group)) {
    p <- p + guides(fill = "none")
  }
  
  # Create custom legends only if group is provided
  if (!is.null(group)) {
    combined_legend_plot <- create_custom_legends(data, cat_c, group, cat_c_colors, group_colors, var_positions, num_vars, dot_size)

    
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
        # Removed preserve = "aspect"
      )
  } else {
    combined_plot <- p
  }
  
  # Dynamic Plot Sizing and Saving
  n_cat_a <- length(unique(plot_data[[cat_a]]))
  n_cat_b <- length(unique(plot_data[[cat_b]]))
  total_width <- max(n_cat_a * base_width_per_cat_a + ifelse(!is.null(group), 3, 6), 4)
  total_height <- max(n_cat_b * base_height_per_cat_b + 3, 4)  
  
  return(combined_plot)
}