utils::globalVariables(c(
  "unique_groups", "x_min", "x_max", "y_min", "y_max","var",
  "x_pos", "y_pos", "x_offset", "y_offset", "x", "y",
  "present", "combined", "value", "count", "gene", "Contrast",
  "Celltype", "avg_log2FC", "p_val_adj", "log_p_val_adj"
))

default_cat_c_colors <- c(
  "Var1" = "#d5cccd",
  "Var2" = "#cb9992",
  "Var3" = "#ad310f",
  "Var4" = "#7e2a20",
  "Var5" = "#FFD700",
  "Var6" = "#FF6622"
)

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
#' @param cat_c_colors A named vector of colors for `cat_c` categories. Defaults to `NULL` (automatically generated).
#' @param group_colors A named vector of colors for the group variable. Defaults to `NULL` (automatically generated).
#' @param custom_theme A ggplot2 theme for customizing the plot's appearance. Defaults to `theme_minimal()`.
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
                      custom_theme = theme_minimal()) {
  
  num_vars <- length(unique(data[[cat_c]]))
  
  # Validate number of variables
  if (num_vars < 1 || num_vars > 9) {
    stop("Unsupported number of categories for cat_c. Must be between 1 and 9.")
  }
  
  # Set default cat_c_colors using ColorBrewer if not provided
  if (is.null(cat_c_colors)) {
    available_colors <- RColorBrewer::brewer.pal.info
    # Choose a palette that can handle the number of categories
    suitable_palettes <- rownames(available_colors[available_colors$category == "qual" & available_colors$maxcolors >= num_vars, ])
    if (length(suitable_palettes) == 0) {
      stop("No suitable ColorBrewer palette found for the number of categories in cat_c.")
    }
    palette_name <- suitable_palettes[1]  # Select the first suitable palette
    default_cat_c_colors <- RColorBrewer::brewer.pal(n = num_vars, name = palette_name)
    names(default_cat_c_colors) <- unique(data[[cat_c]])
    cat_c_colors <- default_cat_c_colors
  }
  
  # Handle group-related logic only if group is provided
  if (!is.null(group)) {
    # Assign default group colors using ColorBrewer if not provided
    if (is.null(group_colors)) {
      unique_groups <- unique(data[[group]])
      num_groups <- length(unique_groups)
      if (num_groups > 9) {
        stop("The number of groups exceeds the maximum colors available in ColorBrewer palettes.")
      }
      # Choose a suitable palette
      available_colors <- RColorBrewer::brewer.pal.info
      suitable_palettes <- rownames(available_colors[available_colors$category == "qual" & available_colors$maxcolors >= num_groups, ])
      if (length(suitable_palettes) == 0) {
        stop("No suitable ColorBrewer palette found for the number of groups.")
      }
      palette_name <- suitable_palettes[2]  # Select a different palette to avoid color overlap
      group_colors_palette <- brewer.pal(n = num_groups, name = palette_name)
      names(group_colors_palette) <- unique_groups
      group_colors <- group_colors_palette
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
  
  # Perform hierarchical clustering
  cat_a_order <- perform_clustering(data, cat_a, cat_b, cat_c)
  
  # Order cat_b based on group and frequency if group is provided
  if (!is.null(group)) {
    cat_b_order <- order_cat_b(data, group, cat_b, group_colors)
  } else {
    cat_b_order <- levels(data[[cat_b]])
  }
  
  # Prepare plot data
  plot_data <- prepare_plot_data(data, cat_a, cat_b, cat_c, group, var_positions, cat_a_order, cat_b_order)
  
  # Prepare box data
  if (!is.null(group)) {
    box_data <- prepare_box_data(data, cat_a, cat_b, group, cat_a_order, cat_b_order)
  }
  
  # Calculate dynamic dot size
  dot_size <- calculate_dot_size(num_vars)
  
  # Start building the main plot
  p <- ggplot()
  
  # Add group-related layers if group is provided
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
    
    # Define legend dimensions
    legend_width <- 0.25  # Adjust as needed
    legend_height <- 0.5  # Adjust based on the number of legend items
    
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
  base_width_per_cat_a <- 0.5  # inches
  base_height_per_cat_b <- 0.3 # inches
  total_width <- max(n_cat_a * base_width_per_cat_a + ifelse(!is.null(group), 3, 6), 4)
  total_height <- max(n_cat_b * base_height_per_cat_b + 3, 4)  
  
  return(combined_plot)
}