utils::globalVariables(c(
  "label_x", "x_min", "x_max", "y_min", "y_max","var",
  "x_pos", "y_pos", "x_offset", "y_offset", "x", "y",
  "present", "combined", "value", "count", "gene", "Contrast",
  "Celltype", "avg_log2FC", "p_val_adj", "log_p_val_adj","adj_logfc","label_y","label_x"
))


#' @title Create Variable Positions
#' @description
#' Generates a data frame containing variable names from `cat_c_colors` and corresponding x and y offsets based on the number of variables.
#' @param cat_c_colors A named vector of colors for variables in category C. The names correspond to variable names.
#' @param num_vars The number of variables. Supported values are "3", "4", "5", or "6".
#' @return A data frame with columns:
#' \describe{
#'   \item{var}{Factor of variable names from `cat_c_colors`.}
#'   \item{x_offset}{Numeric x-axis offset for plotting.}
#'   \item{y_offset}{Numeric y-axis offset for plotting.}
#' }
#' @examples
#' library(dplyr)
#' cat_c_colors <- c("Var1" = "red", "Var2" = "blue", "Var3" = "green")
#' create_var_positions(cat_c_colors, 3)
#' @importFrom dplyr %>% mutate group_by summarise n n_distinct arrange desc left_join distinct pull
#' @importFrom tidyr unite pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom data.table :=
#' @export
create_var_positions <- function(cat_c_colors, num_vars) {
  num_vars <- as.character(num_vars)
  switch(num_vars,
         "6" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(-0.2, 0.2, -0.2, 0.2, -0.2, 0.2),
           y_offset = c(0.2, 0.2, 0, 0, -0.2, -0.2)
         ),
         "5" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(0, -0.2, 0.2, -0.2, 0.2),
           y_offset = c(0, 0.2, 0.2, -0.2, -0.2)
         ),
         "4" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(-0.2, 0.2, -0.2, 0.2),
           y_offset = c(0.2, 0.2, -0.2, -0.2)
         ),
         "3" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(0, -0.2, 0.2),
           y_offset = c(0, 0.2, -0.2)
         ),
         "2" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(-0.2, 0.2),
           y_offset = c(0.2, -0.2)
         ),
         "1" = data.frame(
           var = factor(names(cat_c_colors), levels = names(cat_c_colors)),
           x_offset = c(0),
           y_offset = c(0)
         ),
         stop("Unsupported number of variables for variable positions.")
  )
}

#' @title Perform Hierarchical Clustering on Category A
#' @description
#' Performs hierarchical clustering on category A based on the binary presence of combinations of categories B and C.
#' @param data A data frame containing the variables.
#' @param cat_a The name of the column representing category A.
#' @param cat_b The name of the column representing category B.
#' @param cat_c The name of the column representing category C.
#' @return A vector of category A labels ordered according to the hierarchical clustering.
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' data <- data.frame(
#'   cat_a = rep(letters[1:5], each = 4),
#'   cat_b = rep(LETTERS[1:2], times = 10),
#'   cat_c = sample(c("Var1", "Var2", "Var3"), 20, replace = TRUE)
#' )
#' perform_clustering(data, "cat_a", "cat_b", "cat_c")
#' @export
perform_clustering <- function(data, cat_a, cat_b, cat_c) {
  binary_matrix <- data %>%
    mutate(present = 1) %>%
    group_by(!!sym(cat_a), !!sym(cat_b), !!sym(cat_c)) %>%
    summarise(value = sum(present), .groups = "drop") %>%
    unite("combined", c(!!sym(cat_b), !!sym(cat_c)), sep = "_") %>%
    pivot_wider(id_cols = !!sym(cat_a), 
                names_from = combined, 
                values_from = value, 
                values_fill = 0) %>%
    column_to_rownames(cat_a) %>%
    as.matrix()
  
  cat_a_dist <- dist(binary_matrix, method = "binary")
  cat_a_hclust <- hclust(cat_a_dist, method = "ward.D2")
  cat_a_hclust$labels[cat_a_hclust$order]
}

#' @title Order Category B
#' @description
#' Determines the ordering of category B based on the counts within each group, ordered by group and count.
#' @param data A data frame containing the variables.
#' @param group The name of the column representing the grouping variable.
#' @param cat_b The name of the column representing category B.
#' @param group_colors A named vector of colors for each group. The names correspond to group names.
#' @param reverse_order Reverse the ordering? Default is FALSE.
#' @return A vector of category B labels ordered according to group and count.
#' @examples
#' library(dplyr)
#' data <- data.frame(
#'   group = rep(c("G1", "G2"), each = 5),
#'   cat_b = sample(LETTERS[1:3], 10, replace = TRUE)
#' )
#' group_colors <- c("G1" = "red", "G2" = "blue")
#' order_cat_b(data, "group", "cat_b", group_colors)
#' @importFrom dplyr %>% mutate group_by summarise n n_distinct arrange desc left_join distinct pull
#' @importFrom tidyr unite pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom data.table :=
#' @importFrom rlang sym
#' @export
order_cat_b <- function(data, group, cat_b, group_colors, reverse_order = FALSE) {
  cat_b_order <- data %>%
    mutate(!!sym(group) := factor(!!sym(group), levels = rev(names(group_colors)))) %>%  # Reverse to match legacy code
    group_by(!!sym(group), !!sym(cat_b)) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(!!sym(group), desc(count), !!sym(cat_b)) %>%
    pull(!!sym(cat_b)) %>%
    unique()
  if (reverse_order) {
    cat_b_order <- rev(cat_b_order)
  }
  
  return(cat_b_order)
}

#' @title Prepare Plot Data
#' @description
#' Prepares data for plotting by calculating positions based on provided variable positions and orders.
#' @param data A data frame containing the variables.
#' @param cat_a The name of the column representing category A.
#' @param cat_b The name of the column representing category B.
#' @param cat_c The name of the column representing category C.
#' @param group The name of the column representing the grouping variable.
#' @param var_positions A data frame with variable positions, typically output from `create_var_positions`.
#' @param cat_a_order A vector specifying the order of category A.
#' @param cat_b_order A vector specifying the order of category B.
#' @return A data frame ready for plotting with added x_pos and y_pos columns.
#' @examples
#' library(dplyr)
#' data <- data.frame(
#'   cat_a = rep(letters[1:3], each = 4),
#'   cat_b = rep(LETTERS[1:2], times = 6),
#'   cat_c = rep(c("Var1", "Var2"), times = 6),
#'   group = rep(c("G1", "G2"), times = 6)
#' )
#' var_positions <- data.frame(
#'   var = c("Var1", "Var2"),
#'   x_offset = c(0.1, -0.1),
#'   y_offset = c(0.1, -0.1)
#' )
#' cat_a_order <- c("a", "b", "c")
#' cat_b_order <- c("A", "B")
#' prepare_plot_data(data, "cat_a", "cat_b", "cat_c", "group", var_positions, cat_a_order, cat_b_order)
#' @importFrom dplyr %>% mutate group_by summarise n n_distinct arrange desc left_join distinct pull
#' @importFrom tidyr unite pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom data.table :=
#' @importFrom rlang sym
#' @importFrom stats setNames
#' @export
prepare_plot_data <- function(data, cat_a, cat_b, cat_c, group, var_positions, cat_a_order, cat_b_order) {
  plot_data <- data %>%
    left_join(var_positions, by = setNames("var", cat_c)) %>%
    mutate(
      !!sym(cat_a) := factor(!!sym(cat_a), levels = cat_a_order),
      !!sym(cat_b) := factor(!!sym(cat_b), levels = cat_b_order),
      x_pos = as.numeric(!!sym(cat_a)) + x_offset,
      y_pos = as.numeric(!!sym(cat_b)) + y_offset
    ) %>%
    arrange(!!sym(cat_a), !!sym(group), !!sym(cat_b))
  return(plot_data)
}

#' @title Prepare Box Data
#' @description
#' Prepares data for plotting boxes by calculating box boundaries based on category positions.
#' @param data A data frame containing the variables.
#' @param cat_a The name of the column representing category A.
#' @param cat_b The name of the column representing category B.
#' @param group The name of the column representing the grouping variable.
#' @param cat_a_order A vector specifying the order of category A.
#' @param cat_b_order A vector specifying the order of category B.
#' @return A data frame with box boundaries for plotting.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @examples
#' library(dplyr)
#' data <- data.frame(
#'   cat_a = rep(letters[1:3], each = 2),
#'   cat_b = rep(LETTERS[1:2], times = 3),
#'   group = rep(c("G1", "G2"), times = 3)
#' )
#' cat_a_order <- c("a", "b", "c")
#' cat_b_order <- c("A", "B")
#' prepare_box_data(data, "cat_a", "cat_b", "group", cat_a_order, cat_b_order)
#' @importFrom dplyr %>% mutate group_by summarise n n_distinct arrange desc left_join distinct pull
#' @importFrom tidyr unite pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom data.table :=
#' @importFrom rlang sym
#' @export
prepare_box_data <- function(data, cat_a, cat_b, group, cat_a_order, cat_b_order) {
  box_data <- data %>%
    mutate(
      !!sym(cat_a) := factor(!!sym(cat_a), levels = cat_a_order),
      !!sym(cat_b) := factor(!!sym(cat_b), levels = cat_b_order)
    ) %>%
    distinct(!!sym(cat_a), !!sym(cat_b), !!sym(group)) %>%
    mutate(
      x_min = as.numeric(!!sym(cat_a)) - 0.4,
      x_max = as.numeric(!!sym(cat_a)) + 0.4,
      y_min = as.numeric(!!sym(cat_b)) - 0.4,
      y_max = as.numeric(!!sym(cat_b)) + 0.4
    ) %>%
    arrange(!!sym(cat_a), !!sym(group), !!sym(cat_b))
  return(box_data)
}


#' @title Calculate Dynamic Dot Size
#' @description
#' Calculates the dot size based on the number of variables.
#' @param num_vars Number of variables.
#' @param max_size Maximal dot size for the plot to scale the dot sizes.
#' @param min_size Minimal dot size for the plot to scale the dot sizes.
#' @return A numeric value representing the dot size.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @export
calculate_dot_size <- function(num_vars, max_size, min_size) {
  size <- max(min_size, max_size - (num_vars - 1))
  return(size)
}


#' @title Create Custom Legends
#' @description
#' Creates custom legend plots for `cat_c` and `group`.
#' @param data The original data frame.
#' @param cat_c The name of the `cat_c` variable.
#' @param group The name of the group variable.
#' @param cat_c_colors A named vector of colors for `cat_c`.
#' @param group_colors A named vector of colors for the group variable.
#' @param var_positions Data frame with variable positions.
#' @param num_vars Number of variables in `cat_c`.
#' @param dot_size The size of the dots used in the plot.
#' @return A combined ggplot object of the custom legends.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_manual scale_fill_manual theme_void theme element_text element_blank margin coord_fixed geom_text ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom stats dist hclust
#' @importFrom utils globalVariables
#' @importFrom rlang sym
#' @export
create_custom_legends <- function(data, cat_c, group, cat_c_colors, group_colors, var_positions, num_vars, dot_size) {
  # Create legend_data using var_positions
  legend_data <- var_positions %>%
    mutate(
      x = x_offset + 1,
      y = y_offset + 1
    )
  
  # Adjust label positions using legend_data
  label_offsets <- legend_data %>%
    mutate(
      label_x = x + ifelse(x_offset == 0, 0.25, x_offset * 1.5),  # Adjust multiplier as needed
      label_y = y  # You can adjust y position if needed
    )
  
  # Create the custom legend plot for cat_c
  custom_legend_plot <- ggplot() +
    geom_point(data = legend_data, aes(x = x, y = y, color = var), size = dot_size) +
    geom_point(data = legend_data, aes(x = x, y = y), size = dot_size + 0.5, shape = 1, color = "black") +
    scale_color_manual(values = cat_c_colors, name = cat_c) +
    theme_void() +
    theme(
      legend.position = "none", 
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1  # Enforce square aspect ratio
    ) +
    coord_fixed(ratio = 1, xlim = c(0.5, 2.5), ylim = c(0.5, 1.5), expand = FALSE) +
    geom_text(data = label_offsets, 
              aes(x = label_x, y = label_y, label = var),
              size = 3, 
              color = "black",
              hjust = 0,
              vjust = 0.5) +
    ggtitle("Dice arrangement")
  
  
  # Compute coordinate ranges
  ylim_min <- 0.5
  ylim_max <- length(group_colors) + 0.5
  ylim_range <- ylim_max - ylim_min
  
  xlim_min <- 0.5
  xlim_max <- 1.5  # Keep x-axis narrow to help maintain aspect ratio
  xlim_range <- xlim_max - xlim_min
  
  # Compute aspect ratio
  aspect_ratio <- ylim_range / xlim_range
  
  # Create legend data for group
  legend_data_group <- data.frame(
    group = factor(names(group_colors), levels = names(group_colors)),
    x = 1,
    y = seq(length(group_colors), 1)
  )
  
  # Create the custom legend plot for group
  group_legend_plot <- ggplot() +
    geom_rect(
      data = legend_data_group,
      aes(
        xmin = x - 0.3,
        xmax = x + 0.3,
        ymin = y - 0.3,
        ymax = y + 0.3,
        fill = group
      ),
      color = "grey",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    scale_fill_manual(values = group_colors, name = group) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 50, 5, 5),  # Increase right margin for labels
      aspect.ratio = aspect_ratio  # Set computed aspect ratio
    ) +
    coord_fixed(
      ratio = 1,  # Keep units equal on x and y axes
      xlim = c(xlim_min, xlim_max),
      ylim = c(ylim_min, ylim_max),
      expand = FALSE,
      clip = "off"  # Allow labels to be drawn outside the plotting area
    ) +
    geom_text(
      data = legend_data_group,
      aes(x = x + 0.4, y = y, label = group),
      size = 3,
      color = "black",
      hjust = 0
    )
  # Combine the legend plots vertically
  combined_legend_plot <- cowplot::plot_grid(
    custom_legend_plot, 
    group_legend_plot, 
    ncol = 1, 
    align = 'v', 
    rel_heights = c(2, 1)
  )
  
  return(combined_legend_plot)
}

#' @title Create Custom Domino Legends
#' @description
#' Creates custom legend plots for domino plot, including variable positions within each domino, contrast layout,
#' and explanations for log fold change and p-values.
#' @param contrast_levels A character vector with the two contrast levels.
#' @param var_positions A data frame containing variable positions data.
#' @param var_id A string specifying the variable identifier column name.
#' @param contrast A string specifying the contrast column name.
#' @param logfc_colors A named vector of colors for log fold change.
#' @param logfc_limits A numeric vector of length 2 specifying the limits for the log fold change color scale.
#' @param color_scale_name A string specifying the name of the color scale in the legend.
#' @param size_scale_name A string specifying the name of the size scale in the legend.
#' @param min_dot_size A numeric value indicating the minimum dot size in the plot.
#' @param max_dot_size A numeric value indicating the maximum dot size in the plot.
#' @return A combined ggplot object of the custom legends.
#' @importFrom ggplot2 ggplot aes geom_point geom_rect scale_color_gradient2 scale_size_continuous theme_void theme element_text element_blank margin coord_fixed geom_text element_rect ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom utils globalVariables
#' @export
create_custom_domino_legends <- function(
    contrast_levels,
    var_positions,
    var_id,
    contrast,
    logfc_colors,
    logfc_limits,
    color_scale_name,
    size_scale_name,
    min_dot_size,
    max_dot_size
) {
  # Create a legend showing variable positions inside the domino boxes
  # Split var_positions by contrast
  var_positions_left <- var_positions[var_positions[[contrast]] == contrast_levels[1],]
  var_positions_right <- var_positions[var_positions[[contrast]] == contrast_levels[2],]
  
  # Create a more comprehensive visualization of the domino layout
  domino_data <- data.frame(
    contrast = c("Left", "Right"),
    contrast_label = contrast_levels,
    x = c(1, 2),
    y = c(1, 1)
  )
  
  # Create the var positions legend
  var_legend_data <- data.frame()
  
  # For each var in the left contrast
  for (i in 1:nrow(var_positions_left)) {
    var_entry <- var_positions_left[i,]
    var_legend_data <- rbind(var_legend_data, data.frame(
      var = var_entry[[var_id]],
      x = 1 + var_entry$x_offset - 1, # Adjust so it's centered around x=1
      y = 1 + var_entry$y_offset,
      position = "Left",
      text_x_offset = -0.6, # Text will appear to the left
      hjust_val = 1         # Right-aligned text
    ))
  }
  
  # For each var in the right contrast
  for (i in 1:nrow(var_positions_right)) {
    var_entry <- var_positions_right[i,]
    var_legend_data <- rbind(var_legend_data, data.frame(
      var = var_entry[[var_id]],
      x = 2 + var_entry$x_offset - 2, # Adjust so it's centered around x=2
      y = 1 + var_entry$y_offset,
      position = "Right",
      text_x_offset = 0.6,  # Text will appear to the right
      hjust_val = 0         # Left-aligned text
    ))
  }
  
  # Create the position and variable layout legend plot
  vars_legend_plot <- ggplot() +
    # Draw the domino boxes
    geom_rect(
      data = domino_data,
      aes(
        xmin = x - 0.4,
        xmax = x + 0.4,
        ymin = y - 0.4,
        ymax = y + 0.4
      ),
      fill = "white", color = "darkgrey", alpha = 0.5
    ) +
    # Add dots showing variable positions within each box
    geom_point(
      data = var_legend_data,
      aes(x = x, y = y),
      size = 3,
      color = "darkblue"
    ) +
    # Add variable labels with position-dependent placement
    geom_text(
      data = var_legend_data,
      aes(x = x + text_x_offset, y = y, label = var, hjust = hjust_val),
      size = 3
    ) +
    # Add contrast labels inside boxes
    geom_text(
      data = domino_data,
      aes(x = x, y = y - 0.2, label = contrast_label),
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    # Add position labels above boxes
    geom_text(
      data = domino_data,
      aes(x = x, y = y + 0.7, label = contrast),
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(10, 20, 10, 20),  # Increased margins
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-0.5, 3.5),  # Wider limits to accommodate the text
      ylim = c(0.2, 2.0),
      expand = FALSE,
      clip = "off"  # Allow text outside the plot area
    ) +
    ggtitle("Variable Positions")
  
  # Create data for log fold change legend
  logfc_range <- seq(from = logfc_limits[1], to = logfc_limits[2], length.out = 5)
  logfc_legend_data <- data.frame(
    logfc = logfc_range,
    x = rep(1, length(logfc_range)),
    y = seq(length(logfc_range), 1)
  )
  
  # Create log fold change legend plot
  logfc_legend_plot <- ggplot() +
    geom_point(
      data = logfc_legend_data,
      aes(x = x, y = y, color = logfc),
      size = 3
    ) +
    geom_point(
      data = logfc_legend_data,
      aes(x = x, y = y),
      size = 3.5,
      shape = 1,
      color = "black"
    ) +
    scale_color_gradient2(
      name = color_scale_name,
      low = logfc_colors["low"],
      mid = logfc_colors["mid"],
      high = logfc_colors["high"],
      limits = logfc_limits
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 50, 5, 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(0.5, 1.5),
      ylim = c(0.5, length(logfc_range) + 0.5),
      expand = FALSE,
      clip = "off"
    ) +
    geom_text(
      data = logfc_legend_data,
      aes(x = x + 0.4, y = y, label = sprintf("%.1f", logfc)),
      size = 3,
      color = "black",
      hjust = 0
    ) +
    ggtitle(color_scale_name)
  
  # Create data for p-value legend
  p_val_range <- seq(from = 0, to = 5, length.out = 5)  # -log10(p) from 0 to 5
  p_val_legend_data <- data.frame(
    log_p_val = p_val_range,
    x = rep(1, length(p_val_range)),
    y = seq(length(p_val_range), 1)
  )
  
  # Create p-value legend plot
  p_val_legend_plot <- ggplot() +
    geom_point(
      data = p_val_legend_data,
      aes(x = x, y = y, size = log_p_val),
      color = "black"
    ) +
    geom_point(
      data = p_val_legend_data,
      aes(x = x, y = y, size = log_p_val),
      shape = 1,
      color = "black"
    ) +
    scale_size_continuous(
      name = size_scale_name,
      range = c(min_dot_size, max_dot_size)
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 50, 5, 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(0.5, 1.5),
      ylim = c(0.5, length(p_val_range) + 0.5),
      expand = FALSE,
      clip = "off"
    ) +
    geom_text(
      data = p_val_legend_data,
      aes(x = x + 0.4, y = y, label = sprintf("%.1f", log_p_val)),
      size = 3,
      color = "black",
      hjust = 0
    ) +
    ggtitle(size_scale_name)
  
  # Combine all legend plots
  combined_legend_plot <- cowplot::plot_grid(
    vars_legend_plot,
    logfc_legend_plot,
    p_val_legend_plot,
    ncol = 1,
    align = 'v',
    rel_heights = c(2, 1, 1)  # Give more height to the variable positions legend
  )
  
  return(combined_legend_plot)
}
save_diceplot = function(){
  
}





