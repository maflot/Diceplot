utils::globalVariables(c(
  "label_x", "x_min", "x_max", "y_min", "y_max","var",
  "x_pos", "y_pos", "x_offset", "y_offset", "x", "y",
  "present", "combined", "value", "count", "gene", "Contrast","z", "label", "log_p_val",
  "Celltype", "avg_log2FC", "p_val_adj", "log_p_val_adj","adj_logfc","label_y","label_x","dot_color"
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
  if (is.null(group)) {
    # If group is NULL, order by count only
    cat_b_order <- data %>%
      group_by(!!sym(cat_b)) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(desc(count), !!sym(cat_b)) %>%
      pull(!!sym(cat_b)) %>%
      unique()
  } else {
    # Original logic using group
    cat_b_order <- data %>%
      mutate(!!sym(group) := factor(!!sym(group), levels = rev(names(group_colors)))) %>%
      group_by(!!sym(group), !!sym(cat_b)) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(!!sym(group), desc(count), !!sym(cat_b)) %>%
      pull(!!sym(cat_b)) %>%
      unique()
  }
  
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
    )
  
  # Arrange differently based on whether group is NULL
  if (!is.null(group)) {
    plot_data <- plot_data %>%
      arrange(!!sym(cat_a), !!sym(group), !!sym(cat_b))
  } else {
    plot_data <- plot_data %>%
      arrange(!!sym(cat_a), !!sym(cat_b))
  }
  
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
  # This function should only be called when group is not NULL, but adding a check for safety
  if (is.null(group)) {
    stop("prepare_box_data function cannot be called with group = NULL")
  }
  
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

#' @title Prepare Simple Box Data (no grouping)
#' @description
#' Prepares data for plotting boxes without grouping by calculating box boundaries based on category positions.
#' @param data A data frame containing the variables.
#' @param cat_a The name of the column representing category A.
#' @param cat_b The name of the column representing category B.
#' @param cat_a_order A vector specifying the order of category A.
#' @param cat_b_order A vector specifying the order of category B.
#' @return A data frame with box boundaries for plotting.
#' @importFrom dplyr %>% mutate distinct arrange
#' @importFrom rlang sym
#' @export
prepare_simple_box_data <- function(data, cat_a, cat_b, cat_a_order, cat_b_order) {
  box_data <- data %>%
    mutate(
      !!sym(cat_a) := factor(!!sym(cat_a), levels = cat_a_order),
      !!sym(cat_b) := factor(!!sym(cat_b), levels = cat_b_order)
    ) %>%
    distinct(!!sym(cat_a), !!sym(cat_b)) %>%
    mutate(
      x_min = as.numeric(!!sym(cat_a)) - 0.4,
      x_max = as.numeric(!!sym(cat_a)) + 0.4,
      y_min = as.numeric(!!sym(cat_b)) - 0.4,
      y_max = as.numeric(!!sym(cat_b)) + 0.4
    ) %>%
    arrange(!!sym(cat_a), !!sym(cat_b))
  
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
#' @importFrom ggrepel geom_text_repel
#' @export
create_custom_legends <- function(data, cat_c, group, cat_c_colors, group_colors, var_positions, num_vars, dot_size) {
  # Create legend_data using var_positions
  legend_data <- var_positions %>%
    mutate(
      x = x_offset + 1,
      y = y_offset + 1
    )
  
  custom_legend_plot <- ggplot() +
    geom_point(data = legend_data, aes(x = x, y = y, color = var), size = dot_size) +
    geom_point(data = legend_data, aes(x = x, y = y), size = dot_size + 0.5, shape = 1, color = "black") +
    # Replace geom_text with geom_text_repel
    geom_text_repel(
      data = legend_data,
      aes(x = x, y = y, label = var),
      size = 3.5,
      segment.size = 0.2,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 1,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = cat_c_colors, name = cat_c) +
    theme_void() +
    ggtitle("Dice arrangement")+
    theme(
      legend.position = "none", 
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
  
  
 
  # Create the custom legend plot for group
  if (!is.null(group)){
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
    combined_legend_plot <- cowplot::plot_grid(
      custom_legend_plot, 
      group_legend_plot, 
      ncol = 1, 
      align = 'v', 
      rel_heights = c(2, 1)
    )
    
    return(combined_legend_plot)
  }
  # Combine the legend plots vertically
  print(custom_legend_plot)
  return(custom_legend_plot)
  
}



#' Create custom legends for the domino plot
#'
#' @param contrast_levels A character vector of contrast level names.
#' @param var_positions A data frame containing variable positions.
#' @param var_id A string representing the column name for the variable identifier.
#' @param contrast A string representing the column name for the contrast variable.
#' @param logfc_colors A named vector specifying the colors for the log fold change scale.
#' @param logfc_limits A numeric vector of length 2 specifying the limits for the log fold change color scale.
#' @param color_scale_name A string specifying the name of the color scale in the legend.
#' @param size_scale_name A string specifying the name of the size scale in the legend.
#' @param min_dot_size A numeric value indicating the minimum dot size.
#' @param max_dot_size A numeric value indicating the maximum dot size.
#' @param legend_text_size A numeric value indicating the text size for the legend.
#'
#' @return A ggplot object containing custom legends.
#' @importFrom ggplot2 geom_raster scale_fill_gradient2
#' @keywords internal
create_custom_domino_legends <- function(contrast_levels,
                                         var_positions,
                                         var_id,
                                         contrast,
                                         logfc_colors,
                                         logfc_limits,
                                         color_scale_name,
                                         size_scale_name,
                                         min_dot_size,
                                         max_dot_size,
                                         legend_text_size = 8) {
  # Create a data frame for the color scale legend
  color_legend_data <- data.frame(
    x = rep(1, 5),
    y = 1:5,
    logfc = seq(from = logfc_limits[1],
                to = logfc_limits[2],
                length.out = 5)
  )
  
  # Create a data frame for the size scale legend
  size_legend_data <- data.frame(
    x = rep(1, 5),
    y = 1:5,
    log_p_val = seq(from = 0, to = 5, length.out = 5)
  )
  
  # Get unique variables from the left contrast only (assuming both sides have same variables)
  vars_to_show <- unique(var_positions[var_positions[[contrast]] == contrast_levels[1], var_id])
  
  # Filter var_positions to include only the first contrast
  legend_vars <- var_positions[var_positions[[contrast]] == contrast_levels[1],]
  
  # Create a simplified legend data frame for the positions exactly as in the example
  var_legend_data <- data.frame(
    var = legend_vars[[var_id]],
    x = legend_vars$x_offset,
    y = legend_vars$y_offset
  )
  
  # Create the position legend plot centered - adjust coordinate system to match other legends
  vars_legend_plot <- ggplot() +
    # Add dots showing variable positions using the correct pattern
    geom_point(
      data = var_legend_data,
      aes(x = x, y = y),
      size = 3,
      fill = "white",
      color = "black",
      stroke = 0.8,
      shape = 21
    ) +
    # Add variable labels with proper offset and alignment
    geom_text(
      data = var_legend_data,
      aes(
        # Position text with appropriate offsets
        x = ifelse(x > 0, x + 0.15, 
                   ifelse(x < 0, x - 0.15, 
                          ifelse(x == 0, 0, x))),
        y = ifelse(y > 0, y + 0.15,
                   ifelse(y < 0, y - 0.15, 
                          ifelse(y == 0 & x == 0, 0.15, y))),
        label = var,
        # Adjust text alignment based on position
        hjust = ifelse(x < 0, 1,
                       ifelse(x > 0, 0, 0.5)),
        vjust = ifelse(y > 0, 0,
                       ifelse(y < 0, 1, 
                              ifelse(y == 0 & x == 0, 0, 0.5)))
      ),
      size = legend_text_size/3
    ) +
    # Add title
    ggtitle("Domino Layout") +
    theme_void() +
    theme(
      plot.margin = margin(10, 20, 10, 20),  # Match margins with other legends
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = legend_text_size*1.5, face = "bold", margin = margin(0, 0, 10, 0)),
      # Force plot to use the entire available space
      plot.title.position = "plot"
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-0.75, 0.75),  # Match x-limits with other legends
      ylim = c(-0.5, 0.5),
      expand = FALSE,  # Match with other legends
      clip = "off"
    )
  
  # Create data for log fold change legend - centered
  logfc_range <- seq(from = logfc_limits[1], to = logfc_limits[2], length.out = 5)
  logfc_legend_data <- data.frame(
    logfc = logfc_range,
    x = rep(0, length(logfc_range)),  # Centered at x=0
    y = seq(length(logfc_range), 1)
  )
  
  # Create a gradient data frame for the color bar - centered
  grad_data <- expand.grid(
    x = seq(-0.25, 0.25, length.out = 50),  # Centered around x=0
    y = seq(1, 5, length.out = 50)
  )
  grad_data$z <- rep(seq(logfc_limits[1], logfc_limits[2], length.out = 50), each = 50)
  
  # Create tick positions for the color bar - centered
  tick_positions <- seq(1, 5, length.out = 5)
  tick_labels <- seq(logfc_limits[1], logfc_limits[2], length.out = 5)
  tick_data <- data.frame(
    x = rep(0.35, length(tick_positions)),  # Adjusted for centered bar
    y = tick_positions,
    label = sprintf("%.1f", tick_labels)
  )
  
  logfc_legend_plot <- ggplot() +
    # Add the continuous color bar
    geom_raster(
      data = grad_data,
      aes(x = x, y = y, fill = z)
    ) +
    # Add border around the color bar
    geom_rect(
      aes(xmin = -0.25, xmax = 0.25, ymin = 1, ymax = 5),  # Centered
      fill = NA, color = "black", size = 0.5
    ) +
    # Add scale for the fill
    scale_fill_gradient2(
      name = color_scale_name,
      low = logfc_colors["low"],
      mid = logfc_colors["mid"],
      high = logfc_colors["high"],
      limits = logfc_limits
    ) +
    # Add tick labels
    geom_text(
      data = tick_data,
      aes(x = x, y = y, label = label),
      size = legend_text_size/2.3,
      hjust = 0,
      vjust = 0.5
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 20, 10, 20),  # Equal margins
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = legend_text_size*1.375, face = "bold", margin = margin(0, 0, 10, 0))  # Centered title
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-0.75, 0.75),  # Centered x limits
      ylim = c(0.8, 5.2),
      expand = FALSE,
      clip = "off"
    ) +
    ggtitle(color_scale_name)
  
  # Create data for p-value legend - centered
  p_val_range <- seq(from = 0, to = 5, length.out = 5)  # -log10(p) from 0 to 5
  p_val_legend_data <- data.frame(
    log_p_val = p_val_range,
    x = rep(0, length(p_val_range)),  # Centered at x=0
    y = seq(length(p_val_range), 1)
  )
  
  # Calculate a better size range for dots based on the range of p-values
  size_range <- c(min_dot_size * 1.2, max_dot_size * 0.8)  # Adjusted for better visual scaling
  
  # Create p-value legend plot with centered elements
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
      range = size_range
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 20, 10, 20),  # Equal margins
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = legend_text_size*1.375, face = "bold", margin = margin(0, 0, 10, 0))  # Centered title
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-0.75, 0.75),  # Centered x limits
      ylim = c(0.5, length(p_val_range) + 0.5),
      expand = FALSE,
      clip = "off"
    ) +
    geom_text(
      data = p_val_legend_data,
      aes(x = x + 0.4, y = y, label = sprintf("%.1f", log_p_val)),  # Adjusted position
      size = legend_text_size/2.3,  # Size based on legend_text_size
      color = "black",
      hjust = 0,
      vjust = 0.5  # Center vertically
    ) +
    ggtitle(size_scale_name)
  
  # Combine all legend plots
  combined_legend_plot <- cowplot::plot_grid(
    vars_legend_plot,
    logfc_legend_plot,
    p_val_legend_plot,
    ncol = 1,
    align = 'v',
    axis = 'lr',  # Align plots at their left and right edges
    rel_heights = c(2, 1.5, 1.5)
  )
  
  return(combined_legend_plot)
}



#' Plot Dice Representations on sf Objects
#'
#' @description
#' Creates a ggplot2 layer that places dice representations on spatial features
#' in an sf object. The dice values are determined by a column in the sf object.
#'
#' @param sf_data An sf object containing the spatial features.
#' @param dice_value_col Character. Name of the column in sf_data containing dice values (1-6).
#'   Default is "dice".
#' @param face_color Character vector. Column names in sf_data containing color information
#'   for each dice dot. If NULL (default), all dots are black.
#' @param dice_color Character. Background color of the dice. Default is "white".
#' @param dice_size Numeric. Size of the dice. Default is 3.
#' @param dot_size Numeric. Size of the dots on the dice. If NULL (default), 
#'   it's calculated as 20% of dice_size.
#' @param rectangle_padding Numeric. Padding of the rectangle around the dots, as a 
#'   proportion of dice_size. Default is 0.05.
#' @param ... Additional arguments passed to geom_point for the dots.
#'
#' @return A list of ggplot2 layers (rectangle layer and dots layer).
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(sf)
#' 
#' # Create sample sf data with dice values
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#' nc$dice <- sample(1:6, nrow(nc), replace = TRUE)
#' 
#' # Basic plot with dice
#' ggplot(nc) + 
#'   geom_sf() + 
#'   geom_dice_sf(sf_data = nc)
#'   
#' # Customized dice
#' ggplot(nc) + 
#'   geom_sf() + 
#'   geom_dice_sf(sf_data = nc, dice_color = "lightblue", dice_size = 5)
#' }
#'
#' @importFrom sf st_centroid st_coordinates
#' @importFrom dplyr %>% group_by summarise
#' @importFrom ggplot2 geom_rect geom_point aes
#' @export
geom_dice_sf <- function(sf_data,
                         dice_value_col = "dice",
                         face_color = NULL,
                         dice_color = "white",
                         dice_size = 3,
                         dot_size = NULL,
                         rectangle_padding = 0.05,
                         ...) {
  # Input validation
  if (!inherits(sf_data, "sf")) {
    stop("sf_data must be an sf object.", call. = FALSE)
  }
  
  if (!dice_value_col %in% names(sf_data)) {
    stop(sprintf("Column '%s' not found in sf_data.", dice_value_col), call. = FALSE)
  }
  
  if (!is.null(face_color)) {
    missing_cols <- face_color[!face_color %in% names(sf_data)]
    if (length(missing_cols) > 0) {
      stop(sprintf("Color column(s) not found in sf_data: %s",
                   paste(missing_cols, collapse = ", ")), 
           call. = FALSE)
    }
  }
  
  if (!is.numeric(dice_size) || dice_size <= 0) {
    stop("dice_size must be a positive number.", call. = FALSE)
  }
  
  if (!is.null(dot_size) && (!is.numeric(dot_size) || dot_size <= 0)) {
    stop("dot_size must be a positive number or NULL.", call. = FALSE)
  }
  
  if (!is.numeric(rectangle_padding) || rectangle_padding < 0) {
    stop("rectangle_padding must be a non-negative number.", call. = FALSE)
  }
  
  # Calculate centroids of spatial features
  centroids <- sf::st_centroid(sf_data)
  coords <- sf::st_coordinates(centroids)
  sf_data$centroid_x <- coords[, "X"]
  sf_data$centroid_y <- coords[, "Y"]
  
  # Set scaling parameters
  dot_scale <- if (!is.null(dot_size)) dot_size else dice_size * 0.2
  position_scale <- dice_size * 0.1
  
  # Generate dice dots for each feature
  dice_dots <- do.call(rbind, lapply(seq_len(nrow(sf_data)), function(i) {
    row <- sf_data[i, ]
    n_dots <- as.numeric(as.character(row[[dice_value_col]]))
    
    # Handle custom colors for dots
    if (!is.null(face_color)) {
      if (is.na(n_dots) || n_dots < 1 || n_dots > length(face_color)) {
        warning(sprintf("Dice value %s must be between 1 and %d. Skipping feature %d.",
                        ifelse(is.na(n_dots), "NA", as.character(n_dots)),
                        length(face_color), i),
                call. = FALSE)
        return(NULL)
      }
      palette <- sapply(face_color[seq_len(n_dots)], function(col) row[[col]])
      dot_colors <- palette
    } else {
      if (is.na(n_dots) || n_dots < 1 || n_dots > 6) {
        warning(sprintf("Dice value %s must be between 1 and 6. Skipping feature %d.",
                        ifelse(is.na(n_dots), "NA", as.character(n_dots)), i),
                call. = FALSE)
        return(NULL)
      }
      palette <- rep("black", n_dots)
      dot_colors <- palette
    }
    
    # Create positions for dots based on dice value
    positions <- create_var_positions(palette, n_dots)
    positions$x_offset <- positions$x_offset * position_scale
    positions$y_offset <- positions$y_offset * position_scale
    
    # Center the dots
    center_x <- mean(positions$x_offset)
    center_y <- mean(positions$y_offset)
    positions$x_offset <- positions$x_offset - center_x
    positions$y_offset <- positions$y_offset - center_y
    
    # Calculate final positions
    positions$x <- row$centroid_x + positions$x_offset
    positions$y <- row$centroid_y + positions$y_offset
    positions$feature_id <- i
    positions$dot_color <- dot_colors
    
    positions
  }))
  
  # If no valid dice dots were created, return empty layers
  if (is.null(dice_dots) || nrow(dice_dots) == 0) {
    warning("No valid dice representations could be created.", call. = FALSE)
    return(list(
      ggplot2::geom_blank(),
      ggplot2::geom_blank()
    ))
  }
  
  # Convert to data frame
  dice_dots <- as.data.frame(dice_dots)
  
  # Create rectangles for each dice
  rect_padding <- dice_size * rectangle_padding
  group_rects <- dice_dots %>%
    dplyr::group_by(feature_id) %>%
    dplyr::summarise(
      x_min = min(x) - rect_padding,
      x_max = max(x) + rect_padding,
      y_min = min(y) - rect_padding,
      y_max = max(y) + rect_padding,
      .groups = "drop"
    )
  
  # Create the rectangle layer
  rect_layer <- ggplot2::geom_rect(
    data = group_rects,
    mapping = ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
    fill = dice_color,
    color = "black",
    inherit.aes = FALSE
  )
  
  # Create the dots layer
  dots_layer <- ggplot2::geom_point(
    data = dice_dots,
    mapping = ggplot2::aes(x = x, y = y, color = dot_color),
    size = dot_scale,
    ...
  )
  
  # Return both layers
  list(rect_layer, dots_layer)
}