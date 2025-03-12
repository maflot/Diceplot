

# Change-log v0.1.5

## domino_plot function
- Adopted cleaner naming conventions for parameters:
  - `feature_col` → `x`
  - `celltype_col` → `y`
  - `contrast_col` → `contrast`
  - `logfc_col` → `log_fc`
  - `pval_col` → `p_val`
- Added deprecation warnings for old parameter names that will be removed in a future version >v1.5
- Updated function documentation to reflect new parameter names

## dice_plot function
- Adopted cleaner naming conventions for parameters:
  - `cat_a` → `x`
  - `cat_b` → `y`
  - `cat_c` → `z`
  - `cat_c_colors` → `z_colors`
  - `base_width_per_cat_a` → `base_width_per_x`
  - `base_height_per_cat_b` → `base_height_per_y`
- Added deprecation warnings for old parameter names that will be removed in a future version >v1.5
- Updated function documentation to reflect new parameter names
- Improved handling of factors in both functions
- `cluster_by_row` argument defaults to TRUE, set to FALSE to use factor levels for ordering 
- `cluster_by_column` argument defaults to TRUE, set to FALSE to use factor levels for ordering 
- `show_legend` defaults to TRUE, controls whether to show or omit the legend plot
- `cat_b_order` argument removed, will throw an error in a future version

# Change-log v0.1.4
- handling factors in diceplot
- ```cluster_by_row``` argument defaults to TRUE, false it will use the factor levels for ordering 
- ```cluster_by_column``` argument defaults to TRUE, false it will use the factor levels for ordering 
- ```show_legend``` defaults to TRUE, show or omit the legend plot s
- ```cat_b_order``` argument removed, will throw an error in a future version

# Change-log v0.1.3
- add citation 
- ```cat_c_colors``` allow to pass a colorbrewer palette name for chosing the color palette as an additional option
- ```group_colors``` allow to pass a colorbrewer palette name for chosing the color palette as an additional option
- ```cat_b_order``` allow to pass an explicit row order

# Change-log v0.1.2
  - rename files from dice_plot/ domino_plot to diceplot/dominoplot
  - adapted dice_plot arguments:
    - ```group``` variable is not needed anymore for running the code
    - ```group_color``` is set automatically using RColorBrewer
    - dice_plot is now returning the ggplot object and removed the option to save the plot
      - ```plot_path``` removed
      - ```format```    removed
      - ```output_string``` removed
    - ```max_dot_size``` and ```min_dot_size``` as additional arguments to modify your plot.
    - ```legend_width``` and ```legend_height``` as additional arguments to modify your plots legend.
    - ```base_width_per_cat_a``` and ```base_height_per_cat_b``` as additional arguments to modify your plots legend.
    - ```reverse_ordering``` reverse the clustering order if wanted



