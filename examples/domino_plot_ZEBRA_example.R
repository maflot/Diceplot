
library(dplyr)
library(tidyr)
library(ggplot2)
library(diceplot)

zebra.df = read.csv(file = "data/ZEBRA_sex_degs_set.csv")

genes = c("SPP1","APOE","SERPINA1","PINK1","ANGPT1","ANGPT2","APP","CLU","ABCA7")
zebra.df <- zebra.df %>% filter(gene %in% genes) %>%
  filter(contrast %in% c("MS-CT","AD-CT","ASD-CT","FTD-CT","HD-CT")) %>%
  mutate(cell_type = factor(cell_type, levels = sort(unique(cell_type)))) %>%
  filter(PValue < 0.05)


p <- domino_plot(
  data = zebra.df,
  gene_list = genes,
  var_id = "contrast",
  x = "gene",
  y = "cell_type",
  legend_text_size = 12,
  x_axis_text_size = 12,
  y_axis_text_size = 12,
  contrast = "sex",
  log_fc = "logFC",
  p_val = "FDR",
  logfc_limits = c(min(zebra.df$logFC)-1, max(zebra.df$logFC)-1),
  output_file = "ZEBRA_example1.pdf"
)
print(p)


# First create the base plot without the legend
p <- domino_plot(
  data = zebra.df,
  gene_list = genes,
  var_id = "contrast",
  x = "gene",
  min_dot_size = 0.1,
  max_dot_size = 2,
  legend_text_size = 12,
  x_axis_text_size = 12,
  y_axis_text_size = 12,
  y = "cell_type",
  contrast = "sex",
  log_fc = "logFC",
  p_val = "FDR",
  logfc_limits = c(min(zebra.df$logFC)-2, max(zebra.df$logFC)-2),
)
print(p$domino_plot)

p <- domino_plot(
  data = zebra.df,
  gene_list = genes,
  var_id = "contrast",
  x = "gene",
  min_dot_size = 0.1,
  max_dot_size = 2,
  legend_text_size = 12,
  x_axis_text_size = 12,
  y_axis_text_size = 12,
  y = "cell_type",
  contrast = "sex",
  log_fc = "logFC",
  size_scale_name = "Other scale name",
  color_scale_name = "Other color name",
  p_val = "FDR",
  logfc_limits = c(min(zebra.df$logFC)-2, max(zebra.df$logFC)-2),
)
print(p)
