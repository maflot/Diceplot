library(sf)
library(ggplot2)
library(diceplot)
library(dplyr)
library(cowplot)
library(rnaturalearth)

# Define custom var_positions for a 4-dot dice face
var_positions <- data.frame(
  x_offset = c(-0.3, 0.3, -0.3, 0.3),
  y_offset = c(0.3, 0.3, -0.3, -0.3),
  var = c("log_France", "log_Swiss", "log_Luxembourg", "log_Rheinlandpfalz")
)

germany_states <- ne_states(country = "Germany", returnclass = "sf")
saarland <- germany_states[germany_states$name == "Saarland", ]

cities <- data.frame(
  name = c("Saarbrücken", "Saarlouis", "Homburg", "Britten", "Merzig", "Lebach", "Ottweiler"),
  dice = 4,
  lon = c(6.996, 6.751, 7.339, 6.784, 6.639, 6.913, 7.167),
  lat = c(49.234, 49.315, 49.320, 49.481, 49.442, 49.407, 49.400),
  France = c(14, 12, 38, 27, 18, 27, 36),
  Swiss = c(190, 204, 195, 221, 220, 210, 206),
  Luxembourg = c(51, 31, 67, 23, 17, 35, 52),
  Rheinlandpfalz = c(29, 27, 6, 16, 20, 12, 16)
)

cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
cities_sf$log_France <- log(cities_sf$France)
cities_sf$log_Swiss <- log(cities_sf$Swiss)
cities_sf$log_Luxembourg <- log(cities_sf$Luxembourg)
cities_sf$log_Rheinlandpfalz <- log(cities_sf$Rheinlandpfalz)

create_custom_legends_for_map <- function(var_positions, dot_size, legend_text_size = 9) {
  legend_data <- var_positions %>% mutate(
    x = x_offset + 1,
    y = y_offset + 1
  )
  ggplot() +
    geom_point(data = legend_data, aes(x = x, y = y), size = dot_size, color = "black") +
    geom_point(data = legend_data, aes(x = x, y = y), size = dot_size + 0.5, shape = 1, color = "black") +
    coord_fixed(ratio = 1, xlim = c(0.5, 2.5), ylim = c(0.5, 1.5), expand = FALSE) +
    geom_text(
      data = legend_data,
      aes(
        x = ifelse(x > 0, x + 0.15, ifelse(x < 0, x - 0.15, ifelse(x == 0, 0, x))),
        y = ifelse(y > 0, y + 0.15, ifelse(y < 0, y - 0.15, ifelse(y == 0 & x == 0, 0.15, y))),
        label = var,
        hjust = ifelse(x < 0, 1, ifelse(x > 0, 0, 0.5)),
        vjust = ifelse(y > 0, 0, ifelse(y < 0, 1, ifelse(y == 0 & x == 0, 0, 0.5)))
      ),
      size = legend_text_size/3,
      color = "black"
    ) +
    ggtitle("Dice arrangement") +
    theme_void()
}

legend_plot <- create_custom_legends_for_map(var_positions, dot_size = 3)

main_plot <- ggplot() +
  geom_sf(data = saarland, fill = "lightblue", color = "black") +
  geom_dice_sf(sf_data = cities_sf,
               dice_value_col = "dice",
               face_color = c("log_France", "log_Swiss", "log_Luxembourg", "log_Rheinlandpfalz"),
               dice_size = 0.5,
               dot_size = 3) +
  geom_text(data = cities_sf,
            mapping = aes(x = st_coordinates(cities_sf)[,1],
                          y = st_coordinates(cities_sf)[,2],
                          label = name),
            nudge_y = 0.03,
            size = 3) +
  ggtitle("Saarland with Dice Markers Showing Log-Scaled Distances to Borders") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",
                        name = "Distance (km)",
                        labels = function(x) round(exp(x)),
                        breaks = log(c(10, 25, 50, 100, 200)))

final_plot <- plot_grid(main_plot, legend_plot, ncol = 2, rel_widths = c(4, 1))
ggsave(final_plot,file = "saarland_geom_dice_sf.png",dpi = 150, width = 10, height = 8)






