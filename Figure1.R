library(sf)
library(ggplot2)
library(ggspatial)
library(ggrepel)
library(dplyr)
library(terra)      
library(tidyterra)  
library(patchwork)  

setwd("RRRB/FIGURES")

# ==============================================================================
# 1. COLOR PALETTE & SITES DATA
# ==============================================================================
segment_colors <- c(
  "Upper Basin"                    = "#1b9e77",
  "Middle Basin"                   = "#d95f02",
  "Lower Basin"                    = "#7570b3",
  "Peripheral (Lakes/Rivers/Dams)" = "#e7298a"
)

sites <- data.frame(
  site  = c("Mtr", "Kih", "Kid", "Ifa", "Moh","Kite", "Song", 
            "Ruah", "Rujw","Mans", "Utet", "Rukw", "Ruvu", "Mamb", 
            "Ruhu"),
  lat   = c(-7.128, -8.448, -7.661, -8.191, -7.978, -10.359, -10.625,
            -7.808, -8.708, -7.276, -7.991, -8.397, -6.864, -10.3595, -10.625),
  lon   = c(35.932, 35.195, 36.977, 36.702, 38.979, 39.775, 35.653, 36.897,
            34.388, 39.567, 38.749, 32.902, 37.608, 39.775, 35.653),
  group = c("Upper Basin", "Upper Basin", "Middle Basin", "Middle Basin", 
            "Lower Basin", "Peripheral (Lakes/Rivers/Dams)", "Peripheral (Lakes/Rivers/Dams)", 
            "Upper Basin", "Upper Basin", 
            "Lower Basin", "Lower Basin", "Peripheral (Lakes/Rivers/Dams)",
            "Peripheral (Lakes/Rivers/Dams)", 
            "Peripheral (Lakes/Rivers/Dams)", "Peripheral (Lakes/Rivers/Dams)"),
  stringsAsFactors = FALSE
)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# ==============================================================================
# 2. SPATIAL EXTENT & LAYER CROPPING
# ==============================================================================
xmin <- 32.0; xmax <- 41.5; ymin <- -12.5; ymax <- -5.5
extended_extent <- terra::ext(xmin, xmax, ymin, ymax)

catchments   <- st_read("hybas_lake_af_lev05_v1c.shp")
rivers       <- st_read("HydroRIVERS_v10_af.shp")
lakes        <- st_read("Africa_waterbody.shp")

study_bbox   <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326))
major_rivers <- rivers %>% filter(ORD_FLOW <= 5) %>% st_make_valid() %>% st_crop(study_bbox)

local_master_raster <- terra::rast("TZA_wc2.1_30s_bio.tif")

dem_cropped    <- terra::crop(local_master_raster[[1]], extended_extent)
names(dem_cropped) <- "TZA_elv_srv"

temp_grid      <- terra::crop(local_master_raster[[1]], extended_extent)
names(temp_grid) <- "Temperature"

precip_grid    <- terra::crop(local_master_raster[[12]], extended_extent)
names(precip_grid) <- "Precipitation"

sao_hill_sf    <- st_as_sf(data.frame(lon = 35.25, lat = -8.35, label = "Sao Hill Forest"), coords = c("lon", "lat"), crs = 4326)

# ==============================================================================
# 3. LOCATOR INSET MAP
# ==============================================================================
cat("Assembling Context Inset Map...\n")
world_borders <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>% st_transform(4326)

old_s2_state <- sf_use_s2()
sf_use_s2(FALSE)

tza_context_catchments <- catchments %>% 
  st_make_valid() %>% 
  st_crop(st_bbox(c(xmin = 28.0, xmax = 43.0, ymin = -17.0, ymax = 0.0), crs = 4326))

regional_borders <- world_borders %>%
  st_make_valid() %>%
  st_crop(st_bbox(c(xmin = 28.0, xmax = 43.0, ymin = -17.0, ymax = 0.0), crs = 4326))

sf_use_s2(old_s2_state)

p_inset <- ggplot() +
  geom_sf(data = tza_context_catchments, fill = "grey96", color = "lightblue3", linewidth = 0.15) +
  geom_sf(data = regional_borders, fill = NA, color = "grey40", linewidth = 0.4) +
  geom_sf(data = study_bbox, fill = NA, color = "red", linewidth = 0.8) + 
  geom_sf(data = sites_sf, aes(color = group), size = 1.0) +
  scale_color_manual(values = segment_colors) +
  coord_sf(xlim = c(28.0, 43.0), ylim = c(-17.0, 0.0), expand = FALSE) +
  labs(title = "E") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(2, 2, 2, 2),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6)
  )

# ==============================================================================
# 4. PANELS A, B, C: SPATIAL MAPS
# ==============================================================================
tight_font_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  plot.title = element_text(face = "bold", size = 13),
  axis.title = element_text(face = "bold", size = 10),
  axis.text  = element_text(size = 8),
  legend.title = element_text(face = "bold", size = 10),
  legend.text  = element_text(size = 9),
  plot.margin  = margin(3, 3, 3, 3)
)

# Panel A: Elevation Map
p_topo <- ggplot() +
  geom_spatraster(data = dem_cropped, aes(fill = TZA_elv_srv)) +
  scale_fill_hypso_c(palette = "dem_poster", name = "Elevation (m)") +
  geom_sf(data = major_rivers, color = "royalblue3", linewidth = 0.3, alpha = 0.8) +
  geom_sf(data = lakes, fill = "aliceblue", color = "dodgerblue4", linewidth = 0.25) +
  geom_sf_text(data = sao_hill_sf, aes(label = label), size = 2.8, fontface = "bold.italic", color = "darkgreen") +
  geom_sf(data = sites_sf, aes(color = group), size = 2.2, shape = 16) +
  geom_text_repel(
    data = sites, aes(x = lon, y = lat, label = site),
    size = 2.8, fontface = "bold", box.padding = 0.35, point.padding = 0.2,
    bg.color = "white", bg.r = 0.15, min.segment.length = 0.2
  ) +
  scale_color_manual(values = segment_colors, name = "Basin Segment") +
  guides(color = guide_legend(override.aes = list(size = 3.5))) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  labs(title = "A") +
  tight_font_theme

# Panel B: Temperature Map
p_temp <- ggplot() +
  geom_spatraster(data = temp_grid, aes(fill = Temperature)) +
  scale_fill_gradientn(colors = c("blue", "yellow", "red"), name = "Temp (°C)") +
  geom_sf(data = major_rivers, color = "grey20", linewidth = 0.25, alpha = 0.4) +
  geom_sf(data = sites_sf, aes(color = group), size = 2.2, shape = 16) +
  geom_text_repel(
    data = sites, aes(x = lon, y = lat, label = site),
    size = 2.8, fontface = "bold", box.padding = 0.35, point.padding = 0.2,
    bg.color = "white", bg.r = 0.15, min.segment.length = 0.2
  ) +
  scale_color_manual(values = segment_colors) +
  guides(color = "none") + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "B") +
  tight_font_theme

# Panel C: Precipitation Map
p_climate <- ggplot() +
  geom_spatraster(data = precip_grid, aes(fill = Precipitation)) +
  scale_fill_gradientn(colors = c("#f7fbff", "#9ecae1", "#084594"), name = "Rain (mm)") +
  geom_sf(data = major_rivers, color = "grey20", linewidth = 0.25, alpha = 0.4) +
  geom_sf(data = sites_sf, aes(color = group), size = 2.2, shape = 16) +
  geom_text_repel(
    data = sites, aes(x = lon, y = lat, label = site),
    size = 2.8, fontface = "bold", box.padding = 0.35, point.padding = 0.2,
    bg.color = "white", bg.r = 0.15, min.segment.length = 0.2
  ) +
  scale_color_manual(values = segment_colors) +
  guides(color = "none") + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "C") +
  tight_font_theme

# ==============================================================================
# 5. BIOCLIMATIC PCA
# ==============================================================================
extracted_vals   <- terra::extract(local_master_raster, vect(sites_sf))
extracted_matrix <- as.matrix(extracted_vals[, -1])
rownames(extracted_matrix) <- sites$site
colnames(extracted_matrix) <- paste0("Bio", 1:19)

complete_rows    <- complete.cases(extracted_matrix)
extracted_matrix <- extracted_matrix[complete_rows, ]
sites_pca        <- sites[complete_rows, ]

pca_result <- prcomp(extracted_matrix, center = TRUE, scale. = TRUE)
var_exp    <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

pca_scores <- as.data.frame(pca_result$x[, 1:2]) %>% 
  mutate(pop = rownames(.), group = sites_pca$group)

# Panel D: PCA
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2.8, shape = 16) +
  geom_text_repel(
    aes(label = pop), size = 3.2, fontface = "bold", 
    bg.color = "white", bg.r = 0.15, show.legend = FALSE
  ) +
  scale_color_manual(values = segment_colors) +
  guides(color = "none") + 
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       title = "D") +
  tight_font_theme

# ==============================================================================
# 6. STITCH FIVE-PANEL SYSTEM
# ==============================================================================
cat("Stitching tightly unified 5-panel synthesis...\n")

bottom_row <- (p_climate | p_pca | p_inset) + plot_layout(widths = c(1, 1, 0.55))

stitched_figure <- (p_topo | p_temp) / bottom_row & 
  theme(plot.title = element_text(face = "bold", size = 13))


ggsave("Figure1.jpg", stitched_figure, width = 14, height = 9.5, dpi = 300)

