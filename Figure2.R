# ==============================================================================
# Phylogeny + ADMIXTURE (K=3, K=4)
# ==============================================================================
library(ggtree)
library(treeio)
library(tidyverse)
library(ape)
library(pophelper)
library(reshape2)
library(patchwork)

setwd("RRRB/FIGURES")

# ==============================================================================
# 1. READING TREE DATA & MeTADATA
# ==============================================================================
basin_colors <- c(
  "Upper Basin"                    = "#1b9e77",
  "Middle Basin"                   = "#d95f02",
  "Lower Basin"                    = "#7570b3",
  "Peripheral (Lakes/Rivers/Dams)" = "#e7298a",
  "Outgroup"                       = "#666666"
)

# Shared master palette for all ancestral clusters across K=3 and K=4
cluster_colors <- c(
  "Cluster 1" = "#2b83ba", # Blue
  "Cluster 2" = "#d7191c", # Red
  "Cluster 3" = "#fdae61", # Orange
  "Cluster 4" = "#abdda4"  # Green
)

raw_tree <- read.tree("species_tree_astral.tree")
meta     <- read.csv("phylo_metadata.txt", stringsAsFactors = FALSE)

# Clean strings
raw_tree$tip.label <- str_trim(raw_tree$tip.label)

meta <- meta %>%
  mutate(
    id            = str_trim(id),
    basin_segment = str_trim(basin_segment),
    species       = str_trim(species),
    species_clean = species %>% str_replace_all("['\"]", "") %>% str_squish()
  )

# Outgroup setup
outgroup_taxa <- intersect(c("HapVanKid-01", "HapVanRuh-01", "HapVanRuh-02"), raw_tree$tip.label)

if(length(outgroup_taxa) > 0) {
  rooted_tree <- root(raw_tree, outgroup = outgroup_taxa, resolve.root = TRUE)
} else {
  rooted_tree <- raw_tree
}

if (any(rooted_tree$edge.length == 0)) {
  rooted_tree$edge.length[rooted_tree$edge.length == 0] <- 1e-6
}

# ==============================================================================
# 2. MeTADATA MATCHING
# ==============================================================================meta_tree <- data.frame(label = rooted_tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(meta, by = c("label" = "id")) %>%
  mutate(
    basin_segment = if_else(is.na(basin_segment) | basin_segment == "", "Outgroup", basin_segment),
    species_clean = if_else(is.na(species_clean) | species_clean == "", label, species_clean)
  )

tree_data <- as.treedata(rooted_tree)

# ==============================================================================
# 3. A. PHYLOGENY
# ==============================================================================
p_tree <- ggtree(
  tree_data, 
  layout = "rectangular", 
  branch.length = "branch.length",
  color = "grey20", 
  linewidth = 0.75
) %<+% meta_tree

# Support values
p_tree$data <- p_tree$data %>%
  mutate(
    raw_pp = suppressWarnings(as.numeric(str_extract(label, "^[0-9.]+"))),
    parsed_support = if_else(!isTip & !is.na(raw_pp) & raw_pp >= 0.95,
                             sprintf("%.2f", raw_pp), NA_character_)
  )

p_tree <- p_tree + 
  geom_nodelab(
    aes(label = parsed_support), 
    size = 3.2, hjust = 1.15, vjust = -0.35, color = "grey30"
  ) +
  geom_tippoint(
    aes(fill = basin_segment), 
    shape = 21, color = "black", stroke = 0.5, size = 3.8
  ) +
  geom_tiplab(
    aes(label = species_clean),
    size = 4.3,
    offset = 0.35,
    align = FALSE,
    color = "black",
    fontface = "italic"  # ITALICIZES SPECIES NAMES
  ) +
  scale_fill_manual(values = basin_colors, name = "Basin Segment", na.value = "#666666") +
  theme_tree() + 
  theme(
    legend.title = element_text(face = "bold", size = 12),
    legend.text  = element_text(size = 11),
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.margin  = margin(10, 5, 10, 5)
  ) +
  labs(title = "A")

# Vertical Clade Strips
tip_data <- p_tree$data %>% 
  filter(isTip) %>% 
  arrange(y) %>% 
  select(label, y, species_clean, basin_segment)

max_tree_x <- max(p_tree$data$x, na.rm = TRUE)

clade_blocks <- tip_data %>%
  mutate(block_id = cumsum(species_clean != lag(species_clean, default = first(species_clean)))) %>%
  group_by(block_id, species_clean, basin_segment) %>%
  summarise(
    y_min     = min(y),
    y_max     = max(y),
    tip_start = label[which.min(y)],
    tip_end   = label[which.max(y)],
    count     = n(),
    .groups   = "drop"
  )

for (i in seq_len(nrow(clade_blocks))) {
  cb <- clade_blocks[i, ]
  bar_col <- basin_colors[cb$basin_segment]
  if(is.na(bar_col)) bar_col <- "#666666"
  
  if(cb$count > 1) {
    p_tree <- p_tree + 
      geom_strip(
        taxa1    = cb$tip_start, 
        taxa2    = cb$tip_end, 
        label    = "", 
        color    = bar_col, 
        barsize  = 4.5, 
        offset   = max_tree_x * 0.52  # Calibrated offset to sit cleanly past italic labels
      )
  }
}

p_tree <- p_tree + hexpand(0.60)

# ==============================================================================
# 4. CORRELATE WITH ADMIXTURE DATA
# ==============================================================================
tree_tip_order <- get_taxa_name(p_tree)

q3_raw  <- readQ("/ADMIXTURE/rivercichlids_LDpruned.3.Q")[[1]]
q4_raw  <- readQ("/ADMIXTURE/rivercichlids_LDpruned.4.Q")[[1]]
fam     <- read.table("/ADMIXTURE/rivercichlids_LDpruned.fam", stringsAsFactors = FALSE)
ind_names <- str_trim(fam[, 1])

# Match K=4 columns to K=3 columns
cor_mat <- cor(q3_raw, q4_raw)
best_matches <- apply(cor_mat, 1, which.max)

unmatched_k4 <- setdiff(1:4, best_matches)
new_k4_order <- c(best_matches, unmatched_k4)
q4_aligned   <- q4_raw[, new_k4_order]

# Assign uniform, clean cluster names
colnames(q3_raw)     <- paste("Cluster", 1:3)
colnames(q4_aligned) <- paste("Cluster", 1:4)

format_admix_df <- function(q_mat, ind_names, tip_order) {
  df <- cbind(Sample = ind_names, q_mat)
  df_melt <- melt(df, id.vars = "Sample", variable.name = "Cluster", value.name = "Ancestry")
  
  full_grid <- expand.grid(
    Sample = tip_order, 
    Cluster = paste("Cluster", 1:4),
    stringsAsFactors = FALSE
  )
  
  df_aligned <- full_grid %>%
    left_join(df_melt, by = c("Sample", "Cluster")) %>%
    mutate(
      Sample   = factor(Sample, levels = rev(tip_order)),
      Cluster  = factor(Cluster, levels = paste("Cluster", 1:4)),
      Ancestry = if_else(is.na(Ancestry), 0, Ancestry)
    )
  
  return(df_aligned)
}

q3_df <- format_admix_df(q3_raw, ind_names, tree_tip_order)
q4_df <- format_admix_df(q4_aligned, ind_names, tree_tip_order)

# ==============================================================================
# 5. B. ADMIXTURE PLOTS
# ==============================================================================
p_k3 <- ggplot(q3_df, aes(y = Sample, x = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.95, position = "stack") +
  scale_fill_manual(values = cluster_colors, name = "Ancestry Cluster", drop = FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1.0), labels = c("0", "0.5", "1")) +
  labs(title = "K = 3", x = "Ancestry") +
  theme_minimal() +
  theme(
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.title.x     = element_text(size = 14, face = "bold"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    plot.margin      = margin(10, 6, 10, 2)
  )

p_k4 <- ggplot(q4_df, aes(y = Sample, x = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.95, position = "stack") +
  scale_fill_manual(values = cluster_colors, name = "Ancestry Cluster", drop = FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1.0), labels = c("0", "0.5", "1")) +
  labs(title = "K = 4", x = "Ancestry") +
  theme_minimal() +
  theme(
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.title.x     = element_text(size = 14, face = "bold"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    plot.margin      = margin(10, 2, 10, 6)
  )

p_admix <- (p_k3 | p_k4) + 
  plot_annotation(
    title = "B",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )
    )
  )

# ==============================================================================
# 6. TWO PLOTS INTEGRATION 
# ==============================================================================
combined_figure <- (p_tree | p_admix) + 
  plot_layout(widths = c(2.8, 1.2), guides = "collect") &
  theme(
    legend.position = "top",
    legend.box      = "horizontal"
  )

ggsave("Figure_Phylo_ADMIXTURE_Integrated.png", combined_figure, width = 22, height = 18, dpi = 600)