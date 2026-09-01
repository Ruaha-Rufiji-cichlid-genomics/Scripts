library(bigsnpr)        
library(bigstatsr)
library(adegenet)
library(vegan)          
library(ggplot2)
library(readxl)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(car)
library(scales)     
library(patchwork)   


setwd("RRRB/FIGURES")

# ==============================================================================
# 2. READ PLINK BINARY DATA & FILTER OUT UNMAPPED POPULATIONS
# ==============================================================================
if (!file.exists("./ADMIXTURE/rivercichlids_LDpruned.rds")) {
  snp_readBed("./ADMIXTURE/rivercichlids_LDpruned.bed")
}

bigsnp  <- snp_attach("./ADMIXTURE/rivercichlids_LDpruned.rds")
G_raw   <- bigsnp$genotypes       
fam_raw <- bigsnp$fam            
map     <- bigsnp$map            

plink_ids <- str_trim(as.character(fam_raw$sample.ID))

# Read popmap
popmap <- read.table("popmap.txt", header = FALSE, col.names = c("sample", "pop"), stringsAsFactors = FALSE)
popmap$sample <- str_trim(as.character(popmap$sample))
popmap$pop    <- str_trim(as.character(popmap$pop))

# --- EXCLUDE INSUFICIENT SAMPLE SIZE/UNMAPPED POPULATIONS ---
pops_to_exclude <- c("Unknown", "Wami", "Mansi", "Ruvu", 
                     "Mambi", "Kitele", "Rukwa", "Ruhuhi", "Songea")
cat("Excluding populations without coordinate data:", paste(pops_to_exclude, collapse = ", "), "\n")

popmap_filtered <- popmap %>% 
  filter(!pop %in% pops_to_exclude)

# Match remaining samples to PLINK row indices
retained_inds <- match(popmap_filtered$sample, plink_ids)

if (any(is.na(retained_inds))) {
  plink_fids <- str_trim(as.character(fam_raw$family.ID))
  retained_inds <- match(popmap_filtered$sample, plink_fids)
}

if (any(is.na(retained_inds))) {
  missing_samples <- popmap_filtered$sample[is.na(retained_inds)]
  stop("The following retained samples could not be found in PLINK binary: ", 
       paste(missing_samples, collapse = ", "))
}

cat("Retained", length(retained_inds), "samples across valid populations (Excluded", 
    nrow(fam_raw) - length(retained_inds), "samples total).\n")

n_samples <- length(retained_inds)
n_snps    <- ncol(G_raw)

# ==============================================================================
# 3. HIGH-RESOLUTION ALLELE FREQUENCY FILTERING (VALID SAMPLES ONLY)
# ==============================================================================
cat("Executing Missing Rate Profiling on retained samples...\n")
cnt <- big_counts(G_raw, ind.row = retained_inds)                 
missing_rate <- cnt[4, ] / n_samples   
keep_miss <- which(missing_rate == 0) # Keep strictly complete loci
cat("Retained", length(keep_miss), "SNPs with zero missing data.\n")

cat("Executing Minor Allele Frequency (MAF) Screening...\n")
cnt_keep <- big_counts(G_raw, ind.row = retained_inds, ind.col = keep_miss)   
total_nonmiss <- colSums(cnt_keep[1:3, , drop = FALSE])
af <- (2 * cnt_keep[3, ] + cnt_keep[2, ]) / (2 * total_nonmiss)
maf <- pmin(af, 1 - af)
keep_maf <- which(maf > 0.05)                    

final_keep <- keep_miss[keep_maf]
cat("Final Quality Filter Total:", length(final_keep), "highly polymorphic SNPs.\n")

# ==============================================================================
# 4. COMPUTE POPULATION ALLELE FREQUENCIES
# ==============================================================================
pop <- popmap_filtered$pop
unique_pops <- unique(pop)

freq_matrix <- matrix(NA, nrow = length(unique_pops), ncol = length(final_keep))
rownames(freq_matrix) <- unique_pops

for (i in seq_along(unique_pops)) {
  pop_samples <- popmap_filtered$sample[pop == unique_pops[i]]
  pop_inds    <- match(pop_samples, plink_ids)
  
  stats <- big_colstats(G_raw, ind.row = pop_inds, ind.col = final_keep)
  freq_matrix[i, ] <- stats$sum / (length(pop_inds) * 2)
}

rownames(freq_matrix) <- str_trim(rownames(freq_matrix))
cat("Allele frequency matrix built for", nrow(freq_matrix), "valid populations.\n")

# ==============================================================================
# 5. PREPARE CLIMATE PREDICTORS AND SPATIAL COVARIATES
# ==============================================================================
cat("Assembling Environmental Matrix via Principal Components...\n")
clim_data <- read_excel("extracted_climate_data.xlsx", sheet = 1)

raw_pop_names <- str_trim(as.character(clim_data[[1]]))
clim_matrix   <- as.matrix(clim_data[, -1])
rownames(clim_matrix) <- raw_pop_names

clim_matrix_t <- t(clim_matrix)
rownames(clim_matrix_t) <- str_trim(rownames(clim_matrix_t))

# Verify climate alignment
missing_in_clim <- setdiff(rownames(freq_matrix), rownames(clim_matrix_t))
if (length(missing_in_clim) > 0) {
  stop("Populations missing from climate data: ", paste(missing_in_clim, collapse = ", "))
}

clim_matrix_t <- clim_matrix_t[rownames(freq_matrix), , drop = FALSE]

# Environmental PCA
pca_clim <- prcomp(clim_matrix_t, center = TRUE, scale. = TRUE)
env_pca  <- as.data.frame(scale(pca_clim$x[, 1:2]))
colnames(env_pca) <- c("Env_PC1", "Env_PC2")

cat("Assembling Spatial Predictors via Principal Coordinates (PCoA)...\n")
coords <- read.csv("sites_coords.csv", row.names = 1, sep = "\t")
rownames(coords) <- str_trim(rownames(coords))

# Check for missing site row names
missing_coords <- setdiff(rownames(freq_matrix), rownames(coords))
if (length(missing_coords) > 0) {
  stop("The following populations are missing from sites_coords.csv: ", 
       paste(missing_coords, collapse = ", "))
}

# Align coordinates to frequency matrix
coords_aligned <- coords[rownames(freq_matrix), , drop = FALSE]

# Ensure coordinate columns are numeric
coords_aligned[, 1] <- as.numeric(coords_aligned[, 1])
coords_aligned[, 2] <- as.numeric(coords_aligned[, 2])

# Handle any unexpected NA values in coordinates
if (any(is.na(coords_aligned))) {
  cat("\n!!! WARNING: NA coordinates found! Removing incomplete rows... !!!\n")
  valid_spatial_rows <- complete.cases(coords_aligned)
  
  # Filter all aligned datasets to maintain sample integrity
  coords_aligned <- coords_aligned[valid_spatial_rows, , drop = FALSE]
  freq_matrix    <- freq_matrix[valid_spatial_rows, , drop = FALSE]
  env_pca        <- env_pca[valid_spatial_rows, , drop = FALSE]
  clim_matrix_t  <- clim_matrix_t[valid_spatial_rows, , drop = FALSE]
}

# Compute spatial PCoA eigenvectors (Isolation By Distance control)
geo_dist     <- dist(coords_aligned)
n_sites      <- nrow(coords_aligned)
k_max        <- min(2, n_sites - 1)
pcoa         <- cmdscale(geo_dist, k = k_max, eig = TRUE)

# DEFINE SPATIAL PREDICTORS MATRIX
spatial_pred <- as.data.frame(pcoa$points)
colnames(spatial_pred) <- paste0("Geo_PC", 1:ncol(spatial_pred))

cat("Spatial PCoA successfully assembled across", ncol(spatial_pred), "axes for", 
    nrow(spatial_pred), "populations.\n")

# ==============================================================================
# 6. FORMAL VARIANCE PARTITIONING & MODEL FITTING
# ==============================================================================
cat("\n--- Step 6: Running Variance Partitioning & RDA Models ---\n")

# Scaled environmental predictors
raw_env <- as.data.frame(scale(clim_matrix_t[, c("Bio7", "Bio11", "Bio15", "Bio17")]))
covariates_matrix <- as.matrix(spatial_pred)

# A. Formal Variance Partitioning (Quantify E|S, E∩S, S|E)
var_part <- varpart(freq_matrix, raw_env, covariates_matrix)
cat("\n--- Variance Partitioning Summary ---\n")
print(var_part)

# B. Model 1: Unconditioned RDA (Total Environmental Selection = E|S + E∩S)
rda_unconditioned <- rda(freq_matrix ~ Bio7 + Bio11 + Bio15 + Bio17, data = raw_env)

# C. Model 2: Partial RDA (Pure Micro-Environmental Adaptation = E|S)
rda_conditioned <- rda(freq_matrix ~ Bio7 + Bio11 + Bio15 + Bio17 + Condition(covariates_matrix), data = raw_env)

# D. Verify Variance Inflation Factors (VIF) on Partial RDA
cat("\n--- Variance Inflation Factors (pRDA) ---\n")
print(vif.cca(rda_conditioned))
#print(vif.cca(rda_unconditioned))

# E. Global Significance Tests (ANOVA via Permutation)
cat("\n--- Permutation Test: Unconditioned RDA ---\n")
anova_uncond <- anova.cca(rda_unconditioned, permutations = 999)
print(anova_uncond)

cat("\n--- Permutation Test: Partial RDA ---\n")
anova_cond <- anova.cca(rda_conditioned, permutations = 999)
#print(anova_cond)

# ==============================================================================
# 7. EXTRACT VARIANCE & CANDIDATE ADAPTIVE LOCI
# ==============================================================================
cat("\n--- Step 7: Extracting Loadings & Outlier SNPs ---\n")

# Constrained Eigenvalues & Percentage Variance Explained
eig_u <- eigenvals(rda_unconditioned, model = "constrained")
pct_u <- round((eig_u / sum(eig_u)) * 100, 1)

eig_c <- eigenvals(rda_conditioned, model = "constrained")
pct_c <- round((eig_c / sum(eig_c)) * 100, 1)

# Extract SNP Scores & Outlier Detection from pRDA (>2 SDs on pRDA1 or pRDA2)
snp_scores <- as.data.frame(scores(rda_conditioned, choices = c(1, 2), 
                                   display = "species"))
colnames(snp_scores) <- c("RDA1", "RDA2")

outliers_rda1 <- which(abs(scale(snp_scores$RDA1)) > 2)
outliers_rda2 <- which(abs(scale(snp_scores$RDA2)) > 2)
candidate_idx <- unique(c(outliers_rda1, outliers_rda2))

snp_scores$Type <- "Neutral Loci"
snp_scores$Type[candidate_idx] <- "Candidate Loci"
cat("Identified", length(candidate_idx), "candidate adaptive SNPs (>2 SDs along pRDA axes).\n")

# Helper function to classify Basin Segments matching Figure 1D
assign_basin <- function(df) {
  df %>% mutate(Basin_Segment = case_when(
    str_detect(site, "Ruj|Kih|Mtr|Rujewa|Kihanga|Mtera") ~ "Upper Basin",
    str_detect(site, "Kid|Ifa|Kidatu|Ifakara") ~ "Middle Basin",
    str_detect(site, "Ute|Moh|Utete|Mohoro")   ~ "Lower Basin",
    TRUE                                       ~ "Peripheral"
  ))
}

# Extract Site and Arrow Scores for Model A (Unconditioned)
sites_3a <- as.data.frame(scores(rda_unconditioned, display = "sites", choices = c(1, 2)))
colnames(sites_3a) <- c("RDA1", "RDA2")
sites_3a$site <- rownames(sites_3a)
sites_3a <- assign_basin(sites_3a)

bio_3a <- as.data.frame(scores(rda_unconditioned, display = "bp", choices = c(1, 2)))
colnames(bio_3a) <- c("RDA1", "RDA2")
bio_3a$var <- rownames(bio_3a)
bio_3a <- mutate(bio_3a, RDA1_s = RDA1 * 7.5, RDA2_s = RDA2 * 7.5)

snps_3a <- as.data.frame(scores(rda_unconditioned, display = "species", choices = c(1, 2)))
colnames(snps_3a) <- c("RDA1", "RDA2")

# Extract Site and Arrow Scores for Model B (Conditioned)
sites_3b <- as.data.frame(scores(rda_conditioned, display = "sites", choices = c(1, 2)))
colnames(sites_3b) <- c("RDA1", "RDA2")
sites_3b$site <- rownames(sites_3b)
sites_3b <- assign_basin(sites_3b)

bio_3b <- as.data.frame(scores(rda_conditioned, display = "bp", choices = c(1, 2)))
colnames(bio_3b) <- c("RDA1", "RDA2")
bio_3b$var <- rownames(bio_3b)
bio_3b <- mutate(bio_3b, RDA1_s = RDA1 * 7.5, RDA2_s = RDA2 * 7.5)

# Basin Palette Color Codes matching Figure 1D
basin_colors <- c(
  "Upper Basin"  = "#1b9e77",
  "Middle Basin" = "#d95f02",
  "Lower Basin"  = "#7570b3",
  "Peripheral"   = "#e7298a"
)

# ==============================================================================
# 8. PUBLICATION PATCHWORK TRIPLOT VISUALIZATION (FIGURE 3)
# ==============================================================================
cat("\n--- Step 8: Rendering Patched Publication Figure ---\n")

# PANEL 3A: UNCONDITIONED RDA
p_3a <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(data = snps_3a, aes(x = RDA1, y = RDA2), color = "grey80", size = 0.6, alpha = 0.35) +
  geom_segment(data = bio_3a, aes(x = 0, y = 0, xend = RDA1_s, yend = RDA2_s),
               arrow = arrow(length = unit(0.2, "cm")), color = "#2b83ba", linewidth = 0.85) +
  geom_text(data = bio_3a, aes(x = RDA1_s * 1.15, y = RDA2_s * 1.15, label = var),
            color = "#2b83ba", size = 3.8, fontface = "bold") +
  geom_point(data = sites_3a, aes(x = RDA1, y = RDA2, fill = Basin_Segment),
             shape = 21, color = "black", stroke = 0.8, size = 4.5) +
  geom_text(data = sites_3a, aes(x = RDA1, y = RDA2, label = site),
            size = 3.5, fontface = "bold", vjust = -1.1) +
  scale_fill_manual(values = basin_colors) +
  labs(
    x = paste0("RDA 1 (", pct_u[1], "% Variance)"),
    y = paste0("RDA 2 (", pct_u[2], "% Variance)"),
    title = sprintf("A. Total Environmental Selection (RDA, p = %.3f)", anova_uncond$`Pr(>F)`[1]),
    subtitle = "Includes Spatially Structured Climate Gradients"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11),
    panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.8)
  )

# PANEL 3B: PARTIAL RDA
p_3b <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(data = filter(snp_scores, Type == "Neutral Loci"), 
             aes(x = RDA1, y = RDA2), color = "grey80", size = 0.6, alpha = 0.35) +
  geom_point(data = filter(snp_scores, Type == "Candidate Loci"), 
             aes(x = RDA1, y = RDA2), color = "#d7191c", size = 1.8, alpha = 0.90) +
  geom_segment(data = bio_3b, aes(x = 0, y = 0, xend = RDA1_s, yend = RDA2_s),
               arrow = arrow(length = unit(0.2, "cm")), color = "#2b83ba", linewidth = 0.85) +
  geom_text(data = bio_3b, aes(x = RDA1_s * 1.15, y = RDA2_s * 1.15, label = var),
            color = "#2b83ba", size = 3.8, fontface = "bold") +
  geom_point(data = sites_3b, aes(x = RDA1, y = RDA2, fill = Basin_Segment),
             shape = 21, color = "black", stroke = 0.8, size = 4.5) +
  geom_text(data = sites_3b, aes(x = RDA1, y = RDA2, label = site),
            size = 3.5, fontface = "bold", vjust = -1.1) +
  scale_fill_manual(values = basin_colors) +
  labs(
    x = paste0("pRDA 1 (", pct_c[1], "% Constrained Var)"),
    y = paste0("pRDA 2 (", pct_c[2], "% Constrained Var)"),
    title = sprintf("B. Pure Micro-Environmental Adaptation (pRDA, p = %.3f)", anova_cond$`Pr(>F)`[1]),
    subtitle = "Conditioned on Spatial PCoA (Distance Control)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 11),
    panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.8)
  )

# COMBINE WITH PATCHWORK
figure_3 <- p_3a + p_3b + plot_layout(guides = "collect")

# SAVE 
ggsave("Figure3.png", figure_3, width = 13, height = 6, dpi = 600)
