library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(scales)
library(stringr)
library(readr)
library(tidyr)
library(ggrepel)

setwd("RRRB/FIGURES")

pop_colors <- c(
  "Kihanga" = "#1b9e77", # Upper Basin
  "Mtera"   = "#238b45", # Upper Basin
  "Kidatu"  = "#d95f02", # Middle Basin
  "Ifakara" = "#e6550d", # Middle Basin
  "Mohoro"  = "#7570b3"  # Lower Basin
)

# --- 2. Load and Clean SnpEff Deleterious Mutation Matrix ---
snpeff_raw <- read_tsv(
  "deleterious_mutation_matrix.txt", 
  col_names = TRUE, 
  col_types = cols(.default = col_character()), 
  show_col_types = FALSE, 
  comment = "#"
)

colnames(snpeff_raw)[1:5] <- c("CHROM", "POSITION", "GENE", "EFFECT", "IMPACT")

snpeff_clean <- snpeff_raw %>%
  filter(!is.na(IMPACT) & IMPACT != "") %>%
  separate_rows(IMPACT, sep = ",") %>%
  mutate(
    CHROM    = gsub("^chr", "", as.character(CHROM)),
    POSITION = as.numeric(POSITION)
  )

# --- 3. Load and Harmonize RDA Candidate Data ---
rda_raw <- read.delim("../RDA/rda_candidates_annotated.tsv", header = TRUE, stringsAsFactors = FALSE)

selection_clean <- rda_raw %>%
  mutate(
    CHR       = gsub("^chr", "", as.character(CHR)),
    POSITION  = as.numeric(POS),
    SNP_ID    = paste0(CHR, "_", POSITION),
    RDA_Score = abs(RDA1),
    status    = if_else(is_candidate == TRUE | is_candidate == 1, "RDA_Outlier", "Neutral")
  )

# --- 4. Load iHS Master Files ---
ihs_files <- list(
  "Kihanga" = "ihs_master_Kihanga.csv",
  "Mtera"   = "ihs_master_Mtera.csv",
  "Kidatu"  = "ihs_master_Kidatu.csv",
  "Ifakara" = "ihs_master_Ifakara.csv",
  "Mohoro"  = "ihs_master_Mohoro.csv"
)

ihs_all_pops <- map_dfr(names(ihs_files), function(pop) {
  f <- ihs_files[[pop]]
  if (file.exists(f)) {
    read.csv(f, stringsAsFactors = FALSE) %>%
      mutate(
        CHR        = gsub("^chr", "", as.character(CHR)),
        POSITION   = as.numeric(POSITION),
        SNP_ID     = paste0(CHR, "_", POSITION),
        Population = pop
      )
  }
})

# --- 5. Export HIGH & MODERATE Impact Candidate Genes Summary Table ---
candidate_snp_ids <- unique(c(
  selection_clean %>% filter(status == "RDA_Outlier") %>% pull(SNP_ID),
  ihs_all_pops %>% filter(outlier == TRUE) %>% pull(SNP_ID)
))

safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA_real_)
  max(x_clean)
}

functional_candidates_df <- snpeff_clean %>%
  filter(IMPACT %in% c("HIGH", "MODERATE")) %>%
  mutate(
    CHROM  = gsub("^chr", "", as.character(CHROM)),
    SNP_ID = paste0(CHROM, "_", POSITION)
  ) %>%
  filter(SNP_ID %in% candidate_snp_ids) %>%
  left_join(
    select(selection_clean, SNP_ID, RDA_Score, status),
    by = "SNP_ID"
  ) %>%
  group_by(CHROM, GENE, IMPACT) %>%
  summarise(
    Total_Candidate_SNPs = n_distinct(SNP_ID),
    RDA_Outlier_Count    = sum(status == "RDA_Outlier", na.rm = TRUE),
    Max_RDA_Score        = safe_max(RDA_Score),
    Position_Range       = paste0(min(POSITION), "-", max(POSITION)),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(CHROM), desc(IMPACT), desc(RDA_Outlier_Count))

write.csv(functional_candidates_df, 
          "Astatotilapia_HIGH_MODERATE_candidate_genes.csv", 
          row.names = FALSE)

# --- 6. Build Candidate Gene Annotation Set ---
gene_annotations_high_mod <- selection_clean %>%
  filter(status == "RDA_Outlier") %>%
  inner_join(
    snpeff_clean %>% filter(IMPACT %in% c("HIGH", "MODERATE")), 
    by = c("CHR" = "CHROM", "POSITION" = "POSITION")
  ) %>%
  filter(!is.na(GENE) & GENE != "") %>%
  rename(CHROM = CHR) %>%
  group_by(CHROM, GENE) %>%
  slice_max(order_by = RDA_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

gene_annotations_fallback <- selection_clean %>%
  filter(status == "RDA_Outlier") %>%
  inner_join(
    snpeff_clean, 
    by = c("CHR" = "CHROM", "POSITION" = "POSITION")
  ) %>%
  filter(!is.na(GENE) & GENE != "") %>%
  rename(CHROM = CHR) %>%
  group_by(CHROM) %>%
  slice_max(order_by = RDA_Score, n = 10, with_ties = FALSE) %>%
  ungroup()

all_gene_annotations <- bind_rows(gene_annotations_high_mod, gene_annotations_fallback) %>%
  group_by(CHROM, GENE) %>%
  slice_max(order_by = RDA_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- 7. Load Diversity Tracks (\u03c0) ---
load_pi_track <- function(file_path, pop_name) {
  if (!file.exists(file_path)) return(NULL)
  read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(
      CHROM      = gsub("^chr", "", as.character(CHROM)),
      Population = pop_name
    ) %>%
    select(CHROM, BIN_START, BIN_END, PI, Population)
}

pi_tracks <- bind_rows(
  load_pi_track("AstRuhBKih_pi.windowed.pi", "Kihanga"),
  load_pi_track("AstRuhBMtr_pi.windowed.pi", "Mtera"),
  load_pi_track("AstRuhBKid_pi.windowed.pi", "Kidatu"),
  load_pi_track("AstGiIfa_pi.windowed.pi",   "Ifakara"),
  load_pi_track("AstGiMoh_pi.windowed.pi",   "Mohoro")
)

# --- 8. Compute Windowed Runs of Homozygosity (ROH Density in Mb) ---
load_roh_track <- function(file_path, pop_name, pi_windows) {
  if (!file.exists(file_path) || file.info(file_path)$size == 0) return(NULL)
  
  roh_df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(CHR = gsub("^chr", "", as.character(CHR)))
  
  pop_windows <- pi_windows %>% filter(Population == pop_name)
  if (nrow(pop_windows) == 0) return(NULL)
  
  roh_counts <- sapply(seq_len(nrow(pop_windows)), function(i) {
    chr     <- pop_windows$CHROM[i]
    w_start <- pop_windows$BIN_START[i]
    w_end   <- pop_windows$BIN_END[i]
    
    sub_roh <- roh_df %>% filter(CHR == chr, POS1 <= w_end, POS2 >= w_start)
    if (nrow(sub_roh) == 0) return(0)
    
    overlap_bp <- pmin(sub_roh$POS2, w_end) - pmax(sub_roh$POS1, w_start)
    sum(pmax(0, overlap_bp)) / 1e6
  })
  
  pop_windows$ROH_Mb <- roh_counts
  return(pop_windows)
}

roh_tracks <- bind_rows(
  load_roh_track("AstRuhBKih_roh.hom", "Kihanga", pi_tracks),
  load_roh_track("AstRuhBMtr_roh.hom", "Mtera",   pi_tracks),
  load_roh_track("AstRuhBKid_roh.hom", "Kidatu",   pi_tracks),
  load_roh_track("AstGiIfa_roh.hom",   "Ifakara",  pi_tracks),
  load_roh_track("AstGiMoh_roh.hom",   "Mohoro",   pi_tracks)
)

# --- 9. Vectorized SnpEff Deleterious Load Calculation ---
mut_list <- split(snpeff_clean$POSITION, snpeff_clean$CHROM)

unique_windows <- pi_tracks %>%
  distinct(CHROM, BIN_START, BIN_END) %>%
  arrange(CHROM, BIN_START) %>%
  group_by(CHROM) %>%
  group_modify(~ {
    chr_name   <- .y$CHROM[1]
    chrom_muts <- mut_list[[chr_name]]
    
    if (is.null(chrom_muts) || length(chrom_muts) == 0) {
      .x$Global_Deleterious_Count <- 0
      return(.x)
    }
    
    chrom_muts  <- sort(chrom_muts)
    bin_indices <- findInterval(chrom_muts, .x$BIN_START)
    valid_mask  <- bin_indices > 0 & bin_indices <= nrow(.x)
    valid_mask[valid_mask] <- chrom_muts[valid_mask] <= .x$BIN_END[bin_indices[valid_mask]]
    
    counts <- table(factor(bin_indices[valid_mask], levels = seq_len(nrow(.x))))
    .x$Global_Deleterious_Count <- as.numeric(counts)
    return(.x)
  }) %>%
  ungroup()

master_tracks_5pop <- pi_tracks %>%
  mutate(
    BIN_START = as.numeric(BIN_START),
    BIN_END   = as.numeric(BIN_END)
  ) %>%
  left_join(unique_windows, by = c("CHROM", "BIN_START", "BIN_END")) %>%
  left_join(
    select(roh_tracks, CHROM, BIN_START, BIN_END, Population, ROH_Mb),
    by = c("CHROM", "BIN_START", "BIN_END", "Population")
  ) %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  mutate(
    PI_Weight         = (PI + 1e-6) / (sum(PI, na.rm = TRUE) + 1e-6),
    Deleterious_Count = Global_Deleterious_Count * PI_Weight * 3
  ) %>%
  ungroup()

# --- 10. Multi-Track Generator Loop ---
target_chromosomes <- unique(c(
  selection_clean %>% filter(status == "RDA_Outlier") %>% pull(CHR),
  ihs_all_pops %>% filter(outlier == TRUE) %>% pull(CHR)
))

for (current_chrom in target_chromosomes) {
  
  cat("\n==================================================\n")
  cat(" Generating Synthesis Plot: Chromosome", current_chrom, "\n")
  cat("==================================================\n")
  
  chr_annotations <- all_gene_annotations %>% filter(CHROM == current_chrom)
  selection_chr   <- selection_clean %>% filter(CHR == current_chrom)
  ihs_chr         <- ihs_all_pops %>% filter(CHR == current_chrom)
  
  ihs_outliers_chr <- ihs_chr %>% filter(outlier == TRUE) %>% pull(SNP_ID) %>% unique()
  
  ihs_neutral_bg <- ihs_chr %>% 
    filter(outlier == FALSE) %>%
    distinct(CHR, POSITION, SNP_ID) %>%
    mutate(
      RDA_Score   = 0.0005,
      status      = "Neutral",
      point_class = "Neutral Background"
    )
  
  selection_chr <- selection_chr %>%
    mutate(
      is_rda_outlier = status == "RDA_Outlier",
      is_ihs_outlier = SNP_ID %in% ihs_outliers_chr,
      point_class    = case_when(
        is_rda_outlier & is_ihs_outlier ~ "RDA + iHS Overlap",
        is_rda_outlier                 ~ "RDA Outlier",
        is_ihs_outlier                 ~ "iHS Outlier",
        TRUE                           ~ "Neutral Background"
      )
    )
  
  set.seed(42)
  sampled_bg <- ihs_neutral_bg %>% sample_n(min(2000, nrow(ihs_neutral_bg)))
  
  track1_data <- bind_rows(
    selection_chr,
    filter(sampled_bg, !SNP_ID %in% selection_chr$SNP_ID)
  )
  
  overlap_snps <- selection_chr %>% filter(point_class == "RDA + iHS Overlap")
  has_overlap  <- nrow(overlap_snps) > 0
  
  all_outliers <- selection_chr %>% filter(point_class %in% c("RDA Outlier", "iHS Outlier", "RDA + iHS Overlap"))
  
  if (nrow(all_outliers) > 0) {
    min_pos <- min(all_outliers$POSITION, na.rm = TRUE)
    max_pos <- max(all_outliers$POSITION, na.rm = TRUE)
  } else {
    min_pos <- min(selection_chr$POSITION, na.rm = TRUE)
    max_pos <- max(selection_chr$POSITION, na.rm = TRUE)
  }
  
  buffer_bp  <- 2000000
  zoom_start <- max(0, min_pos - buffer_bp)
  zoom_end   <- max_pos + buffer_bp
  
  track1_zoomed    <- track1_data %>% filter(POSITION >= zoom_start, POSITION <= zoom_end)
  tracks_5pop_zoom <- master_tracks_5pop %>% filter(CHROM == current_chrom, BIN_START >= zoom_start, BIN_END <= zoom_end)
  
  max_y_val <- max(track1_zoomed$RDA_Score, na.rm = TRUE)
  
  shared_theme_zoom <- theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y     = element_text(size = 12, color = "black", margin = margin(r = 8)),
      axis.text.y      = element_text(size = 11, color = "black"),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(face = "bold", size = 15, color = "black"),
      plot.margin      = margin(t = 4, r = 12, b = 4, l = 12)
    )
  
  # --------------------------------------------------------------------------
  # TRACK 1: SELECTION LANDSCAPE (Direct Gene SNP Leader Lines)
  # --------------------------------------------------------------------------
  p1_zoom <- ggplot(track1_zoomed, aes(x = POSITION, y = RDA_Score)) +
    geom_point(data = filter(track1_zoomed, point_class == "Neutral Background"), 
               color = "grey85", alpha = 0.5, size = 1.0) +
    geom_point(data = filter(track1_zoomed, point_class == "RDA Outlier"), 
               color = "#d95f02", size = 2.5, alpha = 0.85) +
    geom_point(data = filter(track1_zoomed, point_class == "iHS Outlier"), 
               color = "#7570b3", size = 2.5, alpha = 0.85) +
    geom_point(data = filter(track1_zoomed, point_class == "RDA + iHS Overlap"), 
               color = "#1b9e77", size = 4.2, shape = 18) +
    geom_text_repel(
      data = filter(chr_annotations, POSITION >= zoom_start & POSITION <= zoom_end),
      aes(x = POSITION, y = RDA_Score, label = GENE),
      fontface = "italic", size = 4.2, color = "darkred",
      box.padding = 0.8, point.padding = 0.5,
      min.segment.length = 0,
      segment.color = "grey30", segment.linewidth = 0.6,
      max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(zoom_start, zoom_end), ylim = c(0, max_y_val * 1.35)) +
    labs(
      y = "RDA Loading", 
      title = paste0("")
    ) +
    shared_theme_zoom
  
  # --------------------------------------------------------------------------
  # TRACK 2: STANDING GENETIC DIVERSITY (\u03c0)
  # --------------------------------------------------------------------------
  p2_zoom <- ggplot(tracks_5pop_zoom, aes(x = BIN_START, y = PI, colour = Population)) + 
    geom_smooth(aes(group = Population), method = "loess", span = 0.25, se = FALSE, linewidth = 1.1) + 
    geom_vline(data = filter(chr_annotations, POSITION >= zoom_start & POSITION <= zoom_end), 
               aes(xintercept = POSITION), colour = "darkred", linetype = "dashed", linewidth = 0.35, alpha = 0.6) +
    scale_color_manual(values = pop_colors) +
    scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
    coord_cartesian(xlim = c(zoom_start, zoom_end)) +
    labs(y = "Diversity (\u03c0)") +
    shared_theme_zoom +
    theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 11))
  
  # --------------------------------------------------------------------------
  # TRACK 3: RUNS OF HOMOZYGOSITY (ROH Density in Mb)
  # --------------------------------------------------------------------------
  p3_zoom <- ggplot(tracks_5pop_zoom, aes(x = BIN_START, y = ROH_Mb, colour = Population)) +
    geom_smooth(aes(group = Population), method = "loess", span = 0.25, se = FALSE, linewidth = 1.1) +
    geom_vline(data = filter(chr_annotations, POSITION >= zoom_start & POSITION <= zoom_end), 
               aes(xintercept = POSITION), colour = "darkred", linetype = "dashed", linewidth = 0.35, alpha = 0.6) +
    scale_color_manual(values = pop_colors) +
    scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
    coord_cartesian(xlim = c(zoom_start, zoom_end)) +
    labs(y = "ROH (Mb)") +
    shared_theme_zoom
  
  # --------------------------------------------------------------------------
  # TRACK 4: DELETERIOUS MUTATION LOAD
  # --------------------------------------------------------------------------
  p4_zoom <- ggplot(tracks_5pop_zoom, aes(x = BIN_START, y = Deleterious_Count, color = Population)) +
    geom_smooth(aes(group = Population), method = "loess", span = 0.25, se = FALSE, linewidth = 1.2) +
    geom_vline(data = filter(chr_annotations, POSITION >= zoom_start & POSITION <= zoom_end), 
               aes(xintercept = POSITION), color = "darkred", linetype = "dashed", linewidth = 0.35, alpha = 0.6) +
    scale_color_manual(values = pop_colors) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
    coord_cartesian(xlim = c(zoom_start, zoom_end)) +
    labs(y = "Deleterious Load", x = paste("Position on Chromosome", current_chrom, "(bp)")) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(), 
      legend.position  = "none",
      axis.title.y     = element_text(size = 12, color = "black", margin = margin(r = 8)),
      axis.title.x     = element_text(size = 12, color = "black", margin = margin(t = 8)),
      axis.text.y      = element_text(size = 11, color = "black"),
      axis.text.x      = element_text(size = 11, color = "black"), 
      axis.ticks       = element_line(color = "black", linewidth = 0.6),
      plot.margin      = margin(t = 2, r = 12, b = 10, l = 12)
    )
  
  comprehensive_synthesis_plot <- p1_zoom / p2_zoom / p3_zoom / p4_zoom + 
    plot_layout(heights = c(1.8, 1.2, 1.2, 1.2))
  
  filename_output <- paste0("RDA5pop_hybrid_focus_chr", current_chrom, ".png")
  ggsave(filename_output, comprehensive_synthesis_plot, width = 13, height = 12, dpi = 300)
}