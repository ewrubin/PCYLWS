# ============================================================
# 03_ordination_permanova.R
# Purpose: Generate Bray-Curtis PCoA ordination plots and run
# PERMANOVA / pairwise PERMANOVA analyses for microbial community
# composition across sample types and sites.

# ============================================================

# =========================
# 1. Load libraries
# =========================
library(phyloseq)
library(RColorBrewer)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(pairwiseAdonis)
library(ggplot2)
library(dplyr)

# =========================
# 2. Define input/output files
# =========================
input_otu <- "output/psmin5_ASV_table.txt"
input_tax <- "output/psmin5_taxa_table.txt"
metadata_file <- "data/metadata/metadata-v2.txt"
output_dir <- "output/ordination_permanova"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 3. Read input data
# =========================
otu <- read.table(input_otu, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
taxon <- read.table(input_tax, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(as.matrix(taxon))
)

print(ps)

# =========================
# 4. Harmonize metadata fields
# =========================
# metadata-v2 uses:
#   type    = coral / eDNA
#   subtype = healthy / diseased / seawater / sediment
# To reproduce Figure 4, we use subtype as the color grouping variable.
sample_data(ps)$month <- as.factor(sample_data(ps)$month)
sample_data(ps)$month.year <- as.factor(sample_data(ps)$month.year)
sample_data(ps)$site <- as.factor(sample_data(ps)$site)
sample_data(ps)$type <- as.factor(sample_data(ps)$type)
sample_data(ps)$subtype <- as.factor(sample_data(ps)$subtype)

# =========================
# 5. Define colors
# =========================
# Colors follow the four major sample categories shown in Figure 4:
# diseased, healthy, seawater, sediment
group_colors <- c(
  "diseased" = "brown1",
  "healthy" = "darkolivegreen3",
  "seawater" = "cadetblue1",
  "sediment" = "darkmagenta"
)

# =========================
# 6. Ordination plot (Bray-Curtis PCoA)
# =========================
# Taxa are kept at ASV level (identity transform) and Bray-Curtis
# dissimilarity is calculated from the abundance table. The plot
# colors samples by subtype and uses point shape for site.
ordination_plot <- ps %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "subtype", shape = "site") +
  theme_classic(base_size = 12) +
  stat_ellipse(aes(color = subtype)) +
  scale_color_manual(values = group_colors, name = "type")

print(ordination_plot)

ggsave(
  filename = file.path(output_dir, "Figure4_pcoa_bray_subtype_site.tiff"),
  plot = ordination_plot,
  units = "in",
  width = 8,
  height = 6,
  dpi = 300,
  compression = "lzw"
)

# =========================
# 7. Bray-Curtis distance matrix
# =========================
bray_dists <- ps %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc("bray")

# =========================
# 8. PERMANOVA: subtype only
# =========================
# Tests whether microbial community composition differs among the
# major sample categories (healthy, diseased, seawater, sediment).
bray_perm_subtype <- bray_dists %>%
  dist_permanova(
    variables = "subtype",
    seed = 1234,
    n_processes = 1,
    n_perms = 999
  )

permanova_subtype_tbl <- perm_get(bray_perm_subtype) %>% as.data.frame()
write.table(
  permanova_subtype_tbl,
  file = file.path(output_dir, "permanova_subtype.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(permanova_subtype_tbl)

# =========================
# 9. PERMANOVA: subtype + site + month
# =========================
# Tests whether composition varies with subtype, site, and month.
bray_perm_multifactor <- bray_dists %>%
  dist_permanova(
    variables = c("subtype", "site", "month"),
    seed = 321,
    n_processes = 1,
    n_perms = 999
  )

permanova_multifactor_tbl <- perm_get(bray_perm_multifactor) %>% as.data.frame()
write.table(
  permanova_multifactor_tbl,
  file = file.path(output_dir, "permanova_subtype_site_month.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(permanova_multifactor_tbl)

# =========================
# 10. Pairwise PERMANOVA by subtype
# =========================
# pairwise.adonis2 compares all subtype pairs using the Bray-Curtis
# distance matrix with Bonferroni-adjusted p-values.
samples_df <- data.frame(sample_data(ps))

ps_dist_matrix <- distance(ps, method = "bray")

pairwise_results <- pairwise.adonis2(
  ps_dist_matrix ~ subtype,
  data = samples_df,
  p.adjust.m = "bonferroni"
)

capture.output(
  print(pairwise_results),
  file = file.path(output_dir, "pairwise_permanova_subtype.txt")
)

print(pairwise_results)

# =========================
# 11. Save R results object
# =========================
saveRDS(
  list(
    phyloseq_object = ps,
    bray_distances = bray_dists,
    permanova_subtype = permanova_subtype_tbl,
    permanova_subtype_site_month = permanova_multifactor_tbl,
    pairwise_permanova_subtype = pairwise_results
  ),
  file = file.path(output_dir, "ordination_permanova_results.rds")
)
