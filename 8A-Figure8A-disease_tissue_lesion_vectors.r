library(phyloseq)
library(vegan)
library(ape)
library(tidyverse)
library(ggrepel)
library(grid)

# --------------------------------------------------
# 1. Read diseased-only files
# --------------------------------------------------
otu_dis <- read.table("disease_otu_table.txt", sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)

tax_dis <- read.table("disease_taxa_table.txt", sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE)

meta_dis <- read.table("disease_metadata.txt", sep = "\t", header = TRUE,
                       row.names = 1, check.names = FALSE)

# Make sure sample order matches
meta_dis <- meta_dis[rownames(otu_dis), , drop = FALSE]

# --------------------------------------------------
# 2. Build phyloseq object
# --------------------------------------------------
tax_mat <- as.matrix(tax_dis)

ps_dis <- phyloseq(
  otu_table(as.matrix(otu_dis), taxa_are_rows = FALSE),
  sample_data(meta_dis),
  tax_table(tax_mat)
)

# --------------------------------------------------
# 3. Agglomerate to Genus
# --------------------------------------------------
ps_dis_genus <- tax_glom(ps_dis, taxrank = "Genus", NArm = FALSE)
ps_dis_genus <- filter_taxa(ps_dis_genus, function(x) sum(x) > 0, TRUE)

# --------------------------------------------------
# 4. Extract OTU matrix
# --------------------------------------------------
otu_gen <- as(otu_table(ps_dis_genus), "matrix")
if (taxa_are_rows(ps_dis_genus)) {
  otu_gen <- t(otu_gen)
}

# --------------------------------------------------
# 5. Metadata
# --------------------------------------------------
meta_gen <- data.frame(sample_data(ps_dis_genus)) %>%
  rownames_to_column("sampleID")

# Keep only samples with lesion data
meta_gen <- meta_gen %>%
  filter(!is.na(lesion.size), !is.na(progr.rate))

# Keep OTU matrix in same order
otu_gen <- otu_gen[meta_gen$sampleID, , drop = FALSE]

# Optional: set month order
meta_gen$month.year <- factor(
  meta_gen$month.year,
  levels = c("Nov.23", "Feb.24", "May.24", "Aug.24")
)

# --------------------------------------------------
# 6. Bray-Curtis + PCoA
# --------------------------------------------------
bray_dis <- vegdist(otu_gen, method = "bray")
pcoa_dis <- ape::pcoa(bray_dis)

ord_dis <- as.data.frame(pcoa_dis$vectors[, 1:2])
ord_dis$sampleID <- rownames(ord_dis)
colnames(ord_dis)[1:2] <- c("Axis.1", "Axis.2")

plot_dis <- ord_dis %>%
  left_join(meta_gen, by = "sampleID")

# --------------------------------------------------
# 7. envfit
# --------------------------------------------------
envfit_dis <- envfit(
  ord_dis[, c("Axis.1", "Axis.2")],
  plot_dis[, c("lesion.size", "progr.rate")],
  permutations = 999
)

print(envfit_dis)

# --------------------------------------------------
# 8. Extract vectors
# --------------------------------------------------
vec_dis <- as.data.frame(scores(envfit_dis, display = "vectors"))
vec_dis$variable <- rownames(vec_dis)

# Scale arrows for plotting
arrow_mult <- 0.35
vec_dis$Axis.1 <- vec_dis$Axis.1 * arrow_mult
vec_dis$Axis.2 <- vec_dis$Axis.2 * arrow_mult

# --------------------------------------------------
# 9. Colors
# --------------------------------------------------
my_colors_dis <- c(
  "Nov.23" = "darkgoldenrod1",
  "Feb.24" = "cornflowerblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4"
)

# --------------------------------------------------
# 10. Plot
# --------------------------------------------------
p_dis <- ggplot(plot_dis,
                aes(x = Axis.1, y = Axis.2,
                    color = month.year, shape = site)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_segment(data = vec_dis,
               aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black",
               linewidth = 0.8) +
  geom_text_repel(data = vec_dis,
                  aes(x = Axis.1, y = Axis.2, label = variable),
                  inherit.aes = FALSE,
                  size = 4) +
  scale_colour_manual(values = my_colors_dis) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("PCoA1 (", round(pcoa_dis$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_dis$values$Relative_eig[2] * 100, 1), "%)"),
    color = "month.year",
    shape = "site"
  )

print(p_dis)

ggsave("PCoA_envfit_diseased_genus_clean.tiff",
       p_dis,
       width = 6, height = 4.5, dpi = 600, compression = "lzw")