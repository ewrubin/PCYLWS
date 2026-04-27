library(phyloseq)
library(vegan)
library(ape)
library(tidyverse)
library(ggrepel)
library(grid)


otu <- read.table("psmin5_ASV_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("psmin5_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadata.txt",sep="\t",header=TRUE,row.names=1)


taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(TAX))

to_remove <- c("LP196","eDNA1")
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

ntaxa(ps) 


chla   <- read.delim("chla-data.txt",   sep = "\t", header = TRUE, check.names = FALSE)
nit    <- read.delim("nitrate-data.txt", sep = "\t", header = TRUE, check.names = FALSE)



# seawater only
ps_sw <- subset_samples(ps, subtype == "seawater")
ps_sw <- prune_samples(!(sample_names(ps_sw) %in% to_remove), ps_sw)
ps_sw <- filter_taxa(ps_sw, function(x) sum(x) > 0, TRUE)

# agglomerate to Genus to match original Figure 5B
ps_sw_genus <- tax_glom(ps_sw, taxrank = "Genus", NArm = FALSE)

# remove taxa that are zero after agglomeration
ps_sw_genus <- filter_taxa(ps_sw_genus, function(x) sum(x) > 0, TRUE)

# metadata
meta_sw <- data.frame(sample_data(ps_sw_genus)) %>%
  rownames_to_column("sampleID") %>%
  left_join(chla_sum, by = c("site", "month.year")) %>%
  left_join(nit_sum,  by = c("site", "month.year"))

# otu matrix
otu_sw <- as(otu_table(ps_sw_genus), "matrix")
if (taxa_are_rows(ps_sw_genus)) {
  otu_sw <- t(otu_sw)
}

# Bray-Curtis + PCoA
bray_sw <- vegdist(otu_sw, method = "bray")
pcoa_sw <- ape::pcoa(bray_sw)

ord_sw <- as.data.frame(pcoa_sw$vectors[, 1:2])
ord_sw$sampleID <- rownames(ord_sw)
colnames(ord_sw)[1:2] <- c("Axis.1", "Axis.2")

plot_sw <- ord_sw %>%
  left_join(meta_sw, by = "sampleID")

# envfit
envfit_sw <- envfit(
  ord_sw[, c("Axis.1", "Axis.2")],
  plot_sw[, c("mean_chla", "mean_nit")],
  permutations = 999
)

print(envfit_sw)	   


# -----------------------------
# envfit
# -----------------------------
envfit_sw <- envfit(
  ord_sw[, c("Axis.1", "Axis.2")],
  plot_sw[, c("mean_chla", "mean_nit")],
  permutations = 999
)

print(envfit_sw)

vec_sw <- as.data.frame(scores(envfit_sw, display = "vectors"))
vec_sw$variable <- rownames(vec_sw)

# -----------------------------
# plot
# -----------------------------
p_sw <- ggplot(plot_sw, aes(x = Axis.1, y = Axis.2,
                            color = month.year, shape = site)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_segment(data = vec_sw,
               aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black",
               linewidth = 0.8) +
  geom_text_repel(data = vec_sw,
                  aes(x = Axis.1, y = Axis.2, label = variable),
                  inherit.aes = FALSE,
                  color = "black",
                  size = 4) +
  scale_colour_manual(values = my_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PCoA1 (", round(pcoa_sw$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_sw$values$Relative_eig[2] * 100, 1), "%)"),
    color = "Month-Year",
    shape = "Site"
  )

print(p_sw)

ggsave("PCoA_envfit_seawater_genus-level.tiff", p_sw,
       width = 6, height = 4.5, dpi = 600, compression = "lzw")
	   
	   
	   
library(phyloseq)
library(vegan)
library(ape)
library(tidyverse)
library(ggrepel)
library(grid)

# --------------------------------------------------
# Sediment only
# --------------------------------------------------
ps_sed <- subset_samples(ps, subtype == "sediment")
ps_sed <- prune_samples(!(sample_names(ps_sed) %in% to_remove), ps_sed)
ps_sed <- filter_taxa(ps_sed, function(x) sum(x) > 0, TRUE)

# --------------------------------------------------
# Agglomerate to Genus to match Figure 5A style
# --------------------------------------------------
ps_sed_genus <- tax_glom(ps_sed, taxrank = "Genus", NArm = FALSE)
ps_sed_genus <- filter_taxa(ps_sed_genus, function(x) sum(x) > 0, TRUE)

# --------------------------------------------------
# Metadata + environmental summaries
# assumes chla_sum and nit_sum already exist
# --------------------------------------------------
meta_sed <- data.frame(sample_data(ps_sed_genus)) %>%
  rownames_to_column("sampleID") %>%
  left_join(chla_sum, by = c("site", "month.year")) %>%
  left_join(nit_sum,  by = c("site", "month.year"))

# optional: put merged metadata back
sample_data(ps_sed_genus) <- sample_data(meta_sed %>% column_to_rownames("sampleID"))

# --------------------------------------------------
# OTU matrix
# --------------------------------------------------
otu_sed <- as(otu_table(ps_sed_genus), "matrix")
if (taxa_are_rows(ps_sed_genus)) {
  otu_sed <- t(otu_sed)
}

# --------------------------------------------------
# Bray-Curtis + PCoA
# --------------------------------------------------
bray_sed <- vegdist(otu_sed, method = "bray")
pcoa_sed <- ape::pcoa(bray_sed)

ord_sed <- as.data.frame(pcoa_sed$vectors[, 1:2])
ord_sed$sampleID <- rownames(ord_sed)
colnames(ord_sed)[1:2] <- c("Axis.1", "Axis.2")

plot_sed <- ord_sed %>%
  left_join(meta_sed, by = "sampleID")

# --------------------------------------------------
# envfit
# --------------------------------------------------
envfit_sed <- envfit(
  ord_sed[, c("Axis.1", "Axis.2")],
  plot_sed[, c("mean_chla", "mean_nit")],
  permutations = 999
)

print(envfit_sed)

# --------------------------------------------------
# Extract vectors
# --------------------------------------------------
vec_sed <- as.data.frame(scores(envfit_sed, display = "vectors"))
vec_sed$variable <- rownames(vec_sed)

# optional: scale arrows a bit for visibility
arrow_mult <- 0.35
vec_sed$Axis.1 <- vec_sed$Axis.1 * arrow_mult
vec_sed$Axis.2 <- vec_sed$Axis.2 * arrow_mult

# --------------------------------------------------
# Month colors to match Figure 5
# --------------------------------------------------
my_colors <- c(
  "Oct.23" = "darkgoldenrod4",
  "Nov.23" = "darkgoldenrod1",
  "Jan.24" = "cornflowerblue",
  "Mar.24" = "darkslateblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4"
)

plot_sed$month.year <- factor(
  plot_sed$month.year,
  levels = c("Oct.23", "Nov.23", "Jan.24", "Mar.24", "May.24", "Aug.24")
)

# --------------------------------------------------
# Plot
# --------------------------------------------------
p_sed_genus <- ggplot(plot_sed, aes(x = Axis.1, y = Axis.2,
                                    color = month.year, shape = site)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_segment(data = vec_sed,
               aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black",
               linewidth = 0.8) +
  geom_text_repel(data = vec_sed,
                  aes(x = Axis.1, y = Axis.2, label = variable),
                  inherit.aes = FALSE,
                  color = "black",
                  size = 4) +
  scale_colour_manual(values = my_colors) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("PCoA1 (", round(pcoa_sed$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_sed$values$Relative_eig[2] * 100, 1), "%)"),
    color = "month.year",
    shape = "site"
  )

print(p_sed_genus)

ggsave("PCoA_envfit_sediment_genus-level.tiff",
       p_sed_genus,
       width = 6, height = 4, dpi = 600, compression = "lzw")	   
	   
	   
	   
> print(envfit_sed)

***VECTORS

            Axis.1   Axis.2     r2 Pr(>r)   
mean_chla -0.38013  0.92493 0.4688  0.002 **
mean_nit   0.00945  0.99996 0.0873  0.219   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

> print(envfit_sw)

***VECTORS

            Axis.1   Axis.2     r2 Pr(>r)    
mean_chla -0.99095 -0.13425 0.7490  0.001 ***
mean_nit  -0.98888  0.14872 0.4898  0.001 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999
	   