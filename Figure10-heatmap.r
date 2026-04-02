library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(cowplot)
library(grid)
library(scales)
library(RColorBrewer)
library(pheatmap)

signgen <- read.table("RA.19.genera.txt", sep = "\t", header = TRUE, row.names = 1)
signgen <- data.matrix(signgen)

# Log-transform relative abundance matrix
# Add a small pseudocount to avoid log10(0) = -Inf
signgen_log <- log10(signgen + 1e-6)

annot <- read.table("Fig10-metadata.txt", sep = "\t", header = TRUE, row.names = 1)
row.names(annot) <- colnames(signgen_log)

ann_colors <- list(
  site = c("Tumon" = "black", "Luminao" = "gray"),
  type = c("healthy" = "limegreen", "diseased" = "brown1"),
  month.year = c(
    "Nov23" = "darkgoldenrod1",
    "Feb24" = "cornflowerblue",
    "May24" = "darkseagreen4",
    "Aug24" = "darkolivegreen4"
  )
)

# TIFF output
tiff("Figure10-heatmap_w_sample_annot.tiff",
     width = 12, height = 8, units = "in", res = 300, compression = "lzw")

pheatmap(
  signgen_log,
  legend = TRUE,
  legend_breaks = c(-6, -4, -2, 0),
  legend_labels = c("-6", "-4", "-2", "0"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D",
  fontsize_row = 12,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annot,
  annotation_colors = ann_colors)

dev.off()

tiff("Figure10-heatmap_no_legend.tiff",
     width = 12, height = 8, units = "in", res = 300, compression = "lzw")

pheatmap(
  signgen_log,
  legend = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D",
  fontsize_row = 12,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annot,
  annotation_colors = ann_colors)

dev.off()





# PDF output
pdf("Figure10-heatmap_w_sample_annot.pdf",
    width = 12, height = 8)

pheatmap(
  signgen_log,
  legend = TRUE,
  legend_breaks = c(-6, -4, -2, 0),
  legend_labels = c("-6", "-4", "-2", "0"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D",
  fontsize_row = 12,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annot,
  annotation_colors = ann_colors)

dev.off()


