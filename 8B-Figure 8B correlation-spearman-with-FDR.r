
library(phyloseq)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(cowplot)
library(grid)
library(scales)
library(RColorBrewer)
library(randomcoloR)
library(microbiome)
library(microViz)
library(tidyr)
library(tibble)

#------------------------------------------
#read in curate diseased samples only data 
#------------------------------------------

otu <- read.table("disease_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("disease_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("disease_metadata.txt",sep="\t",header=TRUE,row.names=1)


taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(TAX))

ntaxa(ps)

#-------------------------------------------
#recode metadata into numerical 
#-------------------------------------------
psQ <- ps %>%
  ps_mutate(
    site = if_else(site == "Tumon", true = 1, false = 0),
    month = recode(month, Nov = 1, Feb = 2, May = 3, Aug = 4),
	colony = recode(colony, A = 1, B = 2, C = 3, D = 4, E=5, F=6, G=7, H=8),
    lesion.size = as.numeric(lesion.size),
    lesion.progr.rate = as.numeric(progr.rate))

#----------------------------------------------    
#Agglomerate to genus level and transforme to RA
#----------------------------------------------

psQG <- tax_agg(psQ, "Genus")
ntaxa(psQG)
ra_ps <- microbiome::transform(psQG, "compositional")
ra_mat <- otu_table(ra_ps) %>% as.data.frame()
head(ra_mat) 


#--------------------------------------------------
# Extract metadata of interest
#--------------------------------------------------
meta <- sample_data(psQG)[, c("lesion.size", "lesion.progr.rate")] %>% as.data.frame()
head(meta)

#---------------------------------------------------
#Calculate spearman correlation 
#---------------------------------------------------

# Step 1: define helper function
get_spearman_stats <- function(x, y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(rho = unname(test$estimate), p = test$p.value)
}

# Step 2: run Spearman tests for lesion variables
results <- lapply(colnames(meta), function(var) {
  res <- t(apply(ra_mat, 2, function(genus) {
    get_spearman_stats(genus, meta[[var]])
  }))
  
  res <- as.data.frame(res)
  res$genus <- rownames(res)
  res$variable <- var
  res
})

results_df <- do.call(rbind, results)

# Step 3: adjust p-values for multiple testing
results_df$p_adj <- p.adjust(results_df$p, method = "BH")

# Look at first rows
head(results_df)

#------------------------------------------------------
#Prepare martix for plotting 
#------------------------------------------------------

cor_mat <- results_df %>%
  select(genus, variable, rho) %>%
  pivot_wider(names_from = variable, values_from = rho) %>%
  column_to_rownames("genus") %>%
  as.matrix()

p_adj_mat <- results_df %>%
  select(genus, variable, p_adj) %>%
  pivot_wider(names_from = variable, values_from = p_adj) %>%
  column_to_rownames("genus") %>%
  as.matrix()

head(cor_mat)
head(p_adj_mat)


sel_taxa <- rownames(cor_mat)[
  apply((abs(cor_mat) > 0.3) & (p_adj_mat < 0.05), 1, any)
]

length(sel_taxa)
sel_taxa

#---------------------------------------------------------
#plot heatmap for lesion size and lesion progression 
#---------------------------------------------------------
set.seed(123)

tiff("spearman-genera-vs-lesion-size-progr-cor-0.3-to-0.3.tiff",
     units="cm", width=16, height=16, res=300)

cor_heatmap(
  data = psQG,
  taxa = sel_taxa,
  cor = "spearman",
  vars = c("lesion.size", "lesion.progr.rate"),
  colors = heat_palette("Blue-Red 2", rev = FALSE, sym = TRUE),
  grid_lwd = 3)

dev.off()



#plot heatmap for other variables 
			
set.seed(123)

tiff("spearman-genera-vs-site-month-colony-cor-0.3-to-0.3.tiff",
     units="cm", width=16, height=16, res=300)

cor_heatmap(
  data = psQG,
  taxa = sel_taxa2,
  cor = "spearman",
  vars = c("site", "month","colony"),
  colors = heat_palette("Blue-Red 2", rev = FALSE, sym = TRUE),
  grid_lwd = 3)


dev.off()

#pdf saving

pdf("spearman-genera-vs-lesion-size-progr.pdf",
    width = 8, height = 8)  # inches (not cm!)

cor_heatmap(
  data = psQG,
  taxa = sel_taxa,
  cor = "spearman",
  vars = c("lesion.size", "lesion.progr.rate"),
  colors = heat_palette("Blue-Red 2", rev = FALSE, sym = TRUE),
  grid_lwd = 3)

dev.off()



pdf("spearman-genera-vs-site-month-colony.pdf",
    width = 8, height = 8)  # inches (not cm!)

cor_heatmap(
  data = psQG,
  taxa = sel_taxa2,
  cor = "spearman",
  vars = c("site", "month","colony"),
  colors = heat_palette("Blue-Red 2", rev = FALSE, sym = TRUE),
  grid_lwd = 3)

dev.off()

#-----------------------------------------------------------------
#alternative plotting options with Complex heat maps 
#-----------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(grid)

# subset matrices to selected taxa
cor_use <- cor_mat[sel_taxa, , drop = FALSE]
p_use   <- p_adj_mat[sel_taxa, , drop = FALSE]

# matrix of significance symbols
sig_mat <- ifelse(p_use < 0.05, "*", "")

# optional: round rho values for display
rho_labels <- matrix(sprintf("%.2f", cor_use),
                     nrow = nrow(cor_use),
                     ncol = ncol(cor_use),
                     dimnames = dimnames(cor_use))

pdf("Complex-spearman-genera-vs-lesion-size-progr-rho-labeled.pdf",
    width = 8, height = 8)

Heatmap(
  cor_use,
  name = "rho",
  col = colorRamp2(c(-1, 0, 1),c("#1B7837", "white", "#2166AC")),

  row_names_gp = gpar(fontsize = 12,fontface = "italic"),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "black", lwd = 1),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)),

  cell_fun = function(j, i, x, y, width, height, fill) {
    # rho value in center
    grid.text(rho_labels[i, j], x, y,
              gp = gpar(fontsize = 9))
    
    # asterisk slightly above center for significant cells
    if (sig_mat[i, j] == "*") {
      grid.text("*", x, y + unit(2.5, "mm"),
                gp = gpar(fontsize = 10, fontface = "bold"))
    }
  }
)

dev.off()

set.seed(123)

tiff("Complex-spearman-genera-vs-lesion-size-progr-rho-labeled.tiff",
     units="in", width=8, height=8, res=300)

Heatmap(
  cor_use,
  name = "rho",
  col = colorRamp2(c(-1, 0, 1),c("#1B7837", "white", "#2166AC")),

  row_names_gp = gpar(fontsize = 12,fontface = "italic"),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "black", lwd = 1),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)),

  cell_fun = function(j, i, x, y, width, height, fill) {
    # rho value in center
    grid.text(rho_labels[i, j], x, y,
              gp = gpar(fontsize = 9))
    
    # asterisk slightly above center for significant cells
    if (sig_mat[i, j] == "*") {
      grid.text("*", x, y + unit(2.5, "mm"),
                gp = gpar(fontsize = 10, fontface = "bold"))
    }
  }
)

dev.off()