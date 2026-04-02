library(phyloseq)
library(RColorBrewer)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(devtools)
library(forcats)

otu <- read.table("psmin5_ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxon <- read.table("psmin5_taxa_table.txt", sep = "\t", header = TRUE, row.names = 1)
samples <- read.table("metadata.txt", sep = "\t", header = TRUE, row.names = 1)

taxon <- as.matrix(taxon)
TAX <- tax_table(taxon)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(TAX))
  ps
to_remove <- c("LP196","eDNA1")
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

#---------------------
#Figure 4 -- ordination by sample type 
#---------------------

mycolors=c("brown1","darkolivegreen3","cadetblue3","darkmagenta") 

ps %>%
  tax_transform("identity", rank = "Genus")%>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "subtype",shape="site") +
  theme_classic(12) +
  stat_ellipse(aes(color = subtype)) +
  scale_color_manual(values = mycolors)
  
ggsave("Figure4-PCoA-bray-Genus-subtype.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw') 
 
 
#------------------------------- 
# healthy samples only
#-------------------------------
ps_h <- subset_samples(ps, subtype == "healthy")
ntaxa(ps_h)
ps_h_rm0 <- filter_taxa(ps_h, function(x) sum(x) > 0, TRUE)
ntaxa(ps_h_rm0)
sample_data(ps_h_rm0)$month.year <- fct_relevel(sample_data(ps_h_rm0)$month.year,
                                       "Nov.23","Feb.24","May.24","Aug.24")

my_colors2 <- c(
  "Nov.23" = "darkgoldenrod1",
  "Feb.24" = "cornflowerblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4")
  
   
ps_h_rm0 %>%
  tax_transform("identity", rank = "Genus")%>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 3, color = "month.year",shape="site") +
  theme_classic(12)+
  scale_colour_manual(values=my_colors2)


ggsave("Fig5C-PCoA-bray-Genus-healthy.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')  
 

#-------------------------------
# diseased samples only
#-------------------------------
ps_d <- subset_samples(ps,subtype == "diseased") 
ntaxa(ps_d)
ps_d_rm0 <-filter_taxa(ps_d, function(x) sum(x) >0, TRUE)
ntaxa(ps_d_rm0) 

sample_data(ps_d_rm0)$month.year <- fct_relevel(sample_data(ps_d_rm0)$month.year,
                                       "Nov.23","Feb.24","May.24","Aug.24")

my_colors2 <- c(
  "Nov.23" = "darkgoldenrod1",
  "Feb.24" = "cornflowerblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4")


ps_d_rm0 %>%
  tax_transform("identity", rank = "Genus")%>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "month.year",shape="site") +
  theme_classic(12) +
  scale_colour_manual(values=my_colors2)


ggsave("Fig5D-PCoA-bray-Genus-diseased.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw') 



#-------------------------------
# seawater samples only
#-------------------------------

ps_w <- subset_samples(ps,subtype == "seawater") 
ntaxa(ps_w)
ps_w_rm0 <-filter_taxa(ps_w, function(x) sum(x) >0, TRUE)
ntaxa(ps_w_rm0) 
my_colors <- c(
  "Oct.23" = "darkgoldenrod4",
  "Nov.23" = "darkgoldenrod1",
  "Jan.24" = "cornflowerblue",
  "Mar.24" = "darkslateblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4")
month.year <- sample_data(ps_w_rm0)$month.year
sample_data(ps_w_rm0)$month.year <- fct_relevel(sample_data(ps_w_rm0)$month.year,
                                       "Oct.23", "Nov.23", "Jan.24","Mar.24","May.24","Aug.24")

ps_w_rm0 %>%
  tax_transform("identity", rank = "Genus")%>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 3, color = "month.year",shape="site") +
  theme_classic(12) +
  scale_colour_manual(values=my_colors)
  
  
ggsave("Fig5B-PCoA-bray-Genus-seawater.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw') 

#-------------------------------
# sediment samples only
#-------------------------------

ps_s <- subset_samples(ps,subtype == "sediment") 
ntaxa(ps_s)
ps_s_rm0 <-filter_taxa(ps_s, function(x) sum(x) >0, TRUE)
ntaxa(ps_s_rm0) 

sample_data(ps_s_rm0)$month.year <- fct_relevel(sample_data(ps_s_rm0)$month.year,
                                               "Oct.23", "Nov.23", "Jan.24","Mar.24","May.24","Aug.24")
my_colors <- c(
  "Oct.23" = "darkgoldenrod4",
  "Nov.23" = "darkgoldenrod1",
  "Jan.24" = "cornflowerblue",
  "Mar.24" = "darkslateblue",
  "May.24" = "darkseagreen2",
  "Aug.24" = "darkolivegreen4")

ps_s_rm0 %>%
  tax_transform("identity", rank = "Genus")%>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 3, color = "month.year",shape="site") +
  theme_classic(12)+
  scale_colour_manual(values=my_colors)

 ggsave("Figure5A-PCoA-bray-Genus-sediment6x4.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')  
 