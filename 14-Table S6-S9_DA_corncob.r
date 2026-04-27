library(microViz)
library(microbiome)
library(devtools)
library(tidyverse)
library(corncob)

otu <- read.table(input_otu, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
taxon <- read.table(input_tax, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(as.matrix(taxon)))

ps

to_remove <- c("LP196","eDNA1")

ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

ps



#differential abundance between sample types: seawater, sediment, diseased tissue and healthy tissue 

ps_order <- ps %>% tax_glom("Order")
sample_data(ps_order)$month.year <- factor(sample_data(ps_order)$month.year)

da_order_type <- differentialTest(formula = ~ subtype+site+month.year,
											phi.formula = ~ 1,
											formula_null = ~ site+month.year,
											phi.formula_null = ~ 1,
											test = "Wald", boot = FALSE,
											data = ps_order,
											fdr_cutoff = 0.05)

plot(da_order_type,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")
ggsave("btw-types-eff-site.tiff", units="in", width=12, height=8, dpi=300, compression = 'lzw')
da_order_type$significant_taxa 



#testing this as an option for differece in sites without spliting the data ------------------------

ps_genus <- ps %>% tax_glom("Genus")
sample_data(ps_genus)$month.year <- factor(sample_data(ps_genus)$month.year)

da_genus_site <- differentialTest(formula = ~ site+type+month.year,
											phi.formula = ~ 1,
											formula_null = ~ type+month.year,
											phi.formula_null = ~ 1,
											test = "Wald", boot = FALSE,
											data = ps_order,
											fdr_cutoff = 0.05)

plot(da_genus_site,level="Genus")+labs(x = "Coefficient and Confidence Interval", y = "Genus")
ggsave("genera-btw-sites-eff-site.tiff", units="in", width=12, height=8, dpi=300, compression = 'lzw')
da_genera_site$significant_taxa 




#--------------------------------------------------------
#Supplementary Table S6
#--------------------------------------------------------

library(phyloseq)
library(dplyr)
library(tidyr)

# Agglomerate at Order level
ps_order <- tax_glom(ps, taxrank = "Order")

# Transform to relative abundance
ps_order_rel <- transform_sample_counts(ps_order, function(x) x / sum(x))

df_order <- psmelt(ps_order_rel)
da_orders <- c("Pseudomonadales", "Rhodobacterales", "Enterobacterales",
               "Cytophagales", "Caulobacterales", "Flavobacteriales",
               "Verrucomicrobiales", "Desulfobacterales", "Cyanobacteriales")

df_da <- df_order %>%
  filter(Order %in% da_orders)

da_orders <- c("Pseudomonadales", "Rhodobacterales", "Enterobacterales",
               "Cytophagales", "Caulobacterales", "Flavobacteriales",
               "Verrucomicrobiales", "Desulfobacterales", "Cyanobacteriales")

df_da <- df_order %>%
  filter(Order %in% da_orders)
  
stats_overall <- df_da %>%
  group_by(Order) %>%
  summarise(
    mean_RA = mean(Abundance),
    min_RA  = min(Abundance),
    max_RA  = max(Abundance),
    sd_RA   = sd(Abundance),
    prevalence = mean(Abundance > 0) * 100
  ) %>%
  arrange(desc(mean_RA))


stats_by_type <- df_da %>%
  group_by(Order, subtype) %>%
  summarise(
    mean_RA = mean(Abundance),
    min_RA  = min(Abundance),
    max_RA  = max(Abundance),
    sd_RA   = sd(Abundance),
    prevalence = mean(Abundance > 0) * 100,
    .groups = "drop"
  )
  
 stats_full <- df_da %>%
  group_by(Order, subtype, site, month.year) %>%
  summarise(
    mean_RA = mean(Abundance),
    min_RA  = min(Abundance),
    max_RA  = max(Abundance),
    prevalence = mean(Abundance > 0) * 100,
    .groups = "drop"
  )


write.csv(stats_by_type, "Table_S6_DA_order_summary.csv", row.names = FALSE)


#-----------------------------------------------------------------------
#Seawater: Differential abundance between sites and time points 
#-----------------------------------------------------------------------
ps_w <- subset_samples(ps,subtype == "seawater") 
ntaxa(ps_w)
ps_w_rm0 <-filter_taxa(ps_w, function(x) sum(x) >0, TRUE)
ntaxa(ps_w_rm0) 
ps_order_w <- ps_w_rm0 %>% tax_glom("Order")
da_order_site_w <- differentialTest(formula = ~ site+month,
											phi.formula = ~ 1,
											formula_null = ~ month,
											phi.formula_null = ~ 1,
											test = "LRT", boot = FALSE,
											data = ps_order_w,
											fdr_cutoff = 0.05)

plot(da_order_site_w,level="Order")
ggsave("seawater-orders-btw-site-eff-month.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')
da_order_site_w$significant_taxa 

da_order_sw_site_FDR <- da_order_site_w$p_fdr
write.table(da_order_sw_site_FDR, "da_order_sw_site_FDR.tsv", sep = "\t", row.names = TRUE)

da_order_month_w <- differentialTest(formula = ~ month+site,
											phi.formula = ~ 1,
											formula_null = ~ site,
											phi.formula_null = ~ 1,
											test = "LRT", boot = FALSE,
											data = ps_order_w,
											fdr_cutoff = 0.05)

plot(da_order_month_w,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")+ 
theme(panel.grid = element_blank ())
ggsave("seawater-orders-month-eff-site.tiff", units="in", width=16, height=8, dpi=300, compression = 'lzw')

da_order_month_w$significant_taxa 
da_order_sw_month_FDR <- da_order_month_w$p_fdr
write.table(da_order_sw_month_FDR, "da_order_sw_month_FDR.tsv", sep = "\t", row.names = TRUE)

#-------------------------------------------------------------------------------
#Sediment: Differentially abundance between sites and across time points 
#-------------------------------------------------------------------------------
ps_s <- subset_samples(ps,subtype == "sediment") 
ntaxa(ps_s)
ps_s_rm0 <-filter_taxa(ps_s, function(x) sum(x) >0, TRUE)
ntaxa(ps_s_rm0) 
ps_order_s <- ps_s_rm0 %>% tax_glom("Order")
ps_order_s_ra <- transform_sample_counts(ps_order_s, function(otu) otu/sum(otu))
write.table(as(otu_table(ps_order_s_ra), "matrix"),"ps_s_order_ra_table.txt",sep="\t",col.names=NA)
da_order_s_s <- differentialTest(formula = ~ site+month,
                                 phi.formula = ~ 1,
                                 formula_null = ~ month,
                                 phi.formula_null = ~ 1,
                                 data = ps_order_s,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)	
da_order_s_s$significant_taxa 
da_order_s_m <- differentialTest(formula = ~ month+site,
                                 phi.formula = ~ 1,
                                 formula_null = ~ site,
                                 phi.formula_null = ~ 1,
                                 data = ps_order_s,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)								 
da_order_s_m$significant_taxa 
da_order_s_s_FDR_ <- da_order_s_s$p_fdr
write.table(da_order_FDR_m, "da_s_s_orders_FDR.tsv", sep = "\t", row.names = TRUE)


da_order_s_m_FDR_ <- da_order_s_m$p_fdr
write.table(da_order_FDR_m, "da_s_m_orders_FDR.tsv", sep = "\t", row.names = TRUE)
						 
#----------------------------------------------------------------------------
#Diseased tissue: Differential abundance between sites and across time points 
#----------------------------------------------------------------------------
ps_d <- subset_samples(ps,subtype == "diseased") 
ntaxa(ps_d)
ps_d_rm0 <-filter_taxa(ps_d, function(x) sum(x) >0, TRUE)
ntaxa(ps_d_rm0) 
ps_order_d <- ps_d_rm0 %>% tax_glom("Order")						 
ps_order_d_ra <- transform_sample_counts(ps_order_d, function(otu) otu/sum(otu))
write.table(as(otu_table(ps_order_d_ra), "matrix"),"ps_d_order_ra_table.txt",sep="\t",col.names=NA)
da_order_d_s <- differentialTest(formula = ~ site+month,
                                 phi.formula = ~ 1,
                                 formula_null = ~ month,
                                 phi.formula_null = ~ 1,
                                 data = ps_order_d,test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)	
plot(da_order_d_s,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")+ 
theme(panel.grid = element_blank ())
da_order_d_s$significant_taxa 	
da_WS_order_FDR_s <- da_order_d_s$p_fdr
write.table(da_WS_order_FDR_s, "da_WS_order_FDR_s.tsv", sep = "\t", row.names = TRUE)
da_order_d_m <- differentialTest(formula = ~ month+site,
                                 phi.formula = ~ 1,
                                 formula_null = ~ site,
                                 phi.formula_null = ~ 1,
                                 data = ps_order_d,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)	
plot(da_order_d_m,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")+ theme(panel.grid = element_blank ())
ggsave("diseased-orders-per-month-e-site.tiff", units="in", width=10, height=6, dpi=300, compression = 'lzw')
da_order_FDR_m <- da_order_d_m$p_fdr
write.table(da_order_FDR_m, "da_WS_month_orders_FDR.tsv", sep = "\t", row.names = TRUE)
