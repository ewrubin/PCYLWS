library(phyloseq)
library(RColorBrewer)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(devtools)
library(pairwiseAdonis)
library(vegan)

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

ps

# Make metadata variables factors
sample_data(ps)$subtype <- factor(sample_data(ps)$subtype)
sample_data(ps)$site <- factor(sample_data(ps)$site)
sample_data(ps)$month.year <- factor(sample_data(ps)$month.year)

# -----------------------------
# Genus-level object
# -----------------------------
ps_genus <- tax_agg(ps, rank = "Genus")

# Bray-Curtis distance matrix at genus level
bray_dist <- phyloseq::distance(ps_genus, method = "bray")

meta <- data.frame(sample_data(ps_genus))

# -----------------------------
# PERMANOVA
# -----------------------------

adonis_full_strata <- adonis2(bray_dist ~ type + site + month.year,data = meta,permutations = 999,strata = meta$colony.name,by = "margin")
print(adonis_full_strata)

# -----------------------------
# Beta-dispersion
# -----------------------------
# Dispersion among sample subtypes
bd_subtype <- betadisper(bray_dist, group = meta$subtype)
anova(bd_subtype)
permutest(bd_subtype, permutations = 999)

# Optional pairwise dispersion comparisons
TukeyHSD(bd_subtype)


# -----------------------------
# Pairwise PERMANOVA
# -----------------------------
pairwise_results <- pairwise.adonis2(
  bray_dist ~ type,
  data = meta,
  p.adjust.m = "bonferroni")
print(pairwise_results)

-------------------------
#------------------------------- 
# healthy samples only
#-------------------------------
ps_h <- subset_samples(ps, subtype == "healthy")

ntaxa(ps_h)

ps_h_rm0 <- filter_taxa(ps_h, function(x) sum(x) > 0, TRUE)

ntaxa(ps_h_rm0)

sample_data(ps_h_rm0)$site <- factor(sample_data(ps_h_rm0)$site)
sample_data(ps_h_rm0)$month.year <- factor(sample_data(ps_h_rm0)$month.year)

ps_genus_h <- tax_agg(ps_h_rm0, rank = "Genus")
bray_dist_h <- phyloseq::distance(ps_genus_h, method = "bray")
meta_h <- data.frame(sample_data(ps_genus_h))

adonis_h_site_time <- adonis2(bray_dist_h ~ site + month.year, data = meta_h, permutations = 999,by = "margin")
print(adonis_h_site)

bd_site_h <- betadisper(bray_dist_h, group = meta_h$site)
anova(bd_site_h)
permutest(bd_site_h, permutations = 999)

bd_time_h <- betadisper(bray_dist_h, group = meta_h$month.year)
anova(bd_time_h)
permutest(bd_time_h, permutations = 999)


pairwise_h_site <- pairwise.adonis2(bray_dist_h ~ site, data = meta_h, p.adjust.m = "bonferroni")
print(pairwise_h_site)
pairwise_h_time <- pairwise.adonis2(bray_dist_h ~ month.year, data = meta_h, p.adjust.m = "bonferroni")
print(pairwise_h_time)



#-------------------------------
# diseased samples only
#-------------------------------

ps_d <- subset_samples(ps, subtype == "diseased")

ntaxa(ps_d)

ps_d_rm0 <- filter_taxa(ps_d, function(x) sum(x) > 0, TRUE)

ntaxa(ps_d_rm0)

sample_data(ps_d_rm0)$site <- factor(sample_data(ps_d_rm0)$site)
sample_data(ps_d_rm0)$month.year <- factor(sample_data(ps_d_rm0)$month.year)
sample_data(ps_d_rm0)$colony.name <- factor(sample_data(ps_d_rm0)$colony.name)

ps_genus_d <- tax_agg(ps_d_rm0, rank = "Genus")
bray_dist_d <- phyloseq::distance(ps_genus_d, method = "bray")
meta_d <- data.frame(sample_data(ps_genus_d))

# make unique colony ID if colony names are site-specific
meta_d$colony_id <- interaction(meta_d$site, meta_d$colony.name, drop = TRUE)

# test site while accounting for repeated colony sampling
adonis_d_site_time <- adonis2(bray_dist_d ~ site + month.year,data = meta_d,permutations = 999,strata = meta_d$colony_id,by = "margin")
print(adonis_d_site_time)

# beta-dispersion by site
bd_site_d <- betadisper(bray_dist_d, group = meta_d$site)
anova(bd_site_d)
permutest(bd_site_d, permutations = 999)

# beta-dispersion by time
bd_time_d <- betadisper(bray_dist_d, group = meta_d$month.year)
anova(bd_time_d)
permutest(bd_time_d, permutations = 999)

pairwise_d_site <- pairwise.adonis2(bray_dist_d ~ site, data = meta_d, p.adjust.m = "bonferroni")
print(pairwise_d_site)
pairwise_d_time <- pairwise.adonis2(bray_dist_d ~ month.year, data = meta_d, p.adjust.m = "bonferroni")
print(pairwise_d_time)


#-------------------------------
# seawater samples only
#-------------------------------

ps_w <- subset_samples(ps,subtype == "seawater") 

ntaxa(ps_w)

ps_w_rm0 <-filter_taxa(ps_w, function(x) sum(x) >0, TRUE)

ntaxa(ps_w_rm0) 

sample_data(ps_w_rm0)$site <- factor(sample_data(ps_w_rm0)$site)
sample_data(ps_w_rm0)$month.year <- factor(sample_data(ps_w_rm0)$month.year)

ps_genus_w <- tax_agg(ps_w_rm0, rank = "Genus")
bray_dist_w <- phyloseq::distance(ps_genus_w, method = "bray")
meta_w <- data.frame(sample_data(ps_genus_w))

adonis_w_site_time <- adonis2(bray_dist_w ~ site + month.year, data = meta_w, permutations = 999,by = "margin")
print(adonis_w_site_time)


bd_site_w <- betadisper(bray_dist_w, group = meta_w$site)
anova(bd_site_w)
permutest(bd_site_w, permutations = 999)


bd_time_w <- betadisper(bray_dist_w, group = meta_w$month.year)
anova(bd_time_w)
permutest(bd_time_w, permutations = 999)


pairwise_w_site <- pairwise.adonis2(bray_dist_w ~ site, data = meta_w, p.adjust.m = "bonferroni")
print(pairwise_w_site)
pairwise_w_time <- pairwise.adonis2(bray_dist_w ~ month.year, data = meta_w, p.adjust.m = "bonferroni")
print(pairwise_w_time)



#-------------------------------
# sediment samples only
#-------------------------------

ps_s <- subset_samples(ps,subtype == "sediment") 

ntaxa(ps_s)

ps_s_rm0 <-filter_taxa(ps_s, function(x) sum(x) >0, TRUE)

ntaxa(ps_s_rm0) 

sample_data(ps_s_rm0)$site <- factor(sample_data(ps_s_rm0)$site)
sample_data(ps_s_rm0)$month.year <- factor(sample_data(ps_s_rm0)$month.year)

ps_genus_s <- tax_agg(ps_s_rm0, rank = "Genus")
bray_dist_s <- phyloseq::distance(ps_genus_s, method = "bray")
meta_s <- data.frame(sample_data(ps_genus_s))

adonis_s_site_time <- adonis2(bray_dist_s ~ site + month.year, data = meta_s, permutations = 999,by = "margin")
print(adonis_s_site_time)



bd_site_s <- betadisper(bray_dist_s, group = meta_s$site)
anova(bd_site_s)
permutest(bd_site_s, permutations = 999)


bd_time_s <- betadisper(bray_dist_s, group = meta_s$month.year)
anova(bd_time_s)
permutest(bd_time_s, permutations = 999)

pairwise_s_site <- pairwise.adonis2(bray_dist_s ~ site, data = meta_s, p.adjust.m = "bonferroni")
print(pairwise_s_site)
pairwise_s_time <- pairwise.adonis2(bray_dist_s ~ month.year, data = meta_s, p.adjust.m = "bonferroni")
print(pairwise_s_time)