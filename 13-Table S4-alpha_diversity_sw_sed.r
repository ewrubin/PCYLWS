

library(phyloseq)
library(microbiome)
library(dplyr)
library(ARTool)
library(emmeans)
library(car)

# ----------------------------
# 1. Read files and make phyloseq object
# ----------------------------
otu <- read.table("psmin5_ASV_table.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
taxon <- read.table("psmin5_taxa_table.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table("metadata.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

taxon <- as.matrix(taxon)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(taxon)
)

# optional: remove same samples as in your other scripts
to_remove <- c("LP196", "eDNA1")
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

# ----------------------------
# 2. Calculate alpha diversity
# ----------------------------
alpha_div <- microbiome::alpha(ps, index = "all")
ps.meta <- meta(ps)

ps.meta$Shannon <- alpha_div$diversity_shannon
ps.meta$InverseSimpson <- alpha_div$diversity_inverse_simpson
ps.meta$Richness <- alpha_div$observed


# ----------------------------
# 3. Define factors
# ----------------------------
ps.meta$site <- factor(ps.meta$site, levels = c("Luminao", "Tumon"))
ps.meta$month.year <- factor(
  ps.meta$month.year,
  levels = c("Oct.23", "Nov.23", "Dec.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24", "Sep.24")
)
ps.meta$type <- factor(ps.meta$type, levels = c("healthy", "diseased", "seawater", "sediment"))

# ----------------------------
# 4. Function to run ART ANOVA
# ----------------------------
run_art_alpha <- function(df, response_var) {
  
  df_sub <- df %>%
    filter(!is.na(.data[[response_var]]),
           !is.na(site),
           !is.na(month.year)) %>%
    mutate(
      site = droplevels(site),
      month.year = droplevels(month.year)
    )
  
  # optional assumption checks on raw data
  print(paste("Checking", response_var))
  print(shapiro.test(df_sub[[response_var]]))
  print(leveneTest(as.formula(paste(response_var, "~ interaction(site, month.year)")), data = df_sub))
  
  # ART ANOVA
  model <- art(
    as.formula(paste(response_var, "~ site * month.year")),
    data = df_sub
  )
  
  anova_tbl <- as.data.frame(anova(model))
  anova_tbl$term <- rownames(anova_tbl)
  rownames(anova_tbl) <- NULL
  
  # post hoc: site differences within each timepoint
  emm_site <- emmeans(artlm(model, "site:month.year"), ~ site | month.year)
  site_contrasts <- as.data.frame(pairs(emm_site, adjust = "holm"))
  
  # post hoc: timepoint differences within each site
  emm_time <- emmeans(artlm(model, "site:month.year"), ~ month.year | site)
  time_contrasts <- as.data.frame(pairs(emm_time, adjust = "holm"))
  
  list(
    model = model,
    anova = anova_tbl,
    site_contrasts = site_contrasts,
    time_contrasts = time_contrasts
  )
}

# ----------------------------
# 5. Subset seawater and sediment
# ----------------------------
seawater_data <- ps.meta %>%
  filter(type == "seawater")

sediment_data <- ps.meta %>%
  filter(type == "sediment")

# ----------------------------
# 6. Run analyses for seawater
# ----------------------------
sw_shannon <- run_art_alpha(seawater_data, "Shannon")
sw_simpson <- run_art_alpha(seawater_data, "InverseSimpson")
sw_richness <- run_art_alpha(seawater_data, "Richness")

# ----------------------------
# 7. Run analyses for sediment
# ----------------------------
sed_shannon <- run_art_alpha(sediment_data, "Shannon")
sed_simpson <- run_art_alpha(sediment_data, "InverseSimpson")
sed_richness <- run_art_alpha(sediment_data, "Richness")

# ----------------------------
# 8. Print results
# ----------------------------
sw_shannon$anova
sw_simpson$anova
sw_richness$anova

sed_shannon$anova
sed_simpson$anova
sed_richness$anova

# Example: inspect post hoc results
sw_shannon$site_contrasts
sw_shannon$time_contrasts

sed_shannon$site_contrasts
sed_shannon$time_contrasts

# ----------------------------
# 9. Optional: save combined ANOVA tables
# ----------------------------
combine_anova <- function(result_list, metric_names, sample_type) {
  out <- data.frame()
  for (i in seq_along(result_list)) {
    tbl <- result_list[[i]]$anova
    tbl$metric <- metric_names[i]
    tbl$sample_type <- sample_type
    out <- rbind(out, tbl)
  }
  out[, c("sample_type", "metric", setdiff(names(out), c("sample_type", "metric")))]
}

sw_anova_all <- combine_anova(
  list(sw_shannon, sw_simpson, sw_richness),
  c("Shannon", "InverseSimpson", "Richness"),
  "seawater"
)

sed_anova_all <- combine_anova(
  list(sed_shannon, sed_simpson, sed_richness),
  c("Shannon", "InverseSimpson", "Richness"),
  "sediment"
)

alpha_art_anova_all <- rbind(sw_anova_all, sed_anova_all)

write.table(
  alpha_art_anova_all,
  file = "alpha_diversity_ART_ANOVA_seawater_sediment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)