# ============================================================
# Purpose:
# 1. Read OTU table, taxonomy table, and metadata from ONE folder
# 2. Calculate alpha diversity metrics
# 3. Run mixed models for coral samples (healthy vs diseased)
# 4. Run linear models for seawater and sediment
#
# Put this script in the SAME folder as:
#   psmin5_ASV_table.txt
#   psmin5_taxa_table.txt
#   metadata.txt
# ============================================================

library(phyloseq)
library(microbiome)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

# ----------------------------
# 1. Read input files
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

# Remove the same samples you removed in the older script
to_remove <- c("LP196", "eDNA1")
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

print(ps)
print(colnames(sample_data(ps)))

# ----------------------------
# 2. Calculate alpha diversity
# ----------------------------
alpha_div <- microbiome::alpha(ps, index = "all")
ps.meta <- meta(ps)

ps.meta$Shannon <- alpha_div$diversity_shannon
ps.meta$InverseSimpson <- alpha_div$diversity_inverse_simpson
ps.meta$Richness <- alpha_div$observed


# ----------------------------
# 3. Define metadata variables
# Updated metadata:
# type = healthy / diseased / seawater / sediment
# ----------------------------
ps.meta$site <- factor(ps.meta$site, levels = c("Luminao", "Tumon"))

time_levels <- c("Oct.23", "Nov.23", "Dec.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24", "Sep.24")
time_levels_present <- time_levels[time_levels %in% unique(ps.meta$month.year)]


ps.meta$timepoint <- factor(ps.meta$month.year, levels = time_levels_present)

ps.meta$colony <- factor(ps.meta$colony.name)
ps.meta$type <- factor(ps.meta$type, levels = c("healthy", "diseased", "seawater", "sediment"))

# Save metadata + alpha diversity values
write.table(
  ps.meta,
  file = "alpha_diversity_metadata_table.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# ----------------------------
# 4. Mixed models for coral samples only
# ----------------------------
coral_data <- ps.meta %>%
  filter(type %in% c("healthy", "diseased")) %>%
  filter(!is.na(site), !is.na(timepoint), !is.na(colony)) %>%
  mutate(
    site = factor(site, levels = c("Luminao", "Tumon")),
    type = factor(type, levels = c("healthy", "diseased")),
    timepoint = factor(timepoint, levels = time_levels_present),
    colony = factor(colony)
  )

run_coral_mixed_model <- function(df, response_var) {

  df_sub <- df %>%
    filter(!is.na(.data[[response_var]]))

  model <- lmer(
    as.formula(paste0(response_var, " ~ site * type + timepoint + (1 | colony)")),
    data = df_sub,
    REML = FALSE
  )

  anova_tbl <- as.data.frame(anova(model))
  anova_tbl$term <- rownames(anova_tbl)
  rownames(anova_tbl) <- NULL

  emm_site <- emmeans(model, pairwise ~ site | type)
  emm_type <- emmeans(model, pairwise ~ type | site)

  list(
    model = model,
    anova = anova_tbl,
    site_emmeans = as.data.frame(emm_site$emmeans),
    site_contrasts = as.data.frame(emm_site$contrasts),
    type_emmeans = as.data.frame(emm_type$emmeans),
    type_contrasts = as.data.frame(emm_type$contrasts)
  )
}

metrics <- c("Shannon", "InverseSimpson", "Richness")
coral_results <- list()

all_coral_anova <- data.frame()
all_coral_site_contrasts <- data.frame()
all_coral_type_contrasts <- data.frame()

for (response_var in metrics) {
  
  coral_results[[response_var]] <- run_coral_mixed_model(coral_data, response_var)
  
  # Add metric name to each output
  anova_tbl <- coral_results[[response_var]]$anova
  anova_tbl$metric <- response_var
  
  site_contrasts_tbl <- coral_results[[response_var]]$site_contrasts
  site_contrasts_tbl$metric <- response_var
  
  type_contrasts_tbl <- coral_results[[response_var]]$type_contrasts
  type_contrasts_tbl$metric <- response_var
  
  # Combine across metrics
  all_coral_anova <- rbind(all_coral_anova, anova_tbl)
  all_coral_site_contrasts <- rbind(all_coral_site_contrasts, site_contrasts_tbl)
  all_coral_type_contrasts <- rbind(all_coral_type_contrasts, type_contrasts_tbl)
}

# Reorder columns so metric comes first
all_coral_anova <- all_coral_anova[, c("metric", setdiff(names(all_coral_anova), "metric"))]
all_coral_site_contrasts <- all_coral_site_contrasts[, c("metric", setdiff(names(all_coral_site_contrasts), "metric"))]
all_coral_type_contrasts <- all_coral_type_contrasts[, c("metric", setdiff(names(all_coral_type_contrasts), "metric"))]

# Display in R
print(all_coral_anova)
print(all_coral_site_contrasts)
print(all_coral_type_contrasts)

# Optional: open in Viewer
View(all_coral_anova)
View(all_coral_site_contrasts)
View(all_coral_type_contrasts)

# Save only 3 files instead of many
write.table(
  all_coral_anova,
  file = "coral_alpha_anova_all.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  all_coral_site_contrasts,
  file = "coral_alpha_site_contrasts_all.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  all_coral_type_contrasts,
  file = "coral_alpha_type_contrasts_all.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

#alpha divesity differences in seawater and sediment 
