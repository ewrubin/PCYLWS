# ============================================================
# 02_alpha_diversity_mixed_models_metadata_v2.R
# Purpose: Calculate alpha diversity metrics from a phyloseq
# object, generate boxplots, and test for site differences using
# linear mixed-effects models for repeated coral-colony samples
# and linear models for seawater and sediment samples.
#
# Updated to match metadata-v2.txt:
#   site
#   month.year
#   type        -> coral / eDNA
#   subtype     -> healthy / diseased / seawater / sediment
#   colony.name -> colony tag number
# ============================================================

library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

input_otu <- "output/psmin5_ASV_table.txt"
input_tax <- "output/psmin5_taxa_table.txt"
metadata_file <- "data/metadata/metadata-v2.txt"
output_dir <- "output/alpha_diversity_mixed_models"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

otu <- read.table(input_otu, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
taxon <- read.table(input_tax, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(as.matrix(taxon))
)

alpha_div <- microbiome::alpha(ps, index = "all")
ps.meta <- meta(ps)

ps.meta$Shannon <- alpha_div$diversity_shannon
ps.meta$InverseSimpson <- alpha_div$diversity_inverse_simpson
ps.meta$Richness <- alpha_div$observed
ps.meta$Chao1 <- alpha_div$chao1

ps.meta$timepoint <- as.factor(ps.meta$month.year)
ps.meta$colony <- as.factor(ps.meta$colony.name)
ps.meta$site <- as.factor(ps.meta$site)
ps.meta$type <- as.factor(ps.meta$type)
ps.meta$subtype <- as.factor(ps.meta$subtype)

cat_levels <- c("L-H", "T-H", "L-D", "T-D", "L-W", "T-W", "L-S", "T-S")
cat_labels <- c(
  "L-H: Luminao-healthy",
  "T-H: Tumon-healthy",
  "L-D: Luminao-diseased",
  "T-D: Tumon-diseased",
  "L-W: Luminao-seawater",
  "T-W: Tumon-seawater",
  "L-S: Luminao-sediment",
  "T-S: Tumon-sediment"
)
mycolors <- c(
  "L-H" = "#1b9e77",
  "T-H" = "#d95f02",
  "L-D" = "#7570b3",
  "T-D" = "#e7298a",
  "L-W" = "#66a61e",
  "T-W" = "#e6ab02",
  "L-S" = "#a6761d",
  "T-S" = "#666666"
)

ps.meta$site_short <- ifelse(ps.meta$site == "Luminao", "L",
                             ifelse(ps.meta$site == "Tumon", "T", as.character(ps.meta$site)))
ps.meta$subtype_short <- dplyr::case_when(
  ps.meta$subtype == "healthy" ~ "H",
  ps.meta$subtype == "diseased" ~ "D",
  ps.meta$subtype == "seawater" ~ "W",
  ps.meta$subtype == "sediment" ~ "S",
  TRUE ~ as.character(ps.meta$subtype)
)
ps.meta$cat <- factor(
  paste(ps.meta$site_short, ps.meta$subtype_short, sep = "-"),
  levels = cat_levels
)

write.table(
  ps.meta,
  file = file.path(output_dir, "alpha_diversity_metadata_table.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)

base_plot_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

plot_shannon <- ggplot(ps.meta, aes(x = cat, y = Shannon, color = cat)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "Shannon diversity") +
  scale_colour_manual(values = mycolors, name = "Site-type", breaks = cat_levels, labels = cat_labels) +
  base_plot_theme
ggsave(file.path(output_dir, "shannon_diversity_boxplot.tiff"), plot_shannon,
       units = "in", width = 6, height = 3.5, dpi = 300, compression = "lzw")

plot_inverse_simpson <- ggplot(ps.meta, aes(x = cat, y = InverseSimpson, color = cat)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "Inverse Simpson diversity") +
  scale_colour_manual(values = mycolors, name = "Site-type", breaks = cat_levels, labels = cat_labels) +
  base_plot_theme
ggsave(file.path(output_dir, "inverse_simpson_diversity_boxplot.tiff"), plot_inverse_simpson,
       units = "in", width = 6, height = 3.5, dpi = 300, compression = "lzw")

plot_richness <- ggplot(ps.meta, aes(x = cat, y = Richness, color = cat)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "ASV richness") +
  scale_colour_manual(values = mycolors, name = "Site-type", breaks = cat_levels, labels = cat_labels) +
  base_plot_theme
ggsave(file.path(output_dir, "richness_boxplot.tiff"), plot_richness,
       units = "in", width = 6, height = 3.5, dpi = 300, compression = "lzw")

coral_data <- ps.meta %>%
  filter(type == "coral") %>%
  filter(subtype %in% c("healthy", "diseased")) %>%
  filter(!is.na(site), !is.na(timepoint), !is.na(colony))

run_coral_mixed_model <- function(df, response_var) {
  df_sub <- df %>% filter(!is.na(.data[[response_var]]))
  if (nrow(df_sub) == 0) return(list(anova = NULL, emmeans = NULL, contrasts = NULL, note = "No rows"))
  if (length(unique(df_sub$site)) < 2) return(list(anova = NULL, emmeans = NULL, contrasts = NULL, note = "One site only"))
  if (length(unique(df_sub$colony)) < 2) return(list(anova = NULL, emmeans = NULL, contrasts = NULL, note = "Not enough colony levels"))

  model <- lmer(
    as.formula(paste0(response_var, " ~ site * subtype + timepoint + (1 | colony)")),
    data = df_sub,
    REML = FALSE
  )

  anova_tbl <- as.data.frame(anova(model))
  anova_tbl$term <- rownames(anova_tbl)
  rownames(anova_tbl) <- NULL

  emm <- emmeans(model, pairwise ~ site | subtype)

  list(
    model = model,
    anova = anova_tbl,
    emmeans = as.data.frame(emm$emmeans),
    contrasts = as.data.frame(emm$contrasts)
  )
}

metrics <- c("Shannon", "InverseSimpson", "Richness")
coral_results <- list()

for (response_var in metrics) {
  coral_results[[response_var]] <- run_coral_mixed_model(coral_data, response_var)
  if (!is.null(coral_results[[response_var]]$anova)) {
    write.table(coral_results[[response_var]]$anova,
                file.path(output_dir, paste0("coral_", response_var, "_anova.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(coral_results[[response_var]]$contrasts)) {
    write.table(coral_results[[response_var]]$contrasts,
                file.path(output_dir, paste0("coral_", response_var, "_site_contrasts.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

run_noncoral_summary <- function(df, response_var, sample_type_name) {
  df_sub <- df %>%
    filter(type == "eDNA", subtype == sample_type_name) %>%
    filter(!is.na(.data[[response_var]]), !is.na(site), !is.na(timepoint))

  if (nrow(df_sub) == 0 || length(unique(df_sub$site)) < 2) return(NULL)

  fit <- lm(as.formula(paste0(response_var, " ~ site + timepoint")), data = df_sub)

  anova_tbl <- as.data.frame(anova(fit))
  anova_tbl$term <- rownames(anova_tbl)
  rownames(anova_tbl) <- NULL

  coeff_tbl <- as.data.frame(summary(fit)$coefficients)
  coeff_tbl$term <- rownames(coeff_tbl)
  rownames(coeff_tbl) <- NULL

  list(model = fit, anova = anova_tbl, coefficients = coeff_tbl)
}

noncoral_types <- c("seawater", "sediment")
noncoral_results <- list()

for (sample_type_name in noncoral_types) {
  for (response_var in metrics) {
    result_name <- paste(sample_type_name, response_var, sep = "_")
    noncoral_results[[result_name]] <- run_noncoral_summary(ps.meta, response_var, sample_type_name)

    if (!is.null(noncoral_results[[result_name]])) {
      write.table(noncoral_results[[result_name]]$anova,
                  file.path(output_dir, paste0(result_name, "_anova.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(noncoral_results[[result_name]]$coefficients,
                  file.path(output_dir, paste0(result_name, "_coefficients.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}

extract_pvalue <- function(tbl, term_name) {
  if (is.null(tbl)) return(NA_real_)
  pcol <- grep("Pr\\(", colnames(tbl), value = TRUE)
  if (length(pcol) == 0) return(NA_real_)
  idx <- which(tbl$term == term_name)
  if (length(idx) == 0) return(NA_real_)
  tbl[[pcol[1]]][idx[1]]
}

synthesis <- data.frame(
  sample_set = character(),
  metric = character(),
  effect = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (metric_name in names(coral_results)) {
  tbl <- coral_results[[metric_name]]$anova
  synthesis <- rbind(
    synthesis,
    data.frame(sample_set = "coral", metric = metric_name, effect = "site", p_value = extract_pvalue(tbl, "site")),
    data.frame(sample_set = "coral", metric = metric_name, effect = "subtype", p_value = extract_pvalue(tbl, "subtype")),
    data.frame(sample_set = "coral", metric = metric_name, effect = "timepoint", p_value = extract_pvalue(tbl, "timepoint")),
    data.frame(sample_set = "coral", metric = metric_name, effect = "site:subtype", p_value = extract_pvalue(tbl, "site:subtype"))
  )
}

for (nm in names(noncoral_results)) {
  obj <- noncoral_results[[nm]]
  if (!is.null(obj)) {
    parts <- strsplit(nm, "_")[[1]]
    synthesis <- rbind(
      synthesis,
      data.frame(sample_set = parts[1], metric = parts[2], effect = "site", p_value = extract_pvalue(obj$anova, "site")),
      data.frame(sample_set = parts[1], metric = parts[2], effect = "timepoint", p_value = extract_pvalue(obj$anova, "timepoint"))
    )
  }
}

write.table(synthesis,
            file.path(output_dir, "alpha_diversity_synthesis_table.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(
  list(
    phyloseq_object = ps,
    metadata_with_alpha = ps.meta,
    coral_results = coral_results,
    noncoral_results = noncoral_results,
    synthesis = synthesis
  ),
  file = file.path(output_dir, "alpha_diversity_model_results.rds")
)

out_path <- "/mnt/data/02_alpha_diversity_mixed_models_metadata_v2.R"
writeLines(readLines(sys.frame(1)$ofile %||% ""), out_path)
