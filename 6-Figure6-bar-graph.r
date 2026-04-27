library(ggplot2)
library(phyloseq)
library(dplyr)
library(tibble)

# ----------------------------
# 1. Load data
# ----------------------------
otu <- read.table("psmin5_ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxon <- read.table("psmin5_taxa_table.txt", sep = "\t", header = TRUE, row.names = 1)
samples <- read.table("metadata.txt", sep = "\t", header = TRUE, row.names = 1)

taxon <- as.matrix(taxon)
TAX <- tax_table(taxon)

ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(TAX)
)

# ----------------------------
# 2. Remove unwanted samples
# ----------------------------
to_remove <- c("LP196", "eDNA1")
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

# ----------------------------
# 3. Agglomerate at Order level
# ----------------------------
ps_order <- tax_glom(ps, taxrank = "Order", NArm = FALSE)

tax_table(ps_order)[, "Order"] <- ifelse(
  is.na(tax_table(ps_order)[, "Order"]) | tax_table(ps_order)[, "Order"] == "",
  "Unclassified",
  as.character(tax_table(ps_order)[, "Order"])
)

# ----------------------------
# 4. Relative abundance
# ----------------------------
ps_order_ra <- transform_sample_counts(ps_order, function(x) x / sum(x))

# ----------------------------
# 5. Melt for ggplot
# ----------------------------
df_order <- psmelt(ps_order_ra)

# ----------------------------
# 6. Select top taxa
# ----------------------------
order_rank <- df_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance))

top_n <- 16

top_orders <- order_rank %>%
  slice_head(n = top_n) %>%
  pull(Order)

df_order_top <- df_order %>%
  mutate(Order_plot = ifelse(Order %in% top_orders, as.character(Order), "Other"))

# ----------------------------
# 7. Create panel variable
# ----------------------------
df_order_top <- df_order_top %>%
  mutate(panel = paste(subtype, site, sep = "\n"))

# ----------------------------
# 8. Set chronological month order and panel order
# ----------------------------
time_levels <- c("Oct.23", "Nov.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24")

df_order_top <- df_order_top %>%
  mutate(
    panel = factor(panel, levels = c(
      "seawater\nLuminao",
      "seawater\nTumon",
      "sediment\nLuminao",
      "sediment\nTumon",
      "healthy\nLuminao",
      "healthy\nTumon",
      "diseased\nLuminao",
      "diseased\nTumon"
    )),
    time_label = factor(`month.year`, levels = time_levels)
  )
# ----------------------------
# 9. Convert date to real Date
# ----------------------------
df_order_top <- df_order_top %>%
  mutate(date2 = as.Date(date, format = "%m/%d/%Y"))

# ----------------------------
# 10. Force "Other" last
# ----------------------------
order_levels <- df_order_top %>%
  filter(Order_plot != "Other") %>%
  group_by(Order_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Order_plot)

order_levels <- c(order_levels, "Other")

df_order_top <- df_order_top %>%
  mutate(Order_plot = factor(Order_plot, levels = order_levels))

# ----------------------------
# 11. Sample ordering + x labels
# ----------------------------
sample_order_df <- df_order_top %>%
  distinct(panel, time_label, date2, Sample) %>%
  arrange(panel, date2, Sample) %>%
  group_by(panel, time_label) %>%
  mutate(
    x_label = ifelse(row_number() == 1, as.character(time_label), ""),
    Sample_plot = paste(panel, Sample, sep = "___")
  ) %>%
  ungroup()

sample_levels <- sample_order_df$Sample_plot

df_order_top <- df_order_top %>%
  mutate(
    Sample_plot = paste(panel, Sample, sep = "___"),
    Sample_plot = factor(Sample_plot, levels = sample_levels)
  )

x_labels <- sample_order_df$x_label
names(x_labels) <- sample_order_df$Sample_plot

# ----------------------------
# 12. Colors
# ----------------------------
color_values <- c(
  "cornflowerblue",
  "darkseagreen1",
  "darkgoldenrod4",
  "darkolivegreen",
  "darkmagenta",
  "antiquewhite2",
  "lightpink1",
  "darkred",
  "darkslateblue",
  "coral1",
  "aquamarine4",
  "burlywood2",
  "darkgoldenrod1",
  "mediumpurple",
  "darkgreen",
  "darkslategray",
  "lightgray"
)

names(color_values) <- levels(df_order_top$Order_plot)

# ----------------------------
# 13. Plot
# ----------------------------
p <- ggplot(df_order_top, aes(x = Sample_plot, y = Abundance, fill = Order_plot)) +
  geom_col(width = 0.95) +
  facet_wrap(~ panel, scales = "free_x", ncol = 2) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = color_values, drop = FALSE, na.value = "black") +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  theme_bw() +
  theme(
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom"
  )

p
# ----------------------------
# 14. Save
# ----------------------------
ggsave(
  "Figure6-top16_orders_ggplot.tiff",
  plot = p,
  units = "in",
  width = 8,
  height = 10,
  dpi = 300,
  compression = "lzw"
)

ggsave("Figure6-top16_orders_ggplot.pdf",
       plot = p,
       width = 8,
       height = 10,
       units = "in")
	   
	   