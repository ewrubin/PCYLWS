library(phyloseq)
library(dplyr)
library(ggplot2)

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
# 3. Agglomerate at Order level and convert to relative abundance
# ----------------------------
ps_order <- tax_glom(ps, taxrank = "Order", NArm = FALSE)

tax_table(ps_order)[, "Order"] <- ifelse(
  is.na(tax_table(ps_order)[, "Order"]) | tax_table(ps_order)[, "Order"] == "",
  "Unclassified",
  as.character(tax_table(ps_order)[, "Order"])
)

ps_order_ra <- transform_sample_counts(ps_order, function(x) x / sum(x))

# ----------------------------
# 4. Melt phyloseq object
# ----------------------------
df_order <- psmelt(ps_order_ra)

# ----------------------------
# 3. Keep only the four focal orders
# ----------------------------
focal_orders <- c("Rhodobacterales",
                  "Enterobacterales",
                  "Flavobacteriales",
                  "Cytophagales")

df_focal <- df_order %>%
  filter(Order %in% focal_orders)

# ----------------------------
# 5. Set factor levels
# ----------------------------
time_levels <- c("Oct.23", "Nov.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24")

df_focal <- df_focal %>%
  mutate(
    type = factor(type, levels = c("seawater", "sediment", "healthy", "diseased")),
    Order = factor(Order, levels = focal_orders),
    site = factor(site, levels = c("Luminao", "Tumon")),
    month.year = factor(`month.year`, levels = time_levels)
  )

# ----------------------------
# 6. Calculate mean relative abundance for each panel group
# ----------------------------
df_mean <- df_focal %>%
  group_by(type, Order, site, month.year) %>%
  summarise(
    mean_RA = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 7. Plot
# ----------------------------
p <- ggplot(df_mean,
            aes(x = month.year, y = mean_RA,
                color = site, linetype = site, group = site)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(type ~ Order, scales = "free_y") +
  scale_color_manual(values = c("Luminao" = "gray60", "Tumon" = "black")) +
  scale_linetype_manual(values = c("Luminao" = "solid", "Tumon" = "dashed")) +
  labs(x = NULL, y = "Mean relative abundance", color = "Site", linetype = "Site") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p

# ----------------------------
# 8. Save
# ----------------------------
ggsave("Fig_copiotrophic_orders_4x4.tiff",
       plot = p,
       width = 8,
       height = 8,
       dpi = 300,
       compression = "lzw")

ggsave("Fig_copiotrophic_orders_4x4.pdf",
       plot = p,
       width = 8,
       height = 8)
	   
	   
	   
library(phyloseq)
library(dplyr)
library(ggplot2)

# ----------------------------
# 1. Agglomerate at Order level and convert to relative abundance
# ----------------------------
ps_order <- tax_glom(ps, taxrank = "Order", NArm = FALSE)

tax_table(ps_order)[, "Order"] <- ifelse(
  is.na(tax_table(ps_order)[, "Order"]) | tax_table(ps_order)[, "Order"] == "",
  "Unclassified",
  as.character(tax_table(ps_order)[, "Order"])
)

ps_order_ra <- transform_sample_counts(ps_order, function(x) x / sum(x))

# ----------------------------
# 2. Melt phyloseq object
# ----------------------------
df_order <- psmelt(ps_order_ra)

# ----------------------------
# 3. Keep only focal orders and tissue samples
# ----------------------------
focal_orders <- c("Rhodobacterales",
                  "Enterobacterales",
                  "Flavobacteriales",
                  "Cytophagales")

df_tissue <- df_order %>%
  filter(Order %in% focal_orders,
         type %in% c("healthy", "diseased"))

# ----------------------------
# 4. Factor levels
# ----------------------------
time_levels <- c("Oct.23", "Nov.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24")

df_tissue <- df_tissue %>%
  mutate(
    Order = factor(Order, levels = focal_orders),
    type = factor(type, levels = c("healthy", "diseased")),
    `month.year` = factor(`month.year`, levels = time_levels),
    cat = factor(cat, levels = c("L-H", "T-H", "L-D", "T-D"))
  )

# ----------------------------
# 5. Mean relative abundance
# ----------------------------
df_mean_tissue <- df_tissue %>%
  group_by(Order, `month.year`, cat) %>%
  summarise(mean_RA = mean(Abundance, na.rm = TRUE), .groups = "drop")

# ----------------------------
# 6. Colors
# ----------------------------
mycolors <- c(
  "L-H" = "darkolivegreen1",
  "T-H" = "darkolivegreen3",
  "L-D" = "brown1",
  "T-D" = "brown3",
  "L-W" = "cadetblue1",
  "T-W" = "cadetblue3",
  "L-S" = "magenta2",
  "T-S" = "magenta4"
)

# use only the 4 tissue colors
tissue_colors <- mycolors[c("L-H", "T-H", "L-D", "T-D")]

# optional linetypes to help distinguish healthy vs diseased
line_types <- c(
  "L-H" = "solid",
  "T-H" = "dashed",
  "L-D" = "solid",
  "T-D" = "dashed"
)

# ----------------------------
# 7. Plot
# ----------------------------
p_tissue <- ggplot(df_mean_tissue,
                   aes(x = `month.year`, y = mean_RA,
                       color = cat, linetype = cat, group = cat)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ Order, scales = "free_y", ncol = 2) +
  scale_color_manual(values = tissue_colors) +
  scale_linetype_manual(values = line_types) +
  labs(x = NULL, y = "Mean relative abundance", color = "Category", linetype = "Category") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_tissue

# ----------------------------
# 8. Save
# ----------------------------
ggsave("Fig_tissue_copiotrophic_orders.tiff",
       plot = p_tissue,
       width = 7,
       height = 5,
       dpi = 300,
       compression = "lzw")

ggsave("Fig_tissue_copiotrophic_orders.pdf",
       plot = p_tissue,
       width = 7,
       height = 5)	

#------------------------------------------------------------------------------------	  
#------------------------------------------------------------------------------------

library(phyloseq)
library(dplyr)
library(ggplot2)

# ----------------------------
# 1. Agglomerate at Order level and convert to relative abundance
# ----------------------------
ps_order <- tax_glom(ps, taxrank = "Order", NArm = FALSE)

tax_table(ps_order)[, "Order"] <- ifelse(
  is.na(tax_table(ps_order)[, "Order"]) | tax_table(ps_order)[, "Order"] == "",
  "Unclassified",
  as.character(tax_table(ps_order)[, "Order"])
)

ps_order_ra <- transform_sample_counts(ps_order, function(x) x / sum(x))

# ----------------------------
# 2. Melt phyloseq object
# ----------------------------
df_order <- psmelt(ps_order_ra)

# ----------------------------
# 3. Keep only focal orders
# ----------------------------
focal_orders <- c("Rhodobacterales",
                  "Enterobacterales",
                  "Flavobacteriales",
                  "Cytophagales")

df_plot <- df_order %>%
  filter(Order %in% focal_orders)

# ----------------------------
# 4. Create habitat group from cat
# ----------------------------
# Assumes cat contains:
# L-H, T-H, L-D, T-D, L-W, T-W, L-S, T-S

df_plot <- df_plot %>%
  mutate(
    habitat_group = case_when(
      cat %in% c("L-W", "T-W", "L-S", "T-S") ~ "Environment",
      cat %in% c("L-H", "T-H", "L-D", "T-D") ~ "Coral tissue",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(habitat_group))

# ----------------------------
# 5. Set factor levels
# ----------------------------
time_levels <- c("Oct.23", "Nov.23", "Jan.24", "Feb.24", "Mar.24", "May.24", "Aug.24")

df_plot <- df_plot %>%
  mutate(
    Order = factor(Order, levels = focal_orders),
    habitat_group = factor(habitat_group, levels = c("Environment", "Coral tissue")),
    `month.year` = factor(`month.year`, levels = time_levels),
    cat = factor(cat, levels = c("L-W", "T-W", "L-S", "T-S", "L-H", "T-H", "L-D", "T-D"))
  )

# ----------------------------
# 6. Calculate mean relative abundance
# ----------------------------
df_mean <- df_plot %>%
  group_by(Order, habitat_group, `month.year`, cat) %>%
  summarise(
    mean_RA = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 7. Colors
# ----------------------------
mycolors <- c(
  "L-H" = "darkolivegreen1",
  "T-H" = "darkolivegreen3",
  "L-D" = "brown1",
  "T-D" = "brown3",
  "L-W" = "cadetblue1",
  "T-W" = "cadetblue3",
  "L-S" = "magenta2",
  "T-S" = "magenta4"
)

# Optional: linetypes to distinguish site
mylinetypes <- c(
  "L-H" = "solid",
  "T-H" = "dashed",
  "L-D" = "solid",
  "T-D" = "dashed",
  "L-W" = "solid",
  "T-W" = "dashed",
  "L-S" = "solid",
  "T-S" = "dashed"
)

# ----------------------------
# 8. Plot
# ----------------------------
p_4x2 <- ggplot(
  df_mean,
  aes(x = `month.year`, y = mean_RA, color = cat, linetype = cat, group = cat)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(Order ~ habitat_group, scales = "free_y") +
  scale_color_manual(values = mycolors) +
  scale_linetype_manual(values = mylinetypes) +
  labs(
    x = NULL,
    y = "Mean relative abundance",
    color = "Category",
    linetype = "Category"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_4x2

p_4x2 <- p_4x2 +
  scale_color_manual(
    values = mycolors,
    name = "Site-type",
    breaks = c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"),
    labels = c(
      "L-H: Luminao-healthy",
      "T-H: Tumon-healthy",
      "L-D: Luminao-diseased",
      "T-D: Tumon-diseased",
      "L-W: Luminao-seawater",
      "T-W: Tumon-seawater",
      "L-S: Luminao-sediment",
      "T-S: Tumon-sediment"
    )
  ) +
  scale_linetype_manual(
    values = mylinetypes,
    name = "Site-type",
    breaks = c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"),
    labels = c(
      "L-H: Luminao-healthy",
      "T-H: Tumon-healthy",
      "L-D: Luminao-diseased",
      "T-D: Tumon-diseased",
      "L-W: Luminao-seawater",
      "T-W: Tumon-seawater",
      "L-S: Luminao-sediment",
      "T-S: Tumon-sediment"
    )
  ) +
  guides(
    linetype = "none",
    color = guide_legend(nrow = 2, byrow = TRUE)
  )
  
p_4x2

# ----------------------------
# 9. Save
# ----------------------------
ggsave(
  "Fig_copiotrophic_orders_4x2.tiff",
  plot = p_4x2,
  width = 8,
  height = 9,
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "Fig_copiotrophic_orders_4x2.pdf",
  plot = p_4x2,
  width = 8,
  height = 9
)