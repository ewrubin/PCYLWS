# =========================================================
# PCYLWS lesion + environmental integration
# =========================================================

library(tidyverse)

# ---------------------------
# 1. Import data
# ---------------------------
lesion <- read.delim("lesion_data.txt", sep = "\t", header = TRUE, check.names = FALSE)
chla   <- read.delim("chla-data.txt",   sep = "\t", header = TRUE, check.names = FALSE)
nit    <- read.delim("nitrate-data.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Quick check
str(lesion)
str(chla)
str(nit)

# ---------------------------
# 2. Standardize month labels
# ---------------------------
# lesion has column "month"
# chla and nitrate have column "month.year"

lesion <- lesion %>%
  rename(month.year = month)

# Make sure site names are consistent
lesion$site <- as.character(lesion$site)
chla$site   <- as.character(chla$site)
nit$site    <- as.character(nit$site)

# Standardize month-year labels if needed
# Here I keep your labels as character strings for simple matching
unique(lesion$month.year)
unique(chla$month.year)
unique(nit$month.year)

# ---------------------------
# 3. Summarize lesion data by site × month
# ---------------------------
lesion_sum <- lesion %>%
  group_by(site, month.year) %>%
  summarise(
    n_lesions = n(),
    mean_progr = mean(progr.rate, na.rm = TRUE),
    sd_progr   = sd(progr.rate, na.rm = TRUE),
    se_progr   = sd_progr / sqrt(n_lesions),
    mean_size  = mean(size, na.rm = TRUE),
    sd_size    = sd(size, na.rm = TRUE),
    .groups = "drop"
  )

lesion_sum

# ---------------------------
# 4. Summarize environmental data by site × month
# ---------------------------
chla_sum <- chla %>%
  group_by(site, month.year) %>%
  summarise(
    n_chla    = n(),
    mean_chla = mean(chla, na.rm = TRUE),
    sd_chla   = sd(chla, na.rm = TRUE),
    se_chla   = sd_chla / sqrt(n_chla),
    .groups = "drop"
  )

nit_sum <- nit %>%
  group_by(site, month.year) %>%
  summarise(
    n_nit    = n(),
    mean_nit = mean(nit, na.rm = TRUE),
    sd_nit   = sd(nit, na.rm = TRUE),
    se_nit   = sd_nit / sqrt(n_nit),
    .groups = "drop"
  )

chla_sum
nit_sum

# ---------------------------
# 5. STRICT overlap dataset
# ---------------------------
# Use only exact overlapping months among lesion, chl a, and nitrate

strict_dat <- lesion_sum %>%
  inner_join(chla_sum, by = c("site", "month.year")) %>%
  inner_join(nit_sum,  by = c("site", "month.year")) %>%
  arrange(site, month.year)

strict_dat

# Save if useful
write.csv(strict_dat, "merged_strict_overlap.csv", row.names = FALSE)

# ---------------------------
# 6. FEBRUARY approximation dataset
# ---------------------------
# Approximate Feb.24 environmental values as mean of Jan.24 and Mar.24 for each site

chla_feb <- chla_sum %>%
  filter(month.year %in% c("Jan.24", "Mar.24")) %>%
  group_by(site) %>%
  summarise(
    month.year = "Feb.24",
    n_chla = sum(n_chla),
    mean_chla = mean(mean_chla, na.rm = TRUE),
    sd_chla = NA_real_,
    se_chla = NA_real_,
    .groups = "drop"
  )

nit_feb <- nit_sum %>%
  filter(month.year %in% c("Jan.24", "Mar.24")) %>%
  group_by(site) %>%
  summarise(
    month.year = "Feb.24",
    n_nit = sum(n_nit),
    mean_nit = mean(mean_nit, na.rm = TRUE),
    sd_nit = NA_real_,
    se_nit = NA_real_,
    .groups = "drop"
  )

# Add approximated Feb to environmental summaries
chla_sum_plus <- bind_rows(chla_sum, chla_feb)
nit_sum_plus  <- bind_rows(nit_sum, nit_feb)

interp_dat <- lesion_sum %>%
  filter(month.year %in% c("Nov.23", "Feb.24", "May.24", "Aug.24")) %>%
  left_join(chla_sum_plus, by = c("site", "month.year")) %>%
  left_join(nit_sum_plus,  by = c("site", "month.year")) %>%
  arrange(site, month.year)

interp_dat

write.csv(interp_dat, "merged_with_feb_approx.csv", row.names = FALSE)

# ---------------------------
# 7. Linear models - STRICT dataset
# ---------------------------
# Very small sample size, so treat results as exploratory

cat("\n====================\n")
cat("STRICT DATASET MODELS\n")
cat("====================\n")

m1_strict <- lm(mean_progr ~ mean_chla, data = strict_dat)
m2_strict <- lm(mean_progr ~ mean_nit,  data = strict_dat)
m3_strict <- lm(mean_progr ~ mean_chla + mean_nit, data = strict_dat)
m4_strict <- lm(mean_progr ~ mean_chla + mean_nit + site, data = strict_dat)

summary(m1_strict)
summary(m2_strict)
summary(m3_strict)
summary(m4_strict)

# ---------------------------
# 8. Linear models - dataset with February approximation
# ---------------------------
cat("\n==============================\n")
cat("INTERPOLATED FEB DATASET MODELS\n")
cat("==============================\n")

m1_interp <- lm(mean_progr ~ mean_chla, data = interp_dat)
m2_interp <- lm(mean_progr ~ mean_nit,  data = interp_dat)
m3_interp <- lm(mean_progr ~ mean_chla + mean_nit, data = interp_dat)
m4_interp <- lm(mean_progr ~ mean_chla + mean_nit + site, data = interp_dat)

summary(m1_interp)
summary(m2_interp)
summary(m3_interp)
summary(m4_interp)

# ---------------------------
# 9. Optional: include lesion size as another response or covariate
# ---------------------------
cat("\n====================================\n")
cat("OPTIONAL MODELS INCLUDING LESION SIZE\n")
cat("====================================\n")

m5_interp <- lm(mean_progr ~ mean_chla + mean_nit + mean_size + site, data = interp_dat)
summary(m5_interp)

m6_interp <- lm(mean_size ~ mean_chla + mean_nit + mean_progr + site, data = interp_dat)
summary(m6_interp)


# ---------------------------
# 10. Quick plots
# ---------------------------

p1 <- ggplot(interp_dat, aes(x = mean_chla, y = mean_progr, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(nudge_y = 0.01, show.legend = FALSE) +
  labs(
    x = "Mean chlorophyll a",
    y = "Mean lesion progression rate",
    title = "Lesion progression vs chlorophyll a"
  ) +
  theme_bw()

p2 <- ggplot(interp_dat, aes(x = mean_nit, y = mean_progr, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(nudge_y = 0.01, show.legend = FALSE) +
  labs(
    x = "Mean nitrate",
    y = "Mean lesion progression rate",
    title = "Lesion progression vs nitrate"
  ) +
  theme_bw()

print(p1)
print(p2)

ggsave("lesion_vs_chla.png", p1, width = 6, height = 5, dpi = 300)
ggsave("lesion_vs_nitrate.png", p2, width = 6, height = 5, dpi = 300)

# ---------------------------
# 11. Correlations (optional)
# ---------------------------
cat("\n====================\n")
cat("CORRELATIONS\n")
cat("====================\n")

cor.test(interp_dat$mean_progr, interp_dat$mean_chla, method = "pearson")
cor.test(interp_dat$mean_progr, interp_dat$mean_nit, method = "pearson")

# If worried about tiny sample size / non-normality:
cor.test(interp_dat$mean_progr, interp_dat$mean_chla, method = "spearman")
cor.test(interp_dat$mean_progr, interp_dat$mean_nit, method = "spearman")


-------------------------------------------------------------------------------
library(tidyverse)

# ---------------------------
# Read merged dataset
# ---------------------------
interp_dat <- read.csv("merged_with_feb_approx.csv")

# Make month order explicit
interp_dat$month.year <- factor(
  interp_dat$month.year,
  levels = c("Nov.23", "Feb.24", "May.24", "Aug.24")
)

# ---------------------------
# Clean plotting theme
# ---------------------------
clean_theme <- theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# ---------------------------
# 1. Lesion progression vs chlorophyll a
# ---------------------------
p_prog_chla <- ggplot(interp_dat,
                      aes(x = mean_chla, y = mean_progr, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text(size = 4, nudge_y = 0.003, show.legend = FALSE) +
  labs(
    title = "Lesion progression vs chlorophyll a",
    x = expression("Mean chlorophyll " * a),
    y = expression("Mean lesion progression rate (cm"^2 * " day"^-1 * ")")
  ) +
  clean_theme

# ---------------------------
# 2. Lesion progression vs nitrate
# ---------------------------
p_prog_nit <- ggplot(interp_dat,
                     aes(x = mean_nit, y = mean_progr, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text(size = 4, nudge_y = 0.003, show.legend = FALSE) +
  labs(
    title = "Lesion progression vs nitrate",
    x = "Mean nitrate (mg/L)",
    y = expression("Mean lesion progression rate (cm"^2 * " day"^-1 * ")")
  ) +
  clean_theme

# ---------------------------
# 3. Lesion size vs chlorophyll a
# ---------------------------
p_size_chla <- ggplot(interp_dat,
                      aes(x = mean_chla, y = mean_size, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text(size = 4, nudge_y = 0.12, show.legend = FALSE) +
  labs(
    title = "Lesion size vs chlorophyll a",
    x = expression("Mean chlorophyll " * a),
    y = expression("Mean lesion size (cm"^2 * ")")
  ) +
  clean_theme

# ---------------------------
# 4. Lesion size vs nitrate
# ---------------------------
p_size_nit <- ggplot(interp_dat,
                     aes(x = mean_nit, y = mean_size, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text(size = 4, nudge_y = 0.12, show.legend = FALSE) +
  labs(
    title = "Lesion size vs nitrate",
    x = "Mean nitrate (mg/L)",
    y = expression("Mean lesion size (cm"^2 * ")")
  ) +
  clean_theme

# ---------------------------
# Print plots
# ---------------------------
print(p_prog_chla)
print(p_prog_nit)
print(p_size_chla)
print(p_size_nit)

# ---------------------------
# Save plots
# ---------------------------
ggsave("lesion_progression_vs_chla_clean.png", p_prog_chla, width = 7, height = 6, dpi = 300)
ggsave("lesion_progression_vs_nitrate_clean.png", p_prog_nit, width = 7, height = 6, dpi = 300)
ggsave("lesion_size_vs_chla_clean.png", p_size_chla, width = 7, height = 6, dpi = 300)
ggsave("lesion_size_vs_nitrate_clean.png", p_size_nit, width = 7, height = 6, dpi = 300)

# ---------------------------
# Optional models for lesion size
# ---------------------------
m_size_chla <- lm(mean_size ~ mean_chla, data = interp_dat)
m_size_nit  <- lm(mean_size ~ mean_nit, data = interp_dat)
m_size_both <- lm(mean_size ~ mean_chla + mean_nit + site, data = interp_dat)

summary(m_size_chla)
summary(m_size_nit)
summary(m_size_both)

# Optional correlations
cor.test(interp_dat$mean_size, interp_dat$mean_chla, method = "spearman")
cor.test(interp_dat$mean_size, interp_dat$mean_nit, method = "spearman")

# ---------------------------
# Save plots as TIFF (publication quality)
# ---------------------------

ggsave("lesion_progression_vs_chla.tiff", p_prog_chla,
       width = 7, height = 6, dpi = 600, compression = "lzw")

ggsave("lesion_progression_vs_nitrate.tiff", p_prog_nit,
       width = 7, height = 6, dpi = 600, compression = "lzw")

ggsave("lesion_size_vs_chla.tiff", p_size_chla,
       width = 7, height = 6, dpi = 600, compression = "lzw")

ggsave("lesion_size_vs_nitrate.tiff", p_size_nit,
       width = 7, height = 6, dpi = 600, compression = "lzw")
	   
---------------------------------------------------------------------------------------	   
	   
library(tidyverse)
library(patchwork)

# =========================================================
# 1. Import raw data
# =========================================================
lesion <- read.delim("lesion_data.txt", sep = "\t", header = TRUE, check.names = FALSE)
chla   <- read.delim("chla-data.txt",   sep = "\t", header = TRUE, check.names = FALSE)
nit    <- read.delim("nitrate-data.txt", sep = "\t", header = TRUE, check.names = FALSE)

# lesion month column name -> month.year for consistency
lesion <- lesion %>%
  rename(month.year = month)

# =========================================================
# 2. Define month order
# =========================================================
env_levels <- c("Oct.23", "Nov.23", "Dec.23", "Jan.24", "Mar.24", "May.24", "Aug.24", "Sep.24")
lesion_levels <- c("Nov.23", "Feb.24", "May.24", "Aug.24")

# =========================================================
# 3. Summarize environmental data
# =========================================================
chla_sum <- chla %>%
  group_by(site, month.year) %>%
  summarise(
    mean_chla = mean(chla, na.rm = TRUE),
    sd_chla   = sd(chla, na.rm = TRUE),
    n_chla    = n(),
    se_chla   = sd_chla / sqrt(n_chla),
    .groups = "drop"
  ) %>%
  mutate(month.year = factor(month.year, levels = env_levels))

nit_sum <- nit %>%
  group_by(site, month.year) %>%
  summarise(
    mean_nit = mean(nit, na.rm = TRUE),
    sd_nit   = sd(nit, na.rm = TRUE),
    n_nit    = n(),
    se_nit   = sd_nit / sqrt(n_nit),
    .groups = "drop"
  ) %>%
  mutate(month.year = factor(month.year, levels = env_levels))

# =========================================================
# 4. Summarize lesion data
# =========================================================
lesion_sum <- lesion %>%
  group_by(site, month.year) %>%
  summarise(
    mean_size  = mean(size, na.rm = TRUE),
    sd_size    = sd(size, na.rm = TRUE),
    mean_progr = mean(progr.rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month.year = factor(month.year, levels = lesion_levels))

# =========================================================
# 5. Approximate February env values as mean of Jan + Mar
# =========================================================
chla_feb <- chla_sum %>%
  filter(month.year %in% c("Jan.24", "Mar.24")) %>%
  group_by(site) %>%
  summarise(
    month.year = factor("Feb.24", levels = lesion_levels),
    mean_chla = mean(mean_chla, na.rm = TRUE),
    se_chla   = NA_real_,
    .groups = "drop"
  )

nit_feb <- nit_sum %>%
  filter(month.year %in% c("Jan.24", "Mar.24")) %>%
  group_by(site) %>%
  summarise(
    month.year = factor("Feb.24", levels = lesion_levels),
    mean_nit = mean(mean_nit, na.rm = TRUE),
    se_nit   = NA_real_,
    .groups = "drop"
  )

# Keep only months overlapping lesion monitoring + Feb approximation
chla_lesion <- chla_sum %>%
  filter(month.year %in% c("Nov.23", "May.24", "Aug.24")) %>%
  select(site, month.year, mean_chla, se_chla) %>%
  mutate(month.year = factor(as.character(month.year), levels = lesion_levels))

nit_lesion <- nit_sum %>%
  filter(month.year %in% c("Nov.23", "May.24", "Aug.24")) %>%
  select(site, month.year, mean_nit, se_nit) %>%
  mutate(month.year = factor(as.character(month.year), levels = lesion_levels))

chla_lesion <- bind_rows(chla_lesion, chla_feb) %>%
  mutate(month.year = factor(as.character(month.year), levels = lesion_levels))

nit_lesion <- bind_rows(nit_lesion, nit_feb) %>%
  mutate(month.year = factor(as.character(month.year), levels = lesion_levels))

# =========================================================
# 6. Merge for integration panel
# =========================================================
plot_dat <- lesion_sum %>%
  left_join(chla_lesion, by = c("site", "month.year")) %>%
  left_join(nit_lesion,  by = c("site", "month.year"))

# =========================================================
# 7. Shared plot settings
# =========================================================
site_cols <- c("Luminao" = "cyan1", "Tumon" = "cyan4")

theme_fig <- theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "plain"),
	axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"))

# =========================================================
# 8. Panel A: chlorophyll a through time
# =========================================================
p9A <- ggplot(chla_sum, aes(x = month.year, y = mean_chla, color = site, group = site)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_chla - se_chla, ymax = mean_chla + se_chla),
                width = 0.15, linewidth = 0.7) +
  scale_color_manual(values = site_cols) +
  labs(
    x = "Month-Year",
    y = expression("Chlorophyll " * a * " (" * mu * "g L"^-1 * ")"),
    color = "site"
  ) +
  theme_fig +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top",legend.direction = "horizontal"
  )
 
# =========================================================
# 9. Panel B: nitrate through time
# =========================================================
p9B <- ggplot(nit_sum, aes(x = month.year, y = mean_nit, color = site, group = site)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_nit - se_nit, ymax = mean_nit + se_nit),
                width = 0.15, linewidth = 0.7) +
  scale_color_manual(values = site_cols) +
  labs(
    x = expression("Month-Year"),
    y = expression("Nitrate (mg L"^-1 * ")"),
    color = "site"
  ) +
  theme_fig +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  theme(
    legend.position = "none"
  )


# =========================================================
# 10. Panel C: lesion size vs chlorophyll a
# =========================================================

  
p9C <- ggplot(plot_dat,
              aes(x = mean_chla, y = mean_size, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = site_cols) +
  labs(
    x = expression("Mean chlorophyll " * a * " (" * mu * "g L"^-1 * ")"),
    y = expression("Mean lesion size (cm"^2 * ")"),
    color = "site"
  ) +
  theme_fig +
  theme(
    legend.position = "none"
  )



# =========================================================
# 11. Panel D: lesion progression vs chlorophyll a
# =========================================================

  
p9D <- ggplot(plot_dat,
              aes(x = mean_chla, y = mean_progr, color = site, label = month.year)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = site_cols) +
  labs(
    x = expression("Mean chlorophyll " * a * " (" * mu * "g L"^-1 * ")"),
    y = expression("Mean lesion progression rate (cm"^2 * " day"^-1 * ")"),
    color = "site"
  ) +
  theme_fig + 
  theme(
    legend.position = "none"
  )




# =========================================================
# 12. Combine into Figure 9
# =========================================================

fig9 <- (p9A | p9B) / (p9C | p9D) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )
	
	
# Display
print(p9A)
print(p9B)
print(p9C)
print(p9D)
print(fig9)


# =========================================================
# 12. Save panels as TIFF
# =========================================================
ggsave("Figure9A_chlorophyll_time.tiff", p9A,
       width = 7, height = 5.5, dpi = 600, compression = "lzw")

ggsave("Figure9B_nitrate_time.tiff", p9B,
       width = 7, height = 5.5, dpi = 600, compression = "lzw")

ggsave("Figure9C_lesion_size_vs_chla.tiff", p9C,
       width = 7, height = 5.5, dpi = 600, compression = "lzw")


ggsave("Figure9D_lesion_progr_vs_chla.tiff", p9D,
       width = 7, height = 5.5, dpi = 600, compression = "lzw")



# =========================================================
# 13. Save combined figure as TIFF
# =========================================================
ggsave("Figure9_revised_v2.tiff", fig9,
       width = 12, height = 10, dpi = 600, compression = "lzw")	   
	   
	   
	   
	   