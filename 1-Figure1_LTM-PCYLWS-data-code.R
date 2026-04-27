# Figure 1. Long-term PCYLWS prevalence and colony abundance by size class
# Input files required in working directory:
#   Lum-Pcyl-WS-sev-sz-cl.txt
#   Tum-Pcyl-WS-sev-sz-cl.txt
# Output:
#   Figure1_PCylWS_prevalence_colony_count.png
#   Figure1_PCylWS_prevalence_colony_count.tiff

library(tidyverse)
library(lubridate)
library(stringr)
library(readr)
library(patchwork)

# ------------------------------------------------------------
# 1. Read and clean data
# ------------------------------------------------------------

lum <- read_tsv("Lum-Pcyl-WS-sev-sz-cl.txt", show_col_types = FALSE) %>%
  mutate(
    DATE = mdy(DATE),
    TRANSECT = str_remove_all(as.character(TRANSECT), "[^0-9]"),
    TRANSECT = parse_integer(TRANSECT),
    SZ_CL = as.integer(SZ_CL),
    WS = na_if(WS, ""),
    WS_SEV = na_if(`WS-SEV`, "")
  ) %>%
  select(-`WS-SEV`)

tum <- read_tsv("Tum-Pcyl-WS-sev-sz-cl.txt", show_col_types = FALSE) %>%
  mutate(
    DATE = mdy(DATE),
    TRANSECT = str_remove_all(as.character(TRANSECT), "[^0-9]"),
    TRANSECT = parse_integer(TRANSECT),
    SZ_CL = as.integer(SZ_CL),
    WS = na_if(WS, ""),
    WS_SEV = na_if(`WS-SEV`, "")
  ) %>%
  select(-`WS-SEV`)

dat <- bind_rows(lum, tum) %>%
  mutate(
    SITE = factor(SITE, levels = c("LUM", "TUM")),
    year = year(DATE)
  )

# Size-class labels used in Figure 1
size_code <- tibble(
  SZ_CL = 1:6,
  size_range = c("< 10 cm", "11 - 30 cm", "31 - 60 cm",
                 "61 - 100 cm", "1 - 2 m", "> 2 m")
)

size_levels <- size_code$size_range
site_cols <- c(LUM = "gray30", TUM = "steelblue")

# ------------------------------------------------------------
# 2. Panel A: mean annual WS prevalence by site and size class
# ------------------------------------------------------------
# Prevalence is first calculated for each transect, year, and size class,
# then averaged across transects/surveys within each year.

ws_size_transect <- dat %>%
  filter(!is.na(SZ_CL)) %>%
  group_by(SITE, year, TRANSECT, SZ_CL) %>%
  summarise(
    total_colonies = n(),
    ws_colonies = sum(WS == "WS", na.rm = TRUE),
    ws_prev_pct = 100 * ws_colonies / total_colonies,
    .groups = "drop"
  )

ws_plot_year <- ws_size_transect %>%
  group_by(SITE, year, SZ_CL) %>%
  summarise(
    mean_ws_prev = mean(ws_prev_pct, na.rm = TRUE),
    sd_ws_prev = sd(ws_prev_pct, na.rm = TRUE),
    n_transects = n(),
    se_ws_prev = sd_ws_prev / sqrt(n_transects),
    .groups = "drop"
  ) %>%
  left_join(size_code, by = "SZ_CL") %>%
  mutate(size_range = factor(size_range, levels = size_levels, ordered = TRUE))

p_prev <- ggplot(ws_plot_year, aes(x = year, y = mean_ws_prev, fill = SITE)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(mean_ws_prev - se_ws_prev, 0),
        ymax = mean_ws_prev + se_ws_prev),
    width = 0.25,
    alpha = 0.6
  ) +
  facet_grid(SITE ~ size_range) +
  scale_fill_manual(values = site_cols) +
  scale_x_continuous(breaks = seq(2011, 2025, by = 2)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Survey year",
    y = "Mean annual WS prevalence (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 3. Panel B: mean annual colony count per transect by site and size class
# ------------------------------------------------------------

col_transect <- dat %>%
  filter(!is.na(SZ_CL)) %>%
  group_by(SITE, year, TRANSECT, SZ_CL) %>%
  summarise(
    n_colonies = n(),
    .groups = "drop"
  )

col_plot_year <- col_transect %>%
  group_by(SITE, year, SZ_CL) %>%
  summarise(
    mean_colonies = mean(n_colonies, na.rm = TRUE),
    sd_colonies = sd(n_colonies, na.rm = TRUE),
    n_transects = n(),
    se_colonies = sd_colonies / sqrt(n_transects),
    .groups = "drop"
  ) %>%
  left_join(size_code, by = "SZ_CL") %>%
  mutate(size_range = factor(size_range, levels = size_levels, ordered = TRUE))

p_col <- ggplot(col_plot_year, aes(x = year, y = mean_colonies, fill = SITE)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(mean_colonies - se_colonies, 0),
        ymax = mean_colonies + se_colonies),
    width = 0.25,
    alpha = 0.6
  ) +
  facet_grid(SITE ~ size_range) +
  scale_fill_manual(values = site_cols) +
  scale_x_continuous(breaks = seq(2011, 2025, by = 2)) +
  labs(
    x = "Survey year",
    y = "Mean annual colony count per transect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 4. Combine and save final Figure 1
# ------------------------------------------------------------

figure1 <- (p_prev / p_col) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

figure1

ggsave(
  filename = "Figure1_PCylWS_prevalence_colony_count.png",
  plot = figure1,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Figure1_PCylWS_prevalence_colony_count.tiff",
  plot = figure1,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300,
  compression = "lzw"
)
