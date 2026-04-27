library(tidyverse)
library(lubridate)

# read data
temp <- read_tsv(
  "Lum-Tum-Aug.2023-Aug.2024-final.txt",
  col_types = cols(
    Date = col_character(),
    Time = col_character(),
    Lum  = col_double(),
    Tum  = col_double()
  )
)

# fix 24:00:00 → 00:00:00 next day
temp <- temp %>%
  mutate(
    Time = if_else(Time == "24:00:00", "00:00:00", Time),
    Date = as.Date(Date, format = "%m/%d/%Y"),
    datetime = as.POSIXct(
      paste(Date, Time),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Pacific/Guam"
    )
  ) %>%
  arrange(datetime)

# reshape to long format
temp_long <- temp %>%
  pivot_longer(
    cols = c(Lum, Tum),
    names_to = "site",
    values_to = "temperature"
  ) %>%
  filter(!is.na(temperature))

--------------------------------------------------
#Hourly temperature time series, first look 

ggplot(temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  scale_color_manual(
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +
  labs(
    x = "Date",
    y = "Temperature (°C)",
    color = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )


ggsave("Temp-first-look.tiff",, units = "in", width = 10, height = 6, dpi = 300, compression = "lzw")
----------------------------------------------------------------------------
#Daily mean 

temp_daily <- temp_long %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(site, date) %>%
  summarise(
    temp_mean = mean(temperature, na.rm = TRUE),
    .groups = "drop")




ggplot(temp_daily, aes(x = date, y = temp_mean, color = site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +
  labs(
    x = "Date",
    y = "Daily mean temperature (°C)",
    color = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave("Temp-daily-mean.tiff",, units = "in", width = 10, height = 6, dpi = 300, compression = "lzw")


-------------------------------------------------------------------------------------
#daily mean better

ggplot(temp_daily, aes(x = date, y = temp_mean, color = site)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(\
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +
  labs(
    x = "Date",
    y = "Daily mean temperature (°C)",
    color = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10)
  )


ggsave("Temp-daily-mean-10x5.tiff",, units = "in", width = 10, height = 5, dpi = 300, compression = "lzw")

---------------------------------------------------------------------------------------------------------
#adding Degree Heating Weeks (DHW-style metrics) and Shaded bleaching thresholds (e.g. 30–31 °C)

#Awesome — here’s a clean DHW-style workflow using your daily mean data, plus shaded bleaching thresholds (30–31 °C).

#What this does

#Shades 30–31 °C (and optionally >31 °C) behind the temperature time series

#Computes a DHW-style metric from daily data:

#Estimate MMM (Maximum Monthly Mean) from your dataset (one-year deployment → best available estimate, but not a true climatology)

#Compute HotSpot = max(T − (MMM + 1), 0)

#Compute DHW ≈ rolling 12-week sum of HotSpot / 7 (°C-weeks)

#One degree Celsius (1 °C) above the MMM is called the "bleaching threshold" temperature. 
--------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(slider)
library(patchwork)

# temp_daily must exist from your earlier step:
# columns: site, date, temp_mean

# 1) estimate MMM per site from this dataset (max of monthly means)
mmm_tbl <- temp_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(site, month) %>%
  summarise(monthly_mean = mean(temp_mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  summarise(MMM = max(monthly_mean, na.rm = TRUE), .groups = "drop")

# 2) compute HotSpot and DHW-style metric (12-week rolling window = 84 days)
temp_dhw <- temp_daily %>%
  left_join(mmm_tbl, by = "site") %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    threshold = MMM + 1,                               # NOAA-style bleaching threshold
    hotspot = pmax(temp_mean - threshold, 0),          # only count above threshold
    dhw = slide_dbl(hotspot, sum, .before = 83, .complete = TRUE) / 7
  ) %>%
  ungroup()
  
  
p_temp <- ggplot(temp_dhw, aes(x = date, y = temp_mean, color = site)) +
  # Shaded bleaching-relevant band(s)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = 31, alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 31, ymax = Inf, alpha = 0.08) +  # optional

  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(30, 31), linetype = "dashed", linewidth = 0.4) +

  scale_color_manual(
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(x = NULL, y = "Daily mean temperature (°C)", color = "Site") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )


p_dhw <- ggplot(temp_dhw, aes(x = date, y = dhw, color = site)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_color_manual(
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(x = "Date", y = "DHW (°C-weeks)", color = "Site") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


p_temp / p_dhw + plot_layout(heights = c(2, 1))


ggsave("Temp-daily-mean-DHW-bleach-threshold.tiff",, units = "in", width = 10, height = 6, dpi = 300, compression = "lzw")


#precipitation data 

library(rnoaa)
library(httr)
library(jsonlite)
library(tidyverse)
library(lubridate)

library(httr)
library(dplyr)
library(tibble)

token <- Sys.getenv("NOAA_KEY")

get_cdo_daily_prcp <- function(stationid, start_date, end_date) {
  resp <- GET(
    "https://www.ncei.noaa.gov/cdo-web/api/v2/data",
    query = list(
      datasetid  = "GHCND",
      stationid  = stationid,
      datatypeid = "PRCP",
      startdate  = start_date,
      enddate    = end_date,
      units      = "metric",
      limit      = 1000
    ),
    add_headers(token = token)
  )
  stop_for_status(resp)
  js <- content(resp, as="parsed", type="application/json")

  if (is.null(js$results) || length(js$results) == 0) {
    return(tibble(date=as.Date(character()), prcp_mm=numeric()))
  }

  bind_rows(js$results) %>%
    transmute(
      date = as.Date(substr(date, 1, 10)),
      prcp_mm = value
    ) %>%
    arrange(date)
}

prcp_df <- get_cdo_daily_prcp(
  stationid  = "GHCND:GQW00041406",
  start_date = "2023-08-01",
  end_date   = "2024-08-31"
)

nrow(prcp_df)
summary(prcp_df$prcp_mm)


write.csv(prcp_df, "Guam-ariport-prep-data.csv",
          row.names = FALSE)
		  
		  
library(ggplot2)

ggplot(prcp_df, aes(x = date, y = prcp_mm)) +
  geom_col() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(x = "Date", y = "Daily precipitation (mm)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
		  
ggsave("Prcp-daily-Guam-airport.tiff",, units = "in", width = 10, height = 5, dpi = 300, compression = "lzw")



p_prcp_w <- ggplot(prcp_week, aes(x = date, y = prcp_7d_mm)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(x = "Date", y = "7-day cumulative precipitation (mm)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())


p_prcp_w

ggsave("Prcp-7day-cumulative-mean.tiff",, units = "in", width = 10, height = 5, dpi = 300, compression = "lzw")

prcp_week <- prcp_df %>%
  arrange(date) %>%
  mutate(prcp_7d_mm = slider::slide_dbl(prcp_mm, sum, .before = 6, .complete = FALSE))


write.csv(prcp_week, "Guam-ariport-prep-7day-mean-data.csv",
          row.names = FALSE)
		  

-------------------------------------------------------------------------------------------------------
#adding sampling annotation with R 

library(tidyverse)
library(lubridate)

# Seasons (final wet ends at last data day)
season_df <- tribble(
  ~season, ~start,        ~end,
  "wet",   "2023-08-01",  "2023-10-31",
  "tran",  "2023-11-01",  "2023-12-31",
  "dry",   "2024-01-01",  "2024-04-30",
  "tran",  "2024-05-01",  "2024-06-30",
  "wet",   "2024-07-01",  "2024-08-31"
) %>%
  mutate(start = as.Date(start), end = as.Date(end))

# Exact sampling dates you provided
events_df <- tribble(
  ~date,         ~label,   ~type,
  "10/24/2023",  "eDNA1",  "eDNA",
  "11/25/2023",  "eDNA2",  "eDNA",
  "01/23/2024",  "eDNA3",  "eDNA",
  "03/20/2024",  "eDNA4",  "eDNA",
  "06/01/2024",  "eDNA5",  "eDNA",
  "08/12/2024",  "eDNA6",  "eDNA",
  "11/28/2023",  "T1",     "tissue",
  "02/28/2024",  "T2",     "tissue",
  "06/03/2024",  "T3",     "tissue",
  "08/26/2024",  "T4",     "tissue"
) %>%
  mutate(date = mdy(date))
library(ggplot2)
library(ggrepel)

ymax <- max(prcp_df$prcp_mm, na.rm = TRUE)

ggplot(prcp_df, aes(x = date, y = prcp_mm)) +

  # season shading behind bars
  geom_rect(
    data = season_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season),
    inherit.aes = FALSE,
    alpha = 0.12
  ) +

  geom_col() +

  # sampling lines (different linetypes for eDNA vs tissue)
  geom_vline(
    data = events_df,
    aes(xintercept = as.numeric(date), linetype = type),
    color = "orange",
    linewidth = 0.5
  ) +

  # top triangles at sampling dates
  geom_point(
    data = events_df,
    aes(x = date, y = ymax * 1.03, shape = type),
    inherit.aes = FALSE,
    size = 2.2,
    color = "orange",
    fill = "orange"
  ) +

  # labels above triangles
  ggrepel::geom_text_repel(
    data = events_df,
    aes(x = date, y = ymax * 1.03, label = label),
    inherit.aes = FALSE,
    nudge_y = ymax * 0.05,
    direction = "x",
    min.segment.length = 0,
    segment.color = "orange",
    size = 3.1
  ) +

  scale_x_date(
    limits = c(as.Date("2023-08-01"), as.Date("2024-08-31")),
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +

  scale_fill_manual(
    values = c(wet = "deepskyblue3", tran = "goldenrod2", dry = "gray60"),
    breaks = c("wet", "tran", "dry"),
    labels = c("Wet", "Transition", "Dry")
  ) +

  scale_linetype_manual(
    values = c(eDNA = "dashed", tissue = "dotted"),
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  scale_shape_manual(
    values = c(eDNA = 25, tissue = 24),  # filled triangles, slightly different
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  labs(x = "Date", y = "Daily precipitation (mm)", fill = NULL, linetype = NULL, shape = NULL) +
  coord_cartesian(xlim = as.Date(c("2023-08-01", "2024-08-31"))) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin = margin(10, 10, 25, 10))


ggsave("Prcp-daily-precipitation-with-sampling.png",, units = "in", width = 10, height = 5, dpi = 300)

library(ggplot2)
library(ggrepel)

ymax_w <- max(prcp_week$prcp_7d_mm, na.rm = TRUE)

p_prcp_w <- ggplot(prcp_week, aes(x = date, y = prcp_7d_mm)) +

  # season shading (behind the line)
  geom_rect(
    data = season_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season),
    inherit.aes = FALSE,
    alpha = 0.12
  ) +

  geom_line(linewidth = 0.7, na.rm = TRUE) +

  # sampling lines
  geom_vline(
    data = events_df,
    aes(xintercept = date, linetype = type),
    color = "orange",
    linewidth = 0.5
  ) +

  # small triangles at top
  geom_point(
    data = events_df,
    aes(x = date, y = ymax_w * 1.03, shape = type),
    inherit.aes = FALSE,
    size = 2.2,
    color = "orange",
    fill = "orange"
  ) +

  # labels
  ggrepel::geom_text_repel(
    data = events_df,
    aes(x = date, y = ymax_w * 1.03, label = label),
    inherit.aes = FALSE,
    nudge_y = ymax_w * 0.06,
    direction = "x",
    min.segment.length = 0,
    segment.color = "orange",
    size = 3.1
  ) +

  scale_x_date(
    limits = c(as.Date("2023-08-01"), as.Date("2024-08-31")),
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +

  scale_fill_manual(
    values = c(wet = "deepskyblue3", tran = "goldenrod2", dry = "gray60"),
    breaks = c("wet", "tran", "dry"),
    labels = c("Wet", "Transition", "Dry")
  ) +

  scale_linetype_manual(
    values = c(eDNA = "dashed", tissue = "dotted"),
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  scale_shape_manual(
    values = c(eDNA = 25, tissue = 24),
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  labs(x = "Date", y = "7-day cumulative precipitation (mm)", fill = NULL, linetype = NULL, shape = NULL) +

  coord_cartesian(ylim = c(0, ymax_w * 1.18), clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.margin = margin(10, 10, 25, 10)
  )

p_prcp_w


ggsave("Prcp-weekly-precipitation-with-sampling.png",, units = "in", width = 10, height = 5, dpi = 300)

------------------------------------------------------------------------
#sampling scheme over the temprature data



library(ggplot2)
library(ggrepel)

ymax_t <- max(temp_daily$temp_mean, na.rm = TRUE)

ggplot(temp_daily, aes(x = date, y = temp_mean, color = site)) +
  geom_line(linewidth = 0.8) +

  # sampling lines (drawn behind triangles/labels)
  geom_vline(
    data = events_df,
    aes(xintercept = date, linetype = type),
    inherit.aes = FALSE,
    color = "orange",
    linewidth = 0.5
  ) +

  # small triangles at top
  geom_point(
    data = events_df,
    aes(x = date, y = ymax_t * 1.01, shape = type),
    inherit.aes = FALSE,
    size = 2.0,
    color = "orange",
    fill = "orange"
  ) +

  # labels
  ggrepel::geom_text_repel(
    data = events_df,
    aes(x = date, y = ymax_t * 1.01, label = label),
    inherit.aes = FALSE,
    nudge_y = 0.12,              # in °C units (adjust if needed)
    direction = "x",
    min.segment.length = 0,
    segment.color = "orange",
    size = 3.1
  ) +

  scale_color_manual(
    values = c(Lum = "gray30", Tum = "steelblue"),
    labels = c("Luminao", "Tumon Bay")
  ) +

  scale_linetype_manual(
    values = c(eDNA = "dashed", tissue = "dotted"),
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  scale_shape_manual(
    values = c(eDNA = 25, tissue = 24),
    breaks = c("eDNA", "tissue"),
    labels = c("eDNA sampling", "Coral tissue sampling")
  ) +

  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +

  labs(
    x = "Date",
    y = "Daily mean temperature (°C)",
    color = "Site",
    linetype = NULL,
    shape = NULL
  ) +

  coord_cartesian(ylim = c(min(temp_daily$temp_mean, na.rm = TRUE), ymax_t * 1.08), clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(10, 10, 20, 10)
  )
