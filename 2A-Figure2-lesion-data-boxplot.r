library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse) 
library(dplyr)
library(ggpubr)
library(car)





data <- read.table("lesion_data.txt",sep="\t",header=TRUE)

data <- mutate(data,size = as.numeric(size))
data <- mutate(data,progr.rate = as.numeric(progr.rate))



data$month <- factor(data$month, levels=c("Nov.23","Feb.24","May.24","Aug.24"))
data$site <- factor(data$site, levels=c("Luminao","Tumon"))
data$colony <- factor(data$colony, levels=c("A","B","C","D","E","F","G","H"))

ggplot(data, aes(x = month, y = progr.rate, fill = month)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 1.0,
    alpha = 0.7,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(drop = TRUE) +
  facet_grid(~ site, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic(base_size = 14) +
  labs(
    x = "Month.year",
    y = expression("Lesion progression rate (" * cm^2 * " day"^-1 * ")")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
  )


ggsave("Fig.2B-lesion-prog-boxplot-site.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')

----------------------------------------------------------------------------------------------------


ggplot(data, aes(x = month, y = size, fill = month)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 0.8,
    alpha = 0.7,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(drop = TRUE) +
  facet_grid(~ site, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 14) +
  labs(
    x = "Month.year",
    y = expression("Lesion size (" * cm^2 * ")")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
  )

ggsave("Fig.2A-lesion-size-boxplot-site.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')



