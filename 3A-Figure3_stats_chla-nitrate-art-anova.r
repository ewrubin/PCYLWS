library(ARTool)
library(emmeans)
library(car)
library(dplyr)

# -------------------------
# Chlorophyll a
# -------------------------
chla_dat <- read.table("chla-data.txt", sep = "\t", header = TRUE)

chla_dat$month.year <- factor(
  chla_dat$month.year,
  levels = c("Oct.23","Nov.23","Dec.23","Jan.24","Mar.24","May.24","Aug.24","Sep.24")
)
chla_dat$site <- factor(chla_dat$site, levels = c("Luminao","Tumon"))

# Optional assumption checks
shapiro.test(chla_dat$chla)
leveneTest(chla ~ interaction(site, month.year), data = chla_dat)

# ART model
chla_art <- art(chla ~ site * month.year, data = chla_dat)
anova(chla_art)

# Post hoc: site differences within each month
chla_emm <- emmeans(artlm(chla_art, "site:month.year"), ~ site | month.year)
pairs(chla_emm, adjust = "holm")

# Optional: month differences within each site
chla_emm2 <- emmeans(artlm(chla_art, "site:month.year"), ~ month.year | site)
pairs(chla_emm2, adjust = "holm")


# -------------------------
# Nitrate
# -------------------------
nit_dat <- read.table("nitrate-data.txt", sep = "\t", header = TRUE)

nit_dat$month.year <- factor(
  nit_dat$month.year,
  levels = c("Oct.23","Nov.23","Dec.23","Jan.24","Mar.24","May.24","Aug.24","Sep.24")
)
nit_dat$site <- factor(nit_dat$site, levels = c("Luminao","Tumon"))

# Optional assumption checks
shapiro.test(nit_dat$nit)
leveneTest(nit ~ interaction(site, month.year), data = nit_dat)

# ART model
nit_art <- art(nit ~ site * month.year, data = nit_dat)
anova(nit_art)

# Post hoc: site differences within each month
nit_emm <- emmeans(artlm(nit_art, "site:month.year"), ~ site | month.year)
pairs(nit_emm, adjust = "holm")

# Optional: month differences within each site
nit_emm2 <- emmeans(artlm(nit_art, "site:month.year"), ~ month.year | site)
pairs(nit_emm2, adjust = "holm")