# ============================================================
# Lesion size and progression analysis with colony nested in site
# ============================================================

# Packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(performance)

# ------------------------------------------------------------
# 1. Read and prepare data
# ------------------------------------------------------------
dat <- read.table("lesion_data.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Clean names if needed
names(dat) <- trimws(names(dat))

# Basic formatting
dat <- dat %>%
  mutate(
    site = factor(site, levels = c("Luminao", "Tumon")),
    month = factor(month, levels = c("Nov.23", "Feb.24", "May.24", "Aug.24")),
    colony = factor(colony),
    
    # nested colony identifier
    site_colony = interaction(site, colony, drop = TRUE),
    
    size = as.numeric(size),
    progr.rate = as.numeric(progr.rate),
    
    # binary activity variable
    active = ifelse(progr.rate > 0, 1, 0),
    active = factor(active, levels = c(0, 1))
  )

# Quick checks
str(dat)
summary(dat)
table(dat$site, dat$month)
table(dat$site, dat$colony)
table(dat$active)
with(dat, tapply(progr.rate, site, summary))

# ------------------------------------------------------------
# 2. Lesion size model
# ------------------------------------------------------------
# Size is positive and right-skewed, so log-transform
dat$log_size <- log(dat$size)

m_size <- lmer(log_size ~ site * month + (1 | site_colony), data = dat)

summary(m_size)
anova(m_size, type = 3)

# Diagnostics
check_model(m_size)

# Residual simulation-style checks
sim_size <- simulateResiduals(m_size)
plot(sim_size)

# Estimated marginal means
emm_size <- emmeans(m_size, ~ site | month)
pairs(emm_size, adjust = "tukey")

emm_size2 <- emmeans(m_size, ~ month | site)
pairs(emm_size2, adjust = "tukey")

# Back-transform model-estimated means to original scale
size_means <- emmeans(m_size, ~ site * month, type = "response")
size_means

# ------------------------------------------------------------
# 3. Lesion activity model (progression > 0 vs 0)
# ------------------------------------------------------------
# This tests whether the probability of an active lesion differs
# by site, month, or their interaction

m_active <- glmmTMB(
  as.numeric(as.character(active)) ~ site * month + (1 | site_colony),
  family = binomial,
  data = dat
)

summary(m_active)
car::Anova(m_active, type = 3)

# Diagnostics
sim_active <- simulateResiduals(m_active)
plot(sim_active)
testDispersion(sim_active)
testZeroInflation(sim_active)

# Estimated marginal means: probability of active lesions
emm_active <- emmeans(m_active, ~ site | month, type = "response")
emm_active

pairs(emm_active, adjust = "tukey")

emm_active2 <- emmeans(m_active, ~ month | site, type = "response")
emm_active2

pairs(emm_active2, adjust = "tukey")

# ------------------------------------------------------------
# 4. Progression rate among active lesions only
# ------------------------------------------------------------
# Restrict to lesions with progression > 0
dat_active <- dat %>%
  filter(progr.rate > 0) %>%
  mutate(log_prog = log(progr.rate))

nrow(dat_active)
table(dat_active$site, dat_active$month)

m_prog <- lmer(log_prog ~ site * month + (1 | site_colony), data = dat_active)

summary(m_prog)
anova(m_prog, type = 3)

# Diagnostics
check_model(m_prog)

sim_prog <- simulateResiduals(m_prog)
plot(sim_prog)

# Estimated marginal means
emm_prog <- emmeans(m_prog, ~ site | month, type = "response")
emm_prog
pairs(emm_prog, adjust = "tukey")

emm_prog2 <- emmeans(m_prog, ~ month | site, type = "response")
emm_prog2
pairs(emm_prog2, adjust = "tukey")

# ------------------------------------------------------------
# 5. Optional: descriptive summaries for plotting/reporting
# ------------------------------------------------------------
desc_size <- dat %>%
  group_by(site, month) %>%
  summarise(
    n = n(),
    mean_size = mean(size, na.rm = TRUE),
    sd_size = sd(size, na.rm = TRUE),
    median_size = median(size, na.rm = TRUE),
    .groups = "drop"
  )

desc_active <- dat %>%
  group_by(site, month) %>%
  summarise(
    n = n(),
    active_n = sum(progr.rate > 0, na.rm = TRUE),
    prop_active = mean(progr.rate > 0, na.rm = TRUE),
    .groups = "drop"
  )

desc_prog <- dat %>%
  filter(progr.rate > 0) %>%
  group_by(site, month) %>%
  summarise(
    n_active = n(),
    mean_prog = mean(progr.rate, na.rm = TRUE),
    sd_prog = sd(progr.rate, na.rm = TRUE),
    median_prog = median(progr.rate, na.rm = TRUE),
    .groups = "drop"
  )

desc_size
desc_active
desc_prog

