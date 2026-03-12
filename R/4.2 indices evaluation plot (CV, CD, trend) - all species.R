library(tidyverse)
source("manuscript/3. stock assessment with spatial abundance change/R/functions/indices evaluation functions.R")

IND_SF <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/summer flounder_ind_eva.csv") %>%
  add_column(SPECIES = "SUMMER FLOUNDER", .before = 1) %>%
  filter(TS == "WHOLE") %>%
  mutate(ITEM = paste0(SCENARIO, "_", q))

IND_SC <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/surfclam_ind_eva.csv") %>%
  add_column(SPECIES = "ATLANTIC SURFCLAM", .before = 1)  %>%
  add_column(SEASON = "ANNUAL", .before = 5) %>%
  filter(TS == "WHOLE") %>%
  mutate(ITEM = paste0(SCENARIO, "_", q))

IND_LS <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/longfin squid_ind_eva.csv") %>%
  add_column(SPECIES = "LONGFIN SQUID", .before = 1)  %>%
  filter(TS == "WHOLE") %>%
  mutate(ITEM = paste0(SCENARIO, "_", q))




# 1. Variability (CV) ----


## 1.1 summer flounder ----

CV_SF_base <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv") %>%
  group_by(SEASON) %>%
  summarise(CV_base = calc_CV(SWAN)) %>%
  ungroup()
  
CV_SF <- IND_SF %>%
  select(ITEM, SCENARIO, q, SEASON, CV.FULL, CV.WEE) %>%
  rename(FULL = CV.FULL, WEE = CV.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CV") %>%
  left_join(CV_SF_base)


CV_SF_PLOT <- ggplot(CV_SF) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CV_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CV, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~SEASON, nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.04,0.9))

png("manuscript/3. stock assessment with spatial abundance change/plot/2 CV_SF.png",  width = 10, height = 3, units = 'in', res = 800)
print(CV_SF_PLOT)
dev.off()

remove(CV_SF_base, CV_SF, CV_SF_PLOT)

## ---------------------------------------------------- ##


## 1.2 surfclam ----

CV_SC_base <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv") %>%
  summarise(CV_base = calc_CV(VALUE))

CV_SC <- IND_SC %>%
  select(ITEM, SCENARIO, q, SEASON, CV.FULL, CV.WEE) %>%
  rename(FULL = CV.FULL, WEE = CV.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CV") %>%
  add_column(CV_base = CV_SC_base$CV_base)

CV_SC_PLOT <- ggplot(CV_SC) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CV_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CV, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~., nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.02,0.6))


png("manuscript/3. stock assessment with spatial abundance change/plot/2 CV_SC.png",  width = 5, height = 3, units = 'in', res = 800)
print(CV_SC_PLOT)
dev.off()

remove(CV_SC_base, CV_SC, CV_SC_PLOT)

## ---------------------------------------------------- ##


## 1.3 squid ----

CV_LS_base <- read.csv("results/indices for assessment/squid/original.indices.csv") %>%
  group_by(SEASON) %>%
  summarise(CV_base = calc_CV(STRATIFIED_MEAN_N)) %>%
  ungroup()

CV_LS <- IND_LS %>%
  select(ITEM, SCENARIO, q, SEASON, CV.FULL, CV.WEE) %>%
  rename(FULL = CV.FULL, WEE = CV.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CV") %>%
  left_join(CV_LS_base)

CV_LS_PLOT <- ggplot(CV_LS) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CV_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CV, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~SEASON, nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.04,0.6))

png("manuscript/3. stock assessment with spatial abundance change/plot/2 CV_LS.png",  width = 10, height = 3, units = 'in', res = 800)
print(CV_LS_PLOT)
dev.off()

remove(CV_LS_base, CV_LS, CV_LS_PLOT)

## ---------------------------------------------------- ##

# ----------------------------------------------------------------------------------------------------------------------------- #


# 2. Consecutive disparity (CD) ----


## 2.1 summer flounder ----

CD_SF_base <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv") %>%
  group_by(SEASON) %>%
  summarise(CD_base = calc_CD(SWAN)) %>%
  ungroup()

CD_SF <- IND_SF %>%
  select(ITEM, SCENARIO, q, SEASON, CD.FULL, CD.WEE) %>%
  rename(FULL = CD.FULL, WEE = CD.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CD") %>%
  left_join(CD_SF_base)

CD_SF_PLOT <- ggplot(CD_SF) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CD_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CD, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~SEASON, nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.05, 1.6))

png("manuscript/3. stock assessment with spatial abundance change/plot/3 CD_SF.png",  width = 10, height = 3, units = 'in', res = 800)
print(CD_SF_PLOT)
dev.off()

remove(CD_SF_base, CD_SF, CD_SF_PLOT)

## ---------------------------------------------------- ##


## 2.2 surfclam ----

CD_SC_base <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv") %>%
  summarise(CD_base = calc_CD(VALUE))

CD_SC <- IND_SC %>%
  select(ITEM, SCENARIO, q, SEASON, CD.FULL, CD.WEE) %>%
  rename(FULL = CD.FULL, WEE = CD.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CD") %>%
  add_column(CD_base = CD_SC_base$CD_base)

CD_SC_PLOT <- ggplot(CD_SC) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CD_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CD, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~., nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.03, 0.75))

png("manuscript/3. stock assessment with spatial abundance change/plot/3 CD_SC.png",  width = 5, height = 3, units = 'in', res = 800)
print(CD_SC_PLOT)
dev.off()

remove(CD_SC_base, CD_SC, CD_SC_PLOT)


## ---------------------------------------------------- ##


## 2.3 squid ----
CD_LS_base <- read.csv("results/indices for assessment/squid/original.indices.csv") %>%
  group_by(SEASON) %>%
  summarise(CD_base = calc_CD(STRATIFIED_MEAN_N)) %>%
  ungroup()

CD_LS <- IND_LS %>%
  select(ITEM, SCENARIO, q, SEASON, CD.FULL, CD.WEE) %>%
  rename(FULL = CD.FULL, WEE = CD.WEE) %>%
  pivot_longer(cols = c(FULL, WEE),
               names_to = "DATASET",
               values_to = "CD") %>%
  left_join(CD_LS_base)

CD_LS_PLOT <- ggplot(CD_LS) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(aes(yintercept = CD_base), color = "seagreen4", linetype = 2) +
  geom_bar(aes(x = factor(q), y = CD, fill = DATASET), position = position_dodge(width = 0.8), width = 0.6, stat = "identity") +
  facet_wrap(SCENARIO~SEASON, nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  # scale_x_discrete(labels = function(x) sub("^S[12]_", "", x)) +
  xlab("q") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.03, 0.75))

png("manuscript/3. stock assessment with spatial abundance change/plot/3 CD_LS.png",  width = 10, height = 3, units = 'in', res = 800)
print(CD_LS_PLOT)
dev.off()

remove(CD_LS_base, CD_LS, CD_LS_PLOT)


## ---------------------------------------------------- ##



# ----------------------------------------------------------------------------------------------------------------------------- #




# 3. Trend consistency (CD) ----


# no figure generated, a table was made manully



# ----------------------------------------------------------------------------------------------------------------------------- #
