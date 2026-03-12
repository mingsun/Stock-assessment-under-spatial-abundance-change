library(tidyverse)
library(ASAPplots)
library(r4ss)
library(RColorBrewer)

results.path <- paste0("manuscript/3. stock assessment with spatial abundance change/results/")


# 1. load results ----

## summer flounder ASAP ----
SF_ASAP <- list.files(paste0(results.path, "stock assessment/summer flounder/"), pattern = "^ASAP_summary.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()

SF_ASAP_base <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_FINAL.csv")


## surfclam SS ----
AS_SS_SSB <- list.files(paste0(results.path, "stock assessment/surfclam/"), pattern = "^SSB_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  mutate(Year = as.numeric(substr(Label, 5, 8)),
         X = "SSB")

AS_SS_SSB_base <- read.csv("results/stock assessment/surfclam/SSB_original.csv") %>%
  mutate(Year = as.numeric(substr(Label, 5, 8)),
         X = "SSB")


AS_SS_F <- list.files(paste0(results.path, "stock assessment/surfclam/"), pattern = "^F_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  rename(Year = Yr)

AS_SS_F_base <- read.csv("results/stock assessment/surfclam/F_original.csv") %>%
  rename(Year = Yr)


AS_SS_ABC <- list.files(paste0(results.path, "stock assessment/surfclam/"), pattern = "^ABC_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()  


## squid index-based ---- 
LS_INDEX <- list.files(paste0(results.path, "stock assessment/squid/"), pattern = "^Summary_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() 

LS_INDEX_base <- read.csv("results/stock assessment/squid/index-based summary_base.csv")


# ------------------------------------------------------------------------ #


pd <- position_dodge(width = 0.8) # gap between scenarios


# 2. SSB ----

## 2.1 summer flounder ----

SSB_SF_df <- SF_ASAP %>%
  select(Dataset, Scenario, q, Year, SSB, SSB_95_lo, SSB_95_hi) %>%
  left_join(SF_ASAP_base %>%
              select(Year, SSB_base = SSB, SSB_base_95_lo = SSB_95_lo, SSB_base_95_hi = SSB_95_hi),
            by = "Year")
 
write.csv(SSB_SF_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/1 SF_SSB.csv", row.names = FALSE)



## 2.2 surfclam ----

SSB_AS_df <- AS_SS_SSB %>%
  filter(Year <= 2023) %>%
  mutate(SSB_95_lo = Value - 1.96 * StdDev, SSB_95_hi = Value + 1.96 * StdDev) %>%
  rename(SSB = Value) %>%
  select(Dataset, Scenario, q, Year, SSB, SSB_95_lo, SSB_95_hi) %>%
  left_join(AS_SS_SSB_base %>%
              select(Year, SSB_base = Value, StdDev_base = StdDev) %>%
              mutate(SSB_base_95_lo = SSB_base - 1.96 * StdDev_base, SSB_base_95_hi = SSB_base + 1.96 * StdDev_base),
            by = "Year")

write.csv(SSB_AS_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/2 AS_SSB.csv", row.names = FALSE)



## 2.3 squid ----

SSB_LS_df <- LS_INDEX %>%
  select(Dataset, Scenario, q, Year, Annualized.B) %>%
  filter(Year != 2020) %>%
  left_join(LS_INDEX_base %>%
              select(Year, Annualized.B_base = Annualized.B),
            by = "Year") %>%
  rename(SSB = Annualized.B,
         SSB_base = Annualized.B_base) %>%
  add_column(SSB_95_lo = NA, SSB_95_hi = NA,
             SSB_base_95_lo = NA, SSB_base_95_hi = NA)
         

write.csv(SSB_LS_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/3 LS_SSB.csv", row.names = FALSE)



## 2.4 plot for all species ----

SSB_ALL_df <- rbind(cbind(SSB_SF_df, Species = "Summer Flounder"), 
                    cbind(select(SSB_AS_df, -StdDev_base), Species = "Atlantic Surfclam"),
                    cbind(SSB_LS_df, Species = "Longfin Squid")) %>%
  mutate(Species = factor(Species,
                          levels = c("Summer Flounder","Longfin Squid", "Atlantic Surfclam"))) %>%
  # terminal year only
  group_by(Species) %>%
  filter(Year == max(Year)) %>%
  ungroup() 

SSB_plot <- ggplot(SSB_ALL_df, aes(x = factor(q))) +
  # base 95% CI band
  geom_rect(aes(ymin = SSB_base_95_lo, ymax = SSB_base_95_hi), xmin = -Inf, xmax = Inf, fill = "grey80", alpha = 0.4) +
  # base median
  geom_hline(aes(yintercept = SSB_base), linetype = 2, color = "seagreen4", linewidth = 0.6) +
  # quantity plot
  geom_errorbar(aes(ymin = SSB_95_lo, ymax = SSB_95_hi, color = Dataset), width = 0.01, position = pd) +
  geom_point(aes(y = SSB, color = Dataset), position = pd, size = 5) +
  facet_grid(Species~Scenario, scale = "free_y") +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  labs(x = "q", y = "SSB[mt]") +
  theme_bw() +
  theme(legend.position = "none") 
  
png("manuscript/3. stock assessment with spatial abundance change/plot/4 SSB_quantity.png",  width = 5, height = 8, units = 'in', res = 800)
print(SSB_plot)
dev.off()


remove(AS_SS_SSB, AS_SS_SSB_base, SSB_SF_df, SSB_AS_df, SSB_LS_df, SSB_ALL_df, SSB_plot)

# ------------------------------------------------------------------------ #




# 3. F ----

## 3.1 summer flounder ----

F_SF_df <- SF_ASAP %>%
  select(Dataset, Scenario, q, Year, F = Freport, F_95_lo = Freport_95_lo, F_95_hi = Freport_95_hi) %>%
  left_join(SF_ASAP_base %>%
              select(Year, F_base = Freport, F_base_95_lo = Freport_95_lo, F_base_95_hi = Freport_95_hi),
            by = "Year") 

write.csv(F_SF_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/4 SF_F.csv", row.names = FALSE)



## 3.2 surfclam ----

F_AS_df <- AS_SS_F %>%
  filter(Year <= 2023) %>%
  mutate(F_95_lo = F_weighted - 1.96 * F_weighted_std, F_95_hi = F_weighted + 1.96 * F_weighted_std) %>%
  rename(F = F_weighted) %>%
  select(Dataset, Scenario, q, Year, F, F_95_lo, F_95_hi) %>%
  left_join(AS_SS_F_base %>%
              select(Year, F_base = F_weighted, F_base_StdDev = F_weighted_std) %>%
              mutate(F_base_95_lo = F_base - 1.96 * F_base_StdDev, F_base_95_hi = F_base + 1.96 * F_base_StdDev),
            by = "Year")

write.csv(F_AS_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/5 AS_F.csv", row.names = FALSE)



## 3.3 squid ----

F_LS_df <- LS_INDEX %>%
  select(Dataset, Scenario, q, Year, Annualized.Exploitation.Indices.final) %>%
  filter(Year != 2020) %>%
  left_join(LS_INDEX_base %>%
              select(Year, Annualized.Exploitation.Indices.final_base = Annualized.Exploitation.Indices.final),
            by = "Year") %>%
  rename(F = Annualized.Exploitation.Indices.final,
         F_base = Annualized.Exploitation.Indices.final_base,) %>%
  add_column(F_95_lo = NA, F_95_hi = NA,
             F_base_95_lo = NA, F_base_95_hi = NA)

write.csv(F_LS_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/6 LS_F.csv", row.names = FALSE)


## 3.4 plot for all species ----

F_ALL_df <- rbind(cbind(F_SF_df, Species = "Summer Flounder"), 
                  cbind(select(F_AS_df, -F_base_StdDev), Species = "Atlantic Surfclam"),
                  cbind(F_LS_df, Species = "Longfin Squid")) %>%
  mutate(Species = factor(Species,
                          levels = c("Summer Flounder","Longfin Squid", "Atlantic Surfclam"))) %>%
  # terminal year only
  group_by(Species) %>%
  filter(Year == max(Year)) %>%
  ungroup() 


F_plot <- ggplot(F_ALL_df, aes(x = factor(q))) +
  # base 95% CI band
  geom_rect(aes(ymin = F_base_95_lo, ymax = F_base_95_hi), xmin = -Inf, xmax = Inf, fill = "grey80", alpha = 0.4) +
  # base median
  geom_hline(aes(yintercept = F_base), linetype = 2, color = "seagreen4", linewidth = 0.6) +
  # quantity plot
  geom_errorbar(aes(ymin = F_95_lo, ymax = F_95_hi, color = Dataset), width = 0.01, position = pd) +
  geom_point(aes(y = F, color = Dataset), position = pd, size = 5) +
  facet_grid(Species~Scenario, scale = "free_y") +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  labs(x = "q", y = "F") +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/5 F_quantity.png",  width = 5, height = 8, units = 'in', res = 800)
print(F_plot)
dev.off()


# ------------------------------------------------------------------------ #

