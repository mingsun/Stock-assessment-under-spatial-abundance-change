library(tidyverse)
library(ASAPplots)
library(r4ss)


q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)


# 1. summer flounder ----

## 1.1 calc retro values ----

assessment.path <- "manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/"

q= 0.2; scenario = "S1"; dataset = "FULL"

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    for (dataset in c("FULL", "WEE")) {
      
      # get the variable needed
      wd <- paste0(assessment.path, dataset, "_", scenario , "_", q, "/")
      asap.name <- paste0("ASAP3_MTA2023_", dataset, "_", scenario ,"_", q ,"_000")  # need to use the suffix 000
      asap <- dget(paste0(wd, asap.name, ".rdat"))
      
      od <- "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/summer flounder/retrospective error/"
      
      # generate retros file
      PlotRetroWrapper(wd, asap.name, asap, save.plots = FALSE, od, plotf = "jpg") # note that the Mohn's rho is a five year average
      
      # extract and modify results file
      raw.results.path = paste0(od, "Retro.rho.values_ASAP3_MTA2023_", dataset, "_", scenario ,"_", q , "_000.csv")
      final.results.path = paste0(od, "ASAP_retro_error_", dataset, "_", scenario ,"_", q, ".csv")
      
      read.csv(raw.results.path) %>%
        add_column(Dataset = dataset, Scenario = scenario, q = q) %>%
        write.csv(final.results.path, row.names = FALSE)
      
      remove(wd, asap.name, asap, od, raw.results.path, final.results.path)

    }
  }
}; remove(q, scenario, dataset)

## ----------------------------------------------------------- ##



## 1.2 combine ----

SF_retro <- list.files(path = "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/summer flounder/retrospective error/",
                       pattern = "^ASAP_retro_error_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  filter(X == "Mohn.rho") %>%
  select(Dataset, Scenario, q, f.rho, ssb.rho) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Summer Flounder")

write.csv(SF_retro, 
          "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/summer flounder/retrospective error/retro_error_combined.csv", row.names = FALSE)

remove(SF_retro)

## ----------------------------------------------------------- ##




# ------------------------------------------------------------------------ #





# 2. squid ----

# For squid with index-based method, we cannot do likelihood component
# we calculate the SNR based on the loess fit

source("R/functions/4.4.2 squid SNR function.R")

LS_retro <- list.files(path = "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/squid/",
                       pattern = "^Summary_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() 


## 2.1 calc SNR for biomass ----

LS_retro_B <- LS_retro %>%
  filter(Year >= 1987 & !is.na(Annualized.B))  %>%
  select(Dataset, Scenario, q, Year, Annualized.B) %>%
  group_by(Dataset, Scenario, q) %>%
  group_modify(~ {
    snr.B <- calculate_snr(df = .x, index_col = "Annualized.B", year_col = "Year",
                         span = 0.5, method = "variance")
    tibble(SNR.B = snr.B)
  }) %>%  
  ungroup()

## ----------------------------------------------------------- ##


## 2.2 calc SNR for exploitation indices ----

LS_retro_F <- LS_retro %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))  %>%
  select(Dataset, Scenario, q, Year, Annualized.Exploitation.Indices.final) %>%
  group_by(Dataset, Scenario, q) %>%
  group_modify(~ {
    snr.F <- calculate_snr(df = .x, index_col = "Annualized.Exploitation.Indices.final", year_col = "Year",
                         span = 0.5, method = "variance")
    tibble(SNR.F = snr.F)
  }) %>%  
  ungroup()


## ----------------------------------------------------------- ##


## 2.3 combine ----

# format to match the other two models
LS_retro <- LS_retro_B %>%
  left_join(LS_retro_F, by = c("Dataset", "Scenario", "q")) %>%
  rename(f.rho = SNR.F, ssb.rho = SNR.B) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Longfin Squid")

write.csv(LS_retro, "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/squid/SNR.csv", row.names = FALSE)

remove(LS_retro, LS_retro_B, LS_retro_F)

## ----------------------------------------------------------- ##




# ------------------------------------------------------------------------ #




# 3. Surfclam ----

## 3.1 calc retro values  ----

assessment.path <- "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/"

q= 0.2; scenario = "S1"; dataset = "FULL"

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    for (dataset in c("FULL", "WEE")) {
      
      # do retro fitting
      retro(dir = paste0(assessment.path, dataset, "_", scenario, "_", q),
            oldsubdir = "",
            newsubdir = "retrospectives",
            subdirstart = "retro",
            years = 0:-5)
      
      
      # evaluate retro results
      retroModels <- SSgetoutput(dirvec = paste0(assessment.path, dataset, "_", scenario, "_", q , "/retrospectives/", paste("retro", 0:-5, sep = "")))
      
      retroSummary <- SSsummarize(retroModels)
      endyrvec <- retroSummary[["endyrs"]] + 0:-5
      Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)
      
      # generate a data frame below in a format consistent with the ASAP output
      Mohnrho_temp <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                                 explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                                 Dataset = dataset, Scenario = scenario, q = q)
      
      write.csv(Mohnrho_temp, paste0("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/surfclam/retrospective error/SS_retro_error_",
                dataset, "_", scenario, "_", q, ".csv"), row.names = FALSE)
      
      print(paste(dataset, scenario, q, "finished", sep = "_") )

      remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_temp)
      
    }
  }
}; remove(q, scenario, dataset)

## ----------------------------------------------------------- ##



## 3.2 combine ----

SC_retro <- list.files(path = "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/surfclam/retrospective error/",
                       pattern = "^SS_retro_error_.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  filter(X == "Mohn.rho") %>%
  select(Dataset, Scenario, q, f.rho, ssb.rho) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Atlantic Surfclam")

write.csv(SC_retro, 
          "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/surfclam/retrospective error/retro_error_combined.csv", row.names = FALSE)

remove(SC_retro)

## ----------------------------------------------------------- ##


# ------------------------------------------------------------------------ #




# 4. combine and plot ----

## 4.1 load base retro error ----

species.folder.list <- c("summer flounder", "squid", "surfclam")
path <- paste0("results/stock assessment/", species.folder.list, "/retrospective error/")

retro_base <- map_df(path, ~ read_csv(file.path(.x, "retro_error_combined.csv"))) %>%
  filter(Scenario == "Base") %>%
  select(Metric, Base_value = Value, Species) %>%
  mutate(Species = factor(Species,
                          levels = c("Summer Flounder", "Longfin Squid", "Surfclam"),
                          labels = c("Summer Flounder", "Longfin Squid", "Atlantic Surfclam"))) 

remove(species.folder.list, path)

## ----------------------------------------------------------- ##


## 4.2 combine with new values

SF_retro <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/summer flounder/retrospective error/retro_error_combined.csv")
LS_retro <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/squid/SNR.csv")
SC_retro <- read.csv("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/surfclam/retrospective error/retro_error_combined.csv")


retro_combine <- rbind(SF_retro, LS_retro, SC_retro) %>%
  left_join(retro_base) %>%
  mutate(Species = factor(Species,
                          levels = c("Summer Flounder", "Longfin Squid", "Atlantic Surfclam")))  %>%
  mutate(Metric = ifelse(Metric == "f.rho", "F", "SSB"))

## ----------------------------------------------------------- ##
  


## 4.3 plot ----

pd <- position_dodge(width = 0.8) # gap between scenarios

### 4.3.1 SF ----

retro_bar_SF <- ggplot(subset(retro_combine, Species == "Summer Flounder"), aes(x = factor(q))) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  geom_hline(aes(yintercept = Base_value), color = "seagreen4", linetype = 2) +
  geom_bar(aes(y = Value, fill = Dataset), position = pd, width = 0.6, stat = "identity", alpha = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  facet_wrap(Metric ~ Scenario, nrow = 1)  +
  labs(x = "q", y = "5-year Mohn's Rho") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_cartesian(ylim = c(-0.12,0.12))

png("manuscript/3. stock assessment with spatial abundance change/plot/7 retro error SF.png",  width = 10, height = 3, units = 'in', res = 800)
print(retro_bar_SF)
dev.off()


### ---------------------------------------------------------- ###



### 4.3.2 LS ----

retro_bar_LS <- ggplot(subset(retro_combine, Species == "Longfin Squid"), aes(x = factor(q))) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  geom_hline(aes(yintercept = Base_value), color = "seagreen4", linetype = 2) +
  geom_bar(aes(y = Value, fill = Dataset), position = pd, width = 0.6, stat = "identity", alpha = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  facet_wrap(Metric ~ Scenario, nrow = 1)  +
  labs(x = "q", y = "5-year Mohn's Rho") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_cartesian(ylim = c(-0.1, 1))

png("manuscript/3. stock assessment with spatial abundance change/plot/7 retro error LS.png",  width = 10, height = 3, units = 'in', res = 800)
print(retro_bar_LS)
dev.off()


### ---------------------------------------------------------- ###


### 4.3.3 SC ----

retro_bar_SC <- ggplot(subset(retro_combine, Species == "Atlantic Surfclam"), aes(x = factor(q))) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  geom_hline(aes(yintercept = Base_value), color = "seagreen4", linetype = 2) +
  geom_bar(aes(y = Value, fill = Dataset), position = pd, width = 0.6, stat = "identity", alpha = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  facet_wrap(Metric ~ Scenario, nrow = 1)  +
  labs(x = "q", y = "5-year Mohn's Rho") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_cartesian(ylim = c(-0.6, 1))

png("manuscript/3. stock assessment with spatial abundance change/plot/7 retro error SC.png",  width = 10, height = 3, units = 'in', res = 800)
print(retro_bar_SC)
dev.off()


### ---------------------------------------------------------- ###


## ----------------------------------------------------------- ##


# ------------------------------------------------------------------------ #


