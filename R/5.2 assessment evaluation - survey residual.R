library(tidyverse)
library(ASAPplots)
library(r4ss)
library(RColorBrewer)


q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)
pd <- position_dodge(width = 0.8)


# 1. summer flounder ----

## 1.1 extract residual ----

assessment.path <- "manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/"

survey_residual_df <- data.frame()

q= 0.2; scenario = "S1"; dataset = "FULL"

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    for (dataset in c("FULL", "WEE")) {
      
      # get the variable needed
      wd <- paste0(assessment.path, dataset, "_", scenario , "_", q, "/")
      asap.name <- paste0("ASAP3_MTA2023_", dataset, "_", scenario ,"_", q ,"_000")
      asap <- dget(paste0(wd, asap.name, ".rdat"))

      # 25 is BTS BIG spring
      residual.spring.temp <- data.frame(Dataset = dataset, Scenario = scenario, q = q,
                                       Season = "SPRING", Year = asap$index.year$ind25, Residual = asap$index.std.resid$ind25)
      
      # 26 is BTS BIG fall
      residual.fall.temp <- data.frame(Dataset = dataset, Scenario = scenario, q = q,
                                       Season = "FALL", Year = asap$index.year$ind26, Residual = asap$index.std.resid$ind26)
      
      
      survey_residual_df <- rbind(survey_residual_df, residual.spring.temp, residual.fall.temp)

    }; remove(wd, asap.name, asap, residual.spring.temp, residual.fall.temp)
  }
}; remove(q, scenario, dataset)


write.csv(survey_residual_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/7 SF_survey_residual.csv", row.names = FALSE)

## ----------------------------------------------------------- ##


## 1.2 combine with original assessment value ----

# get the variable needed
wd_base<- "assessment model/summer flounder/ASAP/base model/"
asap.name_base <-"ASAP3_MTA2023_FINAL"
asap_base <- dget(paste0(wd_base, asap.name_base, ".rdat"))

residual_base <- rbind(data.frame(Season = "SPRING", Year = asap_base$index.year$ind25, Residual_base = asap_base$index.std.resid$ind25),
                       data.frame(Season = "FALL", Year = asap_base$index.year$ind26, Residual_base = asap_base$index.std.resid$ind26)) %>%
  group_by(Season) %>%
  mutate(base_median = median(Residual_base),
         base_lo50 = quantile(Residual_base, 0.25,  names = FALSE),
         base_hi50 = quantile(Residual_base, 0.75,  names = FALSE),
         base_lo95 = quantile(Residual_base, 0.025, names = FALSE),
         base_hi95 = quantile(Residual_base, 0.975, names = FALSE)) %>%
  ungroup()

  
survey_residual_df <- read.csv("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/7 SF_survey_residual.csv") %>%
  left_join(residual_base, by = c("Season", "Year")) 

remove(wd_base, asap.name_base, asap_base, residual_base)

## ----------------------------------------------------------- ##


## 1.3 plot ----

SF_plot <- ggplot(survey_residual_df, aes(x = factor(q))) +
  # base 95% CI band
  geom_rect(aes(ymin = base_lo95, ymax = base_hi95), xmin = -Inf, xmax = Inf, fill = "grey80", alpha = 0.4) +
  # base 50% CI band
  geom_rect(aes(ymin = base_lo50, ymax = base_hi50), xmin = -Inf, xmax = Inf, fill = "grey60", alpha = 0.5) +
  # base median
  geom_hline(aes(yintercept = base_median), linetype = 2, color = "seagreen4", linewidth = 0.6) +
  # box plot of observed values
  geom_boxplot(aes(y = Residual, fill = Dataset), alpha = 0.7, width = 0.6, position = pd) +
  facet_wrap(Scenario~Season, nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  xlab("q") +
  theme_bw() +
  theme(legend.position = "none") 
  

png("manuscript/3. stock assessment with spatial abundance change/plot/6 survey residual SF .png",  width = 10, height = 3, units = 'in', res = 800)
print(SF_plot)
dev.off()

remove(survey_residual_df, SF_plot)

## ----------------------------------------------------------- ##


#  ------------------------------------------------------------------------------------------------ #




# 2. surfclam ----

## 2.1 extract residual ----

assessment.path <- "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/"

survey_residual_df <- data.frame()

q= 0.2; scenario = "S1"; dataset = "FULL"

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    for (dataset in c("FULL", "WEE")) {
      
      # assessment model path
      run_dir <- paste0(assessment.path, dataset, "_", scenario, "_", q)
      
      # create a list of quantities for the outputs
      assessment_results <- SS_output(run_dir, verbose = FALSE)
      
      # calculate survey residual
      residual.temp <- assessment_results$cpue %>%
        select(Fleet_name, Year = Yr, Obs, Exp, SE, Dev) %>%
        filter(Fleet_name == "MCDS") %>%
        mutate(Residual = (log(Obs) - log(Exp))/SE) %>%
        add_column(Dataset = dataset, Scenario = scenario, q = q, .before = 1)
      

      survey_residual_df <- rbind(survey_residual_df, residual.temp)
      
      remove(run_dir, assessment_results, residual.temp)
      
    }
  }
}; remove(q, scenario, dataset)


write.csv(survey_residual_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/8 SC_survey_residual.csv", row.names = FALSE)

## ----------------------------------------------------------- ##


## 2.2 combine with original assessment value ----

run_dir_base <- paste0("assessment model/surfclam/base model")
assessment_results_base <- SS_output(run_dir_base, verbose = FALSE)

residual_base <- assessment_results_base$cpue %>%
  select(Fleet_name, Year = Yr, Obs, Exp, SE, Dev) %>%
  filter(Fleet_name == "MCDS") %>%
  mutate(Residual_base = (log(Obs) - log(Exp))/SE) %>%
  group_by(Fleet_name) %>%
  summarise(base_median = median(Residual_base),
         base_lo50 = quantile(Residual_base, 0.25,  names = FALSE),
         base_hi50 = quantile(Residual_base, 0.75,  names = FALSE),
         base_lo95 = quantile(Residual_base, 0.025, names = FALSE),
         base_hi95 = quantile(Residual_base, 0.975, names = FALSE)) %>%
  ungroup()
  
survey_residual_df <- read.csv("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/8 SC_survey_residual.csv") %>%
  left_join(residual_base, by = c("Fleet_name")) 

remove(run_dir_base, assessment_results_base, residual_base)

## ----------------------------------------------------------- ##


## 2.3 plot ----

SC_plot <- ggplot(survey_residual_df, aes(x = factor(q))) +
  # base 95% CI band
  geom_rect(aes(ymin = base_lo95, ymax = base_hi95), xmin = -Inf, xmax = Inf, fill = "grey80", alpha = 0.4) +
  # base 50% CI band
  geom_rect(aes(ymin = base_lo50, ymax = base_hi50), xmin = -Inf, xmax = Inf, fill = "grey60", alpha = 0.5) +
  # base median
  geom_hline(aes(yintercept = base_median), linetype = 2, color = "seagreen4", linewidth = 0.6) +
  # box plot of observed values
  geom_boxplot(aes(y = Residual, fill = Dataset), alpha = 0.7, width = 0.6, position = pd) +
  facet_wrap(Scenario~., nrow = 1) +
  scale_fill_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  xlab("q") +
  theme_bw() +
  theme(legend.position = "none") 


png("manuscript/3. stock assessment with spatial abundance change/plot/6 survey residual SC .png",  width = 5, height = 3, units = 'in', res = 800)
print(SC_plot)
dev.off()

remove(survey_residual_df, SC_plot)

## ----------------------------------------------------------- ##



# 3. squid ----

# for squid assessed with an index-based method
# here we cannot look at model generated residuals
# instead we look at the nominal residual in random stratified mean values 
# this residual is log scale 


## 3.1 extract new indices ----

indices.path <- "manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/"

survey_residual_df <- data.frame()

q= 0.2; scenario = "S1"; season = "FALL"

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    INDICES_full <- read.csv(paste0(indices.path, scenario, "/FULL_", q, ".csv")) %>%
      add_column(Dataset = "FULL", Scenario = scenario, q = q, .before = 1)

    INDICES_WEE <- read.csv(paste0(indices.path, scenario, "/WEE_", q, ".csv")) %>%
      add_column(Dataset = "WEE", Scenario = scenario, q = q, .before = 1)
    
    residual.temp <- INDICES_full %>%
      bind_rows(INDICES_WEE) %>%
      select(Dataset, Scenario, q, Year = YEAR, Season = SEASON, STRATIFIED_MEAN_N) %>%
      pivot_wider(names_from = Dataset,
                  values_from = c(STRATIFIED_MEAN_N)) %>%
      mutate(Residual = log(WEE) - log(FULL))

      survey_residual_df <- rbind(survey_residual_df, residual.temp)

  }; remove(INDICES_full, INDICES_WEE)
}; remove(q, scenario, season)

write.csv(survey_residual_df, "manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/9 LS_survey_residual.csv", row.names = FALSE)

## ----------------------------------------------------------- ##


## ----------------------------------------------------------- ##


## 3.2 plot ----

survey_residual_df <- read.csv("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/9 LS_survey_residual.csv") %>%
  filter(Year >= 2013) %>%
  group_by(Scenario, Season, q) %>%
  mutate(med_resid = median(Residual, na.rm = TRUE),
         med_sign  = ifelse(med_resid >= 0, "Positive", "Negative")) %>%
  ungroup()

LS_plot <- ggplot(survey_residual_df, aes(x = factor(q))) +
  geom_hline(aes(yintercept = 0), linetype = 2, color = "seagreen4", linewidth = 0.6) +
  # box plot of observed values
  geom_boxplot(aes(y = Residual, color = med_sign), alpha = 0.7, width = 0.6, position = pd) +
  facet_wrap(Scenario~Season, nrow = 1) +
  scale_color_manual(values = c("Negative" =  "steelblue3", "Positive" = "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  xlab("q") +
  theme_bw() +
  theme(legend.position = "none") 


png("manuscript/3. stock assessment with spatial abundance change/plot/6 survey residual LS .png",  width = 10, height = 3, units = 'in', res = 800)
print(LS_plot)
dev.off()

remove(INDICES_base, survey_residual_df, LS_plot)

## ----------------------------------------------------------- ##

#  ------------------------------------------------------------------------------------------------ #





