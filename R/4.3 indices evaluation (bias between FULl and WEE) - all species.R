# in this script, we evaluate the relationship between index bias (Full vs WEE) and survey residual by year


library(tidyverse)
library(broom)

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# 1. summer flounder ----

Bias_df <- data.frame()

q = 0.2; scenario = "S2"; season = "FALL"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    # extract index bias
    SWAN_full <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/", scenario, "/BIG_FULL_", 
                                 q, ".csv")) %>%
      select(Year = YEAR, Season = SEASON, FULL.SWAN = SWAN)
    
    SWAN_WEE <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/", scenario, "/BIG_WEE_", 
                                q, ".csv")) %>%
      select(Year = YEAR, Season = SEASON, WEE.SWAN = SWAN)
    
    # calculate index bias
    SWAN <- merge(SWAN_full, SWAN_WEE) %>%
      group_by(Year, Season) %>%
      mutate(Bias = (WEE.SWAN - FULL.SWAN)/mean(FULL.SWAN)) %>%
      ungroup() %>%
      add_column(Scenario = scenario, q = q, .before = 1)
    
    Bias_df <- rbind(Bias_df, SWAN) 

  }
}; remove(q, scenario, season, SWAN_full, SWAN_WEE, SWAN)

write.csv(Bias_df,
          paste0("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/summer flounder_ind_bias.csv"), row.names = FALSE)


# extract survey residual results and calculate relative change
survey_residual_df <- read.csv("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/7 SF_survey_residual.csv") %>%
  pivot_wider(names_from = Dataset, values_from = Residual) %>%
  group_by(Scenario, q, Year, Season) %>%
  mutate(Residual.change = (WEE - FULL)/mean(FULL)) %>%
  ungroup() %>% 
  left_join(Bias_df) %>%
  select(Scenario, q, Year, Season, Bias, Residual.change)  %>%
  filter(between(Residual.change, -0.5, 0.5)) # remove outlines

# calculate linear relationship between Bias and Residual  
lm_annot  <- survey_residual_df %>%
  left_join(Bias_df) %>%
  group_by(Scenario, q, Season) %>%
  do({
    m <- lm(Residual.change ~ Bias, data = .)
    td <- tidy(m)
    gl <- glance(m)
    
    data.frame(
      slope   = td$estimate[td$term == "Bias"],
      p_value = td$p.value[td$term == "Bias"],
      r2      = gl$r.squared,
      n       = nrow(.)
    )
  }) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "not significant"
    ),
    label = paste0("lm: ", sig, "\n",
                   "p=", signif(p_value, 2), "\n",
                   "R²=", round(r2, 2), ", n=", n)
  ) 
  

  

## ----------------------------------------------------------- ##


##  plot ----

bias_residual_plot <- ggplot(survey_residual_df, aes(x = Bias, y = Residual.change)) +
  geom_point(aes(color = factor(q))) +
  # regression line in each facet (lm fit within facet data automatically)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(Scenario ~ Season ~ q, nrow = 4, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none") +
  # coord_cartesian(ylim = c(-0.5, 0.5)) +
  # panel label (top-left). inherits aes = FALSE so it doesn't need Bias/Residual.change columns
  geom_text(
    data = lm_annot,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.1,
    size = 3,
    inherit.aes = FALSE
  )


png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S7 bias_survey residual SF .png",  width = 8, height = 8, units = 'in', res = 800)
print(bias_residual_plot)
dev.off()

remove(Bias_df, lm_annot, survey_residual_df, bias_residual_plot)

## ----------------------------------------------------------- ##

# ----------------------------------------------------------------------------------------------------------------------------- #



# 2. surfclam ----

Bias_df <- data.frame()

q = 0.2; scenario = "S2"; season = "FALL"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    # extract index bias
    INDICES_full_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/", scenario, "/MCDS_FULL_", 
                                         q, ".csv")) %>%
      select(Year = YEAR, FULL.Value = VALUE)
    
    INDICES_WEE_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/", scenario, "/MCDS_WEE_", 
                                        q, ".csv")) %>%
      select(Year = YEAR, WEE.Value = VALUE)
    
    # calculate index bias
    INDICES_temp <- merge(INDICES_full_temp, INDICES_WEE_temp) %>%
      group_by(Year) %>%
      mutate(Bias = (WEE.Value - FULL.Value)/mean(FULL.Value)) %>%
      ungroup() %>%
      add_column(Scenario = scenario, q = q, .before = 1)
    
    Bias_df <- rbind(Bias_df, INDICES_temp) 
    
  }
}; remove(q, scenario, INDICES_full_temp, INDICES_WEE_temp, INDICES_temp)

write.csv(Bias_df,
          paste0("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/atlantic surfclam_ind_bias.csv"), row.names = FALSE)


# extract survey residual results and calculate relative change
survey_residual_df <- read.csv("manuscript/3. stock assessment with spatial abundance change/plot/SA/tables/8 SC_survey_residual.csv") %>%
  select(-c(Obs, Exp, SE,Dev)) %>%
  pivot_wider(names_from = Dataset, values_from = Residual) %>%
  group_by(Scenario, q, Year) %>%
  mutate(Residual.change = (WEE - FULL)/mean(FULL)) %>%
  ungroup() %>% 
  left_join(Bias_df) %>%
  select(Scenario, q, Year, Bias, Residual.change)  %>%
  filter(between(Residual.change, -10, 10)) # remove outlines

# calculate linear relationship between Bias and Residual  
lm_annot  <- survey_residual_df %>%
  left_join(Bias_df) %>%
  group_by(Scenario, q) %>%
  do({
    m <- lm(Residual.change ~ Bias, data = .)
    td <- tidy(m)
    gl <- glance(m)
    
    data.frame(
      slope   = td$estimate[td$term == "Bias"],
      p_value = td$p.value[td$term == "Bias"],
      r2      = gl$r.squared,
      n       = nrow(.)
    )
  }) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "not significant"
    ),
    label = paste0("lm: ", sig, "\n",
                   "p=", signif(p_value, 2), "\n",
                   "R²=", round(r2, 2), ", n=", n)
  ) 




## ----------------------------------------------------------- ##


##  plot ----

bias_residual_plot <- ggplot(survey_residual_df, aes(x = Bias, y = Residual.change)) +
  geom_point(aes(color = factor(q))) +
  # regression line in each facet (lm fit within facet data automatically)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(Scenario ~ q, nrow = 2, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none") +
  # coord_cartesian(ylim = c(-10, 10)) +
  # panel label (top-left). inherits aes = FALSE so it doesn't need Bias/Residual.change columns
  geom_text(
    data = lm_annot,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.1,
    size = 3,
    inherit.aes = FALSE
  )


png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S7 bias_survey residual SC.png",  width = 8, height = 4, units = 'in', res = 800)
print(bias_residual_plot)
dev.off()

remove(survey_residual_df, SF_plot)

## ----------------------------------------------------------- ##

# ----------------------------------------------------------------------------------------------------------------------------- #

