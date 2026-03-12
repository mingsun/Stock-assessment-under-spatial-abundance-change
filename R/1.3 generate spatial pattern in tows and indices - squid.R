library(tidyverse)
source("manuscript/3. stock assessment with spatial abundance change/R/functions/stratified mean.R")

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# BTS  ----

catch_by_tow <- read.csv("results/indices for assessment/squid/catch_by_tow.csv") %>% # already calibrated
  mutate(NUMBER = CATCH_WT_CAL) %>%
  filter(YEAR != 2023)

# check CV inside and outside
CV_space <- catch_by_tow %>%
  group_by(YEAR, SEASON, OWF)  %>%
  summarise(MEAN = mean(NUMBER, na.rm = TRUE),
            SD =  sd(NUMBER, na.rm = TRUE),
            CV = SD / MEAN) %>%
  ungroup()

remove(CV_space)

 ## S1: inside up, outside down, total stable ----
  # 20%, 40%, 60%, 80%, 100%
  # the amount of increase inside are allocated proportionally to reduce abundance outside

  ### 1.1  generate tow data ----

diff_total <- data.frame()

for (q in q.pattern) {
  
  # adjust the inside abundance
  inside_temp <- catch_by_tow %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- inside_temp %>%
    group_by(YEAR, SEASON) %>%
    summarise(diff_total = sum(NUMBER_adj - NUMBER, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) %>%
    ungroup()
  
  # get the adjustment prop based on the abundance outside
  outside_ratio_temp <- catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    group_by(YEAR, SEASON) %>%
    summarise(TOTAL_ABUNDANCE = sum(NUMBER)) %>%
    left_join(diff_value_temp) %>%
    mutate(PROP = diff_total/TOTAL_ABUNDANCE) %>%
    select(-c(TOTAL_ABUNDANCE, diff_total)) 
  
  # adjust the rate and prop for year with only outside tows
  if ((sum(is.na(outside_ratio_temp$PROP))) != 0) {
    outside_ratio_temp <- outside_ratio_temp %>%
      group_by(YEAR) %>%
      mutate(Rate = ifelse(is.na(Rate), Rate[!is.na(Rate)][1], Rate),
             PROP = ifelse(is.na(PROP) & any(PROP == 0, na.rm = TRUE), 0, PROP)) %>%
      ungroup()
  }
  
  # adjust outside abundance
  outside_temp <- catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    left_join(outside_ratio_temp) %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 - (PROP)))  %>%
    ungroup() %>%
    select(-c(PROP))
  
  # combine the new catch by tow
  catch_by_tow_adjusted <- rbind(outside_temp, inside_temp) %>%
    arrange(ID) %>%
    # below is important: limit the change to the last 10 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NUMBER,
                                  TRUE ~ NUMBER_adj))
  
  write.csv(catch_by_tow_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S1/",
                   q,".csv"), row.names = FALSE)
  
  # append the diff_total for record
  diff_value_temp <- merge(diff_value_temp, outside_ratio_temp)
  diff_total <- rbind(diff_total, diff_value_temp)
  
}; remove(outside_temp, diff_value_temp, inside_temp, catch_by_tow_adjusted, outside_ratio_temp)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S1/abundance_diff.csv"
          , row.names = FALSE)

  ### ------------------------------------------------------------- ###


  ### 1.2 full dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S1/",
                                       q,".csv"))
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  # save by scenario
  write.csv(stratified_mean_N_df,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S1/FULL_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.IND_df <- read.csv("results/indices for assessment/squid/original.indices.csv")
  
  
  SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(SPRING_ratio, FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S1/ratio_FULL_",
                      q,".Rdata"))
  
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, SPRING_ratio, FALL_ratio)

  ### ------------------------------------------------------------- ###


  ### 1.3 WEE dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S1/",
                                       q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  # save by scenario
  write.csv(stratified_mean_N_df,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S1/WEE_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.IND_df <- read.csv("results/indices for assessment/squid/original.indices.csv")
  
  
  SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(SPRING_ratio, FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S1/ratio_WEE_",
                      q,".Rdata"))
  
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, SPRING_ratio, FALL_ratio)

  ### ------------------------------------------------------------- ###

 ## ------------------------------------------------------------- ##



 ## S2: inside up, outside stable, total up ----
  # 20%, 40%, 60%, 80%, 100%

  ### 2.1  generate tow data ----

diff_total <- data.frame()

for (q in q.pattern) {
  
  # adjust the inside abundance
  inside_temp <- catch_by_tow %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- inside_temp %>%
    group_by(YEAR, SEASON) %>%
    summarise(diff_total = sum(NUMBER_adj - NUMBER, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) 
  
  # adjust outside abundance
  outside_temp <- catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    mutate(NUMBER_adj = NUMBER) %>%
    add_column(Rate = q)
  
  # combine the new catch by tow
  catch_by_tow_adjusted <- rbind(outside_temp, inside_temp) %>%
    arrange(ID) %>%
    # below is important: limit the change to the last 10 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NUMBER,
                                  TRUE ~ NUMBER_adj))
  
  # append the diff_total for record
  diff_total <- rbind(diff_total, diff_value_temp)
  
  write.csv(catch_by_tow_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S2/",
                   q,".csv"), row.names = FALSE)
  
}; remove(outside_temp, diff_value_temp, inside_temp, catch_by_tow_adjusted)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S2/abundance_diff.csv"
          , row.names = FALSE)

  ### ------------------------------------------------------------- ###


  ### 2.2 full dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S2/",
                                       q,".csv"))
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  # save by scenario
  write.csv(stratified_mean_N_df,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S2/FULL_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.IND_df <- read.csv("results/indices for assessment/squid/original.indices.csv")
  
  
  SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(SPRING_ratio, FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S2/ratio_FULL_",
                      q,".Rdata"))
  
}; remove(catch_by_tow_temp, stratified_mean_N_df, SPRING_ratio, FALL_ratio)


  ### ------------------------------------------------------------- ###


  ### 2.3 WEE dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/squid/S2/",
                                       q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  # save by scenario
  write.csv(stratified_mean_N_df,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S2/WEE_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.IND_df <- read.csv("results/indices for assessment/squid/original.indices.csv")
  
  
  SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.IND_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(SPRING_ratio, FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/S2/ratio_WEE_",
                      q,".Rdata"))
  
}; remove(catch_by_tow_temp, stratified_mean_N_df, SPRING_ratio, FALL_ratio)


  ### ------------------------------------------------------------- ###


 ## ------------------------------------------------------------- ##



