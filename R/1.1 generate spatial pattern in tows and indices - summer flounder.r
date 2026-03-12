library(tidyverse)
source("manuscript/3. stock assessment with spatial abundance change/R/functions/stratified mean.R")

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# ALB (1982 - 2008) ----

# not needed because it does not cover the most recent 10 years

# ----------------------------------------------------------------------- #



# BIG (2009 - 2022) ----

BIG_catch_by_tow <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv")

# check CV inside and outside
CV_space <- BIG_catch_by_tow %>%
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
  BIG_inside_temp <- BIG_catch_by_tow %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- BIG_inside_temp %>%
    group_by(YEAR, SEASON) %>%
    summarise(diff_total = sum(NUMBER_adj - NUMBER, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) %>%
    ungroup()
  
  # manually add 2017 FALL value because 2017 has no inside tows
  diff_value_temp <- diff_value_temp %>%
    bind_rows(data.frame(YEAR = 2017, SEASON = "FALL", Rate = q, diff_total = 0, PROP = 0))
  
  # get the adjustment prop based on the abundance outside
  BIG_outside_ratio_temp <- BIG_catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    group_by(YEAR, SEASON) %>%
    summarise(TOTAL_ABUNDANCE = sum(NUMBER)) %>%
    left_join(diff_value_temp) %>%
    mutate(PROP = diff_total/TOTAL_ABUNDANCE) %>%
    select(-c(TOTAL_ABUNDANCE, diff_total)) 
  
  # adjust outside abundance
  BIG_outside_temp <- BIG_catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    left_join(BIG_outside_ratio_temp) %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 - (PROP)))  %>%
    ungroup() %>%
    select(-c(PROP))
  
  # combine the new catch by tow
  BIG_catch_by_tow_adjusted <- rbind(BIG_outside_temp, BIG_inside_temp) %>%
    arrange(ID) %>%
    # below is important: limit the change to the last 10 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NUMBER,
                                  TRUE ~ NUMBER_adj))
  
  # save 
  write.csv(BIG_catch_by_tow_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S1/BIG_",
                   q,".csv"), row.names = FALSE)

  
  # append the diff_total for record
  diff_value_temp <- merge(diff_value_temp, BIG_outside_ratio_temp)
  diff_total <- rbind(diff_total, diff_value_temp)
  
}; remove(BIG_outside_temp, diff_value_temp, BIG_inside_temp, BIG_catch_by_tow_adjusted, BIG_outside_ratio_temp)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S1/BIG_abundance_diff.csv"
          , row.names = FALSE)

### ------------------------------------------------------------- ###


  ### 1.2 full dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S1/BIG_",
                                       q,".csv"))
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  ## calculate Swept Area Numbers adjust by swept area 
    # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
    # survey coverage area: spring 27855, fall 17924
  
  AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)
  
  SWAN_temp <- stratified_mean_N_df %>%
    select(YEAR, SEASON, STRATIFIED_MEAN_N, CV) %>%
    arrange(SEASON) %>%
    left_join(AREA_df) %>%
    mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)
  
  
  # save by scenario
  write.csv(SWAN_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S1/BIG_FULL_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")
  

  BIG_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  BIG_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(BIG_SPRING_ratio, BIG_FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S1/BIG_ratio_FULL_",
                      q,".Rdata"))
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, BIG_SPRING_ratio, BIG_FALL_ratio)

### ------------------------------------------------------------- ###



   ### 1.3 WEE dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S1/BIG_",
                                       q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  ## calculate Swept Area Numbers adjust by swept area 
  # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
  # survey coverage area: spring 27855, fall 17924
  
  AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)
  
  SWAN_temp <- stratified_mean_N_df %>%
    select(YEAR, SEASON, STRATIFIED_MEAN_N, CV) %>%
    arrange(SEASON) %>%
    left_join(AREA_df) %>%
    mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)
  
  
  # save by scenario
  write.csv(SWAN_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S1/BIG_WEE_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")
  
  
  BIG_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  BIG_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(BIG_SPRING_ratio, BIG_FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S1/BIG_ratio_WEE_",
                      q,".Rdata"))
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, BIG_SPRING_ratio, BIG_FALL_ratio)


  ### ------------------------------------------------------------- ###


## ------------------------------------------------------------- ##



 ## S2: inside up, outside stable, total up ----
  # 20%, 40%, 60%, 80%, 100%

  ### 2.1  generate tow data ----

diff_total <- data.frame()

for (q in q.pattern) {
  
  # adjust the inside abundance
  BIG_inside_temp <- BIG_catch_by_tow %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR, SEASON) %>%
    mutate(NUMBER_adj = NUMBER * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- BIG_inside_temp %>%
    group_by(YEAR, SEASON) %>%
    summarise(diff_total = sum(NUMBER_adj - NUMBER, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) 
  
  # manually add 2017 FALL value because 2017 has no inside tows
  diff_value_temp <- diff_value_temp %>%
    bind_rows(data.frame(YEAR = 2017, SEASON = "FALL", Rate = q, diff_total = 0, PROP = 0))
  
  # adjust outside abundance
  BIG_outside_temp <- BIG_catch_by_tow %>%
    filter(OWF == "OUTSIDE") %>%
    mutate(NUMBER_adj = NUMBER) %>%
    add_column(Rate = q)
  
  # combine the new catch by tow
  BIG_catch_by_tow_adjusted <- rbind(BIG_outside_temp, BIG_inside_temp) %>%
    arrange(ID) %>%
    # below is important: limit the change to the last 5 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NUMBER,
                                  TRUE ~ NUMBER_adj))
  
  # append the diff_total for record
  diff_total <- rbind(diff_total, diff_value_temp)
  
  write.csv(BIG_catch_by_tow_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S2/BIG_",
                   q,".csv"), row.names = FALSE)
  
}; remove(BIG_outside_temp, diff_value_temp, BIG_inside_temp, BIG_catch_by_tow_adjusted)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S2/BIG_abundance_diff.csv"
          , row.names = FALSE)

  ### ------------------------------------------------------------- ###


  ### 2.2 full dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S2/BIG_",
                                       q,".csv"))
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  ## calculate Swpt Area Numbers adjust by swept area 
  # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
  # survey coverage area: spring 27855, fall 17924
  
  AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)
  
  SWAN_temp <- stratified_mean_N_df %>%
    select(YEAR, SEASON, STRATIFIED_MEAN_N, CV) %>%
    arrange(SEASON) %>%
    left_join(AREA_df) %>%
    mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)
  
  
  # save by scenario
  write.csv(SWAN_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S2/BIG_FULL_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")
  
  
  BIG_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  BIG_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(BIG_SPRING_ratio, BIG_FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S2/BIG_ratio_FULL_",
                      q,".Rdata"))
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, BIG_SPRING_ratio, BIG_FALL_ratio)


  ### ------------------------------------------------------------- ###


  ### 2.3 full dataset abundance indices ----

for (q in q.pattern) {
  
  catch_by_tow_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/summer flounder/S2/BIG_",
                                       q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # calculate stratified mean
  stratified_mean_N_df <- calc_stratified_mean(catch_by_tow_temp)
  
  
  ## calculate Swpt Area Numbers adjust by swept area 
  # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
  # survey coverage area: spring 27855, fall 17924
  
  AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)
  
  SWAN_temp <- stratified_mean_N_df %>%
    select(YEAR, SEASON, STRATIFIED_MEAN_N, CV) %>%
    arrange(SEASON) %>%
    left_join(AREA_df) %>%
    mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)
  
  
  # save by scenario
  write.csv(SWAN_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S2/BIG_WEE_",
                   q,".csv"), row.names = FALSE)
  
  
  # extract indices ratio relative to original
  original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")
  
  
  BIG_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
  BIG_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.BIG_df, SEASON == "FALL")$STRATIFIED_MEAN_N
  
  
  save(BIG_SPRING_ratio, BIG_FALL_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/S2/BIG_ratio_WEE_",
                      q,".Rdata"))
  
}; remove(q, catch_by_tow_temp, stratified_mean_N_df, BIG_SPRING_ratio, BIG_FALL_ratio)


### ------------------------------------------------------------- ###


 ## ------------------------------------------------------------- ##

