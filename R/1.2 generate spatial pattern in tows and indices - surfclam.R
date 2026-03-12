library(tidyverse)
source("manuscript/3. stock assessment with spatial abundance change/R/functions/stratified mean.R")

# RDtrendS: total N stratified mean, 1982-2011, assuming numbers per m2 consistent
# RDscaleS: total N stratified mean times the area size, 1997-2011
# MCDS: similar to RDscaleS but based on sensor tow distance (DISTANCEDETAIL = SENSOR), 2012-2022

stra_area_df <- read.csv("results/stratified.mean.indices/surfclam/strata.area.csv") %>% # strata area as mean
  filter(REGION == "SVAtoSNE") %>%
  mutate(WEIGHT = AREASQNM/sum(AREASQNM))

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)


# 1. RDtrendS ----

# not needed because it does not cover the most recent 10 years

# ------------------------------------------------------ #



# 2. RDscaleS ----

# not needed because it does not cover the most recent 10 years

# ------------------------------------------------------ #





# 3. MCDS ----

tow_MCDS_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")%>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022), REGION == "SVAtoSNE") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # using a new id here because the overlay ID have different stratum coding system
  mutate(OWF = case_when(ID.temp %in% unique(AS_QQ_overlay_df$ID.temp) ~ "INSIDE",
                         !ID.temp %in% unique(AS_QQ_overlay_df$ID.temp) ~ "OUTSIDE")) %>%
  filter(!is.na(NPERTOW))


 ## S1: inside up, outside down, total stable ----
    # 20%, 40%, 60%, 80%, 100%
    # the amount of increase inside are allocated proportionally to reduce abundance outside 


  ### 3.1.1 generate tow data ----

diff_total <- data.frame()

for (q in q.pattern) {
  
  # adjust the inside abundance
  inside_temp <- tow_MCDS_df %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR) %>%
    mutate(NUMBER_adj = NPERTOW * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- inside_temp %>%
    group_by(YEAR) %>%
    summarise(diff_total = sum(NUMBER_adj - NPERTOW, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) %>%
    ungroup()
  
  # get the adjustment prop based on the abundance outside
  outside_ratio_temp <- tow_MCDS_df %>%
    filter(OWF == "OUTSIDE") %>%
    group_by(YEAR) %>%
    summarise(TOTAL_ABUNDANCE = sum(NPERTOW)) %>%
    left_join(diff_value_temp) %>%
    mutate(PROP = diff_total/TOTAL_ABUNDANCE) %>%
    select(-c(TOTAL_ABUNDANCE, diff_total)) 
  
  # adjust outside abundance
  outside_temp <- tow_MCDS_df %>%
    filter(OWF == "OUTSIDE") %>%
    left_join(outside_ratio_temp) %>%
    group_by(YEAR) %>%
    mutate(NUMBER_adj = NPERTOW * (1 - (PROP)))  %>%
    ungroup() %>%
    select(-c(PROP))
  
  # combine the new catch by tow
  tow_MCDS_df_adjusted <- rbind(outside_temp, inside_temp) %>%
    arrange(ID.temp) %>%
    # below is important: limit the change to the last 10 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NPERTOW,
                                  TRUE ~ NUMBER_adj))
  
  write.csv(tow_MCDS_df_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S1/MCDS_tow_",
                   q,".csv"), row.names = FALSE)
  
  # append the diff_total for record
  diff_value_temp <- merge(diff_value_temp, outside_ratio_temp)
  diff_total <- rbind(diff_total, diff_value_temp)
  
}; remove(outside_temp, diff_value_temp, inside_temp, tow_MCDS_df_adjusted, outside_ratio_temp, q)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S1/MCDS_tow_abundance_diff.csv"
          , row.names = FALSE)

remove(diff_total)

### ------------------------------------------------------------- ###



  ### 3.1.2 full dataset abundance indices ----

for (q in q.pattern) { 
  
  tow_MCDS_df_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S1/MCDS_tow_",
                                          q,".csv"))
  
  # mean abundance by stratum
  mean_N_stratum_df <- tow_MCDS_df_temp %>%
    group_by(YEAR, REGION, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(NUMBER_adj)) %>% # mean within a strata
    mutate(VAR_STRATUM = var(NUMBER_adj)) %>% # variance by stratum
    ungroup() %>% 
    select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
    distinct() %>%
    arrange(YEAR, REGION)
  
  # stratified mean
  stratified_mean_N_df <- mean_N_stratum_df %>%
    group_by(YEAR, REGION) %>%
    summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
              VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
              SE = sqrt(VAR)) %>%   # standard deviance
    mutate(CV = SE / STRATIFIED_MEAN_N) %>%
    ungroup()
  
  # adjust to total abundance based on population area
  MCDS_indices_temp <- tow_MCDS_df_temp %>%
    select(YEAR, STRATUM, AREASQNM) %>%
    distinct()  %>%
    group_by(YEAR) %>%
    summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
    ungroup() %>%
    right_join(stratified_mean_N_df) %>%
    mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
           SE = round(total.AREASQM * SE/1000, -1),
           log_SE = sqrt(log(1 + CV^2))) %>%
    select(YEAR, VALUE, VAR, SE, CV, log_SE)
  
  # save by scenario
  write.csv(MCDS_indices_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S1/MCDS_FULL_",
                   q,".csv"), row.names = FALSE)
  
  # extract indices ratio relative to original
  original.MCDS_df <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv")
  
  MCDS_VALUE_ratio <- MCDS_indices_temp$VALUE/original.MCDS_df$VALUE
  MCDS_STDERR_ratio <- MCDS_indices_temp$log_SE/original.MCDS_df$log_SE
  
  save(MCDS_VALUE_ratio, MCDS_STDERR_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S1/MCDS_ratio_FULL_",
                      q,".Rdata"))
  
}; remove(q, tow_MCDS_df_temp, mean_N_stratum_df, stratified_mean_N_df, MCDS_indices_temp, 
          original.MCDS_df, MCDS_VALUE_ratio, MCDS_STDERR_ratio)


  ### ------------------------------------------------------------- ###




  ### 3.1.3 WEE dataset abundance indices ----

for (q in q.pattern) { 
  
  tow_MCDS_df_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S1/MCDS_tow_",
                                      q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # mean abundance by stratum
  mean_N_stratum_df <- tow_MCDS_df_temp %>%
    group_by(YEAR, REGION, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(NUMBER_adj)) %>% # mean within a strata
    mutate(VAR_STRATUM = var(NUMBER_adj)) %>% # variance by stratum
    ungroup() %>% 
    select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
    distinct() %>%
    arrange(YEAR, REGION)
  
  # stratified mean
  stratified_mean_N_df <- mean_N_stratum_df %>%
    group_by(YEAR, REGION) %>%
    summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
              VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
              SE = sqrt(VAR)) %>%   # standard deviance
    mutate(CV = SE / STRATIFIED_MEAN_N) %>%
    ungroup()
  
  # adjust to total abundance based on population area
  MCDS_indices_temp <- tow_MCDS_df_temp %>%
    select(YEAR, STRATUM, AREASQNM) %>%
    distinct()  %>%
    group_by(YEAR) %>%
    summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
    ungroup() %>%
    right_join(stratified_mean_N_df) %>%
    mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
           SE = round(total.AREASQM * SE/1000, -1),
           log_SE = sqrt(log(1 + CV^2))) %>%
    select(YEAR, VALUE, VAR, SE, CV, log_SE)
  
  # save by scenario
  write.csv(MCDS_indices_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S1/MCDS_WEE_",
                   q,".csv"), row.names = FALSE)
  
  # extract indices ratio relative to original
  original.MCDS_df <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv")
  
  MCDS_VALUE_ratio <- MCDS_indices_temp$VALUE/original.MCDS_df$VALUE
  MCDS_STDERR_ratio <- MCDS_indices_temp$log_SE/original.MCDS_df$log_SE
  
  save(MCDS_VALUE_ratio, MCDS_STDERR_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S1/MCDS_ratio_WEE_",
                      q,".Rdata"))
  
}; remove(q, tow_MCDS_df_temp, mean_N_stratum_df, stratified_mean_N_df, MCDS_indices_temp, 
          original.MCDS_df, MCDS_VALUE_ratio, MCDS_STDERR_ratio)


  ### ------------------------------------------------------------- ###

 ## ------------------------------------------------------------ ##



 ## S2: inside up, outside stable, total up ----
    # 20%, 40%, 60%, 80%, 100%

  ### 3.2.1  generate tow data ----

diff_total <- data.frame()

for (q in q.pattern) {
  
  # adjust the inside abundance
  inside_temp <- tow_MCDS_df %>%
    filter(OWF == "INSIDE") %>%
    group_by(YEAR) %>%
    mutate(NUMBER_adj = NPERTOW * (1 + q)) %>%
    ungroup() %>%
    add_column(Rate = q)
  
  # extract difference value
  diff_value_temp <- inside_temp %>%
    group_by(YEAR) %>%
    summarise(diff_total = sum(NUMBER_adj - NPERTOW, na.rm = TRUE)) %>%
    add_column(.before = 'diff_total', Rate = q) %>%
    ungroup()
  
  # adjust outside abundance
  outside_temp <- tow_MCDS_df %>%
    filter(OWF == "OUTSIDE") %>%
    mutate(NUMBER_adj = NPERTOW)  %>%
    add_column(Rate = q)
  
  # combine the new catch by tow
  tow_MCDS_df_adjusted <- rbind(outside_temp, inside_temp) %>%
    arrange(ID.temp) %>%
    # below is important: limit the change to the last 10 years
    mutate(NUMBER_adj = case_when(!YEAR %in% ((max(YEAR) - 9): (max(YEAR))) ~ NPERTOW,
                                  TRUE ~ NUMBER_adj))
  
  write.csv(tow_MCDS_df_adjusted,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S2/MCDS_tow_",
                   q,".csv"), row.names = FALSE)
  
  # append the diff_total for record
  diff_total <- rbind(diff_total, diff_value_temp)
  
}; remove(outside_temp, diff_value_temp, inside_temp, tow_MCDS_df_adjusted, q)


write.csv(diff_total,
          "manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S2/RDtrend_tow_abundance_diff.csv"
          , row.names = FALSE)

remove(diff_total)


  ### ------------------------------------------------------------- ###



  ### 3.2.2 full dataset abundance indices ----

for (q in q.pattern) { 
  
  tow_MCDS_df_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S2/MCDS_tow_",
                                          q,".csv"))
  
  # mean abundance by stratum
  mean_N_stratum_df <- tow_MCDS_df_temp %>%
    group_by(YEAR, REGION, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(NUMBER_adj)) %>% # mean within a strata
    mutate(VAR_STRATUM = var(NUMBER_adj)) %>% # variance by stratum
    ungroup() %>% 
    select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
    distinct() %>%
    arrange(YEAR, REGION)
  
  # stratified mean
  stratified_mean_N_df <- mean_N_stratum_df %>%
    group_by(YEAR, REGION) %>%
    summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
              VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
              SE = sqrt(VAR)) %>%   # standard deviance
    mutate(CV = SE / STRATIFIED_MEAN_N) %>%
    ungroup()
  
  # adjust to total abundance based on population area
  MCDS_indices_temp <- tow_MCDS_df_temp %>%
    select(YEAR, STRATUM, AREASQNM) %>%
    distinct()  %>%
    group_by(YEAR) %>%
    summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
    ungroup() %>%
    right_join(stratified_mean_N_df) %>%
    mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
           SE = round(total.AREASQM * SE/1000, -1),
           log_SE = sqrt(log(1 + CV^2))) %>%
    select(YEAR, VALUE, VAR, SE, CV, log_SE)
  
  # save by scenario
  write.csv(MCDS_indices_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S2/MCDS_FULL_",
                   q,".csv"), row.names = FALSE)
  
  # extract indices ratio relative to original
  original.MCDS_df <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv")
  
  MCDS_VALUE_ratio <- MCDS_indices_temp$VALUE/original.MCDS_df$VALUE
  MCDS_STDERR_ratio <- MCDS_indices_temp$log_SE/original.MCDS_df$log_SE
  
  save(MCDS_VALUE_ratio, MCDS_STDERR_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S2/MCDS_ratio_FULL_",
                      q,".Rdata"))
  
}; remove(q, tow_MCDS_df_temp, mean_N_stratum_df, stratified_mean_N_df, MCDS_indices_temp, 
          original.MCDS_df, MCDS_VALUE_ratio, MCDS_STDERR_ratio)




  ### ------------------------------------------------------------- ###


### 3.2.3 full dataset abundance indices ----

for (q in q.pattern) { 
  
  tow_MCDS_df_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/catch by tow/surfclam/S2/MCDS_tow_",
                                      q,".csv")) %>%
    filter(OWF == "OUTSIDE")
  
  # mean abundance by stratum
  mean_N_stratum_df <- tow_MCDS_df_temp %>%
    group_by(YEAR, REGION, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(NUMBER_adj)) %>% # mean within a strata
    mutate(VAR_STRATUM = var(NUMBER_adj)) %>% # variance by stratum
    ungroup() %>% 
    select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
    distinct() %>%
    arrange(YEAR, REGION)
  
  # stratified mean
  stratified_mean_N_df <- mean_N_stratum_df %>%
    group_by(YEAR, REGION) %>%
    summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
              VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
              SE = sqrt(VAR)) %>%   # standard deviance
    mutate(CV = SE / STRATIFIED_MEAN_N) %>%
    ungroup()
  
  # adjust to total abundance based on population area
  MCDS_indices_temp <- tow_MCDS_df_temp %>%
    select(YEAR, STRATUM, AREASQNM) %>%
    distinct()  %>%
    group_by(YEAR) %>%
    summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
    ungroup() %>%
    right_join(stratified_mean_N_df) %>%
    mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
           SE = round(total.AREASQM * SE/1000, -1),
           log_SE = sqrt(log(1 + CV^2))) %>%
    select(YEAR, VALUE, VAR, SE, CV, log_SE)
  
  # save by scenario
  write.csv(MCDS_indices_temp,
            paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S2/MCDS_WEE_",
                   q,".csv"), row.names = FALSE)
  
  # extract indices ratio relative to original
  original.MCDS_df <- read.csv("results/indices for assessment/surfclam/original.MCDS.csv")
  
  MCDS_VALUE_ratio <- MCDS_indices_temp$VALUE/original.MCDS_df$VALUE
  MCDS_STDERR_ratio <- MCDS_indices_temp$log_SE/original.MCDS_df$log_SE
  
  save(MCDS_VALUE_ratio, MCDS_STDERR_ratio, 
       file =  paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/S2/MCDS_ratio_WEE_",
                      q,".Rdata"))
  
}; remove(q, tow_MCDS_df_temp, mean_N_stratum_df, stratified_mean_N_df, MCDS_indices_temp, 
          original.MCDS_df, MCDS_VALUE_ratio, MCDS_STDERR_ratio)


  ### ------------------------------------------------------------- ###
 ## ------------------------------------------------------------ ##
# ------------------------------------------------------ #


