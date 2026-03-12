library(tidyverse)

source("manuscript/3. stock assessment with spatial abundance change/R/functions/indices evaluation functions.R")

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# 1. summer flounder ----

Indices.Eva <- data.frame()

q = 0.2; scenario = "S2"; season = "FALL"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    
    SWAN_full <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/", scenario, "/BIG_FULL_", 
                                      q, ".csv"))
    
    SWAN_WEE <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/", scenario, "/BIG_WEE_", 
                                      q, ".csv"))

    
    for (season in c("FALL", "SPRING")) {
      
      SWAN_full_temp <- subset(SWAN_full, SEASON == season)
      SWAN_WEE_temp <- subset(SWAN_WEE, SEASON == season)
      
      
      temp.results <- data.frame(SCENARIO = scenario, q = q, TS = "WHOLE", SEASON = season,
                                 # Overall variability (CV)
                                 CV.FULL = calc_CV(SWAN_full_temp$SWAN), CV.WEE = calc_CV(SWAN_WEE_temp$SWAN), 
                                 # Year-to-year variability (Consecutive Disparity, CD)
                                 CD.FULL = calc_CD(SWAN_full_temp$SWAN), CD.WEE = calc_CD(SWAN_WEE_temp$SWAN),
                                 # significant trend difference (Mann-Kendall test based on relative ratio)
                                 MK.P = mk_trend_difference(SWAN_full_temp$SWAN, SWAN_WEE_temp$SWAN)$p_value)
      
      temp.results.last.10 <-  data.frame(SCENARIO = scenario, q = q, TS = "LAST.10.YEARS", SEASON = season,
                                          # Overall variability (CV)
                                          CV.FULL = calc_CV(tail(SWAN_full_temp$SWAN, 10)), CV.WEE = calc_CV(tail(SWAN_WEE_temp$SWAN, 10)), 
                                          # Year-to-year variability (Consecutive Disparity, CD)
                                          CD.FULL = calc_CD(tail(SWAN_full_temp$SWAN, 10)), CD.WEE = calc_CD(tail(SWAN_WEE_temp$SWAN, 10)),
                                          # significant trend difference (Mann-Kendall test based on relative ratio)
                                          MK.P = mk_trend_difference(tail(SWAN_full_temp$SWAN, 10), tail(SWAN_WEE_temp$SWAN, 10))$p_value)
      
      Indices.Eva <- rbind(Indices.Eva, temp.results, temp.results.last.10) 
    }
  }
}; remove(q, scenario, season, SWAN_full, SWAN_WEE, SWAN_full_temp, SWAN_WEE_temp, temp.results, temp.results.last.10)

write.csv(Indices.Eva,
          paste0("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/summer flounder_ind_eva.csv"), row.names = FALSE)


# ----------------------------------------------------------------------------------------------------------------------------- #


# 2. surfclam ----

Indices.Eva <- data.frame()

q = 0.2; scenario = "S1"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    INDICES_full_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/", scenario, "/MCDS_FULL_", 
                                         q, ".csv"))
    
    INDICES_WEE_temp <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/", scenario, "/MCDS_WEE_", 
                                        q, ".csv"))
    
    
    temp.results <- data.frame(SCENARIO = scenario, q = q, TS = "WHOLE", 
                               # Overall variability (CV)
                               CV.FULL = calc_CV(INDICES_full_temp$VALUE), CV.WEE = calc_CV(INDICES_WEE_temp$VALUE), 
                               # Year-to-year variability (Consecutive Disparity, CD)
                               CD.FULL = calc_CD(INDICES_full_temp$VALUE), CD.WEE = calc_CD(INDICES_WEE_temp$VALUE),
                               # significant trend difference (Mann-Kendall test based on relative ratio)
                               MK.P = mk_trend_difference(INDICES_full_temp$VALUE, INDICES_WEE_temp$VALUE)$p_value)
    
    temp.results.last.10 <-  data.frame(SCENARIO = scenario, q = q, TS = "LAST.10.YEARS",
                                        # Overall variability (CV)
                                        CV.FULL = NA, CV.WEE = NA, 
                                        # Year-to-year variability (Consecutive Disparity, CD)
                                        CD.FULL = NA, CD.WEE = NA,
                                        # significant trend difference (Mann-Kendall test based on relative ratio)
                                        MK.P = NA)
    
    Indices.Eva <- rbind(Indices.Eva, temp.results, temp.results.last.10) 
    
  }
}; remove(q, scenario, INDICES_full_temp, INDICES_WEE_temp, temp.results, temp.results.last.10)

write.csv(Indices.Eva,
          paste0("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/surfclam_ind_eva.csv"), row.names = FALSE)



# ----------------------------------------------------------------------------------------------------------------------------- #



# 3. longfin squid ----

Indices.Eva <- data.frame()

q = 0.2; scenario = "S1"; season = "FALL"


for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    INDICES_full <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/", scenario, "/FULL_", 
                                 q, ".csv"))
    
    INDICES_WEE <- read.csv(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/squid/", scenario, "/WEE_", 
                                q, ".csv"))
    
    for (season in c("FALL", "SPRING")) {
      
      INDICES_full_temp <- subset(INDICES_full, SEASON == season)
      INDICES_WEE_temp <- subset(INDICES_WEE, SEASON == season)
      
      
      temp.results <- data.frame(SCENARIO = scenario, q = q, TS = "WHOLE", SEASON = season,
                                 # Overall variability (CV)
                                 CV.FULL = calc_CV(INDICES_full_temp$STRATIFIED_MEAN_N), CV.WEE = calc_CV(INDICES_WEE_temp$STRATIFIED_MEAN_N), 
                                 # Year-to-year variability (Consecutive Disparity, CD)
                                 CD.FULL = calc_CD(INDICES_full_temp$STRATIFIED_MEAN_N), CD.WEE = calc_CD(INDICES_WEE_temp$STRATIFIED_MEAN_N),
                                 # significant trend difference (Mann-Kendall test based on relative ratio)
                                 MK.P = mk_trend_difference(INDICES_full_temp$STRATIFIED_MEAN_N, INDICES_WEE_temp$STRATIFIED_MEAN_N)$p_value)
      
      temp.results.last.10 <-  data.frame(SCENARIO = scenario, q = q, TS = "LAST.10.YEARS", SEASON = season,
                                          # Overall variability (CV)
                                          CV.FULL = calc_CV(tail(INDICES_full_temp$STRATIFIED_MEAN_N, 10)), CV.WEE = calc_CV(tail(INDICES_WEE_temp$STRATIFIED_MEAN_N, 10)), 
                                          # Year-to-year variability (Consecutive Disparity, CD)
                                          CD.FULL = calc_CD(tail(INDICES_full_temp$STRATIFIED_MEAN_N, 10)), CD.WEE = calc_CD(tail(INDICES_WEE_temp$STRATIFIED_MEAN_N, 10)),
                                          # significant trend difference (Mann-Kendall test based on relative ratio)
                                          MK.P = mk_trend_difference(tail(INDICES_full_temp$STRATIFIED_MEAN_N, 10), tail(INDICES_WEE_temp$STRATIFIED_MEAN_N, 10))$p_value)
      
      Indices.Eva <- rbind(Indices.Eva, temp.results, temp.results.last.10) 
      
    }
  }
}; remove(q, scenario, season, INDICES_full, INDICES_WEE, INDICES_full_temp, INDICES_WEE_temp, temp.results, temp.results.last.10)

write.csv(Indices.Eva,
          paste0("manuscript/3. stock assessment with spatial abundance change/results/indices evaluation/longfin squid_ind_eva.csv"), row.names = FALSE)



# ----------------------------------------------------------------------------------------------------------------------------- #










