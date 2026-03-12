# stratified mean and variance function

calc_stratified_mean <- function(obj) { 
  
  mean_N_stratum_df <- obj %>%
    group_by(YEAR, SEASON, STRATUM) %>%
    mutate(MEAN_N_STRATUM = mean(NUMBER_adj)) %>% # mean within a strata
    mutate(VAR_STRATUM = var(NUMBER_adj), na.rm = TRUE) %>% # variance by stratum, same as last line
    ungroup() 
  
  stratified_mean_N_df <- mean_N_stratum_df %>%
    select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
    distinct() %>% # downsize the data frame to a minimal without repetitive rows
    group_by(YEAR, SEASON) %>%
    mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
    summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
              STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
              STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
    mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
           lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
           lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
           CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
    ungroup()
  
  return(stratified_mean_N_df)
  
  }



