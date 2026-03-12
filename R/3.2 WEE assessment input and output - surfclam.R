library(r4ss)
library(tidyverse)

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# 1. modify input survey data  ---------------------------------------------------------------

  # in this section, I applied the ratio between the original_AI and impacted_AI to the original NEFSC indices used in assessment


 ## 1.1 extract the original abundance indices value ----

  # SS assessment input value 

data <- SS_readdat("manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/input data/data.ss")
data$Nsurveys # 6 surveys, RD means research dredge, MCD means modified commercial dredge
data$fleetnames # need to play with RDtrendS, RDscaleS, and MCDs
data$CPUEinfo # units are all numbers according to data.ss, 



 ## 1.2 loop to update indices by scenario ----

q = 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    
    ### a. RDtrendS: index 3, numbers per m2, 1982-2011 ----
    
    # not needed because it does not cover the most recent 10 years
    
    ### ------------------------------------------------------------- ###
    
    
    ## b. RDscaleS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----
    
    # not needed because it does not cover the most recent 10 years
    
    ### ------------------------------------------------------------- ###
    
    
    
    ## c. MCDS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----
    
    load(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/surfclam/",
                scenario, "/MCDS_ratio_WEE_", q,".Rdata"))
    
    MCDS <- data$CPUE %>%
      filter(index == 5) %>%
      mutate(obs = obs * MCDS_VALUE_ratio,
             se_log = se_log * MCDS_STDERR_ratio)
    
    data$CPUE[data$CPUE$index == 5, ] <- MCDS
    
    remove(MCDS_STDERR_ratio, MCDS_VALUE_ratio, MCDS)
    
    
    ### ------------------------------------------------------------- ###
    
    
    ## 1.3 save the modified input into a data.ss ----
    
    # create empty folders
    
    new_folder <- paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/WEE_",
                         scenario, "_", q)
    
    dir.create(new_folder)
    
    
    
    # generate the input data file
    
    SS_writedat(datlist = data, outfile = paste0(new_folder, "/data.ss"), overwrite = TRUE)
    
    
    # copy paste other data
    
    file.copy(from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/input data/control.ss",
              to = new_folder)
    
    file.copy(from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/input data/forecast.ss",
              to = new_folder)
    
    file.copy(from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/input data/starter.ss",
              to = new_folder)
    
    file.copy(from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/input data/ss3.exe",
              to = new_folder)
   
    remove (new_folder) 
  }
}


#  ------------------------------------------------------------------------------------------------ #




# 2. run  model ----

q = 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    
    run_dir <- paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/WEE_",
                          scenario, "_", q)
    
    r4ss::run(dir = run_dir) # skip any folders that already contain a "Report.sso" file
    
    print(paste0(scenario, "_", q, "_finished"))

  }
}


#  ------------------------------------------------------------------------------------------------ #



# 3. extract outputs  ----


q = 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    # assessment model path
    run_dir <- paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/surfclam/WEE_",
                      scenario, "_", q)
    
    # create a list of quantities for the outputs
    assessment_results <- SS_output(run_dir, verbose = FALSE)
    
    output_path <- paste0("manuscript/3. stock assessment with spatial abundance change/results/",
                          "stock assessment/surfclam/")
    
    # see the names of items available
    # sort(names(assessment_results_WEE))

    
    # time-series data and reference points
    
    
    
    
    ## a. SSB ----
    
    SSB_df <- assessment_results$derived_quants %>%
      filter(grepl("^SSB_\\d{4}$", Label))
    
    # SSB-threshold is 1/4 of the maximum historical SSB 
    SSB.MSY <- 1/2 * max(SSB_df$Value)
    SSB.threshold <- 1/2 * SSB.MSY
    
    SSB_df$Ratio <- SSB_df$Value/SSB.threshold
    SSB_df$Ratio.SD <- SSB_df$StdDev/SSB.threshold
    SSB_df$Dataset <- "WEE"
    SSB_df$Scenario <- scenario
    SSB_df$q <- q
    
    write.csv(SSB_df, paste0(output_path, "SSB_WEE_", scenario, "_", q, ".csv"))
    
    ## ------------------------------------------------------ ##
    
    
    
    
    # b. F ----
    
    # F time series, F is combined with number-weighted mean
    F_df <- assessment_results$timeseries %>%
      select(Yr, Area, 
             F1 = `F:_1`, F2 = `F:_2`,
             N1 = `SmryNum_SX:1_GP:1`, 
             N2 = `SmryNum_SX:1_GP:2`) %>%
      group_by(Yr) %>%
      summarize(
        F1 = max(F1, na.rm = TRUE),
        F2 = max(F2, na.rm = TRUE),
        N1 = max(N1, na.rm = TRUE),
        N2 = max(N2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(F_weighted = (F1 * N1 + F2 * N2) / (N1 + N2))
    
    # grab the F std from the derived quants
    Fstd_df <- assessment_results$derived_quants %>%
      filter(grepl("^F_\\d{4}$", Label)) %>%
      mutate(Yr = as.numeric(substr(Label, 3, 6)))  %>%
      select(Yr, Value, StdDev)
    
    
    # combine F with Fstd
    F_df <- F_df %>%
      left_join(Fstd_df) %>%
      mutate(F_weighted_std = F_weighted/Value*StdDev)  %>%
      select(-c(Value, StdDev))
    
    
    # F threshold is defined based on an algorithm provided by Dan
    
    source("R/functions/2.3.1 functions for surfclam stock assessment.r")
    
    GetFref(rlst = assessment_results, rhoF = rho)
    
    F_df$Ratio <- F_df$F_weighted/0.1526799
    F_df$Ratio.SD <- F_df$F_weighted_std/0.1526799
    F_df$Dataset <- "WEE"
    F_df$Scenario <- scenario
    F_df$q <- q
    
    write.csv(F_df, paste0(output_path, "F_WEE_", scenario, "_", q, ".csv"))
    
    
    
    ## ------------------------------------------------------ ##
    
    
    # c. Catch limits ----
    
    ABC_df <- assessment_results$derived_quants %>%
      filter(grepl("ForeCatch_\\d{4}$", Label)) %>%
      add_column(Dataset = "WEE",
                 Scenario = scenario,
                 q = q)
    
    write.csv(ABC_df, paste0(output_path, "ABC_WEE_", scenario, "_", q, ".csv"), row.names = FALSE)
    
    ## ------------------------------------------------------ ##
    
 

  }
}










