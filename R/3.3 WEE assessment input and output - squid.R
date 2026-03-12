library(tidyverse)

# everything is done in spreadsheet
# at E:\Stony Brook job\NYSERDA Offshore wind\manuscript\3. stock assessment with spatial abundance change\stock assessment\squid

# first use the "2023_longfin_index-based method.xlsx" to calculate biomass and indices based on stratified mean
# then use the "index-based summary" to extract the key results for further analyses


# create empty folders

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    new_folder <- paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/squid/WEE_",
                         scenario, "_", q)
    
    dir.create(new_folder, showWarnings = FALSE) 
    
    
    file.copy(
      from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/squid/input data/2023_longfin_index-based method_base.xlsx",
      to   = file.path(new_folder, paste0("Index_script_WEE_", scenario, "_", q, ".xlsx"))
    )
    
    
    file.copy(
      from = "manuscript/3. stock assessment with spatial abundance change/stock assessment/squid/input data/index-based summary_base.xlsx",
      to   = file.path(new_folder, paste0("Summary_WEE_", scenario, "_", q, ".csv"))
    )
    
  }
}


# after finishing the edits in the scripts, paste the summary file to the results folder

q= 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    
    assessment_folder <- paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/squid/WEE_",
                                scenario, "_", q, "/")
    
    results_folder <- paste0("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/squid/")
    
    
    file.copy(
      from = paste0(assessment_folder, "Summary_WEE_", scenario, "_", q, ".csv"),
      to   = results_folder)
    
  }
}



# change column name and add dataset column

q= 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    
    results_folder <- paste0("manuscript/3. stock assessment with spatial abundance change/results/stock assessment/squid/")
    
    
    read.csv(paste0(results_folder, "Summary_WEE_", scenario, "_", q, ".csv")) %>%
      add_column(Dataset = "WEE", .before = 1) %>%
      rename(Scenario = SCENARIO) %>%
      write.csv(paste0(results_folder, "Summary_WEE_", scenario, "_", q, ".csv"), row.names = FALSE)
    
    
  }
}

