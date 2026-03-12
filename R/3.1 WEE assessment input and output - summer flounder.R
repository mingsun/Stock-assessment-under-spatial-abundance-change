library(tidyverse)
library(ASAPplots)

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

# 1. modify input survey data  ---------------------------------------------------------------

  # in this section, I applied the ratio between the original_AI and impacted_AI to the original NEFSC indices used in assessment


 ## 1.1 extract the original abundance indices value ----

  # ASAP assessment input value 

data <- ReadASAP3DatFile("manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/input data/ASAP3_MTA2023_FINAL.DAT")
data$survey.names

  # a test
  # WriteASAP3DatFile("manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/input data/ASAP3_MTA2023_test_output.DAT", data, 
  #                   header.text = "Summer Flounder 2023: test output \n# \n# \n#") 
  # test <- ReadASAP3DatFile("manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/input data/ASAP3_MTA2023_test_output.DAT")

  # 2 is BTS ALB spring, 3 is BTS ALB fall, they are 1982-2008
  # 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
  # col.1 is year 1982-2022, 1982-2008 correspond to row 1-27, 2009-2022 correspond to row 28-41, 
  # col.2 is sum value, col.3 is CV
  # column 4-11 correspond to age 1-8,  col. 11 is sample size


BIG_spring <- data$dat$IAA_mats[[25]][28:41,c(2,4:11)]
BIG_fall <- data$dat$IAA_mats[[26]][28:41,c(2,4:11)]


for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    ## 1.2 extract ratios ----
    
    # BIG indices ratio
    load(paste0("manuscript/3. stock assessment with spatial abundance change/results/abundance indices/summer flounder/", scenario, "/BIG_ratio_WEE_", q,".Rdata"))

    # note: be careful with the round and -999 values, need to be very consistent with the original indices format 
    # manually edit the -999 year in the original data as -999 after adjustment
    BIG_SPRING_ratio[12] = 1 # 2020 not changed
    BIG_FALL_ratio[9] = 1 # 2017 not changed
    BIG_FALL_ratio <- append(BIG_FALL_ratio, 1, after = 11) # 2020 not changed 
    
    
    ## 1.3 apply the ratio to generate scenario-specific output ----
    BIG_spring_temp <- round(BIG_spring * BIG_SPRING_ratio)
    BIG_fall_temp <-   round(BIG_fall *   BIG_FALL_ratio)
    

    ## 1.4 generate ASAP input by scenario ----

    # BIG spring (2009-2022)
    data$dat$IAA_mats[[25]][28:41,c(2,4:11)] <- BIG_spring_temp
    
    # BIG fall (2009-2022)
    data$dat$IAA_mats[[26]][28:41,c(2,4:11)] <- BIG_fall_temp
    
    WriteASAP3DatFile(paste0("manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/input data/ASAP3_MTA2023_WEE_",
                             scenario, "_", q ,".DAT"), data, 
                      header.text = paste0("Summer Flounder 2023: BTS indices ", scenario, "_" , q, " redistribution indices dataset (Sun 2024) \n \n \n "))
    
  }
}

remove(list = ls())

#  ------------------------------------------------------------------------------------------------ #


# 2. run model ----

# in the ASAP program

assessment.path <- "manuscript/3. stock assessment with spatial abundance change/stock assessment/summer flounder/"

q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

q = 0.2; scenario = "S1"
for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) { 
    
    outdir <- paste0(assessment.path, "WEE_", scenario, "_", q)
    
    if (!dir.exists(outdir)) { 
      dir.create(outdir)
      
      file.copy(from = paste0(assessment.path, "input data/ASAP3_MTA2023_WEE_", scenario, "_", q ,".DAT"),
                to   = paste0(outdir, "/ASAP3_MTA2023_WEE_", scenario, "_", q ,".DAT"))
    }
    
  }
}


# in the ASAP program




#  ------------------------------------------------------------------------------------------------ #


# 3. extract quantity ----

for (q in q.pattern) {
  
  for (scenario in c("S1", "S2")) {
    
    ## 3.1 save output ----
    
    # get the variable needed
    wd <- paste0(assessment.path, "WEE_", scenario , "_", q, "/")
    asap.name <- paste0("ASAP3_MTA2023_WEE_", scenario ,"_", q ,"_000")
    asap <- dget(paste0(wd, asap.name, ".rdat"))
    
    gn <- GrabNames(wd, asap.name, asap) # names of the catch fleet and survey
    fleet.names <- gn$fleet.names
    index.names <- gn$index.names
    
    # reads in auxiliary files with values from the *.std, *.par, and *.cor
    a1 <- GrabAuxFiles(wd, asap.name, asap = asap, fleet.names, index.names)
    
    # generate the estimate output, SSB and F
    output.path <- "manuscript/3. stock assessment with spatial abundance change/results/stock assessment/summer flounder/"
    SummarizeASAP(asap, a1, od = output.path)
    
    output <- read.csv(paste0(output.path, "ASAP_summary_ASAP3_MTA2023_WEE_", scenario , "_", q, "_000",".csv")) %>%
      add_column(Dataset = "WEE", Scenario = scenario, q = q) %>%
      write.csv(paste0(output.path, "ASAP_summary_ASAP3_MTA2023_WEE_", scenario , "_", q, ".csv"), row.names = FALSE)
    
    # delete the original output
    file.remove(paste0(output.path, "ASAP_summary_ASAP3_MTA2023_WEE_", scenario , "_", q, "_000",".csv"))
    
  }; remove(wd, asap.name, asap, gn, fleet.names, index.names, a1)
}

#  ------------------------------------------------------------------------------------------------ #







