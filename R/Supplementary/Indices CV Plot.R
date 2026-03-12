library(tidyverse)
source("manuscript/3. stock assessment with spatial abundance change/R/functions/stratified mean.R")


q.pattern <- c(0.2, 0.4, 0.6, 0.8, 1)

indices_path <- "manuscript/3. stock assessment with spatial abundance change/results/abundance indices/"



# 1. summer flounder ----

SF_INDICES <- data.frame()

q = 0.2; scenario = "S1"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    BIG_full_temp <- read.csv(paste0(indices_path, "summer flounder/", scenario, "/BIG_FULL_", q, ".csv")) %>%
      add_column(DATASET = "FULL")
    
    BIG_wee_temp <- read.csv(paste0(indices_path, "summer flounder/", scenario, "/BIG_WEE_", q, ".csv"))%>%
      add_column(DATASET = "WEE")
    
    BIG_temp <- rbind(BIG_full_temp, BIG_wee_temp) %>%
      select(DATASET, YEAR, SEASON, VALUE = CV) %>%
      add_column(SPECIES = "SUMMER FLOUNDER", SCENARIO = scenario, q = q, .before = 1)
    
    SF_INDICES <- rbind(SF_INDICES, BIG_temp)
    
  }
}; remove(q, scenario, BIG_full_temp, BIG_wee_temp, BIG_temp)


SF_plot <- SF_INDICES %>%
  filter(YEAR != 2017, YEAR >= 2013)

# fall plot
SF_FALL <- ggplot(subset(SF_plot, SEASON == "FALL"), aes(x = YEAR, y = VALUE, color = DATASET)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  scale_x_continuous(breaks = c(2013:2022)) +
  ylab("CV") +
  facet_grid(q~SCENARIO) +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S9 SF_FALL_CV.png",  width = 8, height = 5, units = 'in', res = 800)
print(SF_FALL)
dev.off()

# spring plot
SF_SPRING <- ggplot(subset(SF_plot, SEASON == "SPRING"), aes(x = YEAR, y = VALUE, color = DATASET)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  scale_x_continuous(breaks = c(2013:2022)) +
  ylab("STD.INDEX") +
  facet_grid(q~SCENARIO) +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S10 SF_SPRING_CV.png",  width = 8, height = 5, units = 'in', res = 800)
print(SF_SPRING)
dev.off()



# ----------------------------------------------------------------------------------------------------------------------------- #



# 2. squid ----

LS_INDICES <- data.frame()

q = 0.2; scenario = "S1"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    BTS_full_temp <- read.csv(paste0(indices_path, "squid/", scenario, "/FULL_", q, ".csv")) %>%
      add_column(DATASET = "FULL")
    
    BTS_wee_temp <- read.csv(paste0(indices_path, "squid/", scenario, "/WEE_", q, ".csv"))%>%
      add_column(DATASET = "WEE")
    
    BTS_temp <- rbind(BTS_full_temp, BTS_wee_temp) %>%
      select(DATASET, YEAR, SEASON, VALUE = CV) %>%
      add_column(SPECIES = "LONGFIN SQUID", SCENARIO = scenario, q = q, .before = 1)
    
    LS_INDICES <- rbind(LS_INDICES, BTS_temp)
    
  }
}; remove(q, scenario, BTS_full_temp, BTS_wee_temp, BTS_temp)


LS_plot <- LS_INDICES %>%
  filter(YEAR != 2017, YEAR >= 2013) 

# fall plot
LS_FALL <- ggplot(subset(LS_plot, SEASON == "FALL"), aes(x = YEAR, y = VALUE, color = DATASET)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  scale_x_continuous(breaks = c(2013:2022)) +
  ylab("CV") +
  facet_grid(q~SCENARIO) +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S11 LS_FALL_CV.png",  width = 8, height = 5, units = 'in', res = 800)
print(LS_FALL)
dev.off()

# spring plot
LS_SPRING <- ggplot(subset(LS_plot, SEASON == "SPRING"), aes(x = YEAR, y = VALUE, color = DATASET)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  scale_x_continuous(breaks = c(2013:2022)) +
  ylab("CV") +
  facet_grid(q~SCENARIO) +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S12 LS_SPRING_CV.png",  width = 8, height = 5, units = 'in', res = 800)
print(LS_SPRING)
dev.off()

# ----------------------------------------------------------------------------------------------------------------------------- #



# 3. surfclam ----

SC_INDICES <- data.frame()

q = 0.2; scenario = "S1"

for (q in q.pattern) { 
  
  for (scenario in c("S1", "S2")) { 
    
    MCDS_full_temp <- read.csv(paste0(indices_path, "surfclam/", scenario, "/MCDS_FULL_", q, ".csv")) %>%
      add_column(DATASET = "FULL")
    
    MCDS_wee_temp <- read.csv(paste0(indices_path, "surfclam/", scenario, "/MCDS_WEE_", q, ".csv"))%>%
      add_column(DATASET = "WEE")
    
    MCDS_temp <- rbind(MCDS_full_temp, MCDS_wee_temp) %>%
      add_column(SEASON = "ANNUAL") %>%
      select(DATASET, YEAR, SEASON, VALUE = CV) %>%
      add_column(SPECIES = "ATLANTIC SURFCLAM", SCENARIO = scenario, q = q, .before = 1)
    
    SC_INDICES <- rbind(SC_INDICES, MCDS_temp)
    
  }
}; remove(q, scenario, MCDS_full_temp, MCDS_wee_temp, MCDS_temp)


SC_plot <- SC_INDICES %>%
  filter(YEAR != 2017, YEAR >= 2013)

# fall plot
SC_PLOT <- ggplot(SC_plot, aes(x = YEAR, y = VALUE, color = DATASET)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("FULL", "WEE"), name = "DATASET") +
  scale_x_continuous(breaks = c(2013:2022)) +
  ylab("CV") +
  facet_grid(q~SCENARIO) +
  theme_bw() +
  theme(legend.position = "none") 

png("manuscript/3. stock assessment with spatial abundance change/plot/SA/S13 SC_CV.png",  width = 8, height = 5, units = 'in', res = 800)
print(SC_PLOT)
dev.off()



# ----------------------------------------------------------------------------------------------------------------------------- #



