library(dplyr)
library(readr)

remotes::install_local("../inc2prev-main")
library(inc2prev)

setwd("/Users/lukeburton/Desktop/PhD_Code/Year_1/evalbehtransmission")

# 1. Process ONS prevalence/incidence data -------------------------------------
ons_estimates <- read_csv("data/ONS/estimates_national.csv")

# 2. Process Google Mobility data ----------------------------------------------
google_mobility_raw <- read_csv("data/Global_Mobility_Report.csv")
google_mobility_UK <- google_mobility_raw |>
  filter(country_region == "United Kingdom",
         is.na(sub_region_1),
         is.na(sub_region_2))
write.csv(google_mobility_UK, "data/google_mobility_UK.csv")

# 3. Process CoMix data --------------------------------------------------------