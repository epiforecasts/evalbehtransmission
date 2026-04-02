# Reads raw Google Mobility data from data-raw/, filters to UK national level,
# and saves to data-processed/.
#
# Google Mobility categories (% change from Jan 3 – Feb 6 2020 baseline):
#   retail_and_recreation, grocery_and_pharmacy, parks,
#   transit_stations, workplaces, residential

library(dplyr)
library(readr)

# Read in mobility data. WARNING: large file without filtering (1GB)
google_mobility_raw <- read_csv("data-raw/Global_Mobility_Report.csv")

google_mobility_UK <- google_mobility_raw |>
  filter(country_region == "United Kingdom",
         is.na(sub_region_1),
         is.na(sub_region_2))

write_csv(google_mobility_UK, "data-processed/google_mobility_UK.csv")
