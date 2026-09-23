# Reads raw Google Mobility data, filters to UK national level, and saves to
# data-processed/.
#
# Google Mobility categories (% change from Jan 3 – Feb 6 2020 baseline):
#   retail_and_recreation, grocery_and_pharmacy, parks,
#   transit_stations, workplaces, residential

library(dplyr)
library(readr)

## Config ----------------------------------------------------------------------

mobility_config <- list(

  # FALSE reads the global report from data-raw/, which is 1GB and must be downloaded first
  use_remote   = TRUE,

  # Google publishes one file per country per year, together far smaller than the global report
  remote_url   = "https://www.gstatic.com/covid19/mobility/%d_GB_Region_Mobility_Report.csv",
  remote_years = 2020:2022,

  local_path   = "data-raw/Mobility/Global_Mobility_Report.csv",
  output_path  = "data-processed/google_mobility_UK.csv"
)

## Load ------------------------------------------------------------------------

# Both routes carry the same columns and the same national rows
read_mobility <- function(config = mobility_config) {
  if (config$use_remote) {
    lapply(config$remote_years, \(year) {
      read_csv(sprintf(config$remote_url, year), show_col_types = FALSE)
    }) |> bind_rows()
  } else {
    read_csv(config$local_path, show_col_types = FALSE) |>
      filter(country_region == "United Kingdom")
  }
}

# Rows with no sub-region are the national series
google_mobility_UK <- read_mobility() |>
  filter(is.na(sub_region_1), is.na(sub_region_2)) |>
  arrange(date)

write_csv(google_mobility_UK, mobility_config$output_path)
message(sprintf("Saved %d rows to %s", nrow(google_mobility_UK),
                mobility_config$output_path))
