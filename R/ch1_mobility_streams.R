# Google Mobility categories and subset used in modelling
# This is sourced by ch1_data.R, ch1_covariates.R and ch1_descriptive.R

# Names are short forms used downstream for readability. Values are raw Google columns
mobility_categories <- c(
  retail_recreation = "retail_and_recreation_percent_change_from_baseline",
  grocery_pharmacy  = "grocery_and_pharmacy_percent_change_from_baseline",
  parks             = "parks_percent_change_from_baseline",
  transit           = "transit_stations_percent_change_from_baseline",
  workplaces        = "workplaces_percent_change_from_baseline",
  residential       = "residential_percent_change_from_baseline"
)

# Parks and residential excluded, following Nouvellet et al. (2021)
# Tomori et al. 2021 exclude parks as "this is expected to vary considerably during seasons"
mobility_retained <- c("retail_recreation", "grocery_pharmacy",
                       "transit", "workplaces")

# Raw Google column names with the suffix stripped, for plotting against the unprocessed file
mobility_short_names <- function(streams = mobility_retained) {
  sub("_percent_change_from_baseline", "", mobility_categories[streams])
}