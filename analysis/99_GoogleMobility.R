library(tidyverse)
library(janitor) # For cleaning names
setwd("Desktop/PhD_Code/Year_1/evalbehtransmission/R")

# Read in mobility data. WARNING: large file without filtering (1GB)
mobility_data <- read_csv("~/Desktop/Global_Mobility_Report.csv")
mobility_data <- read_csv("~/Desktop/PhD_Code/Year_1/evalbehtransmission/data/Global_Mobility_Report.csv")

# Filter for England national data
# In this case, UK is the country and England the sub_region_1
# For this reason, we also want null sub_region_2
UK_data <- mobility_data |>
  filter(country_region == "United Kingdom",
         is.na(sub_region_1),
         is.na(sub_region_2))

UK <- mobility_data |> filter(country_region == "United Kingdom")
unique(UK$sub_region_1)
unique(UK$sub_region_2)

# Reshape for plotting
UK_long <- UK_data |>
  pivot_longer(
    cols = ends_with("percent_change_from_baseline"),
    names_to = "location_type",
    values_to = "pct_change"
  ) |>
  # Clean up the names (e.g., "retail_and_recreation_..." -> "Retail And Recreation")
  mutate(location_type = str_remove(location_type, "_percent_change_from_baseline") |>
           str_replace_all("_", " ") |>
           str_to_title())


# Plot time series per location
ggplot(UK_long, aes(x = date, y = pct_change, color = location_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + # Baseline
  geom_line(alpha = 0.4) +
  facet_wrap(~location_type) +
  labs(
    title = "Google Mobility, UK",
    subtitle = "Percentage change relative to Jan-Feb 2020 baseline",
    x = "Date",
    y = "% Change from baseline"
  ) +
  theme_classic() +
  theme(legend.position = "none")


# Every country parks time series
all_country_park <- mobility_data |>
  filter(is.na(sub_region_1) & is.na(sub_region_2)) |>
  select(country = country_region, date, parks = parks_percent_change_from_baseline)

ggplot(all_country_park |> filter(str_starts(country, "M")),
       aes(x = date, y = parks)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(alpha = 0.4) +
  facet_wrap(~country) +
  labs(
    title = "Google Mobility, Park setting",
    subtitle = "Percentage change relative to Jan-Feb 2020 baseline",
    x = "Date",
    y = "% Change from baseline"
  ) +
  theme_classic() +
  theme(legend.position = "none")



M_countries <- unique(mobility_data$country_region) |>
  str_subset("^M")

M_countries <- mobility_data |>
  filter(str_starts(country_region, "M")) |>
  distinct(country_region)
