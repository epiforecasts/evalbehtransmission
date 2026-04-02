# Reads raw ONS/inc2prev estimates from data-raw/ONS/, filters to England national level,
# and saves processed incidence and Rt outputs to data-processed/.
#
# inc2prev output columns:
#   name     - variable name ("infections", "R", "r", "est_prev", ...)
#   q05–q95  - posterior quantiles (absolute counts for infections; unitless for R)
#   variable - geography (e.g. "England")
#   date     - daily date

library(dplyr)
library(readr)

ons_estimates <- read_csv("data-raw/ONS/estimates_national.csv") |>
  filter(variable == "England")

# Daily infection counts (median + 90% credible interval)
incidence_national <- ons_estimates |>
  filter(name == "infections") |>
  select(date, median = q50, q05 = q5, q95) |>
  arrange(date)

# Rt estimates from inc2prev
rt_national <- ons_estimates |>
  filter(name == "R") |>
  select(date, median = q50, q05 = q5, q95) |>
  arrange(date)

write_csv(incidence_national, "data-processed/incidence_national.csv")
write_csv(rt_national,        "data-processed/rt_national.csv")
