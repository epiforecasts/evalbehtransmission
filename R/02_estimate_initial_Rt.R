library(dplyr)
library(readr)
library(EpiEstim)

# Load incidence ---
ons_estimates <- read_csv("data/ONS/estimates_national.csv") |>
  filter(variable == "England")
ons_incidence <- ons_estimates |>
  filter(name == "infections") |>
  select(date, median = q50) |>
  arrange(date)

# EpiEstim requires data.frame with 'I' column with no NAs or non-integers
incidence_df <- ons_incidence |>
  mutate(I = round(median)) |>
  filter(!is.na(I), I >= 0)

# Serial interval for SARS-CoV-2 based on Nishiura et al. 2020
si_mean <- 4.7
si_sd   <- 2.9

# Estimate Rt using EpiEstim (parametric SI, 7-day sliding window)
rt_config <- make_config(
  mean_si   = si_mean,
  std_si    = si_sd,
  t_start   = 2:(nrow(incidence_df) - 6),   # 7-day window
  t_end     = 8:nrow(incidence_df)
)

rt_fit <- estimate_R(
  incid  = incidence_df$I,
  method = "parametric_si",
  config = rt_config
)

# Extract results and attach dates for alignment
rt_estimates <- rt_fit$R |>
  as_tibble() |>
  mutate(
    date = incidence_df$date[t_end], # Date = last day of window
    .before = 1
  ) |>
  select(
    date,
    Rt_mean   = `Mean(R)`,
    Rt_median = `Median(R)`,
    Rt_q025   = `Quantile.0.025(R)`,
    Rt_q975   = `Quantile.0.975(R)`
  )

# Save into data file
write_csv(rt_estimates, "data/ONS/rt_estimates_England.csv")

