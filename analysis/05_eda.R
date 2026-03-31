library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)

# ---------------------------------------------------------------------------
# 1. Load data ---------------------------------------------------------------
# ---------------------------------------------------------------------------

ons_estimates <- read_csv("data-raw/ONS/estimates_national.csv") |>
  filter(variable == "England")

ons_incidence <- ons_estimates |>
  filter(name == "infections") |>
  select(date, median = q50)

mobility_UK <- read_csv("data-processed/google_mobility_UK.csv")
mobility_long <- mobility_UK |>
  pivot_longer(
    cols = ends_with("percent_change_from_baseline"),
    names_to = "location_type",
    values_to = "pct_change"
  ) |>
  mutate(location_type = str_remove(location_type, "_percent_change_from_baseline") |>
           str_replace_all("_", " ") |>
           str_to_title())

# Requires 01_process_ons.R to have been run first; alternatively use EpiEstim Rt from 06_estimate_rt.R: read_csv("data-processed/rt_estimates_England.csv")
rt_estimates <- read_csv("data-processed/rt_national.csv") |>
  rename(Rt_median = median, Rt_q05 = q05, Rt_q95 = q95) |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-12-31")))

# ---------------------------------------------------------------------------
# 2. Summary statistics ------------------------------------------------------
# ---------------------------------------------------------------------------

summary(ons_incidence)
summary(mobility_UK |> select(date, ends_with("percent_change_from_baseline")))
summary(rt_estimates)

# ---------------------------------------------------------------------------
# 3. Time series plots -------------------------------------------------------
# ---------------------------------------------------------------------------

## Compare Google Mobility against ONS incidence for each location type -------

mobility_long_inc <- mobility_long |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-12-31")))

ons_incidence_inc <- ons_incidence |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-12-31")))

mobility_range_inc <- range(mobility_long_inc$pct_change, na.rm = TRUE)
ons_range_inc      <- range(ons_incidence_inc$median,     na.rm = TRUE)
scale_factor_inc   <- diff(mobility_range_inc) / diff(ons_range_inc)
offset_inc         <- mobility_range_inc[1] - ons_range_inc[1] * scale_factor_inc

ons_scaled <- ons_incidence_inc |>
  mutate(scaled_inc = median * scale_factor_inc + offset_inc)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = mobility_long_inc,
            aes(x = date, y = pct_change, color = location_type),
            alpha = 0.5) +
  geom_line(data = ons_scaled |>
              crossing(location_type = unique(mobility_long_inc$location_type)),
            aes(x = date, y = scaled_inc),
            color = "black", linewidth = 0.8, alpha = 0.7) +
  scale_y_continuous(
    name = "% Change from baseline",
    sec.axis = sec_axis(
      transform = ~ (. - offset_inc) / scale_factor_inc,
      name = "ONS daily incidence"
    )
  ) +
  facet_wrap(~location_type) +
  labs(
    title    = "Google Mobility & ONS incidence, UK",
    subtitle = "Black line: ONS incidence (median daily estimate). Coloured = Mobility % change from baseline. 2020-2021.",
    x = "Date"
  ) +
  theme_classic() +
  theme(legend.position = "none")


## Compare Google Mobility against Rt ----------------------------------------

mobility_long_rt <- mobility_long |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-12-31")))

mobility_range_rt <- range(mobility_long_rt$pct_change, na.rm = TRUE)
rt_range          <- range(rt_estimates$Rt_median, na.rm = TRUE)
scale_factor_rt   <- diff(mobility_range_rt) / diff(rt_range)
offset_rt         <- mobility_range_rt[1] - rt_range[1] * scale_factor_rt

rt_scaled <- rt_estimates |>
  mutate(
    scaled_rt  = Rt_median * scale_factor_rt + offset_rt,
    scaled_q05  = Rt_q05     * scale_factor_rt + offset_rt,
    scaled_q95 = Rt_q95    * scale_factor_rt + offset_rt
  )

rt1_scaled <- 1 * scale_factor_rt + offset_rt

ggplot() +
  geom_hline(yintercept = 0,         linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = rt1_scaled, linetype = "dotted", color = "black", alpha = 0.5) +
  geom_line(data = mobility_long_rt,
            aes(x = date, y = pct_change, color = location_type),
            alpha = 0.5) +
  geom_ribbon(data = rt_scaled |>
                crossing(location_type = unique(mobility_long_rt$location_type)),
              aes(x = date, ymin = scaled_q05, ymax = scaled_q95),
              fill = "black", alpha = 0.15) +
  geom_line(data = rt_scaled |>
              crossing(location_type = unique(mobility_long_rt$location_type)),
            aes(x = date, y = scaled_rt),
            color = "black", linewidth = 0.8, alpha = 0.7) +
  scale_y_continuous(
    name     = "% Change from baseline",
    sec.axis = sec_axis(
      transform = ~ (. - offset_rt) / scale_factor_rt,
      name      = "Estimated Rt"
    )
  ) +
  facet_wrap(~location_type) +
  labs(
    title    = "Google Mobility & Estimated Rt, England",
    subtitle = "Black line: EpiEstim Rt",
    x        = "Date"
  ) +
  theme_classic() +
  theme(legend.position = "none")


## Summary plots --------------------------------------------------------------

# Plot 1: ONS Daily Incidence
p1 <- ggplot(ons_incidence |> filter(between(date, as.Date("2020-01-01"), as.Date("2021-12-31"))),
             aes(x = date, y = median)) +
  geom_line() +
  labs(title    = "ONS Daily Incidence, England",
       subtitle = "inc2prev median estimate",
       x = "Date", y = "Daily infections") +
  theme_classic()

# Plot 2: Rt with 95% CrI ribbon
p2 <- ggplot(rt_estimates |> filter(between(date, as.Date("2020-01-01"), as.Date("2021-12-31"))),
             aes(x = date)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = Rt_q05, ymax = Rt_q95), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = Rt_median)) +
  labs(title    = "Estimated Rt, England",
       subtitle = "inc2prev",
       x = "Date", y = "Rt") +
  theme_classic()

print(p1)
print(p2)
