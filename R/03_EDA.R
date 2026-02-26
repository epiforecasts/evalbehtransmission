library(ggplot2)
library(tidyr)
library(stringr)


## Compare Google Mobility against ONS for each location -----------------------

ons_estimates <- read_csv("data/ONS/estimates_national.csv") |>
  filter(variable == "England")
ons_incidence <- ons_estimates |>
  filter(name == "infections") |>
  select(date, median = q50)

mobility_UK <- read_csv("data/google_mobility_UK.csv")
mobility_long <- mobility_UK |>
  pivot_longer(
    cols = ends_with("percent_change_from_baseline"),
    names_to = "location_type",
    values_to = "pct_change"
  ) |>
  mutate(location_type = str_remove(location_type, "_percent_change_from_baseline") |>
           str_replace_all("_", " ") |>
           str_to_title())

mobility_range <- range(mobility_long$pct_change, na.rm = TRUE)
ons_range <- range(ons_incidence$median, na.rm = TRUE)
scale_factor <- diff(mobility_range) / diff(ons_range)
offset <- mobility_range[1] - ons_range[1] * scale_factor

ons_scaled <- ons_incidence |>
  mutate(scaled_inc = median * scale_factor + offset)

# Plot scaled Google Mobility against ONS median prevalence by location
ggplot() +
  # Mobility line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = mobility_long,
            aes(x = date, y = pct_change, color = location_type),
            alpha = 0.5) +
  # ONS incidence overlaid (same axis, scaled)
  geom_line(data = ons_scaled |>
              crossing(location_type = unique(mobility_long$location_type)),
            aes(x = date, y = scaled_inc),
            color = "black", linewidth = 0.8, alpha = 0.7) +
  # Secondary axis label
  scale_y_continuous(
    name = "% Change from baseline",
    sec.axis = sec_axis(
      trans = ~ (. - offset) / scale_factor,
      name = "ONS daily incidence"
    )
  ) +
  facet_wrap(~location_type) +
  labs(
    title = "Google Mobility & ONS incidence, UK",
    subtitle = "Black line: ONS incidence (median daily estimate). Coloured = Mobility % change from baseline",
    x = "Date"
  ) +
  theme_classic() +
  theme(legend.position = "none")











## Same but with Rt estimated from EpiEstim and ONS-inferred median incidence
rt_estimates <- read_csv("data/ONS/rt_estimates_England.csv") |>
  mutate(date = as.Date(date)) |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-03-01")))

mobility_UK <- read_csv("data/google_mobility_UK.csv") |>
  mutate(date = as.Date(date)) |>
  filter(between(date, as.Date("2020-03-01"), as.Date("2021-03-01")))
mobility_long <- mobility_UK |>
  pivot_longer(
    cols      = ends_with("percent_change_from_baseline"),
    names_to  = "location_type",
    values_to = "pct_change"
  ) |>
  mutate(location_type = str_remove(location_type, "_percent_change_from_baseline") |>
           str_replace_all("_", " ") |>
           str_to_title())

# Scale Rt to mobility axis
mobility_range <- range(mobility_long$pct_change, na.rm = TRUE)
rt_range       <- range(rt_estimates$Rt_median,   na.rm = TRUE)
scale_factor   <- diff(mobility_range) / diff(rt_range)
offset         <- mobility_range[1] - rt_range[1] * scale_factor

rt_scaled <- rt_estimates |>
  mutate(
    scaled_rt    = Rt_median * scale_factor + offset,
    scaled_q025  = Rt_q025   * scale_factor + offset,
    scaled_q975  = Rt_q975   * scale_factor + offset
  )

# Rt = 1 line scaled to mobility axis
rt1_scaled <- 1 * scale_factor + offset

ggplot() +
  geom_hline(yintercept = 0,         linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = rt1_scaled, linetype = "dotted", color = "black", alpha = 0.5) +
  # Mobility
  geom_line(data = mobility_long,
            aes(x = date, y = pct_change, color = location_type),
            alpha = 0.5) +
  # Rt ribbon + line
  geom_ribbon(data = rt_scaled |>
                crossing(location_type = unique(mobility_long$location_type)),
              aes(x = date, ymin = scaled_q025, ymax = scaled_q975),
              fill = "black", alpha = 0.15) +
  geom_line(data = rt_scaled |>
              crossing(location_type = unique(mobility_long$location_type)),
            aes(x = date, y = scaled_rt),
            color = "black", linewidth = 0.8, alpha = 0.7) +
  scale_y_continuous(
    name     = "% Change from baseline",
    sec.axis = sec_axis(
      trans = ~ (. - offset) / scale_factor,
      name  = "Estimated Rt"
    )
  ) +
  facet_wrap(~location_type) +
  labs(
    title    = "Google Mobility & Estimated Rt, England",
    subtitle = "Black line: EpiEstim Rt (median + 95% CrI). Coloured = Mobility % change from baseline",
    x        = "Date"
  ) +
  theme_classic() +
  theme(legend.position = "none")
