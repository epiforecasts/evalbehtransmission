# Chapter 1, issue #46: forecast fan plots against observed incidence.
# One representative origin per pandemic period, all four models.

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)

## Config ----------------------------------------------------------------------

ch1_fan_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  output_dir    = "outputs/ch1",

  probs         = c(0.05, 0.25, 0.5, 0.75, 0.95),

  # Single model for the rolling overlay: showing all four would stack 280 fans
  rolling_model = "combined"
)

## Load ------------------------------------------------------------------------

periods <- read_csv(ch1_fan_config$periods_path, show_col_types = FALSE) |>
  select(origin = date, period)

forecasts <- read_csv(ch1_fan_config$forecast_path, show_col_types = FALSE) |>
  left_join(periods, by = "origin")

## Choose origins --------------------------------------------------------------
# The middle origin of each period, so the window sits inside its own regime
# rather than straddling a boundary.

chosen_origins <- forecasts |>
  distinct(origin, period) |>
  arrange(origin) |>
  group_by(period) |>
  slice(ceiling(n() / 2)) |>
  ungroup()

cat("\n--- Origins shown ---\n")
print(as.data.frame(chosen_origins |> arrange(origin)), row.names = FALSE)

## Quantiles -------------------------------------------------------------------
# Via scoringutils rather than hand-rolled quantiles, then pivoted for plotting.

fan <- forecasts |>
  semi_join(chosen_origins, by = c("origin", "period")) |>
  as_forecast_sample(forecast_unit = c("model", "horizon", "origin", "period")) |>
  as_forecast_quantile(probs = ch1_fan_config$probs) |>
  as_tibble() |>
  tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                     names_prefix = "q") |>
  mutate(model  = factor(model, levels = names(ch1_models)),
         period = factor(period, levels = chosen_origins$period[order(chosen_origins$origin)]),
         target_date = origin + horizon * 7)

## Fan plot --------------------------------------------------------------------

dir.create(ch1_fan_config$output_dir, recursive = TRUE, showWarnings = FALSE)

fig_1_8 <- ggplot(fan, aes(x = target_date)) +
  geom_ribbon(aes(ymin = q0.05, ymax = q0.95, fill = model), alpha = 0.25) +
  geom_ribbon(aes(ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.45) +
  geom_line(aes(y = q0.5, colour = model), linewidth = 0.6) +
  geom_point(aes(y = observed), colour = "grey20", size = 1) +
  # Both axes free: each panel spans only its own four target weeks, and
  # incidence levels differ by an order of magnitude across periods.
  facet_grid(period ~ model, scales = "free") +
  labs(title = "Weekly forecasts at 1-4 weeks, by model and pandemic period",
       subtitle = "Ribbons are 50% and 90% prediction intervals; points are observed weekly incidence",
       x = "Target week", y = "Weekly infections") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text.y = element_text(angle = 0))

ggsave(file.path(ch1_fan_config$output_dir, "fig_1_8_forecast_fans.png"),
       fig_1_8, width = 12, height = 14, dpi = 300, bg = "white")

## Rolling overlay -------------------------------------------------------------
# Every origin's fan laid along the observed series, for one model. Origins step
# weekly so all target weeks share a boundary and the observed totals form one
# consistent series.

rolling_quantiles <- forecasts |>
  filter(model == ch1_fan_config$rolling_model) |>
  as_forecast_sample(forecast_unit = c("model", "horizon", "origin", "period")) |>
  as_forecast_quantile(probs = ch1_fan_config$probs) |>
  as_tibble() |>
  tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                     names_prefix = "q") |>
  mutate(target_date = origin + horizon * 7)

observed_weekly <- forecasts |>
  distinct(target_date, observed) |>
  arrange(target_date)

fig_1_8_rolling <- ggplot() +
  geom_ribbon(data = rolling_quantiles,
              aes(x = target_date, ymin = q0.05, ymax = q0.95, group = origin),
              fill = "firebrick", alpha = 0.05) +
  geom_ribbon(data = rolling_quantiles,
              aes(x = target_date, ymin = q0.25, ymax = q0.75, group = origin),
              fill = "firebrick", alpha = 0.10) +
  geom_line(data = rolling_quantiles,
            aes(x = target_date, y = q0.5, group = origin),
            colour = "firebrick", linewidth = 0.3, alpha = 0.5) +
  geom_line(data = observed_weekly, aes(x = target_date, y = observed),
            colour = "grey30", linewidth = 0.6) +
  labs(title = paste0("Rolling-origin forecasts, ", ch1_fan_config$rolling_model,
                      " model"),
       subtitle = "Red = forecast fans (50% and 90% PI) | Grey = observed weekly incidence",
       x = "Target week", y = "Weekly infections") +
  theme_minimal()

ggsave(file.path(ch1_fan_config$output_dir, "fig_1_8_rolling_forecasts.png"),
       fig_1_8_rolling, width = 10, height = 4.5, dpi = 300, bg = "white")

## Coverage of the shown origins -----------------------------------------------
# A quick check that the panels match the coverage reported in scoring.

cat("\n--- Observed inside the intervals, for the origins shown ---\n")
print(as.data.frame(fan |>
  group_by(model) |>
  summarise(within_50 = round(mean(observed >= q0.25 & observed <= q0.75), 2),
            within_90 = round(mean(observed >= q0.05 & observed <= q0.95), 2),
            .groups = "drop")), row.names = FALSE)

message("Saved fan plot to ", ch1_fan_config$output_dir)
