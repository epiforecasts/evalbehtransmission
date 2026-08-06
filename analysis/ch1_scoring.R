# Chapter 1, issue #45: score the rolling-origin forecasts with {scoringutils}.
# CRPS and its over/under/dispersion decomposition come from the sample
# forecasts; interval coverage from the same forecasts as quantiles.

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)

## Config ----------------------------------------------------------------------

ch1_scoring_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  output_dir    = "outputs/ch1",

  # Offset added before logging so zero counts stay defined
  log_offset    = 1,

  # Quantiles needed for the 50% and 90% interval coverage metrics
  probs         = c(0.05, 0.25, 0.5, 0.75, 0.95)
)

print_table <- function(x, digits = 3) {
  print(as.data.frame(mutate(x, across(where(is.numeric), \(v) round(v, digits)))),
        row.names = FALSE)
}

## Load ------------------------------------------------------------------------
# Periods join on origin, so a forecast is labelled by the regime it was made in
# rather than the one it lands in.

periods <- read_csv(ch1_scoring_config$periods_path, show_col_types = FALSE) |>
  select(origin = date, period)

forecasts <- read_csv(ch1_scoring_config$forecast_path, show_col_types = FALSE) |>
  left_join(periods, by = "origin")

forecast_unit <- c("model", "horizon", "origin", "period")

## Score -----------------------------------------------------------------------
# Log transform before scoring so low- and high-incidence periods contribute
# comparably, rather than peaks dominating.

fc_log <- forecasts |>
  as_forecast_sample(forecast_unit = forecast_unit) |>
  transform_forecasts(fun = log_shift, offset = ch1_scoring_config$log_offset,
                      append = FALSE)

sample_scores <- score(fc_log)

coverage_scores <- fc_log |>
  as_forecast_quantile(probs = ch1_scoring_config$probs) |>
  score()

## Table 1.5: model x horizon --------------------------------------------------

table_1_5 <- summarise_scores(sample_scores, by = c("model", "horizon")) |>
  as_tibble() |>
  select(model, horizon, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(
    summarise_scores(coverage_scores, by = c("model", "horizon")) |>
      as_tibble() |>
      select(model, horizon, interval_coverage_50, interval_coverage_90),
    by = c("model", "horizon")
  ) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model, horizon)

cat("\n--- Table 1.5: log-CRPS by model and horizon ---\n")
print_table(table_1_5)

## Table 1.6: model x period ---------------------------------------------------

table_1_6 <- summarise_scores(sample_scores, by = c("model", "period")) |>
  as_tibble() |>
  select(model, period, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(forecasts |> distinct(origin, period) |> count(period, name = "n_origins"),
            by = "period") |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(period, model)

cat("\n--- Table 1.6: log-CRPS by model and period ---\n")
print_table(table_1_6)

## Relative to the incidence-only baseline -------------------------------------

relative_crps <- table_1_5 |>
  select(model, horizon, crps) |>
  group_by(horizon) |>
  mutate(rel_crps = crps / crps[model == "baseline"]) |>
  ungroup()

cat("\n--- CRPS relative to baseline (< 1 is better) ---\n")
print_table(relative_crps |>
  select(model, horizon, rel_crps) |>
  tidyr::pivot_wider(names_from = horizon, values_from = rel_crps,
                     names_prefix = "wk_"))

## Figures ---------------------------------------------------------------------

dir.create(ch1_scoring_config$output_dir, recursive = TRUE, showWarnings = FALSE)

crps_by_origin <- summarise_scores(sample_scores, by = c("model", "origin")) |>
  as_tibble() |>
  mutate(model = factor(model, levels = names(ch1_models)))

fig_1_7 <- ggplot(crps_by_origin, aes(x = origin, y = crps, colour = model)) +
  geom_line(linewidth = 0.5) +
  labs(title = "Log-CRPS over time, averaged across horizons",
       x = "Forecast origin", y = "log-CRPS", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_1_7_crps_over_time.png"),
       fig_1_7, width = 10, height = 4.5, dpi = 300, bg = "white")

fig_1_9 <- table_1_5 |>
  select(model, horizon, overprediction, underprediction, dispersion) |>
  tidyr::pivot_longer(-c(model, horizon), names_to = "component") |>
  ggplot(aes(x = factor(horizon), y = value, fill = component)) +
  geom_col() +
  facet_wrap(~model, nrow = 1) +
  labs(title = "CRPS decomposition by model and horizon",
       x = "Horizon (weeks)", y = "log-CRPS", fill = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_1_9_crps_decomposition.png"),
       fig_1_9, width = 10, height = 4.5, dpi = 300, bg = "white")

## Save tables -----------------------------------------------------------------

write_csv(table_1_5, file.path(ch1_scoring_config$output_dir,
                               "table_1_5_crps_horizon.csv"))
write_csv(table_1_6, file.path(ch1_scoring_config$output_dir,
                               "table_1_6_crps_period.csv"))

message("Saved tables and figures to ", ch1_scoring_config$output_dir)
