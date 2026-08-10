# Chapter 1, issue #45: score the rolling-origin forecasts with {scoringutils}
# CRPS and its over/under/dispersion decomposition come from the sample forecasts
# Interval coverage comes from the same forecasts, converted to quantiles

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)

## Config ----------------------------------------------------------------------

ch1_scoring_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  output_dir    = "outputs/ch1",

  # Ensures log transforms are well-defined
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
# This follows epiforecasts/CovidAgeGroupForecast, which also cuts on forecast_date

periods <- read_csv(ch1_scoring_config$periods_path, show_col_types = FALSE) |>
  select(origin = date, period)

# Place in chronological order to avoid sorting alphabetically
period_levels <- periods |> arrange(origin) |> distinct(period) |> pull(period)
 
forecasts <- read_csv(ch1_scoring_config$forecast_path, show_col_types = FALSE) |>
  left_join(periods, by = "origin")

forecast_unit <- c("model", "horizon", "origin", "period")

## Score -----------------------------------------------------------------------
# Log transform before scoring so low- and high-incidence periods contribute comparably
# Otherwise CRPS scales with incidence and peaks dominate the average

# Appending keeps the natural scale alongside the log scale, in a scale column
forecasts_both_scales <- forecasts |>
  as_forecast_sample(forecast_unit = forecast_unit) |>
  transform_forecasts(fun = log_shift, offset = ch1_scoring_config$log_offset)

all_scores <- score(forecasts_both_scales)
log_scores <- all_scores[all_scores$scale == "log", ]

# Interval coverage is invariant to monotone transforms, so the log scale gives the same values
coverage_scores <- forecasts_both_scales[forecasts_both_scales$scale == "log", ] |>
  as_forecast_quantile(probs = ch1_scoring_config$probs) |>
  score()

## Table 1.5: model x horizon --------------------------------------------------

table_1_5 <- summarise_scores(log_scores, by = c("model", "horizon")) |>
  as_tibble() |>
  select(model, horizon, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(
    summarise_scores(coverage_scores, by = c("model", "horizon")) |>
      as_tibble() |>
      select(model, horizon, interval_coverage_50, interval_coverage_90),
    by = c("model", "horizon")
  ) |>
  left_join(
    # Standard deviation of CRPS across origins
    summarise_scores(log_scores, by = c("model", "horizon"), fun = sd) |>
      as_tibble() |>
      select(model, horizon, crps_sd = crps),
    by = c("model", "horizon")
  ) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model, horizon)

cat("\n--- Table 1.5: log-CRPS by model and horizon, averaged over all origins ---\n")
print_table(table_1_5)

## Table 1.6: model x period ---------------------------------------------------

table_1_6 <- summarise_scores(log_scores, by = c("model", "period")) |>
  as_tibble() |>
  select(model, period, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(forecasts |> distinct(origin, period) |> count(period, name = "n_origins"),
            by = "period") |>
  mutate(model = factor(model, levels = names(ch1_models)),
         period = factor(period, levels = period_levels)) |> # Arrange chronologically, not alphabetically
  arrange(period, model)

cat("\n--- Table 1.6: log-CRPS by model and period, averaged over that period's origins ---\n")
print_table(table_1_6)

## Relative to the incidence-only baseline -------------------------------------
# scoringutils pairs each model against the baseline on shared forecasts
# adj_pval comes from a permutation test, adjusted for multiple comparisons
# It tests whether the difference in score between two models is statistically significant

relative_crps <- get_pairwise_comparisons(log_scores, compare = "model",
                                          by = "horizon", metric = "crps",
                                          baseline = "baseline") |>
  as_tibble() |>
  filter(compare_against == "baseline", model != "baseline") |>
  select(model, horizon, mean_scores_ratio, adj_pval) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model, horizon)

cat("\n--- CRPS relative to baseline (< 1 is better) ---\n")
print_table(relative_crps)

## Figures ---------------------------------------------------------------------

dir.create(ch1_scoring_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Shortest and longest horizons only, as CRPS grows roughly fourfold between them
# Averaging across all four would hide the 1-week horizon
crps_by_origin <- summarise_scores(log_scores, by = c("model", "origin", "horizon")) |>
  as_tibble() |>
  filter(horizon %in% c(1, 4)) |>
  mutate(model = factor(model, levels = names(ch1_models)),
         horizon = factor(horizon, levels = c(1, 4),
                          labels = c("1 week ahead", "4 weeks ahead")))

fig_1_7 <- ggplot(crps_by_origin, aes(x = origin, y = crps, colour = model)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~horizon, ncol = 1, scales = "free_y") +
  labs(title = "Log-CRPS over time, by horizon",
       x = "Forecast origin", y = "log-CRPS", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_1_7_crps_over_time.png"),
       fig_1_7, width = 10, height = 6.5, dpi = 300, bg = "white")

## Save tables -----------------------------------------------------------------

write_csv(table_1_5, file.path(ch1_scoring_config$output_dir,
                               "table_1_5_crps_horizon.csv"))
write_csv(table_1_6, file.path(ch1_scoring_config$output_dir,
                               "table_1_6_crps_period.csv"))

message("Saved tables and figures to ", ch1_scoring_config$output_dir)
