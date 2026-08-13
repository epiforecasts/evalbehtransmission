# Chapter 1, issue #45: score the rolling-origin forecasts with {scoringutils}
# CRPS and its over/under/dispersion decomposition come from the sample forecasts
# Interval coverage comes from the same forecasts, converted to quantiles
# Tables run from overall, to horizon, to period, then decomposition

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)
library(patchwork)

## Config ----------------------------------------------------------------------

ch1_scoring_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  scores_path   = "data-processed/ch1_scores.csv",
  output_dir    = "outputs/ch1",

  # Ensures log transforms are well-defined
  log_offset    = 1,

  # Quantiles needed for the 50% and 90% interval coverage metrics
  probs         = c(0.05, 0.25, 0.5, 0.75, 0.95),

  # The decomposition-over-time figure shows one horizon at a time, in weeks ahead
  decomposition_horizon_weeks = 2
)

print_table <- function(x, digits = 3) {
  print(as.data.frame(mutate(x, across(where(is.numeric), \(v) round(v, digits)))),
        row.names = FALSE)
}

# Ratio to the baseline within the current grouping, as in epiforecasts/CovidAgeGroupForecast
# Only meaningful for positive, lower-is-better metrics, so CRPS but not the signed bias
crps_ratio_to_baseline <- function(score, model) {
  stopifnot("baseline" %in% model) # Without it the ratio would silently return Inf
  score / score[model == "baseline"]
}

## Load ------------------------------------------------------------------------
# Periods join on origin, so a forecast is labelled by the regime it was made in
# This follows epiforecasts/CovidAgeGroupForecast, which also cuts on forecast_date

periods <- read_csv(ch1_scoring_config$periods_path, show_col_types = FALSE)

# Place in chronological order to avoid sorting alphabetically
period_levels <- periods |> arrange(date) |> distinct(period) |> pull(period)

forecasts <- read_csv(ch1_scoring_config$forecast_path, show_col_types = FALSE) |>
  left_join(periods |> select(origin = date, period), by = "origin")

forecast_unit <- c("model", "horizon", "origin", "period")

# Observed weekly infection totals, drawn above the score figures for context
observed_weekly <- forecasts |>
  distinct(target_date, observed) |>
  arrange(target_date)

# Check each target week carries one observed total, not one per origin reaching it
stopifnot(nrow(observed_weekly) == n_distinct(forecasts$target_date))

# Date bounds shared by every figure, from the first origin to the last week scored
# Padded so the bars centred on the end origins are not clipped
score_date_limits <- c(min(forecasts$origin) - 5, max(forecasts$target_date) + 5)

## Score -----------------------------------------------------------------------
# Log transform before scoring so low- and high-incidence periods contribute comparably
# Otherwise CRPS scales with incidence and peaks dominate the average

# transform_forecasts() appends, leaving the natural scale untouched and scoring each separately
forecasts_both_scales <- forecasts |>
  as_forecast_sample(forecast_unit = forecast_unit) |>
  transform_forecasts(fun = log_shift, offset = ch1_scoring_config$log_offset)

all_scores     <- score(forecasts_both_scales)
log_scores     <- all_scores[all_scores$scale == "log", ]
natural_scores <- all_scores[all_scores$scale == "natural", ]

# Coverage needs interval bounds, so the same forecasts are converted to quantiles
# It is 0 or 1 per forecast, averaging to the proportion covered, and only these two columns are used
# Invariant to monotone transforms, so the log scale gives the same values
coverage_scores <- forecasts_both_scales[forecasts_both_scales$scale == "log", ] |>
  as_forecast_quantile(probs = ch1_scoring_config$probs) |>
  score()

## Per-origin scores -----------------------------------------------------------
# One row per model, horizon, origin and scale, so figures read scores rather than recompute them

ch1_scores <- all_scores |>
  as_tibble() |>
  select(model, horizon, origin, period, scale,
         crps, overprediction, underprediction, dispersion, bias,
         median_abs_error = ae_median) |>
  left_join(coverage_scores |>
              as_tibble() |>
              select(model, horizon, origin,
                     interval_coverage_50, interval_coverage_90),
            by = c("model", "horizon", "origin")) |>
  mutate(model  = factor(model, levels = names(ch1_models)),
         period = factor(period, levels = period_levels)) |>
  arrange(model, origin, horizon, scale)

write_csv(ch1_scores, ch1_scoring_config$scores_path)

# The time series figures below all work from the log scale, so filter once here
log_scores_per_origin <- filter(ch1_scores, scale == "log")

## Study periods ---------------------------------------------------------------
# Periods are counted by the origin they contain, matching how forecasts are labelled above

period_summary_table <- periods |>
  group_by(period) |>
  summarise(start = min(date), end = max(date), days = n(),
            mean_stringency = round(mean(stringency), 1), .groups = "drop") |>
  left_join(forecasts |> distinct(origin, period) |> count(period, name = "n_origins"),
            by = "period") |>
  mutate(n_origins = coalesce(n_origins, 0L),
         period    = factor(period, levels = period_levels)) |>
  arrange(period)

cat("\n--- Study periods, with the forecast origins falling in each ---\n")
print_table(period_summary_table, 1)

## Overall, across all horizons ------------------------------------------------
# One row per model, averaging every origin and horizon
# rel_crps is the ratio to baseline, so below 1 means the covariates helped

crps_overall_table <- summarise_scores(log_scores, by = "model") |>
  as_tibble() |>
  select(model, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(summarise_scores(coverage_scores, by = "model") |>
              as_tibble() |>
              select(model, interval_coverage_50, interval_coverage_90),
            by = "model") |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model) |>
  mutate(rel_crps = crps_ratio_to_baseline(crps, model), .after = crps)

cat("\n--- log-CRPS by model, averaged over all origins and horizons ---\n")
print_table(crps_overall_table)

## By horizon ------------------------------------------------------------------
# Averaged over origins, so rel_crps compares models within each horizon

crps_horizon_table <- summarise_scores(log_scores, by = c("model", "horizon")) |>
  as_tibble() |>
  select(model, horizon, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(
    summarise_scores(coverage_scores, by = c("model", "horizon")) |>
      as_tibble() |>
      select(model, horizon, interval_coverage_50, interval_coverage_90),
    by = c("model", "horizon")
  ) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model, horizon) |>
  group_by(horizon) |>
  mutate(rel_crps = crps_ratio_to_baseline(crps, model), .after = crps) |>
  ungroup()

cat("\n--- log-CRPS by model and horizon ---\n")
print_table(crps_horizon_table)

## By period -------------------------------------------------------------------
# Averaged over that period's origins and all horizons, so cells with few origins are thin

crps_period_table <- summarise_scores(log_scores, by = c("model", "period")) |>
  as_tibble() |>
  select(model, period, crps, overprediction, underprediction, dispersion, bias) |>
  left_join(period_summary_table |> select(period, n_origins), by = "period") |>
  mutate(model  = factor(model, levels = names(ch1_models)),
         period = factor(period, levels = period_levels)) |>
  arrange(period, model) |>
  group_by(period) |>
  mutate(rel_crps = crps_ratio_to_baseline(crps, model), .after = crps) |>
  ungroup()

cat("\n--- log-CRPS by model and period ---\n")
print_table(crps_period_table)

## Natural scale, for the supplement -------------------------------------------
# Kept separate as CRPS on this scale is dominated by the highest-incidence weeks

crps_natural_table <- summarise_scores(natural_scores, by = c("model", "horizon")) |>
  as_tibble() |>
  select(model, horizon, crps, overprediction, underprediction, dispersion, bias) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model, horizon) |>
  group_by(horizon) |>
  mutate(rel_crps = crps_ratio_to_baseline(crps, model), .after = crps) |>
  ungroup() |>
  # Scores here are counts of infections, so round them rather than print decimals
  mutate(across(c(crps, overprediction, underprediction, dispersion), round))

cat("\n--- Natural-scale CRPS by model and horizon (supplementary) ---\n")
print_table(crps_natural_table)

## Shared plot elements --------------------------------------------------------
# Period shading and the incidence panel are reused by the score figures below

dir.create(ch1_scoring_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Bands are clamped to the plotted window, as ggplot drops any rect reaching outside it
period_bands <- periods |>
  group_by(period) |>
  summarise(start = min(date), end = max(date), .groups = "drop") |>
  arrange(start) |>
  mutate(period = factor(period, levels = period)) |>
  filter(end >= score_date_limits[1], start <= score_date_limits[2]) |>
  mutate(start = pmax(start, score_date_limits[1]),
         end   = pmin(end,   score_date_limits[2]))

# Monthly ticks throughout, so every panel in the chapter reads on the same axis
score_date_scale <- function() {
  scale_x_date(limits = score_date_limits, expand = c(0, 0),
               date_breaks = "1 month", date_labels = "%d-%b")
}

# Shades each period and fixes the date axis, so stacked panels line up
shade_periods <- function() {
  list(
    geom_rect(data = period_bands, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = period),
              alpha = 0.10),
    score_date_scale()
  )
}

# Scores are placed at their origin, so the incidence panel runs on the same date axis
# It carries the period legend, leaving the score panels to label their own components
incidence_panel <- ggplot(observed_weekly, aes(x = target_date, y = observed)) +
  shade_periods() +
  geom_col(fill = "lightblue3", width = 6) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  labs(y = "Weekly infections", x = NULL, fill = NULL) +
  # order = 2 puts the period key below the model key when both are collected
  guides(fill = guide_legend(nrow = 1, order = 2, override.aes = list(alpha = 0.5))) +
  theme_minimal()

component_colours <- c(overprediction  = "#C0392B",
                       underprediction = "#27AE60",
                       dispersion      = "#2E86C1")

# Outlines each facet, so stacked rows read as separate panels
panel_outline <- theme(panel.border  = element_rect(colour = "grey80", fill = NA),
                       panel.spacing = unit(0.5, "lines"))

## Fig 1.7: score over time ----------------------------------------------------
# Averaged over horizons 1-4, so the longer horizons dominate the level

crps_by_origin <- log_scores_per_origin |>
  group_by(model, origin) |>
  summarise(crps = mean(crps), .groups = "drop")

crps_over_time_panel <- ggplot(crps_by_origin, aes(x = origin, y = crps, colour = model)) +
  shade_periods() +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  # Model key above the period key, and no second period key from this panel
  guides(colour = guide_legend(order = 1), fill = "none") +
  labs(x = "Forecast origin", y = "log-CRPS", colour = NULL, fill = NULL) +
  theme_minimal()

fig_1_7 <- incidence_panel / crps_over_time_panel +
  plot_layout(heights = c(1, 2), guides = "collect") +
  plot_annotation(
    title = "Observed incidence and forecast score over time",
    subtitle = "log-CRPS averaged over horizons 1-4, plotted at the forecast origin"
  ) &
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_1_7_crps_over_time.png"),
       fig_1_7, width = 10, height = 7, dpi = 300, bg = "white")

## Several metrics over time ---------------------------------------------------
# The same origins as above, read across three metrics rather than log-CRPS alone
# Each is averaged over horizons 1-4, so one point per model and origin

metric_labels <- c(crps             = "log-CRPS",
                   median_abs_error = "Median absolute error",
                   bias             = "Bias")

metrics_by_origin <- log_scores_per_origin |>
  group_by(model, origin) |>
  summarise(across(all_of(names(metric_labels)), mean), .groups = "drop") |>
  tidyr::pivot_longer(all_of(names(metric_labels)),
                      names_to = "metric", values_to = "score") |>
  mutate(metric = factor(metric, levels = names(metric_labels),
                         labels = metric_labels))

metrics_panel <- ggplot(metrics_by_origin, aes(x = origin, y = score, colour = model)) +
  # Bias is signed, so only that panel needs a zero reference
  geom_hline(data = tibble(metric = factor(metric_labels[["bias"]], levels = metric_labels)),
             aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  panel_outline +
  score_date_scale() +
  labs(x = "Forecast origin", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.placement = "outside")

fig_metrics_over_time <- (incidence_panel + guides(fill = "none")) / metrics_panel + # Period key dropped, as this figure has no shading below
  plot_layout(heights = c(1, 3)) +
  plot_annotation(
    title = "Observed incidence and forecast performance over time",
    subtitle = "Each metric averaged over horizons 1-4, plotted at the forecast origin"
  )

ggsave(file.path(ch1_scoring_config$output_dir, "fig_metrics_over_time.png"),
       fig_metrics_over_time, width = 10, height = 9, dpi = 300, bg = "white")

## CRPS decomposition, overall -------------------------------------------------
# Components sum to CRPS, so stacking shows what each model's score is made of

decomposition_by_horizon <- crps_horizon_table |>
  select(model, horizon, overprediction, underprediction, dispersion) |>
  tidyr::pivot_longer(-c(model, horizon), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = names(component_colours)))

fig_crps_decomposition <- ggplot(decomposition_by_horizon,
                                 aes(x = model, y = score, fill = component)) +
  geom_col(width = 0.7) +
  facet_wrap(~horizon, nrow = 1,
             labeller = labeller(horizon = \(h) paste(h, "weeks ahead"))) +
  scale_fill_manual(values = component_colours,
                    labels = c("Overprediction", "Underprediction", "Dispersion")) +
  labs(title = "log-CRPS decomposition by model and horizon",
       subtitle = "Components sum to log-CRPS, averaged over all forecast origins",
       x = NULL, y = "Contribution to log-CRPS", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_crps_decomposition.png"),
       fig_crps_decomposition, width = 10, height = 4.5, dpi = 300, bg = "white")

## CRPS decomposition over time ------------------------------------------------
# Incidence sits above the bars so poor scores line up with peaks and dips

decomposition_over_time <- log_scores_per_origin |>
  filter(horizon == ch1_scoring_config$decomposition_horizon_weeks) |>
  select(model, origin, overprediction, underprediction, dispersion) |>
  tidyr::pivot_longer(-c(model, origin), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = names(component_colours)))

crps_decomposition_panel <- ggplot(decomposition_over_time,
                               aes(x = origin, y = score, fill = component)) +
  geom_col(width = 6) +
  facet_grid(model ~ .) +
  scale_fill_manual(values = component_colours) +
  score_date_scale() +
  labs(x = "Forecast origin", y = "log-CRPS contribution", fill = NULL) +
  theme_minimal() +
  panel_outline

fig_crps_decomposition_time <- incidence_panel / crps_decomposition_panel +
  plot_layout(heights = c(1, 4), guides = "collect") +
  plot_annotation(
    title = "Observed incidence and log-CRPS decomposition over time",
    subtitle = sprintf("Horizon %d weeks ahead, plotted at the forecast origin",
                       ch1_scoring_config$decomposition_horizon_weeks)
  ) &
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_crps_decomposition_time.png"),
       fig_crps_decomposition_time, width = 10, height = 9, dpi = 300, bg = "white")

## Coverage, for the supplement ------------------------------------------------
# Reference lines mark nominal coverage

coverage_by_horizon <- crps_horizon_table |>
  select(model, horizon, `50` = interval_coverage_50, `90` = interval_coverage_90) |>
  tidyr::pivot_longer(-c(model, horizon), names_to = "interval", values_to = "coverage") |>
  mutate(nominal  = as.numeric(interval) / 100,
         interval = paste0(interval, "% prediction interval"))

fig_coverage <- ggplot(coverage_by_horizon,
                       aes(x = horizon, y = coverage, colour = model)) +
  geom_hline(aes(yintercept = nominal), linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~interval) +
  scale_x_continuous(breaks = 1:4) +
  ylim(0, 1) +
  labs(title = "Empirical interval coverage against nominal, by horizon",
       subtitle = "Dashed lines mark nominal coverage",
       x = "Horizon (weeks ahead)", y = "Proportion of forecasts covered", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(ch1_scoring_config$output_dir, "fig_coverage_by_horizon.png"),
       fig_coverage, width = 9, height = 4.5, dpi = 300, bg = "white")

## Save tables -----------------------------------------------------------------

write_csv(period_summary_table, file.path(ch1_scoring_config$output_dir, "table_period_summary.csv"))
write_csv(crps_overall_table,   file.path(ch1_scoring_config$output_dir, "table_crps_overall.csv"))
write_csv(crps_horizon_table,   file.path(ch1_scoring_config$output_dir, "table_crps_horizon.csv"))
write_csv(crps_period_table,    file.path(ch1_scoring_config$output_dir, "table_crps_period.csv"))
write_csv(crps_natural_table,   file.path(ch1_scoring_config$output_dir, "table_crps_natural.csv"))

message("Saved tables and figures to ", ch1_scoring_config$output_dir)
