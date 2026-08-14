# Chapter 1, issue #46: forecast fan plots against observed incidence
# One representative origin per pandemic period, all four models

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)

## Config ----------------------------------------------------------------------

ch1_fan_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  output_dir    = "outputs/ch1",

  probs         = c(0.05, 0.25, 0.5, 0.75, 0.95),

  scores_path = "data-processed/ch1_scores.csv",

  # Illustrative origin, chosen after scoring as the one where models disagree most
  zoom_origin = as.Date("2020-11-18"),

  # Weeks of observed incidence drawn before the origin, for context
  zoom_history_weeks = 6,

  # Upper limit on the rolling figure, clipping only the most extreme forecast tails
  rolling_y_max = 2.5e6
)

## Load ------------------------------------------------------------------------

periods <- read_csv(ch1_fan_config$periods_path, show_col_types = FALSE) |>
  select(origin = date, period)

forecasts <- read_csv(ch1_fan_config$forecast_path, show_col_types = FALSE) |>
  left_join(periods, by = "origin")

# An origin outside the OxCGRT date range would become an unlabelled panel
stopifnot(!any(is.na(forecasts$period)))

## Choose origins --------------------------------------------------------------
# The middle origin of each period
# The training window is trailing, so for shorter periods most of it precedes the regime

chosen_origins <- forecasts |>
  distinct(origin, period) |>
  arrange(origin) |>
  group_by(period) |>
  slice(ceiling(n() / 2)) |>
  ungroup()

cat("\n--- Origins shown ---\n")
print(as.data.frame(chosen_origins |> arrange(origin)), row.names = FALSE)

# Highlight any periods with no origins, whether missed by the weekly grid or outside its range
cat("Periods with no origins:",
    paste(setdiff(unique(periods$period), unique(forecasts$period)), collapse = ", "), "\n")

## Quantiles at the chosen origins ---------------------------------------------
# Quantiles from the samples via as_forecast_quantile(), pivoted for plotting
# semi_join keeps only the one origin per period selected above

chosen_origin_quantiles <- forecasts |>
  semi_join(chosen_origins, by = c("origin", "period")) |>
  # Identifies one forecast; unlisted columns are dropped, so period rides along as a label
  as_forecast_sample(forecast_unit = c("model", "horizon", "origin", "period")) |>
  as_forecast_quantile(probs = ch1_fan_config$probs) |> # Collapses each forecast's samples to the quantiles in probs
  as_tibble() |>
  tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                     names_prefix = "q") |>
  mutate(model = factor(model, levels = names(ch1_models)),
         # Row label gives the period and the origin date forecast from
         period_origin_label = paste0(period, "\n", format(origin)),
         period_origin_label = factor(period_origin_label,
                                      levels = unique(period_origin_label[order(origin)])))

## Fan plot --------------------------------------------------------------------

dir.create(ch1_fan_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Borders and tinted strips so each row and column reads as a separate cell
facet_cells <- theme(panel.border     = element_rect(colour = "grey80", fill = NA),
                     panel.spacing    = unit(0.5, "lines"),
                     strip.background = element_rect(fill = "grey94", colour = "grey80"))

# Horizon rather than date on the x-axis, so each panel covers the same 1-4 weeks
# Incidence levels differ across periods, hence free_y
# facet_grid frees the y-axis by row, so the four models share a scale within a period
fig_1_8 <- ggplot(chosen_origin_quantiles, aes(x = horizon)) +
  geom_ribbon(aes(ymin = q0.05, ymax = q0.95, fill = model), alpha = 0.25) +
  geom_ribbon(aes(ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.45) +
  geom_line(aes(y = q0.5, colour = model), linewidth = 0.6) +
  geom_point(aes(y = observed), colour = "grey20", size = 1.2) +
  facet_grid(period_origin_label ~ model, scales = "free_y") +
  scale_x_continuous(breaks = 1:4) +
  labs(title = "Weekly forecasts at 1-4 weeks, by model and pandemic period",
       subtitle = paste("One forecast origin per period, taken from its midpoint;",
                        "ribbons are 50% and 90% prediction intervals,",
                        "points = observed weekly incidence"),
       x = "Horizon (weeks ahead)", y = "Weekly infections") +
  theme_minimal() +
  facet_cells +
  theme(legend.position = "none",
        strip.text.y = element_text(angle = 0))

ggsave(file.path(ch1_fan_config$output_dir, "fig_1_8_forecast_fans.png"),
       fig_1_8, width = 12, height = 14, dpi = 300, bg = "white")

## Rolling forecasts, every origin ---------------------------------------------
# One panel per model, each holding the fan from every origin against the observed series
# Origins step weekly, so all target weeks fall on the same weekday (Wednesday)

rolling_quantiles <- forecasts |>
  as_forecast_sample(forecast_unit = c("model", "horizon", "origin", "period")) |>
  as_forecast_quantile(probs = ch1_fan_config$probs) |>
  as_tibble() |>
  tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                     names_prefix = "q") |>
  mutate(model       = factor(model, levels = names(ch1_models)),
         target_date = origin + horizon * 7)

observed_weekly <- forecasts |>
  distinct(target_date, observed) |>
  arrange(target_date)

# Check each target week carries one observed total, not one per origin reaching it
stopifnot(nrow(observed_weekly) == n_distinct(forecasts$target_date))

# Facet labels name the model, so the colour carries no extra information and needs no legend
rolling_forecast_plot <- ggplot(rolling_quantiles, aes(x = target_date, group = origin)) +
  geom_ribbon(aes(ymin = q0.05, ymax = q0.95, fill = model), alpha = 0.06) +
  geom_ribbon(aes(ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.12) +
  geom_line(aes(y = q0.5, colour = model), linewidth = 0.3, alpha = 0.6) +
  geom_line(data = observed_weekly, inherit.aes = FALSE,
            aes(x = target_date, y = observed), colour = "black", linewidth = 0.5) +
  geom_point(data = observed_weekly, inherit.aes = FALSE,
             aes(x = target_date, y = observed), colour = "black", size = 0.9) +
  facet_wrap(~model, ncol = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d-%b") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  # Shared limit across panels, clipping the few upper tails that run far above any observed week
  coord_cartesian(ylim = c(0, ch1_fan_config$rolling_y_max)) +
  labs(title = "Rolling-origin forecasts against observed incidence, by model",
       subtitle = "One fan per forecast origin, 50% and 90% prediction intervals; black = observed weekly incidence",
       x = "Target week", y = "Weekly infections") +
  theme_minimal() +
  facet_cells +
  theme(legend.position = "none")

ggsave(file.path(ch1_fan_config$output_dir, "fig_rolling_forecasts_by_model.png"),
       rolling_forecast_plot, width = 10, height = 10, dpi = 300, bg = "white")

## Single origin, zoomed --------------------------------------------------------
# A single origin, zoomed in so the models can be compared side by side
# Each panel shows that model's mean log-CRPS here

zoom_fan <- rolling_quantiles |> filter(origin == ch1_fan_config$zoom_origin)

# Observed incidence either side of the origin, so the forecast is read against its run-up
zoom_observed <- observed_weekly |>
  filter(target_date >= ch1_fan_config$zoom_origin - ch1_fan_config$zoom_history_weeks * 7,
         target_date <= max(zoom_fan$target_date))

# Mean over horizons 1-4, matching how this origin is scored elsewhere
zoom_scores <- read_csv(ch1_fan_config$scores_path, show_col_types = FALSE) |>
  filter(scale == "log", origin == ch1_fan_config$zoom_origin) |>
  group_by(model) |>
  summarise(crps = mean(crps), .groups = "drop") |>
  mutate(model = factor(model, levels = names(ch1_models)),
         label = sprintf("log-CRPS %.2f", crps))

stopifnot(nrow(zoom_fan) > 0, nrow(zoom_scores) == length(ch1_models))


fig_forecast_zoom <- ggplot(zoom_fan, aes(x = target_date)) +
  geom_vline(xintercept = ch1_fan_config$zoom_origin, linetype = "dashed",
             colour = "grey50") +
  geom_ribbon(aes(ymin = q0.05, ymax = q0.95, fill = model), alpha = 0.25) +
  geom_ribbon(aes(ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.45) +
  geom_line(aes(y = q0.5, colour = model), linewidth = 0.7) +
  geom_line(data = zoom_observed, inherit.aes = FALSE,
            aes(x = target_date, y = observed), colour = "grey20", linewidth = 0.6) +
  geom_point(data = zoom_observed, inherit.aes = FALSE,
             aes(x = target_date, y = observed), colour = "grey20", size = 1.2) +
  geom_text(data = zoom_scores, inherit.aes = FALSE,
            aes(x = min(zoom_observed$target_date), y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 3.2, colour = "grey20") +
  facet_wrap(~model, ncol = 2) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  labs(title = paste("Forecasts from", format(ch1_fan_config$zoom_origin),
                     "by model"),
       subtitle = paste("Dashed line marks the forecast origin; ribbons are 50% and 90%",
                        "prediction intervals, black = observed weekly incidence"),
       x = "Target week", y = "Weekly infections") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(ch1_fan_config$output_dir, "fig_forecast_zoom.png"),
       fig_forecast_zoom, width = 10, height = 6, dpi = 300, bg = "white")

## Coverage of the shown origins -----------------------------------------------
# All four models, pooling horizons 1-4, over the one origin per period plotted above
# n is therefore small and this is not the headline figure
# ch1_scoring.R reports coverage over every origin, split by horizon

cat("\n--- Observed inside the intervals, origins shown, horizons 1-4 pooled ---\n")
print(as.data.frame(chosen_origin_quantiles |>
  group_by(model) |>
  summarise(n         = n(),
            within_50 = round(mean(observed >= q0.25 & observed <= q0.75), 2),
            within_90 = round(mean(observed >= q0.05 & observed <= q0.95), 2),
            .groups = "drop")), row.names = FALSE)

message("Saved fan and rolling plots to ", ch1_fan_config$output_dir)
