# Appendix A.1.5: model fit and forecast at a single origin, for all four models
# Everything is in weekly totals, so the fit and the forecast share one axis and one unit,
# and the fit is shown on the same quantity the CRPS scores
# Forecast quantiles are read from data-processed/ch1_forecasts.csv, so the panel shows
# the forecasts the chapter scores. Only the training-window fit is recomputed, which is
# deterministic and is checked against ch1_window_coefficients.csv below

source("analysis/ch1_gam.R")

library(scoringutils)
library(ggplot2)

## Config ----------------------------------------------------------------------

a1_5_config <- list(
  forecast_path = "data-processed/ch1_forecasts.csv",
  scores_path   = "data-processed/ch1_scores.csv",
  coef_path     = "data-processed/ch1_window_coefficients.csv",
  output_dir    = "stashed/appendix/outputs",

  # Matches ch1_rolling.R, so the refitted window is the one that produced the forecasts
  window_weeks = 8,

  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),

  # A second origin from Summer relaxation, as an alternative to the headline choice
  alternative_period = "Summer relaxation"
)

dir.create(a1_5_config$output_dir, recursive = TRUE, showWarnings = FALSE)

dat       <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
forecasts <- read_csv(a1_5_config$forecast_path, show_col_types = FALSE)
scores    <- read_csv(a1_5_config$scores_path, show_col_types = FALSE)
stored_coefficients <- read_csv(a1_5_config$coef_path, show_col_types = FALSE)

## Choose the origins to show --------------------------------------------------
# The chapter discusses where the models diverge, so origins are picked on the spread of
# log-CRPS across the four models rather than by hand

origin_spread <- scores |>
  filter(scale == "log") |>
  group_by(origin, period, model) |>
  summarise(crps = mean(crps), .groups = "drop") |> # Mean over horizons 1-4
  group_by(origin, period) |>
  summarise(spread = max(crps) - min(crps), .groups = "drop") |>
  arrange(desc(spread))

cat("\n--- Origins by spread in mean log-CRPS across the four models ---\n")
print(as.data.frame(head(origin_spread, 5) |>
        mutate(spread = round(spread, 3))), row.names = FALSE)

chosen_origins <- c(
  headline    = origin_spread$origin[1],
  alternative = origin_spread |>
    filter(period == a1_5_config$alternative_period) |>
    slice(1) |>
    pull(origin)
)

cat("\nShowing", format(chosen_origins[["headline"]]), "(largest spread) and",
    format(chosen_origins[["alternative"]]), "(largest spread within",
    paste0(a1_5_config$alternative_period, ")\n"))

## Fitted values over the training window --------------------------------------
# Refits the same no-smooth window the forecast came from, which is deterministic
# Days are summed into 7-day blocks ending at the origin, so the training weeks fall on the
# same grid as the forecast weeks and the two sit on one continuous axis

window_fitted_weekly <- function(this_origin) {
  window_start <- this_origin - a1_5_config$window_weeks * 7 + 1

  weekly <- lapply(names(ch1_models), function(model_name) {
    model_fit <- fit_renewal_gam(dat, ch1_models[[model_name]], use_smooth = FALSE,
                                 family = ch1_family(),
                                 fit_from = window_start, fit_to = this_origin)

    # The forecast freezes covariates at the last fitted row, which must be the origin
    stopifnot(max(model_fit$model_data$date) == this_origin)

    # Refitting must reproduce ch1_rolling.R exactly, or the fit shown is not the fit scored
    stored <- stored_coefficients |>
      filter(origin == this_origin, model == model_name, !used_smooth) |>
      arrange(term)
    refitted <- model_fit$coefficients |> arrange(term)
    stopifnot(identical(stored$term, refitted$term),
              max(abs(stored$estimate - refitted$estimate)) < 1e-6)

    tibble(date     = model_fit$model_data$date,
           observed = model_fit$model_data$incidence,
           fitted   = as.numeric(fitted(model_fit$fit))) |>
      # Block 1 is the seven days ending at the origin, block 2 the seven before it
      mutate(days_before_origin = as.numeric(this_origin - date),
             week_index         = ceiling((days_before_origin + 1) / 7),
             week_end           = this_origin - (week_index - 1) * 7) |>
      group_by(week_end) |>
      summarise(observed = sum(observed), fitted = sum(fitted), n_days = n(),
                .groups = "drop") |>
      mutate(model = model_name)
  }) |> bind_rows()

  # Partial weeks would put a short total next to full ones and read as a drop
  stopifnot(all(weekly$n_days == 7))

  weekly |> mutate(model = factor(model, levels = names(ch1_models)))
}

## Forecast quantiles at one origin --------------------------------------------
# Same route as ch1_forecast_plots.R, so the ribbons match the main-body figures

origin_quantiles <- function(this_origin) {
  forecasts |>
    filter(origin == this_origin) |>
    as_forecast_sample(forecast_unit = c("model", "horizon", "origin")) |>
    as_forecast_quantile(probs = a1_5_config$probs) |>
    as_tibble() |>
    tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                       names_prefix = "q") |>
    mutate(model       = factor(model, levels = names(ch1_models)),
           target_date = this_origin + horizon * 7)
}

## Panel --------------------------------------------------------------------------
# One row of four panels, each a continuous weekly series running through the origin

build_panel <- function(this_origin) {

  fitted_weekly <- window_fitted_weekly(this_origin)
  fan           <- origin_quantiles(this_origin)

  period_label <- scores |> filter(origin == this_origin) |> slice(1) |> pull(period)

  # Observed weekly totals either side of the origin, as one series
  observed_weekly <- bind_rows(
    fitted_weekly |> distinct(week_end, observed) |> rename(date = week_end),
    fan |> distinct(target_date, observed) |> rename(date = target_date)
  ) |> distinct(date, observed)

  # Mean over horizons 1-4, matching how this origin is scored elsewhere
  origin_scores <- scores |>
    filter(scale == "log", origin == this_origin) |>
    group_by(model) |>
    summarise(crps = mean(crps), .groups = "drop") |>
    mutate(model = factor(model, levels = names(ch1_models)),
           label = sprintf("log-CRPS %.2f", crps))

  stopifnot(nrow(origin_scores) == length(ch1_models))

  ggplot(mapping = aes(x = date)) +
    geom_vline(xintercept = this_origin, linetype = "dashed", colour = "grey50") +
    geom_ribbon(data = fan, aes(x = target_date, ymin = q0.05, ymax = q0.95, fill = model),
                alpha = 0.25) +
    geom_ribbon(data = fan, aes(x = target_date, ymin = q0.25, ymax = q0.75, fill = model),
                alpha = 0.45) +
    geom_line(data = fan, aes(x = target_date, y = q0.5, colour = model), linewidth = 0.7) +
    geom_line(data = fitted_weekly, aes(x = week_end, y = fitted, colour = model),
              linewidth = 0.7) +
    geom_point(data = observed_weekly, aes(y = observed), colour = "grey20", size = 1.1) +
    geom_text(data = origin_scores, inherit.aes = FALSE,
              aes(x = min(observed_weekly$date), y = Inf, label = label),
              hjust = 0, vjust = 1.5, size = 3, colour = "grey20") +
    facet_wrap(~model, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    # Suffixes chosen per value, as weekly totals differ by an order of magnitude across origins
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k",
                                                     big.mark = ",")) +
    labs(title = sprintf("Model fit and forecast at %s (%s)",
                         format(this_origin, "%d %b %Y"), period_label),
         subtitle = paste("Weekly totals; fitted left of the origin, forecast right with",
                          "50% and 90% intervals, points observed"),
         x = NULL, y = "Weekly infections") +
    theme_minimal() +
    theme(legend.position = "none")
}

for (which_origin in names(chosen_origins)) {
  origin_date <- chosen_origins[[which_origin]]
  ggsave(file.path(a1_5_config$output_dir,
                   sprintf("fig_a1_5_%s_%s.png", which_origin, format(origin_date))),
         build_panel(origin_date), width = 12, height = 3.6, dpi = 300, bg = "white")
}

message("A.1.5 written to ", a1_5_config$output_dir)
