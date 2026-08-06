# Chapter 1, issues #43 and #44: rolling-origin window loop.
# Each window is fitted once and yields both explanatory outputs (coefficients,
# deviance) and forecasts. s(t) is fitted per window for explanation only, since
# a spline cannot be extrapolated past the training range.

# MASS is used via MASS::mvrnorm() only, since attaching it masks dplyr::select.
source("analysis/ch1_gam.R")

## Config ----------------------------------------------------------------------

ch1_rolling_config <- list(
  window_weeks  = 12,
  step_days     = 7,
  first_origin  = as.Date("2020-07-01"),
  last_origin   = as.Date("2021-10-31"),
  horizon_weeks = 4,
  n_sim         = 200,

  # Stored uncalled so each fit gets a fresh family object: nb() carries its
  # estimated theta.
  family        = nb,

  # Far smaller than the full-period fit: k = 20 over an 84-day window is close
  # to interpolation and the NB fit fails to converge.
  smooth_k      = 5,

  coef_path     = "data-processed/ch1_window_coefficients.csv",
  deviance_path = "data-processed/ch1_window_deviance.csv",
  forecast_path = "data-processed/ch1_forecasts.csv"
)

## Covariate projection --------------------------------------------------------
# The only place future covariate values are set: every covariate is held at its
# origin value for the whole horizon, so log(Rt) is constant across the forecast
# and one design row serves every horizon day. Lambda_t is the only quantity that
# moves, and it is recomputed in the renewal loop rather than predicted here.

project_covariates <- function(model_data, covariates) {
  origin_row <- tail(model_data, 1)
  out <- data.frame(log_Lambda = 0, t = origin_row$t)
  for (covariate in covariates) out[[covariate]] <- origin_row[[covariate]]
  out
}

## Forecast one window ---------------------------------------------------------
# Parameter uncertainty from the coefficient covariance, observation noise from
# the fitted family. Incidence is projected through the renewal equation.

simulate_forecast <- function(fitted_model, covariates, incidence_history,
                              horizon_days, n_sim, gi_weights, theta) {

  design_row <- predict(fitted_model$fit,
                        newdata = project_covariates(fitted_model$model_data, covariates),
                        type = "lpmatrix")

  coef_draws <- MASS::mvrnorm(n_sim, mu = coef(fitted_model$fit),
                              Sigma = vcov(fitted_model$fit))
  if (n_sim == 1) coef_draws <- matrix(coef_draws, nrow = 1)

  lapply(seq_len(n_sim), function(sim) {
    log_rt   <- as.numeric(design_row %*% coef_draws[sim, ])
    inc_path <- incidence_history

    for (day in seq_len(horizon_days)) {
      lambda   <- compute_lambda_last(inc_path, gi_weights)
      expected <- exp(log_rt) * pmax(lambda, 1e-8)
      inc_path <- c(inc_path, rnbinom(1, mu = expected, size = theta))
    }

    tibble(sample_id = sim,
           day       = seq_len(horizon_days),
           predicted = tail(inc_path, horizon_days))
  }) |> bind_rows()
}

## One origin ------------------------------------------------------------------

run_window <- function(dat, origin, config = ch1_rolling_config,
                       gi_weights = make_gi_weights()) {

  window_start <- origin - config$window_weeks * 7 + 1
  horizon_days <- config$horizon_weeks * 7

  observed <- dat |> filter(date > origin) |> head(horizon_days)
  if (nrow(observed) < horizon_days) return(NULL)

  incidence_history <- dat |> filter(date <= origin) |> pull(incidence)

  results <- lapply(names(ch1_models), function(model_name) {
    covariates <- ch1_models[[model_name]]

    window_config <- modifyList(ch1_gam_config, list(smooth_k = config$smooth_k))

    fits <- list(
      no_smooth = fit_renewal_gam(dat, covariates, family = config$family(),
                                  fit_from = window_start, fit_to = origin),
      smooth    = fit_renewal_gam(dat, covariates, use_smooth = TRUE,
                                  family = config$family(),
                                  fit_from = window_start, fit_to = origin,
                                  config = window_config)
    )

    coefficients <- bind_rows(lapply(names(fits), function(s) {
      fits[[s]]$coefficients |> mutate(used_smooth = s == "smooth")
    })) |> mutate(origin = origin, model = model_name)

    deviance <- tibble(
      origin      = origin,
      model       = model_name,
      used_smooth = c(FALSE, TRUE),
      dev_expl    = c(fits$no_smooth$dev_expl, fits$smooth$dev_expl)
    )

    # Forecast from the no-smooth fit only
    theta <- fits$no_smooth$fit$family$getTheta(TRUE)
    forecasts <- simulate_forecast(fits$no_smooth, covariates, incidence_history,
                                   horizon_days, config$n_sim, gi_weights, theta) |>
      mutate(origin = origin, model = model_name,
             target_date = observed$date[day],
             observed = observed$incidence[day])

    list(coefficients = coefficients, deviance = deviance, forecasts = forecasts)
  })

  list(
    coefficients = bind_rows(lapply(results, \(r) r$coefficients)),
    deviance     = bind_rows(lapply(results, \(r) r$deviance)),
    forecasts    = bind_rows(lapply(results, \(r) r$forecasts))
  )
}

## Run all origins -------------------------------------------------------------

dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)

origins <- seq(ch1_rolling_config$first_origin, ch1_rolling_config$last_origin,
               by = ch1_rolling_config$step_days)

message("Running ", length(origins), " origins x ", length(ch1_models), " models")

windows <- lapply(seq_along(origins), function(i) {
  if (i %% 10 == 0) message("  origin ", i, "/", length(origins))
  run_window(dat, origins[i])
})
windows <- Filter(Negate(is.null), windows)

## Assemble and save -----------------------------------------------------------

window_coefficients <- bind_rows(lapply(windows, \(w) w$coefficients)) |>
  select(origin, model, term, used_smooth, estimate, se, lower, upper)

# Increment over the baseline for the same origin and smooth setting
window_deviance <- bind_rows(lapply(windows, \(w) w$deviance)) |>
  group_by(origin, used_smooth) |>
  mutate(dev_expl_increment = dev_expl - dev_expl[model == "baseline"]) |>
  ungroup()

# Daily samples aggregated to weekly totals, matching the 1-4 week horizons
forecasts <- bind_rows(lapply(windows, \(w) w$forecasts)) |>
  mutate(horizon = ceiling(day / 7)) |>
  group_by(origin, model, sample_id, horizon) |>
  summarise(predicted   = sum(predicted),
            observed    = sum(observed),
            target_date = max(target_date),
            .groups = "drop")

write_csv(window_coefficients, ch1_rolling_config$coef_path)
write_csv(window_deviance, ch1_rolling_config$deviance_path)
write_csv(forecasts, ch1_rolling_config$forecast_path)

cat("\norigins run:", length(windows),
    "| coefficient rows:", nrow(window_coefficients),
    "| deviance rows:", nrow(window_deviance),
    "| forecast rows:", nrow(forecasts), "\n")

cat("\nMean deviance explained by model (no smooth):\n")
print(as.data.frame(window_deviance |> filter(!used_smooth) |>
  group_by(model) |>
  summarise(mean_dev_expl  = round(100 * mean(dev_expl), 1),
            mean_increment = round(100 * mean(dev_expl_increment), 1),
            .groups = "drop")), row.names = FALSE)
