# Chapter 1, issues #43 and #44: rolling-origin window loop
# Each window is fitted once, giving both explanatory outputs and forecasts
# s(t) is fitted per window for explanation only, as a spline cannot be extrapolated

# MASS is used via MASS::mvrnorm() only, since attaching it masks dplyr::select.
source("analysis/ch1_gam.R")

## Config ----------------------------------------------------------------------

ch1_rolling_config <- list(
  window_weeks  = 8,
  step_days     = 7,
  first_origin  = as.Date("2020-07-01"), # Early enough that every window is full, given CIS start + gi_max
  # Pre-vaccine cut-off, so the last origin on the weekly grid is the final one at or before this date
  last_origin   = as.Date("2021-01-01"),
  horizon_weeks = 4,
  n_sim         = 200, # Number of simulations (samples) from each forecast origin

  # Forecasts are simulated, so fix the seed to make a run reproducible
  seed          = 42,
  family        = ch1_family, # Negative binomial, set in ch1_gam.R

  # Far smaller than the full-period fit
  # k = 20 over a window this short is close to interpolation, and the NB fit fails to converge
  smooth_k      = 5,

  coef_path     = "data-processed/ch1_window_coefficients.csv",
  deviance_path = "data-processed/ch1_window_deviance.csv",
  forecast_path = "data-processed/ch1_forecasts.csv",

  deviance_table_path = "outputs/ch1/table_deviance_explained.csv"
)

# ch1_gam_config, but with smaller k for shorter window
window_config <- modifyList(ch1_gam_config, list(smooth_k = ch1_rolling_config$smooth_k))

## Covariates at the forecast origin --------------------------------------------
# Every covariate is carried forward from its origin value across the horizon (LOCF)
# log(Rt) is constant across the forecast, so one design row serves every horizon day
# Lambda_t is not projected here, but recomputed each day inside the renewal loop
# Returns a single row holding one value per covariate, taken from the origin date

covariates_at_origin <- function(model_data, covariates) {
  origin_row <- tail(model_data, 1) # Last fitted day of the window, being the origin
  # predict() requires every variable in the formula, but lpmatrix drops the offset
  # log_Lambda is a filler (actual Lambda_t calculated in later loop)
  newdata <- data.frame(log_Lambda = 0, t = origin_row$t)
  for (covariate in covariates) newdata[[covariate]] <- origin_row[[covariate]]
  newdata
}

## Forecast one window ---------------------------------------------------------
# Parameter uncertainty from the coefficient covariance, observation noise from the fitted family
# Incidence is projected forward through the renewal equation
# Returns predicted incidence for each simulation and horizon day

simulate_forecast <- function(fitted_model, covariates, incidence_history,
                              horizon_days, n_sim, gi_weights, theta) {

  # One row of the model matrix at the origin covariate values (X_last in model_rtglm.R)
  # Multiplying it by a coefficient draw gives log(Rt) for that draw
  design_row <- predict(fitted_model$fit,
                        newdata = covariates_at_origin(fitted_model$model_data, covariates),
                        type = "lpmatrix")

  # n_sim draws from the asymptotic multivariate normal of the fitted coefficients
  # The intercept and betas are the estimated terms, so they carry the parameter uncertainty
  coef_draws <- MASS::mvrnorm(n_sim, mu = coef(fitted_model$fit),
                              Sigma = vcov(fitted_model$fit))
  if (n_sim == 1) coef_draws <- matrix(coef_draws, nrow = 1) # MASS::mvrnorm() returns numeric vector otherwise

  # One simulated trajectory per coefficient draw
  lapply(seq_len(n_sim), function(sim) {
    log_rt   <- as.numeric(design_row %*% coef_draws[sim, ])
    inc_path <- incidence_history # Clearer naming for trajectory

    # Step forward one day at a time across the horizon, recomputing Lambda_t as the path grows
    for (day in seq_len(horizon_days)) {
      Lambda_t <- compute_lambda_last(inc_path, gi_weights)
      expected <- exp(log_rt) * pmax(Lambda_t, 1e-8) # Floor keeps Lambda_t positive if the path hits zero
      inc_path <- c(inc_path, rnbinom(1, mu = expected, size = theta)) # NB observation noise sampled per day per trajectory, fed forward into renewal
    }

    tibble(sample_id = sim,
           day       = seq_len(horizon_days),
           predicted = tail(inc_path, horizon_days)) # Drops the history, keeping just the forecast
  }) |> bind_rows()
}

## One origin ------------------------------------------------------------------
# Fits every model on the window ending at the origin, then forecasts past it
# Returns coefficients, deviance explained and forecasts, with rows for all four models

run_window <- function(dat, origin, rolling_config = ch1_rolling_config,
                       gi_weights = make_gi_weights()) {

  window_start <- origin - rolling_config$window_weeks * 7 + 1
  horizon_days <- rolling_config$horizon_weeks * 7

  # True observed incidence over the horizon, for scoring
  observed_horizon <- dat |> filter(date > origin) |> head(horizon_days)
  if (nrow(observed_horizon) < horizon_days) {
    message("  origin ", format(origin), ": forecast horizon longer than observed horizon, skipping")
    return(NULL)
  }

  incidence_history <- dat |> filter(date <= origin) |> pull(incidence)

  results <- lapply(names(ch1_models), function(model_name) {
    covariates <- ch1_models[[model_name]]

    fits <- list(
      no_smooth = fit_renewal_gam(dat, covariates, family = rolling_config$family(),
                                  fit_from = window_start, fit_to = origin,
                                  config = window_config),
      smooth    = fit_renewal_gam(dat, covariates, use_smooth = TRUE,
                                  family = rolling_config$family(),
                                  fit_from = window_start, fit_to = origin,
                                  config = window_config)
    )

    # The forecast freezes covariates at the last fitted row, which must be the origin
    stopifnot(max(fits$no_smooth$model_data$date) == origin)

    # One row per parametric term, split by smoothing or not
    coefficient_rows <- lapply(names(fits), function(smooth_setting) {
      fits[[smooth_setting]]$coefficients |>
        mutate(used_smooth = (smooth_setting == "smooth"))
    }) |>
      bind_rows() |>
      mutate(origin = origin, model = model_name)

    # Deviance explained by each fit, split by smoothing or not
    deviance_rows <- tibble(
      origin      = origin,
      model       = model_name,
      used_smooth = c(FALSE, TRUE),
      dev_expl    = c(fits$no_smooth$dev_expl, fits$smooth$dev_expl)
    )

    # Forecast from the no-smooth fit only, where Var = mu + mu^2/theta
    # theta is a point estimate on the raw scale, and small values mean overdispersion
    # Poisson has no theta, and Inf recovers the Poisson limit so the family can be swapped
    theta <- if (is.null(fits$no_smooth$fit$family$getTheta)) Inf else fits$no_smooth$fit$family$getTheta(TRUE)
    forecasts <- simulate_forecast(fits$no_smooth, covariates, incidence_history,
                                   horizon_days, rolling_config$n_sim, gi_weights, theta) |>
      mutate(origin = origin, model = model_name,
             target_date = observed_horizon$date[day],
             observed = observed_horizon$incidence[day])

    list(coefficients = coefficient_rows, deviance = deviance_rows, forecasts = forecasts)
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

set.seed(ch1_rolling_config$seed)

message("Running ", length(origins), " origins x ", length(ch1_models), " models")

windows <- lapply(seq_along(origins), function(i) {
  if (i %% 10 == 0) message("  origin ", i, "/", length(origins)) # Progress, every 10th origin
  run_window(dat, origins[i])
})
windows <- windows[!sapply(windows, is.null)] # Drops any skipped origins

## Assemble and save -----------------------------------------------------------

# Each element of windows is one origin, already covering all four models
# This puts every origin into one table
window_coefficients <- bind_rows(lapply(windows, \(w) w$coefficients)) |>
  select(origin, model, term, used_smooth, estimate, se, lower, upper)

# Increment over the baseline for the same origin and smooth setting
window_deviance <- bind_rows(lapply(windows, \(w) w$deviance)) |>
  group_by(origin, used_smooth) |>
  mutate(dev_expl_increment = dev_expl - dev_expl[model == "baseline"]) |>
  ungroup()

# Aggregate daily forecasts to weekly forecasts (1-4 week horizons)
# Sum across individual trajectories (indexed by sample_id), preserving day-to-day autocorrelation
forecasts <- bind_rows(lapply(windows, \(w) w$forecasts)) |>
  mutate(horizon = ceiling(day / 7)) |> # Week of horizon (1-4)
  group_by(origin, model, sample_id, horizon) |>
  summarise(predicted   = sum(predicted),
            observed    = sum(observed),
            target_date = max(target_date), # Last day of the horizon week
            .groups = "drop")

write_csv(window_coefficients, ch1_rolling_config$coef_path)
write_csv(window_deviance, ch1_rolling_config$deviance_path)
write_csv(forecasts, ch1_rolling_config$forecast_path)

cat("\norigins run:", length(windows),
    "| coefficient rows:", nrow(window_coefficients),
    "| deviance rows:", nrow(window_deviance),
    "| forecast rows:", nrow(forecasts), "\n")

# Averaged over windows, so baseline against a covariate model with and without s(t)
# dev_expl_increment is dropped here, as baseline explains nothing without s(t) and the
# increment then equals the model's own deviance explained
deviance_summary <- window_deviance |>
  group_by(model, used_smooth) |>
  summarise(mean_dev_expl = round(100 * mean(dev_expl), 1), .groups = "drop") |>
  tidyr::pivot_wider(names_from = used_smooth, values_from = mean_dev_expl,
                     names_prefix = "smooth_") |>
  rename(no_smooth = smooth_FALSE, with_smooth = smooth_TRUE) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(model)

cat("\nMean deviance explained by model, with and without s(t):\n")
print(as.data.frame(deviance_summary), row.names = FALSE)

write_csv(deviance_summary, ch1_rolling_config$deviance_table_path)
