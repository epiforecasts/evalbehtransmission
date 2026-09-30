# Producing synthetic outbreaks to check whether s(t) absorbs signal belonging to covariates
# Addresses issue #55
# This script defines functions only, which are called from ch1_synthetic_experiment.qmd
# Daily time step throughout
# Paths are relative to the repo root, as in Chapter 1

source("analysis/ch1_gam.R") # GI weights, renewal formula, fit_renewal_gam(), fitted_rt(), plus mgcv and dplyr

## Config ----------------------------------------------------------------------

synthetic_config <- list(

  n_days = 180, # length of simulated outbreak, to revise as replicates need only a burn-in plus one 8-week window

  seed_days = 7, # days of constant incidence before the renewal step takes over
  seed_size = 50,

  # log(Rt) intercept, determining whether outbreak grows or shrinks
  # Covariates are z-scored, so this is log(Rt) at average study-period behaviour, not log(R0)
  # Slightly below 0 keeps the real mobility outbreak from growing throughout
  beta0 = -0.05,

  # 0 checks null behaviour (false positives)
  beta1 = c(0, 0.25),

  # mgcv nb() parameterisation, variance = mu + mu^2 / theta
  # Smaller theta means more overdispersion
  # For large counts the SD is about mu / sqrt(theta), so 10 gives about 30% day-to-day noise
  theta = 10, # held constant for now

  smooth_k = c(3, 5, 10), # basis dimension for s(t), Chapter 1 windows use 5

  covariate_noise_sd = 0.5, # noise added to the piecewise series

  # c(t) scale, and the period of the periodic form, one cycle per 8-week window
  trend_amplitude = 0.2,
  trend_period = 56,

  n_draws = 100, # 1000 needed for final results
  seed = 42,

  start_date = as.Date("2020-04-01") # arbitrary, fit_renewal_gam() needs a date column
)

## Covariates ------------------------------------------------------------------

# Piecewise-constant covariate
# s(t) may approximate this closely
# Pick levels so Rt moves either side of 1, avoiding blowup or elimination
make_piecewise_covariate <- function(n_days = synthetic_config$n_days, block_days = 14) {
  levels <- c(0.6, 1.0, 0.4, -0.5, -1.0, -0.6, 0.2, 0.8, 0.1, -0.4, -0.9, -0.3, 0.3)
  x <- rep(levels, each = block_days, length.out = n_days)
  (x - mean(x)) / sd(x)
}

# Piecewise series plus noise, so the two differ only by the noise
# Fluctuates more, so s(t) cannot immediately approximate it
make_noisy_covariate <- function(piecewise = make_piecewise_covariate(),
                                 noise_sd = synthetic_config$covariate_noise_sd,
                                 seed = synthetic_config$seed) {
  set.seed(seed)
  x <- piecewise + rnorm(length(piecewise), sd = noise_sd)
  (x - mean(x)) / sd(x)
}

# z-scored real Chapter 1 series
load_real_covariate <- function(stream = "mobility", path = ch1_gam_config$input_path) {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  dat[[stream]][!is.na(dat[[stream]])]
}

## Temporal trend c(t) ----------------------------------------------------------

# To make: c(t) for each form in the settings grid
# none = 0, periodic = a * sin(2 * pi * t / P), concurvity = a * smoothed covariate
make_temporal_trend <- function() {

}

## Simulation ------------------------------------------------------------------

# simulate the outbreak using the renewal formula
# feed each Lambda_t into the next day, to preserve autocorrelation and accumulate variation
# noise = FALSE returns expected counts, noise = TRUE draws them from nb(theta)
# temporal_trend is c(t), the unmeasured temporal driver, 0 by default
simulate_outbreak <- function(covariate,
                              beta1,
                              temporal_trend = 0,
                              noise = TRUE,
                              config = synthetic_config,
                              gi_weights = make_gi_weights()) {

  n_days <- length(covariate)
  rt <- exp(config$beta0 + beta1 * covariate + temporal_trend)

  incidence <- numeric(n_days)
  incidence[seq_len(config$seed_days)] <- config$seed_size

  for (t in (config$seed_days + 1):n_days) {
    # incidence[t] is still 0 here, and lags start at 1, so this is Lambda_t
    lambda <- compute_lambda_last(incidence[seq_len(t)], gi_weights)
    expected <- rt[t] * lambda
    # rounded when noiseless, as nb() needs integer counts
    incidence[t] <- if (noise) rnbinom(1, mu = expected, size = config$theta) else round(expected)
  }

  tibble(date = config$start_date + seq_len(n_days) - 1,
         covariate = covariate, temporal_trend = temporal_trend, rt = rt, incidence = incidence)
}

## Results ---------------------------------------------------------------------

# One row per fit: one model on one simulated outbreak
# Estimate, standard error, Wald interval, whether it covers beta1
# Also the edf of s(t)
# No default beta1, as each settings row passes its own true value
summarise_beta1_recovery <- function(fit, beta1) {

  coefficients <- summary(fit)$p.table
  estimate <- coefficients["covariate", "Estimate"]
  se       <- coefficients["covariate", "Std. Error"]

  tibble(estimate = estimate,
         se       = se,
         lower    = estimate - 1.96 * se,
         upper    = estimate + 1.96 * se,
         covers   = lower <= beta1 & beta1 <= upper,
         edf      = if (length(fit$smooth)) sum(summary(fit)$edf) else NA_real_)
}

## Settings --------------------------------------------------------------------

# Simulated truth (covariate, beta1, c(t), noise) crossed with fitted model (X_t, s(t), k)
# Deterministic outbreaks only require one draw
synthetic_settings <- expand.grid(covariate_type = c("piecewise", "noisy", "real"),
                                  beta1          = synthetic_config$beta1,
                                  temporal_trend = c("none", "periodic", "concurvity"),
                                  noise          = c(FALSE, TRUE),
                                  with_covariate = c(FALSE, TRUE),
                                  use_smooth     = c(FALSE, TRUE),
                                  smooth_k       = synthetic_config$smooth_k,
                                  stringsAsFactors = FALSE)

## Replicates ------------------------------------------------------------------

# To make: for each settings row, simulate a burn-in plus one 8-week window,
# fit the four models once, and repeat over n_draws noise draws
run_synthetic_replicates <- function() {

}

## Rolling windows -------------------------------------------------------------

# Adapted from run_window() in ch1_rolling.R
# Fits the four nested models on each window of one simulated outbreak
# Forecasts s(t) models by holding s(t_max) fixed, i.e. constant R(t)
# The log R(t) analogue of Mellor et al.'s constant growth rate
run_synthetic_windows <- function() {

}

# True log R(t) against fitted log R(t) per day, fitted from fitted_rt() in ch1_gam.R
# Optionally split into intercept, covariate effect and s(t)
# Terms are centred by mgcv, so s(t) carries variation and not level
log_rt_components <- function() {

}

## Run -------------------------------------------------------------------------

# To add: run replicates and rolling windows, and save results as .rds for the qmd to read
