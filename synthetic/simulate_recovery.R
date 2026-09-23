# Producing synthetic outbreaks to check whether s(t) absorbs signal belonging to covariates
# Addresses issue #55
# This script defines functions only, which are called from synthetic_check.qmd
# Daily time step throughout

library(mgcv)
library(dplyr)
library(EpiEstim)

source("../R/compute_lambda.R") # Retrieve function to compute total infectiousness

## Config ----------------------------------------------------------------------

synthetic_config <- list(

  n_days = 180, # length of simulated outbreak

  seed_days = 7, # days of constant incidence before the renewal step takes over
  seed_size = 50,

  # log(Rt) intercept, determining whether outbreak grows or shrinks
  # Covariates are z-scored, so 0 puts Rt at 1 when the covariate sits at its mean
  beta0 = 0,

  # Initially set coefficient to 0 to check null behaviour (false positives)
  beta1 = 0,

  # mgcv nb() parameterisation, variance = mu + mu^2 / theta
  # Smaller theta means more overdispersion, 2 initially limits noise
  theta = 2,

  smooth_k = 20, # basis dimension for s(t)

  # Matches ch1_gam_config
  gi_mean = 5.5,
  gi_sd   = 2.1,
  gi_max  = 21
)

## Generation interval ---------------------------------------------------------

make_gi_weights <- function(config = synthetic_config) {
  si <- discr_si(k = 0:config$gi_max, mu = config$gi_mean, sigma = config$gi_sd)
  si <- si / sum(si)
  si[-1] # drops day 0, so I_t does not enter its own Lambda_t
}

## Covariates ------------------------------------------------------------------

# Piecewise-constant covariate
# s(t) may approximate this closely
# Pick levels so Rt moves either side of 1, avoiding blowup or elimination
make_piecewise_covariate <- function(n_days = synthetic_config$n_days, block_days = 14) {
  levels <- c(0.6, 1.0, 0.4, -0.5, -1.0, -0.6, 0.2, 0.8, 0.1, -0.4, -0.9, -0.3, 0.3)
  x <- rep(levels, each = block_days, length.out = n_days)
  (x - mean(x)) / sd(x)
}

# make the covariate fluctuate more, so that s(t) cannot immediately approximate it
make_rough_covariate <- function(n_days = synthetic_config$n_days, block_days = 3) {
  x <- rep(rnorm(ceiling(n_days / block_days)), each = block_days)[seq_len(n_days)]
  (x - mean(x)) / sd(x)
}

## Simulation ------------------------------------------------------------------

# simulate the outbreak using the renewal formula
# feed each Lambda_t into the next day, to preserve autocorrelation and accumulate variation
# noise = FALSE returns expected counts, noise = TRUE draws them from nb(theta)
simulate_outbreak <- function(covariate,
                              beta1 = synthetic_config$beta1,
                              noise = TRUE,
                              config = synthetic_config,
                              gi_weights = make_gi_weights(config)) {

  n_days <- length(covariate)
  rt <- exp(config$beta0 + beta1 * covariate)

  incidence <- numeric(n_days)
  incidence[seq_len(config$seed_days)] <- config$seed_size

  for (t in (config$seed_days + 1):n_days) {
    # incidence[t] is still 0 here, and lags start at 1, so this is Lambda_t
    lambda <- compute_lambda_last(incidence[seq_len(t)], gi_weights)
    expected <- rt[t] * lambda
    # rounded when noiseless, as nb() needs integer counts
    incidence[t] <- if (noise) rnbinom(1, mu = expected, size = config$theta) else round(expected)
  }

  tibble(t = seq_len(n_days), covariate = covariate, rt = rt, incidence = incidence)
}

## Fitting ---------------------------------------------------------------------

# Fits log(Rt) = beta0 [+ s(t, k)] + beta1 * X, with log_Lambda as a fixed offset
# sim is one outbreak from simulate_outbreak, holding t, covariate, rt and incidence
fit_synthetic_renewal <- function(sim,
                                  use_smooth = FALSE,
                                  config = synthetic_config,
                                  gi_weights = make_gi_weights(config)) {

  frame <- sim |>
    mutate(Lambda_t   = compute_lambda(incidence, gi_weights),
           log_Lambda = log(Lambda_t)) |>
    filter(is.finite(log_Lambda))

  smooth <- if (use_smooth) sprintf("s(t, k = %d)", config$smooth_k)
  model_formula <- reformulate(c("1", smooth, "covariate", "offset(log_Lambda)"),
                               response = "incidence")

  gam(model_formula, family = nb(), data = frame, method = "REML")
}

## Results ---------------------------------------------------------------------

# One row per fit: one model on one simulated outbreak
# Estimate, standard error, Wald interval, whether it covers beta1
# Also the edf of s(t)
summarise_beta1_recovery <- function(fit, beta1 = synthetic_config$beta1) {

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

# Repeats simulate-and-fit from the same known beta1, stacking the rows above
# Varies covariate type, s(t) on/off and noise on/off
replicate_beta1_recovery <- function(n_draws = 100,
                                     beta1 = synthetic_config$beta1,
                                     config = synthetic_config) {

  gi_weights <- make_gi_weights(config)

  settings <- expand.grid(covariate_type = c("piecewise", "rough"),
                          use_smooth     = c(FALSE, TRUE),
                          noise          = c(FALSE, TRUE),
                          stringsAsFactors = FALSE)

  
  lapply(seq_len(nrow(settings)), function(setting_row) {
    setting <- settings[setting_row, ]

    lapply(seq_len(n_draws), function(draw) {
      covariate <- if (setting$covariate_type == "piecewise") {
        make_piecewise_covariate(config$n_days)
      } else {
        make_rough_covariate(config$n_days)
      }

      sim <- simulate_outbreak(covariate, beta1, setting$noise, config, gi_weights)
      fit <- fit_synthetic_renewal(sim, setting$use_smooth, config, gi_weights)

      summarise_beta1_recovery(fit, beta1) |>
        mutate(covariate_type = setting$covariate_type,
               use_smooth     = setting$use_smooth,
               noise          = setting$noise,
               draw           = draw)
    }) |> bind_rows()
  }) |> bind_rows()
}

# Splits fitted log Rt per day into intercept, covariate effect and s(t)
# Returns the true covariate effect alongside them, as a per-day series
# Terms are centred by mgcv, so s(t) carries variation and not level
log_rt_components <- function() {

}
