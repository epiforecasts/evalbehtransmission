# Chapter 1, issue #41: renewal GAM engine.
# Fits log(Rt) = b0 [+ s(t)] + sum(bk * Xk) with log(Lambda_t) as a fixed offset.
# This script defines functions only, for use in other scripts
# Includes building model data frame, formula, fit GAM (or GLM), coefficient tables, Rt etc

library(mgcv)
library(dplyr)
library(readr)
library(EpiEstim)

source("R/compute_lambda.R")
source("R/inc2prev_path.R")

## Config ----------------------------------------------------------------------

ch1_gam_config <- list(
  
  # Single generation interval used throughout
  # To do: add Delta switch to mean 4.7 around May 2021 + check inc2prev formulation/epiparameter
  gi_mean  = 5.5,
  
  gi_sd    = 2.1,
  gi_max   = 21, # Overly conservative given the proportion of infections that happen before 21d

  smooth_k = 20,

  input_path = "data-processed/ch1_covariates.csv"
)

# Define list of models to compare, based on covariates included
ch1_models <- list(
  baseline = character(0),
  contacts = "contacts", # CoMix-augmented model using the corresponding contact covariate
  mobility = "mobility",
  combined = c("contacts", "mobility")
)

# Set distribution family to negative binomial, after the diagnostics work in ch1_diagnostics.R
# revealed how overly narrow Poisson SEs were and poorly captured dispersion. NB involves theta parameter
ch1_family <- nb

## Generation interval ---------------------------------------------------------

# Drop day-0 elements, avoids circular counting and needing I_t to estimate I_t
make_gi_weights <- function(config = ch1_gam_config) {
  si <- discr_si(k = 0:config$gi_max, mu = config$gi_mean, sigma = config$gi_sd) # discr_si=0 at k=0 by default
  si <- si / sum(si)
  si[-1]
}

## Model frame -----------------------------------------------------------------

# Build data frame with infectiousness (Lambda_t) and lagged covariates for use in subsequent models
# Use the full incidence series when calculating Lambda_t as this needs gi_max days ideally
# Otherwise you end up with extreme Rt estimates when only very recent counts are included in convolution
# After building this entire frame, can restrict to just [fit_from, fit_to] period
build_model_frame <- function(dat, covariates, lag = 0, gi_weights,
                              fit_from = NULL, fit_to = NULL) {
  frame <- dat |>
    arrange(date) |>
    mutate(
      Lambda_t   = compute_lambda(incidence, gi_weights),
      Lambda_t   = if_else(Lambda_t <= 0, NA_real_, Lambda_t),
      log_Lambda = log(Lambda_t),
      t          = as.numeric(date - min(date))
    )

  # Apply lag to variables
  for (covariate in covariates) frame[[covariate]] <- dplyr::lag(frame[[covariate]], lag)

  if (!is.null(fit_from)) frame <- filter(frame, date >= fit_from)
  if (!is.null(fit_to))   frame <- filter(frame, date <= fit_to)

  frame <- frame |>
    filter(!is.na(log_Lambda), !is.na(incidence), # Trivially exclude days with missing counts
           if_all(all_of(covariates), \(x) !is.na(x))) # Keep rows where every covariate in the model is present

  # Lambda_t needs gi_max days of prior incidence, so the earliest windows are shorter than requested
  # Warn rather than let a truncated window pass unnoticed
  if (!is.null(fit_from) && !is.null(fit_to)) {
    requested <- as.numeric(fit_to - fit_from) + 1
    if (nrow(frame) < requested) {
      warning(sprintf("Window %s to %s: %d of %d days usable",
                      format(fit_from), format(fit_to), nrow(frame), requested),
              call. = FALSE)
    }
  }

  frame
}

## Formula ---------------------------------------------------------------------

# Function to assemble the renewal formula based on given covariates
# Can toggle smoothing term on/off, with parameter k
# Covariates are z-scored, so intercept = log(Rt) at mean covariate values, not log(R0)
# s(t) can't be extrapolated, so use this for explanation and not forecasting
renewal_formula <- function(covariates, use_smooth = FALSE, k) {
  smooth <- if (use_smooth) sprintf("s(t, k = %d)", k)
  reformulate(c("1", smooth, covariates, "offset(log_Lambda)"), response = "incidence")
}

# CAVEAT: covariates have different ranges in each window, given global z-scoring
# Cannot compare intercepts exactly equally across all windows as they extrapolate different severities

## Fit -------------------------------------------------------------------------

# Formula to fit GLM/GAM model and return model, formula, data, coefficients, deviance explained etc
# n_obs indicates number of data points used in a given model fit - affects performance and comparison
# A mobility model with complete data may trivially differ from a contact model with incomplete data

# Defaults: no behavioural covariates, lag, s(t) term, Poisson log link, and full date range
fit_renewal_gam <- function(data, covariates = character(0), lag = 0,
                            use_smooth = FALSE, family = poisson(link = "log"),
                            fit_from = NULL, fit_to = NULL,
                            config = ch1_gam_config) {

  model_data    <- build_model_frame(data, covariates, lag, make_gi_weights(config),
                                     fit_from, fit_to)
  model_formula <- renewal_formula(covariates, use_smooth, config$smooth_k)
  fit           <- gam(model_formula, family = family, data = model_data,
                       method = "REML")

  list(
    fit          = fit,
    formula      = model_formula,
    model_data   = model_data,
    coefficients = coefficient_table(fit),
    dev_expl     = summary(fit)$dev.expl,
    fitted_rt    = fitted_rt(fit, model_data),
    n_obs        = nrow(model_data)
  )
}

## Results from a fit ----------------------------------------------------------

# All parametric terms for GLM model, except s(t) with no single coefficient
# Intercept = log(Rt) at mean covariate values
coefficient_table <- function(fit) {
  tab <- summary(fit)$p.table

  tibble(
    term     = rownames(tab),
    estimate = tab[, "Estimate"],
    se       = tab[, "Std. Error"],
    lower    = estimate - 1.96 * se, # Gives uncertainty per coefficient, not joint mvrnorm uncertainty in forecasts
    upper    = estimate + 1.96 * se
  )
}

# exp(predict(fit)) gives exponential of the full linear predictor
# Subtract the offset before exponentiating for Rt estimate
fitted_rt <- function(fit, model_data) {
  tibble(
    date = model_data$date,
    # as.numeric() drops the matrix attributes predict() returns
    Rt   = as.numeric(exp(predict(fit, type = "link") - model_data$log_Lambda))
  )
}

## Toy check -------------------------------------------------------------------

# Runs only when this file is executed directly, not when other scripts use the functions above
# Pass full series and restrict fit to example five month period Sep 2020-Jan 2021

# NOTE: uses Poisson default which gives overly confident estimates
# SEs drop to 0.000, contacts flip to -0.001 in combined model, s(t) with k=20 hits 100% deviance explained
if (sys.nframe() == 0) {

  dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
  from <- as.Date("2020-09-01")
  to   <- as.Date("2021-01-31")

  for (model_name in names(ch1_models)) {
    fitted_model <- fit_renewal_gam(dat, covariates = ch1_models[[model_name]],
                                    fit_from = from, fit_to = to)

    cat("\n---", model_name, "---\n")
    cat(deparse1(fitted_model$formula), "\n")
    cat("n =", fitted_model$n_obs,
        "| deviance explained =", round(100 * fitted_model$dev_expl, 2), "%", # Fraction of null (intercept-only) deviance explained by the model
        "| Rt range:", round(range(fitted_model$fitted_rt$Rt), 2), "\n")
    print(as.data.frame(fitted_model$coefficients |>
            mutate(across(where(is.numeric), \(x) round(x, 3)))),
          row.names = FALSE)
  }

  # Introduce smoothing term to explore whether how much residual temporal variation this absorbs
  smooth_fit <- fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE,
                              fit_from = from, fit_to = to)
  cat("\n--- combined + s(t) ---\n")
  cat(deparse1(smooth_fit$formula), "\n")
  cat("deviance explained =", round(100 * smooth_fit$dev_expl, 2), "%\n")
}
