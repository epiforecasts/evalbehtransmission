# Chapter 1, issue #41: renewal GAM engine.
# Fits log(Rt) = b0 [+ s(t)] + sum(bk * Xk) with log(Lambda_t) as a fixed offset.
# Defines functions only — the other ch1 model scripts source this file.

library(mgcv)
library(dplyr)
library(readr)
library(EpiEstim)

source("R/compute_lambda.R")
source("R/inc2prev_path.R")

## Config ----------------------------------------------------------------------

ch1_gam_config <- list(
  # Single generation interval throughout. The Delta switch to mean 4.7 from
  # ~May 2021 is not implemented, and sigma is an unsourced placeholder.
  gi_mean  = 5.5,
  gi_sd    = 2.1,
  gi_max   = 21,

  smooth_k = 20,

  input_path = "data-processed/ch1_covariates.csv"
)

ch1_models <- list(
  baseline = character(0),
  contacts = "contacts",
  mobility = "mobility",
  combined = c("contacts", "mobility")
)

# Observation family, chosen in ch1_diagnostics.R and used by every script after
# it. Stored uncalled so each fit gets a fresh object: nb() carries its theta.
ch1_family <- nb

## Generation interval ---------------------------------------------------------

# Drops the day-0 element, so the interval starts at day 1 and assumes no
# same-day transmission.
make_gi_weights <- function(config = ch1_gam_config) {
  si <- discr_si(k = 0:config$gi_max, mu = config$gi_mean, sigma = config$gi_sd)
  si <- si / sum(si)
  si[-1]
}

## Model frame -----------------------------------------------------------------

# Lambda_t and the lagged covariates are built here, before any model is
# specified, so the covariate set and lag are settled separately from the fit.
# Always pass the full incidence series: Lambda_t needs gi_max days of prior
# incidence, so it is truncated for the first days of any pre-filtered subset.
# Use fit_from/fit_to to restrict which rows are fitted, after Lambda_t is built.
build_model_frame <- function(dat, covariates, lag = 0, gi_weights,
                              fit_from = NULL, fit_to = NULL) {
  out <- dat |>
    arrange(date) |>
    mutate(
      Lambda_t   = compute_lambda(incidence, gi_weights),
      Lambda_t   = if_else(Lambda_t <= 0, NA_real_, Lambda_t),
      log_Lambda = log(Lambda_t),
      t          = as.numeric(date - min(date))
    )

  for (v in covariates) out[[v]] <- dplyr::lag(out[[v]], lag)

  if (!is.null(fit_from)) out <- filter(out, date >= fit_from)
  if (!is.null(fit_to))   out <- filter(out, date <= fit_to)

  out |>
    filter(!is.na(log_Lambda), !is.na(incidence),
           if_all(all_of(covariates), \(x) !is.na(x)))
}

## Formula ---------------------------------------------------------------------

# The only place the formula is assembled. The intercept is log(R0). s(t) cannot
# be extrapolated, so the smooth model is for explanation, not forecasting.
renewal_formula <- function(covariates, use_smooth = FALSE, k = 20) {
  smooth <- if (use_smooth) sprintf("s(t, k = %d)", k)
  reformulate(c("1", smooth, covariates, "offset(log_Lambda)"), response = "incidence")
}

## Fit -------------------------------------------------------------------------

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

# All parametric terms, including the intercept, which is log(R0). s(t) has no
# single coefficient so does not appear.
coefficient_table <- function(fit) {
  tab <- summary(fit)$p.table

  tibble(
    term     = rownames(tab),
    estimate = tab[, "Estimate"],
    se       = tab[, "Std. Error"],
    lower    = estimate - 1.96 * se,
    upper    = estimate + 1.96 * se
  )
}

# Rt is the linear predictor with the offset removed.
fitted_rt <- function(fit, model_data) {
  tibble(
    date = model_data$date,
    Rt   = exp(predict(fit, type = "link") - model_data$log_Lambda)
  )
}

## Toy check -------------------------------------------------------------------
# Runs only when this file is executed directly, not when another script sources
# it. Passes the full series and restricts the fit to five months, so Lambda_t
# is built with complete history.

if (sys.nframe() == 0) {

  dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
  from <- as.Date("2020-09-01")
  to   <- as.Date("2021-01-31")

  for (v in names(ch1_models)) {
    m <- fit_renewal_gam(dat, covariates = ch1_models[[v]],
                         fit_from = from, fit_to = to)

    cat("\n---", v, "---\n")
    cat(deparse1(m$formula), "\n")
    cat("n =", m$n_obs,
        "| deviance explained =", round(100 * m$dev_expl, 2), "%",
        "| Rt range:", round(range(m$fitted_rt$Rt), 2), "\n")
    print(as.data.frame(m$coefficients |>
            mutate(across(where(is.numeric), \(x) round(x, 3)))),
          row.names = FALSE)
  }

  m_smooth <- fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE,
                              fit_from = from, fit_to = to)
  cat("\n--- combined + s(t) ---\n")
  cat(deparse1(m_smooth$formula), "\n")
  cat("deviance explained =", round(100 * m_smooth$dev_expl, 2), "%\n")
}
