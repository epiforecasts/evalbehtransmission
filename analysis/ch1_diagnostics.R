# Chapter 1, issue #42: run baseline diagnostics on model fit and suitability of observation family
# i.e. me stress-testing whether Poisson is suitable, whether residuals are autocorrelated
# and whether the renewal model is even appropriate for inc2prev data

source("analysis/ch1_gam.R")

library(ggplot2)
library(patchwork)

## Config ----------------------------------------------------------------------

ch1_diag_config <- list(
  periods_path = "data-processed/ch1_periods.csv",
  plot_dir     = "outputs/ch1",
  smooth_k     = c(5, 10, 20, 40), # Test how readily a spline term will just absorb everything

  # Matches window_weeks in ch1_rolling.R, for the window-level checks below
  window_weeks = 8,

  # Matches the rolling origin range, so the window checks cover the same span
  first_origin = as.Date("2020-07-01"),
  last_origin  = as.Date("2021-01-01")
)

# Chosen from the family comparison below, and used for everything after it.
chosen_family <- "nb"   # ch1_family in ch1_gam.R must match

make_family <- function(family_name) {
  if (family_name == "poisson") poisson(link = "log") else nb()
}

# Prints numeric columns rounded, 3dp unless overridden
print_table <- function(x, digits = 3) {
  print(as.data.frame(mutate(x, across(where(is.numeric), \(v) round(v, digits)))),
        row.names = FALSE)
}

## Fit every model under both families, once -----------------------------------

dat     <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
periods <- read_csv(ch1_diag_config$periods_path, show_col_types = FALSE) |>
  select(date, period)

# Ensures periods are in chronological order, otherwise character handled alphabetically
period_levels <- periods |> arrange(date) |> distinct(period) |> pull(period)

# Fit each model with Poisson and NegBin specification
# Done across the entire period, with the window-level dispersion check below as the counterpart
fits <- lapply(c(poisson = "poisson", nb = "nb"), function(family_name) {
  lapply(ch1_models, \(covariates) {
    fit_renewal_gam(dat, covariates, family = make_family(family_name))
  })
})

## Family choice ---------------------------------------------------------------

# Pearson dispersion measures whether variance in model matches expected theoretical variance
# A value of 1 indicates about as variable as family assumes - good
# > 1 = overdispersion, possibly due to unobserved heterogeneity
# < 1 = underdispersion
dispersion <- function(fitted_model) {
  sum(residuals(fitted_model$fit, type = "pearson")^2) / df.residual(fitted_model$fit)
}

# Poisson dispersion far above 1 shows its variance assumption is much too tight
# NB lands near 1 by construction, as theta is fitted, so the Poisson value carries the evidence
family_comparison <- tibble(
  family     = c("poisson", "nb"),
  dispersion = c(dispersion(fits$poisson$baseline), dispersion(fits$nb$baseline)),
  aic        = c(AIC(fits$poisson$baseline$fit), AIC(fits$nb$baseline$fit)),
  intercept  = c(fits$poisson$baseline$coefficients$estimate,
                 fits$nb$baseline$coefficients$estimate),
  se         = c(fits$poisson$baseline$coefficients$se,
                 fits$nb$baseline$coefficients$se)
)

cat("\n--- Baseline under each family ---\n")
print_table(family_comparison, 4)

# Poisson gives incredibly confident coefficient estimates compared to NB
# This implies forecasts under a Poisson specification will be very narrow
# Generally, intercept absorbs some of the level of transmission when contacts/mobility aren't at their average
covariate_estimates_by_family <- lapply(names(fits), function(family_name) {
  lapply(names(ch1_models), function(model_name) {
    fits[[family_name]][[model_name]]$coefficients |>
      filter(term != "(Intercept)") |> # Removes baseline model, as only interested in coefficients
      mutate(model = model_name, family = family_name)
  }) |> bind_rows()
}) |> bind_rows() |>
  select(model, term, family, estimate, se) |>
  tidyr::pivot_wider(names_from = family, values_from = c(estimate, se)) |>
  mutate(se_ratio = se_nb / se_poisson)

cat("\n--- Covariate estimates and SEs by family ---\n")
print_table(covariate_estimates_by_family, 4)

# The fits above span the whole period, so they absorb non-stationarity the models never see
# Repeat on non-overlapping windows, reflecting the model training window
# Poisson dispersion stays between 36 and 1754, so it is not an artefact of the longer fit
window_origins <- seq(ch1_diag_config$first_origin, ch1_diag_config$last_origin,
                      by = ch1_diag_config$window_weeks * 7)

window_dispersion <- lapply(window_origins, function(origin) {
  window_start <- origin - ch1_diag_config$window_weeks * 7 + 1
  window_fits <- lapply(c(poisson = "poisson", nb = "nb"), function(family_name) {
    fit_renewal_gam(dat, ch1_models$baseline, family = make_family(family_name),
                    fit_from = window_start, fit_to = origin)
  })
  tibble(origin  = origin,
         poisson = dispersion(window_fits$poisson),
         nb      = dispersion(window_fits$nb))
}) |> bind_rows()

cat("\n--- Baseline dispersion on non-overlapping training windows ---\n")
print_table(window_dispersion, 2)

## Residuals by period ---------------------------------------------------------

residuals_table <- lapply(names(ch1_models), function(model_name) {
  fitted_model <- fits[[chosen_family]][[model_name]]
  tibble(date     = fitted_model$model_data$date,
         residual = residuals(fitted_model$fit, type = "deviance"),
         model    = model_name) |>
    left_join(periods, by = "date")
}) |> bind_rows() |>
  mutate(model  = factor(model, levels = names(ch1_models)),
         period = factor(period, levels = period_levels))

# Residuals are grouped by period after fitting, rather than models being fit per period (which biases results)
# Covariates in this pooled fit help in some periods, improving Lockdown 3 from -0.94 to 0.52 while Winter tiers worsens
# n is reported as periods differ in length, and Lockdown 1 only covers 8 days here
residuals_by_period <- residuals_table |>
  group_by(period, model) |>
  summarise(mean_resid = mean(residual), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = mean_resid) |>
  left_join(residuals_table |> filter(model == "combined") |> count(period, name = "n"),
            by = "period")

cat("\n--- Mean deviance residual by period ---\n")
print_table(residuals_by_period, 2)

## Autocorrelation -------------------------------------------------------------
# Correlated residuals break the independence assumption behind the SEs
# If adding covariates reduce autocorrelation, they may capture structural/temporal variation beyond the baseline model
# No day-of-week spike expected at lag 7, as the mobility and incidence series are both smoothed
# Contacts cut lag-14 from 0.48 to 0.23, so covariates absorb some of the temporal structure

# acf() assumes evenly spaced observations, and rows with missing covariates are dropped
# Confirm each model frame is contiguous before reading anything into the lags
stopifnot(vapply(ch1_models, function(covariates) {
  all(diff(build_model_frame(dat, covariates, 0, make_gi_weights())$date) == 1) # TRUE if every retained date is one day apart
}, logical(1)))

acfs <- lapply(names(ch1_models), function(model_name) {
  acf_out <- acf(residuals(fits[[chosen_family]][[model_name]]$fit, type = "deviance"),
                 lag.max = 28, plot = FALSE)
  tibble(lag = as.numeric(acf_out$lag), acf = as.numeric(acf_out$acf),
         model = model_name)
}) |> bind_rows()

pooled_acf <- acfs |>
  filter(lag %in% c(1, 7, 14)) |>
  tidyr::pivot_wider(names_from = lag, values_from = acf, names_prefix = "lag_")

cat("\n--- Residual autocorrelation ---\n")
print_table(pooled_acf)

# The ACF above uses one fit across the whole period
# Systematic over- or under-prediction across long stretches produces autocorrelation on its own
# Refitting per window removes that source, and lag 1 still stays at 0.72-0.96 throughout
# Covariates cut lag 7 and 14 in most windows, with lag 14 negative in some fits
# Baseline decays almost linearly where Rt moved monotonically
# Covariates decay faster in those same windows
window_acf <- lapply(window_origins, function(origin) {
  window_start <- origin - ch1_diag_config$window_weeks * 7 + 1
  lapply(names(ch1_models), function(model_name) {
    window_fit <- fit_renewal_gam(dat, ch1_models[[model_name]],
                                  family = make_family(chosen_family),
                                  fit_from = window_start, fit_to = origin)
    acf_out <- acf(residuals(window_fit$fit, type = "deviance"), lag.max = 28, plot = FALSE)
    tibble(origin = origin,
           model  = model_name,
           lag    = as.numeric(acf_out$lag),
           acf    = as.numeric(acf_out$acf))
  }) |> bind_rows()
}) |> bind_rows()

window_acf_summary <- window_acf |>
  filter(lag %in% c(1, 7, 14)) |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  arrange(origin, model) |>
  tidyr::pivot_wider(names_from = lag, values_from = acf, names_prefix = "lag_")

cat("\n--- Autocorrelation on non-overlapping training windows ---\n")
print_table(window_acf_summary)

## Validation against inc2prev Rt ----------------------------------------------

# https://github.com/epiforecasts/inc2prev
# Not a forwards renewal model: infections are a Gaussian process convolved with a PCR detection curve
# Rt is derived afterwards as R[s] = infections[s + seeding_time] / infectiousness[s] (stan/functions/rt.stan)
# R is then averaged over a 3-day centred window, via smooth = 1 (stan/inc2prev.stan:209)
# Their generation time is redrawn per posterior sample, mean 3.64 and sd 3.08 truncated at 15 (R/model.R:20-23)
# Ours is fixed at mean 5.5 and sd 2.1 truncated at 21, from EpiEstim::discr_si
# Only medians are compared, taken from each series after the fact

combined_fit <- fits[[chosen_family]]$combined # Uses all covariates, unlike naive

# Copy their discretisation, being plain gamma CDF differences over 1..max, renormalised
inc2prev_gi <- function(mu, sd, max_val) {
  pmf <- diff(pgamma(1:(max_val + 1), shape = (mu / sd)^2, rate = mu / sd^2))
  pmf / sum(pmf)
}

rt_compare <- bind_rows(
  read_csv(inc2prev_path("outputs/estimates_national.csv"), show_col_types = FALSE) |>
    filter(variable == "England", name == "R") |>
    transmute(date, Rt = q50, source = "inc2prev"), # Median estimate for Rt
  combined_fit$model_data |>
    transmute(date, Rt = incidence / Lambda_t, source = "naive"), # Point estimate from our initial renewal
  combined_fit$fitted_rt |> mutate(source = "fitted"),
  build_model_frame(dat, character(0), 0, inc2prev_gi(3.64, 3.08, 15)) |> # Calculate naive Rt copying their GI
    transmute(date, Rt = incidence / Lambda_t, source = "naive, inc2prev GI")
) |>
  filter(date >= min(combined_fit$model_data$date),
         date <= max(combined_fit$model_data$date))

rt_wide <- rt_compare |> tidyr::pivot_wider(names_from = source, values_from = Rt)

# Compare Rt estimates from inc2prev, our original renewal + covariate fit, and ours copying their GI
rt_agreement <- lapply(c("naive", "fitted", "naive, inc2prev GI"), function(series_name) {
  # Compare on the dates both series cover, as each starts at a different point
  common_dates <- complete.cases(rt_wide$inc2prev, rt_wide[[series_name]])
  tibble(series    = series_name,
         n         = sum(common_dates),
         corr      = cor(rt_wide$inc2prev[common_dates],
                         rt_wide[[series_name]][common_dates]),
         mean_diff = mean(rt_wide[[series_name]][common_dates] -
                            rt_wide$inc2prev[common_dates]),
         # Above 1 means wider swings than inc2prev
         sd_ratio  = sd(rt_wide[[series_name]][common_dates]) /
                       sd(rt_wide$inc2prev[common_dates]))
}) |> bind_rows()

cat("\n--- Rt agreement with inc2prev ---\n")
print_table(rt_agreement)
# Copying their GI reduces mean difference and sd_ratio compared to our original GI

dir.create(ch1_diag_config$plot_dir, recursive = TRUE, showWarnings = FALSE)

# Compare inc2prev Rt with the naive renewal estimate, under each generation interval
# Matching their interval accounts for nearly all the amplitude gap, leaving timing agreement unchanged
p_rt <- rt_compare |>
  filter(source != "fitted") |> # Model output rather than a check on the renewal implementation
  ggplot(aes(x = date, y = Rt, colour = source)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  # inc2prev in grey as the reference, with the two naive variants coloured
  scale_colour_manual(values = c("inc2prev"           = "grey40",
                                 "naive"              = "indianred",
                                 "naive, inc2prev GI" = "steelblue")) +
  labs(title = "Renewal Rt against inc2prev, under each generation interval",
       subtitle = "Naive (shifted gamma, mean 5.5, sd 2.1, max 21); matched (gamma, mean 3.64, sd 3.08, max 15)",
       x = "Date", y = expression(R[t]), colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "rt_validation.png"), p_rt,
       width = 10, height = 4, dpi = 300, bg = "white")


## Smooth sensitivity to k -----------------------------------------------------

# See whether s(t) saturates with increased flexibility
# At k=10, edf=8.7, and at k=40 this is 38.7, suggesting s(t) continues to absorb variability
# Deviance explained rises from 91.7% to 100.0%, the latter being a degenerate fit
# Good argument for excluding s(t) in the primary analysis - it takes over what covariates could explain

k_check <- lapply(ch1_diag_config$smooth_k, function(k) {
  fitted_model <- fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE,
                                  family = make_family(chosen_family),
                                  config = modifyList(ch1_gam_config,
                                                      list(smooth_k = k)))
  tibble(k = k,
         edf = sum(summary(fitted_model$fit)$edf),
         dev_expl = 100 * fitted_model$dev_expl)
}) |> bind_rows()

cat("\n--- s(t) with combined covariates, varying k ---\n")
print_table(k_check, 1)

## Plots -----------------------------------------------------------------------

dir.create(ch1_diag_config$plot_dir, recursive = TRUE, showWarnings = FALSE)

p_resid <- ggplot(residuals_table, aes(x = date, y = residual, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~model, ncol = 1) +
  labs(title = "Deviance residuals", x = NULL, y = "Residual") +
  theme_minimal() + theme(legend.position = "none")

p_period <- ggplot(residuals_table, aes(x = reorder(period, date), y = residual, fill = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_boxplot(outlier.size = 0.3) +
  labs(title = "Residuals by pandemic period", x = NULL, y = "Residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

p_acf <- ggplot(acfs, aes(x = lag, y = acf, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  labs(title = "Residual autocorrelation",
       subtitle = "One pooled fit across the whole study period",
       x = "Lag (days)", y = "ACF") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "baseline_diagnostics.png"),
       p_resid / p_period / p_acf, width = 10, height = 11, dpi = 300, bg = "white")

# Same layout as the pooled ACF above, one panel per window, with baseline drawn heavier
p_window_acf <- window_acf |>
  mutate(model = factor(model, levels = names(ch1_models))) |>
  ggplot(aes(x = lag, y = acf, colour = model, linewidth = model == "baseline")) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1), guide = "none") +
  facet_wrap(~origin, ncol = 2) +
  labs(title = "Residual autocorrelation, by training window",
       x = "Lag (days)", y = "Autocorrelation", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "window_autocorrelation.png"),
       p_window_acf, width = 10, height = 8, dpi = 300, bg = "white")

# Same residuals under each family, as a visual check on the family choice
# Free y scales, since Poisson residuals sit an order of magnitude wider
p_family <- lapply(names(fits), function(family_name) {
  tibble(date     = fits[[family_name]]$combined$model_data$date,
         residual = residuals(fits[[family_name]]$combined$fit, type = "deviance"),
         family   = family_name)
}) |> bind_rows() |>
  ggplot(aes(x = date, y = residual)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~family, ncol = 1, scales = "free_y") +
  labs(title = "Deviance residuals by observation family",
       subtitle = "Combined model, fitted once over the full study period",
       x = NULL, y = "Deviance residual") +
  theme_minimal()

ggsave(file.path(ch1_diag_config$plot_dir, "family_residuals.png"), p_family,
       width = 10, height = 5, dpi = 300, bg = "white")

## Save tables -----------------------------------------------------------------

write_csv(family_comparison,             file.path(ch1_diag_config$plot_dir, "table_family_comparison.csv"))
write_csv(covariate_estimates_by_family, file.path(ch1_diag_config$plot_dir, "table_covariates_by_family.csv"))
write_csv(window_dispersion,             file.path(ch1_diag_config$plot_dir, "table_window_dispersion.csv"))
write_csv(residuals_by_period,           file.path(ch1_diag_config$plot_dir, "table_residuals_by_period.csv"))
write_csv(pooled_acf,                    file.path(ch1_diag_config$plot_dir, "table_residual_acf.csv"))
write_csv(window_acf_summary,            file.path(ch1_diag_config$plot_dir, "table_window_acf.csv"))
write_csv(rt_agreement,                  file.path(ch1_diag_config$plot_dir, "table_rt_agreement.csv"))
write_csv(k_check,                       file.path(ch1_diag_config$plot_dir, "table_smooth_k_sensitivity.csv"))

message("Saved diagnostics to ", ch1_diag_config$plot_dir)
