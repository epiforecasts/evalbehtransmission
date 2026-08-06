# Chapter 1, issue #42: baseline diagnostics and choice of observation family.

source("analysis/ch1_gam.R")

library(ggplot2)
library(patchwork)

## Config ----------------------------------------------------------------------

ch1_diag_config <- list(
  periods_path = "data-processed/ch1_periods.csv",
  plot_dir     = "outputs/ch1",
  smooth_k     = c(5, 10, 20, 40)
)

# Chosen from the family comparison below, and used for everything after it.
chosen_family <- "nb"   # ch1_family in ch1_gam.R must match

make_family <- function(family_name) {
  if (family_name == "poisson") poisson(link = "log") else nb()
}

print_table <- function(x, digits = 3) {
  print(as.data.frame(mutate(x, across(where(is.numeric), \(v) round(v, digits)))),
        row.names = FALSE)
}

## Fit every model under both families, once -----------------------------------

dat     <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
periods <- read_csv(ch1_diag_config$periods_path, show_col_types = FALSE) |>
  select(date, period)

fits <- lapply(c(poisson = "poisson", nb = "nb"), function(family_name) {
  lapply(ch1_models, \(covariates) {
    fit_renewal_gam(dat, covariates, family = make_family(family_name))
  })
})

## Family choice ---------------------------------------------------------------

# Pearson dispersion > 1 means more variability than the family allows.
dispersion <- function(fitted_model) {
  sum(residuals(fitted_model$fit, type = "pearson")^2) / df.residual(fitted_model$fit)
}

cat("\n--- Baseline under each family ---\n")
print_table(tibble(
  family     = c("poisson", "nb"),
  dispersion = c(dispersion(fits$poisson$baseline), dispersion(fits$nb$baseline)),
  aic        = c(AIC(fits$poisson$baseline$fit), AIC(fits$nb$baseline$fit)),
  intercept  = c(fits$poisson$baseline$coefficients$estimate,
                 fits$nb$baseline$coefficients$estimate),
  se         = c(fits$poisson$baseline$coefficients$se,
                 fits$nb$baseline$coefficients$se)
), 4)

# The practical consequence for the covariate coefficients.
cat("\n--- Covariate estimates and SEs by family ---\n")
print_table(lapply(names(fits), function(family_name) {
  lapply(setdiff(names(ch1_models), "baseline"), function(model_name) {
    fits[[family_name]][[model_name]]$coefficients |>
      filter(term != "(Intercept)") |>
      mutate(model = model_name, family = family_name)
  }) |> bind_rows()
}) |> bind_rows() |>
  select(model, term, family, estimate, se) |>
  tidyr::pivot_wider(names_from = family, values_from = c(estimate, se)) |>
  mutate(se_ratio = se_nb / se_poisson), 4)

## Residuals by period ---------------------------------------------------------

residuals <- lapply(names(ch1_models), function(model_name) {
  fitted_model <- fits[[chosen_family]][[model_name]]
  tibble(date     = fitted_model$model_data$date,
         residual = residuals(fitted_model$fit, type = "deviance"),
         model    = model_name) |>
    left_join(periods, by = "date")
}) |> bind_rows() |>
  mutate(model = factor(model, levels = names(ch1_models)))

cat("\n--- Mean deviance residual by period ---\n")
print_table(residuals |> group_by(period, model) |>
              summarise(mean_resid = mean(residual), .groups = "drop") |>
              tidyr::pivot_wider(names_from = model, values_from = mean_resid), 2)

## Autocorrelation -------------------------------------------------------------
# Correlated residuals break the independence assumption behind the SEs. The
# question is whether covariates reduce it, not whether the baseline has it.

acfs <- lapply(names(ch1_models), function(model_name) {
  acf_out <- acf(residuals(fits[[chosen_family]][[model_name]]$fit, type = "deviance"),
                 lag.max = 28, plot = FALSE)
  tibble(lag = as.numeric(acf_out$lag), acf = as.numeric(acf_out$acf),
         model = model_name)
}) |> bind_rows()

cat("\n--- Residual autocorrelation ---\n")
print_table(acfs |> filter(lag %in% c(1, 7, 14)) |>
              tidyr::pivot_wider(names_from = lag, values_from = acf,
                                 names_prefix = "lag_"))

## Validation against inc2prev Rt ----------------------------------------------
# inc2prev derives its own Rt from the same infections series, Cori et al. style
# with a discretised gamma generation time.

combined_fit <- fits[[chosen_family]]$combined

rt_compare <- bind_rows(
  read_csv(inc2prev_path("outputs/estimates_national.csv"), show_col_types = FALSE) |>
    filter(variable == "England", name == "R") |>
    transmute(date, Rt = q50, source = "inc2prev"),
  combined_fit$model_data |>
    transmute(date, Rt = incidence / Lambda_t, source = "naive"),
  combined_fit$fitted_rt |> mutate(source = "fitted")
) |>
  filter(date >= min(combined_fit$model_data$date),
         date <= max(combined_fit$model_data$date))

rt_wide <- rt_compare |> tidyr::pivot_wider(names_from = source, values_from = Rt)

cat("\n--- Rt agreement with inc2prev ---\n")
print_table(lapply(c("naive", "fitted"), function(series_name) {
  # inc2prev R starts after its infections series, and the fitted Rt only spans
  # the model rows, so compare on the dates all three cover.
  common_dates <- complete.cases(rt_wide$inc2prev, rt_wide[[series_name]])
  tibble(series    = series_name,
         n         = sum(common_dates),
         corr      = cor(rt_wide$inc2prev[common_dates],
                         rt_wide[[series_name]][common_dates]),
         mean_diff = mean(rt_wide[[series_name]][common_dates] -
                            rt_wide$inc2prev[common_dates]))
}) |> bind_rows())

## Smooth sensitivity to k -----------------------------------------------------
# Whether s(t) saturates regardless of how flexible it is allowed to be.

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

p_resid <- ggplot(residuals, aes(x = date, y = residual, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~model, ncol = 1) +
  labs(title = "Deviance residuals", x = NULL, y = "Residual") +
  theme_minimal() + theme(legend.position = "none")

p_period <- ggplot(residuals, aes(x = reorder(period, date), y = residual, fill = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_boxplot(outlier.size = 0.3) +
  labs(title = "Residuals by pandemic period", x = NULL, y = "Residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

p_acf <- ggplot(acfs, aes(x = lag, y = acf, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  labs(title = "Residual autocorrelation", x = "Lag (days)", y = "ACF") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "baseline_diagnostics.png"),
       p_resid / p_period / p_acf, width = 10, height = 11, dpi = 300, bg = "white")

# Only the naive ratio and inc2prev: both apply the renewal equation to the same
# infections series, so any gap is generation interval and smoothing. A fitted
# Rt would be a model result, not a check on the renewal implementation.
p_rt <- rt_compare |>
  filter(source %in% c("naive", "inc2prev")) |>
  ggplot(aes(x = date, y = Rt, colour = source)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  labs(title = "Renewal Rt against inc2prev",
       x = "Date", y = expression(R[t]), colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "rt_validation.png"), p_rt,
       width = 10, height = 4, dpi = 300, bg = "white")

message("Saved diagnostics to ", ch1_diag_config$plot_dir)
