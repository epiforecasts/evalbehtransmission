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

# Chosen from the baseline family comparison below, and used for every check
# after it.
ch1_family <- nb()

## Fit baseline under each family ----------------------------------------------

dat     <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)
periods <- read_csv(ch1_diag_config$periods_path, show_col_types = FALSE) |>
  select(date, period)

fit_poisson <- fit_renewal_gam(dat, family = poisson(link = "log"))
fit_nb      <- fit_renewal_gam(dat, family = nb())

## Dispersion ------------------------------------------------------------------

# Pearson dispersion > 1 indicates more variability than the family allows.
dispersion <- function(m) {
  r <- residuals(m$fit, type = "pearson")
  sum(r^2) / df.residual(m$fit)
}

cat("\n--- Observation family ---\n")
cat("Poisson: dispersion =", round(dispersion(fit_poisson), 1),
    "| AIC =", round(AIC(fit_poisson$fit)), "\n")
cat("NB:      dispersion =", round(dispersion(fit_nb), 1),
    "| AIC =", round(AIC(fit_nb$fit)),
    "| theta =", round(fit_nb$fit$family$getTheta(TRUE), 2), "\n")

cat("\nIntercept (log R0) and its SE under each family:\n")
print(as.data.frame(bind_rows(
  fit_poisson$coefficients |> mutate(family = "poisson"),
  fit_nb$coefficients      |> mutate(family = "nb")
) |> select(family, estimate, se, lower, upper) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))), row.names = FALSE)

## Residuals -------------------------------------------------------------------

residual_frame <- function(m, label) {
  tibble(
    date     = m$model_data$date,
    residual = residuals(m$fit, type = "deviance"),
    model  = label
  ) |>
    left_join(periods, by = "date")
}

resids <- lapply(names(ch1_models), function(v) {
  residual_frame(fit_renewal_gam(dat, ch1_models[[v]], family = ch1_family), v)
}) |> bind_rows() |>
  mutate(model = factor(model, levels = names(ch1_models)))

cat("\n--- Mean deviance residual by period and model ---\n")
print(as.data.frame(resids |>
  group_by(period, model) |>
  summarise(mean_resid = round(mean(residual), 2),
            sd_resid   = round(sd(residual), 2), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = c(mean_resid, sd_resid))),
  row.names = FALSE)

## Autocorrelation -------------------------------------------------------------
# Correlated residuals break the independence assumption behind the standard
# errors. Run across models: the question is whether covariates reduce it,
# not whether the baseline has it, which is near-guaranteed under constant Rt.

acf_frame <- function(m, label, lag_max = 28) {
  a <- acf(residuals(m$fit, type = "deviance"), lag.max = lag_max, plot = FALSE)
  tibble(lag = as.numeric(a$lag), acf = as.numeric(a$acf), model = label)
}

acfs <- lapply(names(ch1_models), function(v) {
  acf_frame(fit_renewal_gam(dat, ch1_models[[v]], family = nb()), v)
}) |> bind_rows()

cat("\n--- Residual autocorrelation by model (NB) ---\n")
print(as.data.frame(acfs |> filter(lag %in% c(1, 7, 14)) |>
  mutate(acf = round(acf, 3)) |>
  tidyr::pivot_wider(names_from = lag, values_from = acf,
                     names_prefix = "lag_")), row.names = FALSE)

## Smooth sensitivity to k -----------------------------------------------------
# Whether s(t) saturates regardless of how flexible it is allowed to be.

k_check <- lapply(ch1_diag_config$smooth_k, function(k) {
  cfg <- modifyList(ch1_gam_config, list(smooth_k = k))
  m   <- fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE, config = cfg)
  tibble(k = k, edf = sum(summary(m$fit)$edf), dev_expl = m$dev_expl)
}) |> bind_rows()

cat("\n--- s(t) with combined covariates, varying k ---\n")
print(as.data.frame(k_check |>
  mutate(dev_expl = round(100 * dev_expl, 2), edf = round(edf, 1))), row.names = FALSE)

## Plots -----------------------------------------------------------------------

dir.create(ch1_diag_config$plot_dir, recursive = TRUE, showWarnings = FALSE)

p_resid <- ggplot(resids, aes(x = date, y = residual, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~model, ncol = 1) +
  labs(title = "Deviance residuals by model", x = NULL, y = "Residual") +
  theme_minimal() + theme(legend.position = "none")

p_period <- resids |>
  ggplot(aes(x = reorder(period, date), y = residual, fill = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_boxplot(outlier.size = 0.3) +
  labs(title = "Residuals by pandemic period", x = NULL, y = "Residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

p_acf <- ggplot(acfs, aes(x = lag, y = acf, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  labs(title = "Residual autocorrelation by model", x = "Lag (days)", y = "ACF") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_diag_config$plot_dir, "baseline_diagnostics.png"),
       p_resid / p_period / p_acf, width = 10, height = 11, dpi = 300, bg = "white")

message("Saved diagnostics to ", ch1_diag_config$plot_dir)
