# Appendix A.1.4: in-sample fit with and without s(t), and residual autocorrelation from the full-period fits
# This script refits models, as the Chapter 1 pipeline saves no fitted R(t) and ACF at three lags only

source("analysis/ch1_gam.R")

library(ggplot2)

## Config ----------------------------------------------------------------------

a1_4_config <- list(
  output_dir = "analysis/appendix/outputs",

  lag_max = 28
)

dir.create(a1_4_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Daily incidence and covariates from the Chapter 1 pipeline
dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)

## Full-period combined fit, with and without s(t) -----------------------------
# Use same k as the other full-period fits (20)

combined_fits <- list(
  no_smooth = fit_renewal_gam(dat, ch1_models$combined, use_smooth = FALSE,
                              family = ch1_family()),
  smooth    = fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE,
                              family = ch1_family())
)

# Both fits must use the same days for comparability
stopifnot(identical(combined_fits$no_smooth$model_data$date,
                    combined_fits$smooth$model_data$date))

cat("Combined model fitted on", combined_fits$no_smooth$n_obs, "days |",
    "deviance explained", round(100 * combined_fits$no_smooth$dev_expl, 1),
    "% without s(t),", round(100 * combined_fits$smooth$dev_expl, 1), "% with\n")

# R(t) is never observed, so the naive renewal estimate (incidence / Λt) is the reference
naive_rt <- combined_fits$no_smooth$model_data |>
  transmute(date, Rt = incidence / Lambda_t, series = "Renewal estimate")

# Dataframe of naive and fitted R(t) by date, with and without s(t)
fitted_rt_frame <- bind_rows(
  naive_rt,
  combined_fits$no_smooth$fitted_rt |> mutate(series = "Without s(t)"),
  combined_fits$smooth$fitted_rt    |> mutate(series = "With s(t)")
) |>
  mutate(series = factor(series, levels = c("Renewal estimate", "Without s(t)",
                                            "With s(t)")))

p_fitted_rt <- ggplot(fitted_rt_frame, aes(x = date, y = Rt, colour = series)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  scale_colour_manual(values = c("Renewal estimate" = "grey60",
                                 "Without s(t)"     = "steelblue",
                                 "With s(t)"        = "firebrick")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "In-sample fit of the combined model, with and without s(t)",
       subtitle = sprintf("One fit over the whole study period, s(t) with k = %d",
                          ch1_gam_config$smooth_k),
       x = NULL, y = "R(t)", colour = NULL) +
  theme_minimal() +
  theme(legend.position   = "bottom",
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(a1_4_config$output_dir, "fig_a1_4a_smooth_absorption.png"),
       p_fitted_rt, width = 10, height = 4.5, dpi = 300, bg = "white")

## Residual autocorrelation ----------------------------------------------------
# Full-period fits, matching the ACF computed in ch1_diagnostics.R
# Assumes days are equally spaced, which is checked in analysis/ch1_diagnostics.R

# Dataframe of residual autocorrelation by lag, for each model
acf_frame <- lapply(names(ch1_models), function(model_name) {
  model_fit <- fit_renewal_gam(dat, ch1_models[[model_name]], family = ch1_family())
  acf_out   <- acf(residuals(model_fit$fit, type = "deviance"),
                   lag.max = a1_4_config$lag_max, plot = FALSE)
  tibble(model = model_name,
         lag   = as.numeric(acf_out$lag),
         acf   = as.numeric(acf_out$acf))
}) |> bind_rows() |>
  mutate(model = factor(model, levels = names(ch1_models)))

cat("\n--- Residual autocorrelation, full-period fits ---\n")
print(as.data.frame(acf_frame |>
        filter(lag %in% c(1, 7, 14, 21, 28)) |>
        select(model, lag, acf) |>
        tidyr::pivot_wider(names_from = lag, values_from = acf, names_prefix = "lag_") |>
        mutate(across(where(is.double), \(x) round(x, 3)))), row.names = FALSE)

# All four models in one panel, baseline drawn heavier as the reference
p_acf <- ggplot(acf_frame, aes(x = lag, y = acf, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(aes(linewidth = model == "baseline")) +
  scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1), guide = "none") +
  labs(title = "Residual autocorrelation",
       subtitle = "Deviance residuals, one pooled fit across the whole study period",
       x = "Lag (days)", y = "Autocorrelation", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(a1_4_config$output_dir, "fig_a1_4b_residual_acf.png"),
       p_acf, width = 10, height = 4.5, dpi = 300, bg = "white")

message("A.1.4 written to ", a1_4_config$output_dir)
