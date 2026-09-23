# Appendix A.1.4: what s(t) absorbs, and residual autocorrelation from the full-period fits
# The coefficient comparison reads data-processed/ch1_window_coefficients.csv, which
# ch1_rolling.R already writes with both smooth settings, so no window is refitted here
# Deviance explained is in the main body, so no table is repeated here

source("analysis/ch1_gam.R")

library(ggplot2)

## Config ----------------------------------------------------------------------

a1_4_config <- list(
  coef_path  = "data-processed/ch1_window_coefficients.csv",
  output_dir = "stashed/appendix/outputs",

  lag_max = 28
)

dir.create(a1_4_config$output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)

## Full-period combined fit, with and without s(t) -----------------------------
# Same k as everywhere else outside the rolling windows, from ch1_gam_config

combined_fits <- list(
  no_smooth = fit_renewal_gam(dat, ch1_models$combined, use_smooth = FALSE,
                              family = ch1_family()),
  smooth    = fit_renewal_gam(dat, ch1_models$combined, use_smooth = TRUE,
                              family = ch1_family())
)

# The two fits must cover the same days for the panel below to be a like-for-like comparison
stopifnot(identical(combined_fits$no_smooth$model_data$date,
                    combined_fits$smooth$model_data$date))

cat("Combined model fitted on", combined_fits$no_smooth$n_obs, "days |",
    "deviance explained", round(100 * combined_fits$no_smooth$dev_expl, 1),
    "% without s(t),", round(100 * combined_fits$smooth$dev_expl, 1), "% with\n")

# Rt is never observed, so incidence / Lambda_t is the reference: the renewal estimate taken
# straight from the data with no model in between, which is what a fit has to reproduce
naive_rt <- combined_fits$no_smooth$model_data |>
  transmute(date, Rt = incidence / Lambda_t, series = "Renewal estimate")

fitted_rt_frame <- bind_rows(
  naive_rt,
  combined_fits$no_smooth$fitted_rt |> mutate(series = "Without s(t)"),
  combined_fits$smooth$fitted_rt    |> mutate(series = "With s(t)")
) |>
  mutate(series = factor(series, levels = c("Renewal estimate", "Without s(t)",
                                            "With s(t)")))

# With s(t) the fit tracks the renewal estimate closely, leaving the covariates little to
# explain; without it the fitted Rt is a function of the two covariates alone
p_fit <- ggplot(fitted_rt_frame, aes(x = date, y = Rt, colour = series)) +
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
  # Horizontal major lines only, so the grid does not read as a box behind the series
  theme(legend.position   = "bottom",
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(a1_4_config$output_dir, "fig_a1_4a_smooth_absorption.png"),
       p_fit, width = 10, height = 4.5, dpi = 300, bg = "white")

## Coefficient attenuation across the rolling windows --------------------------
# Read rather than refitted, since ch1_rolling.R already fits every window both ways

window_coefficients <- read_csv(a1_4_config$coef_path, show_col_types = FALSE) |>
  filter(term %in% c("contacts", "mobility"), model == "combined") |>
  select(origin, term, used_smooth, estimate) |>
  tidyr::pivot_wider(names_from = used_smooth, values_from = estimate,
                     names_prefix = "smooth_") |>
  rename(without_smooth = smooth_FALSE, with_smooth = smooth_TRUE)

# One row per origin per term, or the pivot has silently collapsed something
stopifnot(nrow(window_coefficients) == 2 * n_distinct(window_coefficients$origin),
          !anyNA(window_coefficients))

smooth_attenuation <- window_coefficients |>
  group_by(term) |>
  summarise(n_origins             = n(),
            median_without_smooth = median(without_smooth),
            median_with_smooth    = median(with_smooth),
            n_shrunk              = sum(abs(with_smooth) < abs(without_smooth)),
            .groups = "drop") |>
  mutate(median_ratio = median_with_smooth / median_without_smooth)

cat("\n--- Combined-model coefficients across origins, with and without s(t) ---\n")
print(as.data.frame(smooth_attenuation |>
        mutate(across(where(is.double), \(x) round(x, 3)))), row.names = FALSE)

## Residual autocorrelation ----------------------------------------------------
# Full-period fits, matching the ACF computed in ch1_diagnostics.R, which also asserts
# that each model frame is contiguous before reading anything into the lags

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

# All four overlaid in one panel, as in ch1_diagnostics.R, so the models compare directly
# Baseline drawn heavier as the reference, matching outputs/ch1/window_autocorrelation.png
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

write_csv(smooth_attenuation,
          file.path(a1_4_config$output_dir, "table_a1_4_smooth_attenuation.csv"))

message("A.1.4 written to ", a1_4_config$output_dir)
