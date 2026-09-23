# Appendix A.1.1: residuals from the incidence-only baseline under each observation family
# Same figure as outputs/ch1/family_residuals.png, which covers the combined model, but
# fitted on the baseline, which is what the A.1.1 text refers to
# Dispersion ratios and AIC are already in outputs/ch1/table_family_comparison.csv

source("analysis/ch1_gam.R")

library(ggplot2)

## Config ----------------------------------------------------------------------

a1_1_config <- list(
  output_dir = "stashed/appendix/outputs",

  # Existing non-overlapping check, used below to verify the per-origin dispersion
  window_dispersion_path = "outputs/ch1/table_window_dispersion.csv",

  # Matches ch1_rolling.R, so dispersion is measured on the windows actually forecast from
  window_weeks = 8,
  first_origin = as.Date("2020-07-01"),
  last_origin  = as.Date("2021-01-01"),
  step_days    = 7
)

dir.create(a1_1_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# As in ch1_diagnostics.R, so a family is named once and constructed where it is used
make_family <- function(family_name) {
  if (family_name == "poisson") poisson(link = "log") else ch1_family()
}

## Fit the baseline under each family ------------------------------------------

dat <- read_csv(ch1_gam_config$input_path, show_col_types = FALSE)

# Offset-only, so the only difference between the two fits is the variance assumption
baseline <- lapply(c(poisson = "poisson", nb = "nb"), function(family_name) {
  fit_renewal_gam(dat, ch1_models$baseline, family = make_family(family_name))
})

# Both families must use the same rows, or the two panels are not comparable
stopifnot(identical(baseline$poisson$model_data$date, baseline$nb$model_data$date))

cat("Baseline fitted on", nrow(baseline$poisson$model_data), "days,",
    format(min(baseline$poisson$model_data$date)), "to",
    format(max(baseline$poisson$model_data$date)), "\n")

## Residuals -------------------------------------------------------------------

residual_frame <- lapply(names(baseline), function(family_name) {
  tibble(date     = baseline[[family_name]]$model_data$date,
         residual = residuals(baseline[[family_name]]$fit, type = "deviance"),
         family   = family_name)
}) |> bind_rows() |>
  # Poisson first, so the panel showing the problem precedes the one showing the fix
  mutate(family = factor(family, levels = c("poisson", "nb"),
                         labels = c("Poisson", "Negative binomial")))

residual_range <- residual_frame |>
  group_by(family) |>
  summarise(min_residual = min(residual), max_residual = max(residual), .groups = "drop")

cat("\n--- Deviance residual range by family ---\n")
print(as.data.frame(residual_range |>
        mutate(across(where(is.numeric), \(x) round(x, 1)))), row.names = FALSE)

# Free y, since the Poisson residuals sit an order of magnitude wider
p_family <- ggplot(residual_frame, aes(x = date, y = residual)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~family, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Deviance residuals under Poisson and negative binomial observation models",
       subtitle = "Incidence-only baseline (intercept and log(Λt) offset), full study period",
       x = NULL, y = "Deviance residual") +
  theme_minimal()

ggsave(file.path(a1_1_config$output_dir, "fig_a1_1_family_residuals.png"), p_family,
       width = 10, height = 5, dpi = 300, bg = "white")

## Dispersion by training window -----------------------------------------------
# The fit above spans the whole period, so it absorbs non-stationarity no training window
# ever sees. Refitting the baseline on each rolling window shows the same result holds

# Pearson chi-square over residual degrees of freedom, as in ch1_diagnostics.R
# One indicates the family's variance assumption matches the data
dispersion <- function(fitted_model) {
  sum(residuals(fitted_model$fit, type = "pearson")^2) / df.residual(fitted_model$fit)
}

origins <- seq(a1_1_config$first_origin, a1_1_config$last_origin,
               by = a1_1_config$step_days)

window_dispersion <- lapply(origins, function(origin) {
  window_start <- origin - a1_1_config$window_weeks * 7 + 1
  window_fits <- lapply(c(poisson = "poisson", nb = "nb"), function(family_name) {
    fit_renewal_gam(dat, ch1_models$baseline, family = make_family(family_name),
                    fit_from = window_start, fit_to = origin)
  })
  tibble(origin  = origin,
         poisson = dispersion(window_fits$poisson),
         nb      = dispersion(window_fits$nb))
}) |> bind_rows()

# The four non-overlapping origins ch1_diagnostics.R already reports must come out the same,
# or these windows are not the ones the existing table describes
existing_dispersion <- read_csv(a1_1_config$window_dispersion_path,
                                show_col_types = FALSE) |> arrange(origin)
refitted_check <- window_dispersion |>
  filter(origin %in% existing_dispersion$origin) |> arrange(origin)
stopifnot(nrow(refitted_check) == nrow(existing_dispersion),
          max(abs(refitted_check$poisson - existing_dispersion$poisson)) < 1e-6,
          max(abs(refitted_check$nb - existing_dispersion$nb)) < 1e-6)

cat("Poisson dispersion across", nrow(window_dispersion), "windows:",
    round(min(window_dispersion$poisson)), "to", round(max(window_dispersion$poisson)),
    "| negative binomial:", round(min(window_dispersion$nb), 2), "to",
    round(max(window_dispersion$nb), 2), "\n")

# Log scale, as the Poisson values span two orders of magnitude and the negative binomial
# values sit at one. The line at one is where the variance assumption would be correct
p_dispersion <- window_dispersion |>
  tidyr::pivot_longer(-origin, names_to = "family", values_to = "dispersion") |>
  mutate(family = factor(family, levels = c("poisson", "nb"),
                         labels = c("Poisson", "Negative binomial"))) |>
  ggplot(aes(x = origin, y = dispersion, colour = family)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.6) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_colour_manual(values = c(Poisson = "firebrick",
                                 `Negative binomial` = "steelblue")) +
  labs(title = "Pearson dispersion of the incidence-only baseline at each forecast origin",
       subtitle = sprintf("Intercept and log(Λt) offset only, refitted on each %d-week training window",
                          a1_1_config$window_weeks),
       x = "Forecast origin", y = "Pearson dispersion", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(a1_1_config$output_dir, "fig_a1_1_window_dispersion.png"), p_dispersion,
       width = 10, height = 4.5, dpi = 300, bg = "white")

message("A.1.1 written to ", a1_1_config$output_dir)
