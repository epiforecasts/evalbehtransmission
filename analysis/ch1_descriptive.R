# Chapter 1, issue #11: descriptive figures.
# Fig 1.1 stacks Rt, the standardised covariates and stringency over the study
# window. Rt is taken from inc2prev rather than re-estimated, so these stay
# independent of the model.

source("analysis/ch1_gam.R")

library(ggplot2)
library(patchwork)

## Config ----------------------------------------------------------------------

ch1_desc_config <- list(
  covariates_path = "data-processed/ch1_covariates.csv",
  periods_path    = "data-processed/ch1_periods.csv",
  data_path       = "data-processed/ch1_data.csv",
  output_dir      = "outputs/ch1",

  # Raw file, so the mobility panel can show the collapse before the study window
  mobility_path     = "data-processed/google_mobility_UK.csv",

  # Mirrors ch1_covariates.R, which this script does not source
  mobility_retained = c("retail_and_recreation", "grocery_and_pharmacy",
                        "transit_stations", "workplaces")
)

## Load ------------------------------------------------------------------------

covariates <- read_csv(ch1_desc_config$covariates_path, show_col_types = FALSE)
periods    <- read_csv(ch1_desc_config$periods_path, show_col_types = FALSE)
daily      <- read_csv(ch1_desc_config$data_path, show_col_types = FALSE)

window <- range(covariates$date)

inc2prev_rt <- read_csv(inc2prev_path("outputs/estimates_national.csv"),
                        show_col_types = FALSE) |>
  filter(variable == "England", name == "R") |>
  select(date, Rt = q50, Rt_lower = q5, Rt_upper = q95) |>
  filter(date >= window[1], date <= window[2])

## Period bands ----------------------------------------------------------------
# Shared across the stacked panels so the same shading lines up.

bands <- periods |>
  group_by(period) |>
  summarise(start = min(date), end = max(date), .groups = "drop") |>
  arrange(start) |>
  mutate(period = factor(period, levels = period))

period_layer <- function() {
  list(
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = period), alpha = 0.10),
    scale_x_date(limits = window, expand = c(0, 0))
  )
}

## Fig 1.1 ---------------------------------------------------------------------

p_rt <- ggplot(inc2prev_rt, aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), alpha = 0.3, fill = "grey20") +
  geom_line(aes(y = Rt), linewidth = 0.6) +
  labs(y = expression(R[t]), x = NULL, fill = NULL) +
  theme_minimal() + theme(legend.position = "none")

p_covariates <- covariates |>
  select(date, contacts, mobility) |>
  tidyr::pivot_longer(-date, names_to = "covariate") |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = value, colour = covariate), linewidth = 0.6, na.rm = TRUE) +
  scale_colour_manual(values = c(contacts = "firebrick", mobility = "steelblue")) +
  labs(y = "z-score", x = NULL, colour = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(fill = "none")

# Targeted stretches drawn separately: OxCGRT records the strictest policy in
# the jurisdiction, so stringency there is an upper bound rather than an average.
p_stringency <- ggplot(periods, aes(x = date)) +
  period_layer() +
  geom_line(aes(y = stringency), linewidth = 0.6, colour = "grey20") +
  geom_point(data = periods |> filter(stay_at_home_flag == 0),
             aes(y = stringency), colour = "darkorange", size = 0.6) +
  labs(y = "Stringency", x = "Date", fill = NULL,
       caption = "Orange marks days when stay-at-home applied to part of England only") +
  theme_minimal() +
  theme(legend.position = "bottom")

fig_1_1 <- p_rt / p_covariates / p_stringency +
  plot_annotation(
    title = "Rt, behavioural covariates and policy stringency, England",
    subtitle = "Shading marks pandemic periods; Rt is taken directly from inc2prev"
  )

ggsave(file.path(ch1_desc_config$output_dir, "fig_1_1_rt_covariates_stringency.png"),
       fig_1_1, width = 11, height = 9, dpi = 300, bg = "white")

## Descriptive methods figures -------------------------------------------------

cis_prevalence <- read_csv(inc2prev_path("outputs/estimates_national.csv"),
                           show_col_types = FALSE) |>
  filter(variable == "England", name == "est_prev") |>
  select(date, prevalence = q50, lower = q5, upper = q95) |>
  filter(date >= window[1], date <= window[2])

p_prevalence <- ggplot(cis_prevalence, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = prevalence), colour = "steelblue", linewidth = 0.6) +
  labs(title = "ONS CIS prevalence, England",
       x = "Date", y = "Proportion PCR positive") +
  theme_minimal()

ggsave(file.path(ch1_desc_config$output_dir, "fig_cis_prevalence.png"),
       p_prevalence, width = 10, height = 4, dpi = 300, bg = "white")

# All six categories over the full series, so the March 2020 collapse is visible
# even though it precedes the study window.
p_mobility <- read_csv(ch1_desc_config$mobility_path, show_col_types = FALSE) |>
  select(date, ends_with("_percent_change_from_baseline")) |>
  tidyr::pivot_longer(-date, names_to = "category") |>
  mutate(category = sub("_percent_change_from_baseline", "", category),
         retained = if_else(category %in% ch1_desc_config$mobility_retained,
                            "retained", "excluded")) |>
  ggplot(aes(x = date, y = value, colour = retained)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = window, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.4) +
  facet_wrap(~category, ncol = 2) +
  scale_colour_manual(values = c(retained = "steelblue", excluded = "grey60")) +
  labs(title = "Google Mobility, six categories",
       subtitle = "Blue enter the composite, grey are excluded; dashed lines mark the study window",
       x = "Date", y = "% change from baseline", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_desc_config$output_dir, "fig_mobility_six_categories.png"),
       p_mobility, width = 10, height = 7, dpi = 300, bg = "white")

message("Saved descriptive figures to ", ch1_desc_config$output_dir)
