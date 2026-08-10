# Chapter 1, issue #11: descriptive figures
# Fig 1.1 stacks Rt, the standardised covariates and stringency over the study window
# Rt is taken from inc2prev rather than re-estimated, keeping it independent of the model

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("R/inc2prev_path.R")
source("R/ch1_mobility_streams.R")

## Config ----------------------------------------------------------------------

ch1_desc_config <- list(
  covariates_path = "data-processed/ch1_covariates.csv",
  periods_path    = "data-processed/ch1_periods.csv",
  output_dir      = "outputs/ch1",

  # Raw file, so the mobility panel can show the collapse before the study window
  mobility_path   = "data-processed/google_mobility_UK.csv",

  # Marks days when stay-at-home applied to part of England only
  show_targeted   = TRUE
)

## Load ------------------------------------------------------------------------

covariates <- read_csv(ch1_desc_config$covariates_path, show_col_types = FALSE)
periods    <- read_csv(ch1_desc_config$periods_path, show_col_types = FALSE)

window <- range(covariates$date)

# Read once, then filter by name for each series used below
inc2prev_national <- read_csv(inc2prev_path("outputs/estimates_national.csv"),
                              show_col_types = FALSE) |>
  filter(variable == "England", date >= window[1], date <= window[2]) |>
  select(date, name, median = q50, lower = q5, upper = q95)

## Period bands ----------------------------------------------------------------
# Shared across the stacked panels so the same shading lines up.

period_bands <- periods |>
  group_by(period) |>
  summarise(start = min(date), end = max(date), .groups = "drop") |>
  arrange(start) |>
  mutate(period = factor(period, levels = period))

period_layer <- function() {
  list(
    geom_rect(data = period_bands, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = period), alpha = 0.10),
    scale_x_date(limits = window, expand = c(0, 0))
  )
}

## Fig 1.1 ---------------------------------------------------------------------

dir.create(ch1_desc_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Top panel: inc2prev Rt with its 90% credible interval as a reference series
p_rt <- inc2prev_national |>
  filter(name == "R") |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey20") +
  geom_line(aes(y = median), linewidth = 0.6) +
  labs(y = expression(R[t]), x = NULL, fill = NULL) +
  theme_minimal()

# Middle panel: z-scored covariates from CoMix and Google Mobility used in models
p_covariates <- covariates |>
  select(date, contacts, mobility) |>
  tidyr::pivot_longer(-date, names_to = "covariate") |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = value, colour = covariate), linewidth = 0.6, na.rm = TRUE) +
  scale_colour_manual(values = c(contacts = "firebrick", mobility = "steelblue")) +
  labs(y = "z-score", x = NULL, colour = NULL, fill = NULL) +
  theme_minimal()

# Bottom panel: OxCGRT Stringency Index
p_stringency <- ggplot(periods, aes(x = date)) +
  period_layer() +
  geom_line(aes(y = stringency), linewidth = 0.6, colour = "grey20") +
  labs(y = "Stringency Index", x = "Date", fill = NULL) +
  theme_minimal()


# Optional overlay if showing targeted interventions such as stay-at-home being targeted, not national
# Otherwise, OxCGRT Stringency Index records strictest policy, giving above-average impression
if (ch1_desc_config$show_targeted) {
  p_stringency <- p_stringency +
    geom_point(data = periods |> filter(stay_at_home_flag == 0),
               aes(y = stringency), colour = "darkorange", size = 0.6) +
    labs(caption = "Orange marks days when stay-at-home applied to part of England only")
}

# Keeps the panels the same width so the x-axes align
fig_1_1 <- p_rt / p_covariates / p_stringency +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Rt, behavioural covariates and OxCGRT Stringency Index, England",
    subtitle = "Shading marks pandemic periods; Rt is taken directly from inc2prev"
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(ch1_desc_config$output_dir, "fig_1_1_rt_covariates_stringency.png"),
       fig_1_1, width = 11, height = 9, dpi = 300, bg = "white")

## Descriptive methods figures -------------------------------------------------

# Prevalence is what CIS measures directly, and what inc2prev deconvolves into incidence
# Included to show the measured series behind the modelled outcome, not as a model input
p_prevalence <- inc2prev_national |>
  filter(name == "est_prev") |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey40") +
  geom_line(aes(y = median), colour = "grey20", linewidth = 0.6) +
  labs(title = "ONS CIS prevalence, England",
       subtitle = "inc2prev median with 90% credible interval",
       x = "Date", y = "Proportion PCR positive") +
  theme_minimal()

ggsave(file.path(ch1_desc_config$output_dir, "fig_cis_prevalence.png"),
       p_prevalence, width = 10, height = 4, dpi = 300, bg = "white")

# Infection incidence, the outcome the models are fitted to
p_incidence <- inc2prev_national |>
  filter(name == "infections") |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = median), colour = "steelblue", linewidth = 0.6) +
  labs(title = "Daily infection incidence, England",
       subtitle = "inc2prev median with 90% credible interval",
       x = "Date", y = "Daily infections") +
  theme_minimal()

ggsave(file.path(ch1_desc_config$output_dir, "fig_incidence.png"),
       p_incidence, width = 10, height = 4, dpi = 300, bg = "white")

# All six categories over the full series, so the March 2020 decline is visible
# even though it precedes the study window
# Shows which categories enter the composite and which are dropped
p_mobility <- read_csv(ch1_desc_config$mobility_path, show_col_types = FALSE) |>
  select(date, ends_with("_percent_change_from_baseline")) |>
  tidyr::pivot_longer(-date, names_to = "category") |>
  mutate(category = sub("_percent_change_from_baseline", "", category),
         retained = if_else(category %in% mobility_short_names(),
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
