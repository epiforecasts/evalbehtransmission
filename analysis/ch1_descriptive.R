# Chapter 1, issue #11: descriptive figures
# Fig 1.1 stacks Rt, the standardised covariates and stringency over the study window
# Rt is taken from inc2prev rather than re-estimated, keeping it independent of the model

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("R/inc2prev_path.R")
source("R/ch1_mobility_streams.R")
source("R/ch1_contact_covariate.R") # Names the contact panel after whichever series is used
source("R/ch1_period_scale.R") # Shared period colours, so shading matches the other figures

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
  mutate(period = factor(period, levels = ch1_period_levels))

period_layer <- function() {
  list(
    geom_rect(data = period_bands, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = period), alpha = 0.06),
    scale_fill_period(),
    scale_x_date(limits = window, expand = c(0, 0),
                 date_breaks = "1 month", date_labels = "%b %Y")
  )
}

## Fig 1.1 ---------------------------------------------------------------------

dir.create(ch1_desc_config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Bands are labelled on the top panel only, so the stack needs no period legend
period_labels <- period_bands |>
  mutate(midpoint = start + (end - start) / 2)

# Top panel: inc2prev Rt with its 90% credible interval as a reference series
p_rt <- inc2prev_national |>
  filter(name == "R") |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey20") +
  geom_line(aes(y = median), linewidth = 0.6) +
  geom_text(data = period_labels, inherit.aes = FALSE,
            aes(x = midpoint, y = Inf, label = period),
            angle = 90, hjust = 1.05, vjust = 0.4, size = 2.6, colour = "grey30") +
  labs(title = "Reproduction number, England: inc2prev posterior median and 90% credible interval",
       y = expression(R[t]), x = NULL, fill = NULL) +
  theme_minimal()

# Contacts and mobility are drawn separately, as each is a distinct behavioural stream
# Both are the model-ready covariates, z-scored over the study period
p_contacts <- covariates |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = contacts), colour = "firebrick", linewidth = 0.6, na.rm = TRUE) +
  labs(title = paste0(contact_covariate_label(), ", England"),
       y = "z-score", x = NULL, fill = NULL) +
  theme_minimal()

p_mobility_covariate <- covariates |>
  ggplot(aes(x = date)) +
  period_layer() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = mobility), colour = "steelblue", linewidth = 0.6, na.rm = TRUE) +
  labs(title = "Google Mobility composite, UK: trailing 7-day mean of four categories",
       y = "z-score", x = "Date", fill = NULL) +
  theme_minimal()

# Keeps the panels the same width so the x-axes align
fig_1_1 <- p_rt / p_contacts / p_mobility_covariate +
  plot_annotation(
    title = "Reproduction number and behavioural covariates",
    subtitle = paste("Covariates are z-scored over the study period, as entered in the model;",
                     "mobility is UK-wide, the other series are England")
  ) &
  guides(fill = "none") & # Drops the period key, as the bands are labelled directly
  theme(plot.title = element_text(size = 10))

ggsave(file.path(ch1_desc_config$output_dir, "fig_1_1_rt_covariates.png"),
       fig_1_1, width = 11, height = 9, dpi = 300, bg = "white")

## Descriptive methods figures -------------------------------------------------

# Prevalence is what CIS measures directly, and what inc2prev deconvolves into incidence
# Included to show the measured series behind the modelled outcome, not as a model input
# The q columns are on the percentage scale, unlike the per-capita mean and median
p_prevalence <- inc2prev_national |>
  filter(name == "est_prev") |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey40") +
  geom_line(aes(y = median), colour = "grey20", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Modelled SARS-CoV-2 PCR positivity, England",
       subtitle = "inc2prev posterior median and 90% credible interval, fitted to ONS CIS",
       x = "Date", y = "Testing positive (%)") +
  theme_minimal()

ggsave(file.path(ch1_desc_config$output_dir, "fig_cis_prevalence.png"),
       p_prevalence, width = 10, height = 4, dpi = 300, bg = "white")

# Infection incidence, the outcome the models are fitted to
# Only the median enters the models, so the interval here is context rather than an input
p_incidence <- inc2prev_national |>
  filter(name == "infections") |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = median), colour = "steelblue", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  labs(title = "Estimated daily SARS-CoV-2 infection incidence, England",
       subtitle = "inc2prev posterior median and 90% credible interval; the median is the model outcome",
       x = "Date", y = "Daily new infections") +
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
  labs(title = "Google Mobility categories, United Kingdom",
       subtitle = "Dashed lines mark the study window; percentages are relative to the Jan-Feb 2020 baseline",
       x = "Date", y = "% change from baseline", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(ch1_desc_config$output_dir, "fig_mobility_categories.png"),
       p_mobility, width = 10, height = 7, dpi = 300, bg = "white")

message("Saved descriptive figures to ", ch1_desc_config$output_dir)
