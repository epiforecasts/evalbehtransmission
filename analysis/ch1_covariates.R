# Chapter 1, issue #38: build model-ready covariates from the analysis dataset.

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

## Config ----------------------------------------------------------------------

ch1_cov_config <- list(
  contact_covariate = "eigenvalue",   # or "mean_contacts"

  # Equal-weighted composite; parks and residential excluded
  mobility_streams = c("retail_recreation", "grocery_pharmacy",
                       "transit", "workplaces"),

  # Trailing window: t-6..t, so no future information at the forecast origin
  smooth_window = 7,

  lag = 0,

  input_path  = "data-processed/ch1_data.csv",
  output_path = "data-processed/ch1_covariates.csv",
  plot_dir    = "outputs/ch1"
)

## Transformations -------------------------------------------------------------

select_contact_covariate <- function(dat, which) {
  switch(which,
    eigenvalue    = dat$comix_eigen,
    mean_contacts = stop("mean_contacts has not been computed yet"),
    stop("Unknown contact covariate: ", which)
  )
}

# Trailing mean over t-k+1..t. sides = 1 uses only current and past values, so
# no future information enters at the forecast origin. Any 7 consecutive days
# contains each weekday once, so a 7-day window removes day-of-week exactly.
# Note dplyr masks stats::filter, hence the explicit namespace.
trailing_mean <- function(x, k) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
}

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

apply_lag <- function(x, l) if (l == 0) x else dplyr::lag(x, l)

## Build -----------------------------------------------------------------------

build_covariates <- function(config = ch1_cov_config) {
  dat <- read_csv(config$input_path, show_col_types = FALSE)

  # Averaging and smoothing are both linear, so combining first is equivalent
  # to smoothing each stream separately.
  composite <- rowMeans(dat[, config$mobility_streams])

  dat |>
    transmute(
      date,
      incidence,
      contacts_raw    = select_contact_covariate(dat, config$contact_covariate),
      mobility_raw    = composite,
      mobility_smooth = trailing_mean(composite, config$smooth_window)
    ) |>
    mutate(
      contacts = apply_lag(zscore(contacts_raw), config$lag),
      mobility = apply_lag(zscore(mobility_smooth), config$lag)
    )
}

## Plots -----------------------------------------------------------------------

plot_covariates <- function(cov, config = ch1_cov_config) {
  dir.create(config$plot_dir, recursive = TRUE, showWarnings = FALSE)

  p_mobility <- ggplot(cov, aes(x = date)) +
    geom_line(aes(y = mobility_raw), colour = "grey70", linewidth = 0.3) +
    geom_line(aes(y = mobility_smooth), colour = "steelblue", linewidth = 0.7,
              na.rm = TRUE) +
    labs(title = "Mobility composite: raw and trailing 7-day mean",
         x = NULL, y = "% change from baseline") +
    theme_minimal()

  p_contacts <- ggplot(cov, aes(x = date, y = contacts_raw)) +
    geom_step(colour = "firebrick", linewidth = 0.6, na.rm = TRUE) +
    labs(title = "CoMix contact matrix dominant eigenvalue",
         x = NULL, y = expression(rho(C))) +
    theme_minimal()

  p_z <- cov |>
    select(date, contacts, mobility) |>
    pivot_longer(-date, names_to = "covariate") |>
    ggplot(aes(x = date, y = value, colour = covariate)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title = "Standardised covariates", x = "Date", y = "z-score") +
    theme_minimal() +
    theme(legend.position = "bottom")

  combined <- p_mobility / p_contacts / p_z
  ggsave(file.path(config$plot_dir, "covariates.png"), combined,
         width = 10, height = 9, dpi = 300, bg = "white")

  combined
}

## Summary ---------------------------------------------------------------------

report_covariates <- function(cov) {
  vars <- c("contacts_raw", "mobility_raw", "mobility_smooth",
            "contacts", "mobility")

  out <- tibble(
    covariate = vars,
    n_missing = sapply(cov[vars], \(x) sum(is.na(x))),
    mean      = round(sapply(cov[vars], mean, na.rm = TRUE), 3),
    sd        = round(sapply(cov[vars], sd, na.rm = TRUE), 3),
    min       = round(sapply(cov[vars], min, na.rm = TRUE), 2),
    max       = round(sapply(cov[vars], max, na.rm = TRUE), 2)
  )
  print(as.data.frame(out), row.names = FALSE)
  invisible(out)
}

## Run -------------------------------------------------------------------------

ch1_covariates <- build_covariates()

report_covariates(ch1_covariates)
plot_covariates(ch1_covariates)

write_csv(ch1_covariates, ch1_cov_config$output_path)
message(sprintf("Saved to %s", ch1_cov_config$output_path))
