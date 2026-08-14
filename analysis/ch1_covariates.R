# Chapter 1, issue #38: build model-ready covariates from the analysis dataset.

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("R/ch1_mobility_streams.R") # Stream selection moved here for common source across data/covariate/model scripts
source("R/ch1_contact_covariate.R") # Which contact series is used, shared with the descriptive script

## Config ----------------------------------------------------------------------

ch1_cov_config <- list(
  
  contact_covariate = contact_covariate,   # set in R/ch1_contact_covariate.R

  # Smooth over a 7-day window (as per Davies et al.)
  # Use trailing window: t-6 to t so no future information used at forecast origin
  smooth_window = 7,

  # Covariates are left unlagged here, with model lags set in ch1_gam.R

  input_path  = "data-processed/ch1_data.csv",
  output_path = "data-processed/ch1_covariates.csv",
  plot_dir    = "outputs/ch1"
)


## Transformations -------------------------------------------------------------

# Allows switching between the contact matrix eigenvalue and mean contacts per participant
# Mean contacts are age-standardised, as children join in May 2020 and report far more
# contacts, so a plain sample mean would step up on recruitment alone
select_contact_covariate <- function(dat, which) {
  switch(which,
    eigenvalue    = dat$comix_eigen,
    mean_contacts = dat$mean_contacts_standardised,
    stop("Unknown contact covariate: ", which)
  )
}

# Take the trailing mean over a window of length k i.e. mean over t-k+1 to t
# sides = 1 looks at past values, rather than centering on today
# Any 7 days contain each weekday once, so 7-day avg removes day-of-week effect exactly
# dplyr masks stats::filter if package not specified
trailing_mean <- function(x, k) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
}

# z-score over entire study period
# Uses full-period mean and sd, but a global affine rescale is only a reparameterisation to the GAM
# This gives no leakage into fitted values or forecasts
# Global keeps beta comparable across windows, and avoids huge skew during stable periods
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

## Build -----------------------------------------------------------------------

build_covariates <- function(config = ch1_cov_config) {
  dat <- read_csv(config$input_path, show_col_types = FALSE)

  # Average the mobility streams before smoothing, with equal weighting (Davies et al. has simplex alternative with weights)
  # Both transformations are linear, so this is equivalent to smoothing then averaging, just tidier 
  # Each stream has different variance, so averaging then z-scoring != z-scoring then averaging. Nouvellet does not consider this.
  
  # Equal-weighted composite; no worse than individual streams (Nouvellet et al. 2021, at least for early 2020)
  composite <- rowMeans(dat[, mobility_retained], na.rm = TRUE)
  composite[!is.finite(composite)] <- NA

  dat |>
    transmute(
      date,
      incidence,
      # contacts_raw is whichever contact series the config selects
      # Both are carried through as well, so either can be plotted or compared
      contacts_raw    = select_contact_covariate(dat, config$contact_covariate),
      comix_eigen     = dat$comix_eigen,
      mean_contacts   = dat$mean_contacts_standardised,
      mobility_raw    = composite,
      mobility_smooth = trailing_mean(composite, config$smooth_window)
    ) |>
    mutate(
      contacts = zscore(contacts_raw),
      mobility = zscore(mobility_smooth)
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

  # Stepwise, as one value per survey round is carried forward to daily
  p_eigenvalue <- ggplot(cov, aes(x = date, y = comix_eigen)) +
    geom_step(colour = "firebrick", linewidth = 0.6, na.rm = TRUE) +
    labs(title = "CoMix contact matrix dominant eigenvalue",
         x = NULL, y = expression(rho(C))) +
    theme_minimal()

  p_mean_contacts <- ggplot(cov, aes(x = date, y = mean_contacts)) +
    geom_line(colour = "firebrick", linewidth = 0.6, na.rm = TRUE) +
    labs(title = "CoMix mean contacts per participant, trailing 14-day, age-standardised",
         x = NULL, y = "Contacts per day") +
    theme_minimal()

  p_zscored <- cov |>
    select(date, contacts, mobility) |>
    pivot_longer(-date, names_to = "covariate") |>
    ggplot(aes(x = date, y = value, colour = covariate)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title = "Standardised covariates", x = "Date", y = "z-score") +
    theme_minimal() +
    theme(legend.position = "bottom")

  annotation <- patchwork::plot_annotation(
    title = "Behavioural covariate construction",
    subtitle = "Raw streams and the standardised series entering the model")

  # Both contact series, for comparing them
  combined <- p_mobility / p_eigenvalue / p_mean_contacts / p_zscored + annotation
  ggsave(file.path(config$plot_dir, "covariates.png"), combined,
         width = 10, height = 11, dpi = 300, bg = "white")

  # Mean contacts only, matching the series actually used
  without_eigenvalue <- p_mobility / p_mean_contacts / p_zscored + annotation
  ggsave(file.path(config$plot_dir, "covariates_mean_contacts.png"), without_eigenvalue,
         width = 10, height = 9, dpi = 300, bg = "white")

  invisible(combined) # Auto-printing at top level would open a device and write Rplots.pdf
}

## Summary ---------------------------------------------------------------------

report_covariates <- function(cov) {
  vars <- c("contacts_raw", "comix_eigen", "mean_contacts",
            "mobility_raw", "mobility_smooth", "contacts", "mobility")

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
