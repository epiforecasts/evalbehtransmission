# Chapter 1, issue #37: combine incidence, contacts and mobility into a single
# daily dataset covering the study period.

library(dplyr)
library(readr)
library(tidyr)

source("R/inc2prev_path.R")
source("R/ch1_mobility_streams.R") # Source Mobility streams used - unified source of names and subsets

## Config ----------------------------------------------------------------------

ch1_data_config <- list(
  study_start   = as.Date("2020-04-01"), # inc2prev start until Omicron period
  study_end     = as.Date("2021-11-30"), # UPDATE if vaccination period removed
  use_remote    = TRUE,   # FALSE reads inc2prev from data-raw/inc2prev-main/ - requires local copy
  contacts_path = "data-processed/comix_eigenvalues.csv", # Step-wise raw contact matrix eigenvalues from process_comix.R
  mobility_path = "data-processed/google_mobility_UK.csv", # No England series available. England = 84% of UK, intervention timings varied
  output_path   = "data-processed/ch1_data.csv"
)

## Loaders ---------------------------------------------------------------------

load_incidence <- function() {
  read_csv(inc2prev_path("outputs/estimates_national.csv", ch1_data_config$use_remote),
           show_col_types = FALSE) |>
    filter(variable == "England", name == "infections") |>
    select(date, incidence = q50) |> # Posterior median for now. Full uncertainty reserved as extension
    arrange(date)
}

# comix_eigen is rho(C), the contact matrix eigenvalue, not rho(K): no susceptibility scaling has been applied yet
# To do: either create antibody-scaled NGM, or calculate raw mean contacts (avoids circularity using rho(C) as covariate of log(Rt))
# After calculating mean contact series, save into data-processed/, read in here, and add comix_raw_contacts = (name)
load_contacts <- function() {
  read_csv(ch1_data_config$contacts_path, show_col_types = FALSE) |>
    filter(!is.na(lambda1)) |>
    select(date, comix_eigen = lambda1) |> # date is median day filled out for that survey wave, -1 as they reported for previous day
    arrange(date)
}

load_mobility <- function() {
  read_csv(ch1_data_config$mobility_path, show_col_types = FALSE) |>
    select(date, all_of(mobility_categories)) |> # Simultaneously select and rename mobility columns
    arrange(date)
}

## Weekly to daily -------------------------------------------------------------

# Convert daily series to weekly, carrying last observed value forwards
# Results in step-wise constant time series, avoiding leaking future information not available
expand_to_daily <- function(weekly, value_col, dates) {
  tibble(date = dates) |> # Create daily series
    left_join(weekly, by = "date") |> # Anchor each CoMix survey round to it's date in the daily series
    arrange(date) |>
    fill(all_of(value_col), .direction = "down") # Carry each value forwards until it hits the next wave - gives NAs until first CoMix wave
}

## Assemble --------------------------------------------------------------------

build_ch1_data <- function() {
  
  dates <- seq(ch1_data_config$study_start, ch1_data_config$study_end, by = "day") # Complete dates for analysis

  contacts_daily <- expand_to_daily(load_contacts(), "comix_eigen", dates) # Only series not daily by default

  tibble(date = dates) |>
    left_join(load_incidence(), by = "date") |>
    left_join(contacts_daily, by = "date") |>
    left_join(load_mobility(), by = "date") |>
    arrange(date) # For insurance
}

## Coverage check --------------------------------------------------------------

# Checks data has been loaded properly, in raw format, for the full analysis period
report_coverage <- function(dat) {
  cat("\nRows:", nrow(dat),
      "| Range:", format(min(dat$date)), "to", format(max(dat$date)), "\n\n")

  out <- tibble(
    column      = names(dat)[-1],
    n_missing   = sapply(dat[-1], \(x) sum(is.na(x))),
    pct_present = round(100 * sapply(dat[-1], \(x) mean(!is.na(x))), 1),
    min         = round(sapply(dat[-1], min, na.rm = TRUE), 2), 
    max         = round(sapply(dat[-1], max, na.rm = TRUE), 2)
  )
  print(as.data.frame(out), row.names = FALSE)

  # Print the first fortnight of non-zero CoMix eigenvalues to illustrate stepwise processing
  first <- which(!is.na(dat$comix_eigen))[1]
  cat("\nContacts held constant between waves:\n")
  print(as.data.frame(dat[first:(first + 13), c("date", "comix_eigen")]), row.names = FALSE)

  invisible(out) # Avoid printing table twice
}

## Run -------------------------------------------------------------------------

ch1_data <- build_ch1_data()
report_coverage(ch1_data)

write_csv(ch1_data, ch1_data_config$output_path)
message(sprintf("Saved to %s", ch1_data_config$output_path))
