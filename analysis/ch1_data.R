# Chapter 1, issue #37: combine incidence, contacts and mobility into a single
# daily dataset covering the study period.

library(dplyr)
library(readr)
library(tidyr)

## Config ----------------------------------------------------------------------

ch1_config <- list(
  study_start   = as.Date("2020-04-01"),
  study_end     = as.Date("2021-11-30"),
  use_remote    = TRUE,   # FALSE reads inc2prev from data-raw/inc2prev-main/
  contacts_path = "data-processed/comix_eigenvalues.csv",
  mobility_path = "data-processed/google_mobility_UK.csv",
  output_path   = "data-processed/ch1_analysis_data.csv"
)

# All six categories carried through; retention is decided later.
mobility_categories <- c(
  retail_recreation = "retail_and_recreation_percent_change_from_baseline",
  grocery_pharmacy  = "grocery_and_pharmacy_percent_change_from_baseline",
  parks             = "parks_percent_change_from_baseline",
  transit           = "transit_stations_percent_change_from_baseline",
  workplaces        = "workplaces_percent_change_from_baseline",
  residential       = "residential_percent_change_from_baseline"
)

inc2prev_path <- function(path) {
  if (ch1_config$use_remote) {
    paste0("https://raw.githubusercontent.com/epiforecasts/inc2prev/refs/heads/main/", path)
  } else {
    file.path("data-raw/inc2prev-main", path)
  }
}

## Loaders ---------------------------------------------------------------------

load_incidence <- function() {
  read_csv(inc2prev_path("outputs/estimates_national.csv"), show_col_types = FALSE) |>
    filter(variable == "England", name == "infections") |>
    select(date, incidence = q50) |>
    arrange(date)
}

# comix_eigen is rho(C), the contact matrix eigenvalue, not rho(K): no
# susceptibility scaling has been applied.
load_contacts <- function() {
  read_csv(ch1_config$contacts_path, show_col_types = FALSE) |>
    filter(!is.na(lambda1)) |>
    select(date, comix_eigen = lambda1) |>
    arrange(date)
}

load_mobility <- function() {
  read_csv(ch1_config$mobility_path, show_col_types = FALSE) |>
    select(date, all_of(unname(mobility_categories))) |>
    rename(!!!mobility_categories) |>
    arrange(date)
}

## Weekly to daily -------------------------------------------------------------

# Last observation carried forward from the survey date: piecewise constant, and
# uses no future information since these feed forecasts.
expand_to_daily <- function(weekly, value_col, dates) {
  tibble(date = dates) |>
    left_join(weekly, by = "date") |>
    arrange(date) |>
    fill(all_of(value_col), .direction = "down")
}

## Assemble --------------------------------------------------------------------

build_ch1_data <- function() {
  dates <- seq(ch1_config$study_start, ch1_config$study_end, by = "day")

  contacts_daily <- expand_to_daily(load_contacts(), "comix_eigen", dates)

  tibble(date = dates) |>
    left_join(load_incidence(), by = "date") |>
    left_join(contacts_daily, by = "date") |>
    left_join(load_mobility(), by = "date") |>
    arrange(date)
}

## Coverage check --------------------------------------------------------------

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

  # Show the first fortnight of non-missing contacts, where the steps are visible
  first <- which(!is.na(dat$comix_eigen))[1]
  cat("\nContacts held constant between waves:\n")
  print(as.data.frame(dat[first:(first + 13), c("date", "comix_eigen")]), row.names = FALSE)

  invisible(out)
}

## Run -------------------------------------------------------------------------

ch1_data <- build_ch1_data()
report_coverage(ch1_data)

write_csv(ch1_data, ch1_config$output_path)
message(sprintf("Saved to %s", ch1_config$output_path))
