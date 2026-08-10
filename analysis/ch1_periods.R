# Chapter 1, issue #39: define pandemic periods from OxCGRT and provide the
# Stringency Index for annotating figures. Descriptive only, never a covariate.

library(dplyr)
library(readr)
library(ggplot2)

# I've switched to 'compact', not 'simplified' version of OxCGRT
# This gives targeted vs national flag, whereas Stringency Index alone reflects the most stringent application

## Config ----------------------------------------------------------------------

ch1_period_config <- list(
  study_start = as.Date("2020-04-01"),
  study_end   = as.Date("2021-11-30"),

  # oxcgrt_path = "data-raw/OxCGRT/OxCGRT_simplified_v1.csv",
  oxcgrt_path = "data-raw/OxCGRT/OxCGRT_compact_subnational_v1.csv",
  output_path = "data-processed/ch1_periods.csv",
  plot_dir    = "outputs/ch1"
)

# Named policy regimes based on relevant policies in place
# C1 (school closing) and C6 (stay-at-home) are two of eight containment indicators
# Used illustratively. Restrictions largely ended after 19/07/2021 until Omicron plans
period_starts <- tibble::tribble(
  ~start,                    ~period,
  as.Date("2020-04-01"),     "Lockdown 1",
  as.Date("2020-05-13"),     "Summer relaxation",
  as.Date("2020-10-12"),     "Autumn tiers",
  as.Date("2020-11-05"),     "Lockdown 2",
  as.Date("2020-12-03"),     "Winter tiers", # C6 dropped to 1 on 3/12, despite winter tiers announced 2/12
  as.Date("2021-01-05"),     "Lockdown 3",
  as.Date("2021-03-08"),     "Staged reopening",
  as.Date("2021-07-19"),     "Post-restrictions"
)

## Load ------------------------------------------------------------------------

load_oxcgrt_england <- function(config = ch1_period_config) {
  read_csv(config$oxcgrt_path, show_col_types = FALSE) |>
    filter(RegionCode == "UK_ENG") |>
    mutate(date = as.Date(as.character(Date), "%Y%m%d")) |>
    # Simplified file:
    # select(date,
    #        school_closing = C1M_combined_numeric,
    #        stay_at_home   = C6M_combined_numeric,
    #        stringency     = StringencyIndex_Average) |>
    # Compact file, which adds flags: 1 = policy applied across England,
    # 0 = targeted to part of it only.
    select(date,
           school_closing      = `C1M_School.closing`,
           school_closing_flag = C1M_Flag,
           stay_at_home        = `C6M_Stay.at.home.requirements`,
           stay_at_home_flag   = C6M_Flag,
           stringency          = StringencyIndex_Average) |>
    filter(date >= config$study_start, date <= config$study_end) |>
    arrange(date)
}

## Assign periods --------------------------------------------------------------

# Splits the full study period into bins, with 'period' label for the relevant OxCGRT period
assign_periods <- function(dat, starts = period_starts) {
  dat |>
    mutate(
      period = cut(date,
                   breaks = c(starts$start, max(date) + 1), # Add a 1-day buffer to cover the final window
                   labels = starts$period,
                   right  = FALSE)) # Give left-closed intervals, as periods use start date
}

## Checks ----------------------------------------------------------------------

# For each period, report the start and end dates, including duration and mean stringency
# Gives idea of how many forecasts originate in each period too
report_periods <- function(dat) {
  cat("Rows:", nrow(dat), "| unlabelled dates:", sum(is.na(dat$period)), "\n\n")

  out <- dat |>
    group_by(period) |>
    summarise(start = min(date), end = max(date), days = n(),
              mean_stringency = round(mean(stringency), 1),
              .groups = "drop") |>
    arrange(start) # Undo the alphabetical grouping

  print(as.data.frame(out), row.names = FALSE)
  invisible(out)
}

## Plot ------------------------------------------------------------------------

plot_stringency <- function(dat, config = ch1_period_config) {
  dir.create(config$plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Create rectangles for each period alongside the Stringency line plot
  bands <- dat |>
    group_by(period) |>
    summarise(start = min(date), end = max(date), .groups = "drop")

  p <- ggplot(dat, aes(x = date)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = period), alpha = 0.25) +
    geom_line(aes(y = stringency), linewidth = 0.7) +
    labs(title = "OxCGRT Stringency Index, England, with pandemic periods",
         x = "Date", y = "Stringency Index (0-100)", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(config$plot_dir, "periods_stringency.png"), p,
         width = 11, height = 5, dpi = 300, bg = "white")
  p
}

## Run -------------------------------------------------------------------------

ch1_periods <- load_oxcgrt_england() |> assign_periods()

report_periods(ch1_periods)
plot_stringency(ch1_periods)

write_csv(ch1_periods, ch1_period_config$output_path)
message(sprintf("Saved to %s", ch1_period_config$output_path))
