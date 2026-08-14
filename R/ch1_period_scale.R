# Shared fill scale for the OxCGRT period bands, so shading matches across every figure
# Sourced by ch1_descriptive.R, ch1_scoring.R and ch1_window_plots.R

# Chronological, matching the period_starts table in ch1_periods.R
ch1_period_levels <- c("Lockdown 1", "Summer relaxation", "Autumn tiers",
                       "Lockdown 2", "Winter tiers", "Lockdown 3")

# ggplot's default hue palette over all six periods, so fig_1_1 is unchanged
ch1_period_colours <- setNames(scales::hue_pal()(length(ch1_period_levels)),
                               ch1_period_levels)

# drop = FALSE holds each period's colour when a figure's date range excludes one
scale_fill_period <- function(...) {
  ggplot2::scale_fill_manual(values = ch1_period_colours, drop = FALSE, ...)
}
