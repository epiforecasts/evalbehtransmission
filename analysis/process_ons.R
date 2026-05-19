# Reads raw ONS/inc2prev estimates from data-raw/ONS/, filters to England national level,
# and saves processed incidence and Rt outputs to data-processed/.
#
# inc2prev output columns:
#   name     - variable name ("infections", "R", "r", "est_prev", ...)
#   q05–q95  - posterior quantiles (absolute counts for infections; unitless for R)
#   variable - geography (e.g. "England")
#   date     - daily date

library(dplyr)
library(readr)
library(ggplot2)

ons_estimates <- read_csv("data-raw/ONS/estimates_national.csv") |>
  filter(variable == "England")

# Daily infection counts (median + 90% credible interval)
incidence_national <- ons_estimates |>
  filter(name == "infections") |>
  select(date, median = q50, q05 = q5, q95) |>
  arrange(date)

# Rt estimates from inc2prev
rt_national <- ons_estimates |>
  filter(name == "R") |>
  select(date, median = q50, q05 = q5, q95) |>
  arrange(date)

write_csv(incidence_national, "data-processed/incidence_national.csv")
write_csv(rt_national,        "data-processed/rt_national.csv")

# Age-stratified prevalence from ONS CIS (proportion positive within each age group)
# Each date has multiple estimates based on different sources? Must fix to avoid overcounting
# UPDATE: There's clearly data from different publication dates. Filter to the latest one.
cis_age <- read_csv("data-raw/inc2prev-main/data-processed/cis_age.csv") |>
  filter(geography == "England", level == "age_school") |>
  slice_max(publication_date, by = c(start_date, end_date, lower_age_limit)) |>
  select(start_date, end_date, lower_age_limit,
         prev = proportion_pos,
         prev_low_95 = proportion_pos_low_95,
         prev_high_95 = proportion_pos_high_95)

age_populations <- read_csv("data-raw/inc2prev-main/data-processed/populations.csv") |>
  filter(geography == "England", level == "age_school") |>
  select(lower_age_limit, age_pop = population)

total_pop <- sum(age_populations$age_pop)

# Create age group labels for plotting e.g. lower=2 becomes 2-10
age_limits <- sort(unique(age_populations$lower_age_limit))
age_labels <- c(paste0(head(age_limits, -1), "-", tail(age_limits, -1) -1),
                paste0(tail(age_limits, 1), "+"))
names(age_labels) <- age_limits


# Join CIS by age to population sizes
cis_age_pop <- cis_age |>
  mutate(mid_date = start_date + (end_date - start_date) / 2) |>
  left_join(age_populations, by = "lower_age_limit") |>
  mutate(
    age_group = factor(age_labels[as.character(lower_age_limit)], levels = age_labels),
    prev_in_total_pop = prev * age_pop / total_pop
    )

population_prev <- cis_age_pop |> summarise(prev_avg = sum(prev_in_total_pop), .by = mid_date)

## Plot 1: Prevalence per capita in total population
prev_totalpop_plot <- ggplot(population_prev, aes(x = mid_date, y = prev_avg)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "Prevalence per capita (total population)",
    title = "ONS CIS: COVID-19 population prevalence in England"
  ) +
  theme_classic()

ggsave("outputs/prevalence_total_pop.png", prev_totalpop_plot, width = 10, height = 5)

## Plot 2: Prevalence per capita within each age group, with population average overlay
prev_byage_plot <- ggplot(cis_age_pop, aes(x = mid_date, y = prev,
                                           colour = age_group, fill = age_group)) +
  #geom_ribbon() +
  geom_line() +
  geom_line(data = population_prev, aes(x = mid_date, y = prev_avg),
            colour = "black", linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "Prevalence within age group",
    title = "ONS CIS: COVID-19 prevalence by age group in England",
    caption = "Dashed line: population average"
  ) +
  theme_classic()

ggsave("outputs/prevalence_by_age.png", prev_byage_plot, width = 10, height = 5)
