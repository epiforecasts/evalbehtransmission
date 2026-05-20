library(dplyr)
library(readr)
library(ggplot2)

## Age groups, populations and labels --------

age_populations_labelled <- read_csv("data-raw/inc2prev-main/data-processed/populations.csv") |>
  filter(geography == "England", level == "age_school") |>
  select(lower_age_limit, age_pop = population)

age_limits <- sort(unique(age_populations_labelled$lower_age_limit))
age_labels <- c(paste0(head(age_limits, -1), "-", tail(age_limits, -1) - 1),
                paste0(tail(age_limits, 1), "+"))
names(age_labels) <- age_limits

age_populations_labelled <- age_populations_labelled |>
  mutate(age_group = factor(age_labels[as.character(lower_age_limit)], levels = age_labels))

total_pop <- sum(age_populations_labelled$age_pop)


## PREVALENCE plots -----------

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

# Join CIS by age to population sizes
cis_age_pop <- cis_age |>
  mutate(mid_date = start_date + (end_date - start_date) / 2) |>
  left_join(age_populations_labelled, by = "lower_age_limit") |>
  mutate(prev_in_total_pop = prev * age_pop / total_pop)

population_prev <- cis_age_pop |> summarise(prev_avg = sum(prev_in_total_pop), .by = mid_date)

# Plot 1: Prevalence per capita in total population
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

# Plot 2: Prevalence per capita within each age group, with population average overlay
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


## INCIDENCE plots -----------
incidence_age <- read_csv("data-raw/inc2prev-main/outputs/estimates_age.csv") |>
  filter(level == "age_school", name == "infections") |> # Implicitly, is this data for only England?
  select(date, age_group = variable, inc = q50, inc_q05 = q5, inc_q95 = q95) |>
  mutate(age_group = factor(age_group, levels = levels(age_populations_labelled$age_group)))

incidence_age_pop <- incidence_age |>
  left_join(age_populations_labelled, by = "age_group") |>
  mutate(inc_per_capita = inc / age_pop,
         inc_in_total_pop = inc / total_pop)

population_inc <- incidence_age_pop |> summarise(inc_avg = sum(inc_in_total_pop), .by = date)

# Plot 3: Incidence per capita in total population
inc_totalpop_plot <- ggplot(population_inc, aes(x = date, y = inc_avg)) +
  geom_line() +
  labs(
    x = NULL,
    y = "Daily incidence per capita (total population)",
    title = "ONS CIS: COVID-19 population incidence in England"
  ) +
  theme_classic()

ggsave("outputs/incidence_total_pop.png", inc_totalpop_plot, width = 10, height = 5)

# Plot 4: Incidence per capita within each age group, with population average overlay
inc_byage_plot <- ggplot(incidence_age_pop, aes(x = date, y = inc_per_capita, colour = age_group)) +
  geom_line() +
  geom_line(data = population_inc, aes(x = date, y = inc_avg),
            colour = "black", linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
  labs(
    x = NULL,
    y = "Daily incidence per capita (within age group)",
    title = "ONS CIS: COVID-19 incidence by age group in England",
    caption = "Dashed line: population average"
  ) +
  theme_classic()

ggsave("outputs/incidence_by_age.png", inc_byage_plot, width = 10, height = 5)
