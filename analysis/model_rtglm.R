# Renewal equation embedded in a Poisson GLM following Nouvellet et al.
# Piecewise-constant Rt via 7-day time windows as categorical variables.
# Fit up to a single origin date and compare Rt estimates to naive incidence/lambda_t.

library(EpiEstim) # discr_si()
library(tidyverse)
library(patchwork)

source("R/compute_lambda.R")

si_distr   <- discr_si(k = 0:21, mu = 5.5, sigma = 2.1)
si_distr   <- si_distr / sum(si_distr)
gi_weights <- si_distr[-1]

# DATA PATH thing - will tidy this to avoid link and package repetition across scripts (or use submodules)
use_remote <- TRUE

inc2prev_path <- function(path) {
  if (use_remote) {
    paste0("https://raw.githubusercontent.com/epiforecasts/inc2prev/refs/heads/main/", path)
  } else {
    file.path("data-raw/inc2prev-main", path)
  }
}

dat <- read_csv(inc2prev_path("outputs/estimates_national.csv"),
                show_col_types = FALSE) |>
  filter(variable == "England", name == "infections") |>
  dplyr::select(date, incidence = q50) |>
  arrange(date) |>
  mutate(
    Lambda_t   = compute_lambda(incidence, gi_weights),
    Lambda_t   = if_else(Lambda_t == 0, NA_real_, Lambda_t),
    log_Lambda = log(Lambda_t),
    t          = seq_len(n())
  )


## Fit for a SINGLE origin date -------

origin_date <- as.Date("2021-10-11")

train  <- dat |> filter(date <= origin_date)
fit_df <- train |>
  filter(!is.na(Lambda_t)) |>
  mutate(
    days_since_start = as.integer(date - min(date)),
    time_window      = factor(days_since_start %/% 7)
  )

# log(E[I(t)]) = beta_w(t) + log(lambda_t), where w(t) is the 7-day window index
# exp(beta_w) = Rt for window w — no centering, directly interpretable
glm_fit <- glm(incidence ~ 0 + time_window + offset(log_Lambda),
               family = poisson(link = "log"),
               data   = fit_df)

# Extract Rt per window
se <- summary(glm_fit)$coefficients[, "Std. Error"]

window_rt <- tibble(
  time_window = levels(fit_df$time_window),
  Rt          = exp(coef(glm_fit)),
  Rt_lower    = exp(coef(glm_fit) - 1.96 * se),
  Rt_upper    = exp(coef(glm_fit) + 1.96 * se)
) |>
  mutate(time_window = as.integer(time_window))

# Map window Rt back onto daily training data
rt_glm <- fit_df |>
  mutate(
    time_window = as.integer(as.character(time_window)),
    Rt_naive    = incidence / Lambda_t
  ) |>
  left_join(window_rt, by = "time_window")


## Visualise -------

p_rt_glm <- rt_glm |>
  filter(date >= min(date) + lubridate::days(7)) |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "steelblue", alpha = 0.2) +
  geom_step(aes(y = Rt), colour = "steelblue", linewidth = 0.8) +
  geom_line(aes(y = Rt_naive), colour = "firebrick", linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  labs(title    = "GLM piecewise-constant Rt (7-day windows)",
       subtitle = expression("Blue = GLM estimate. Red dashed = incidence / " * lambda[t]),
       x = "Date", y = expression(R[t])) +
  theme_minimal()

p_rt_glm
