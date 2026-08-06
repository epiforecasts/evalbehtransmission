# Renewal equation embedded in a Poisson GLM following Nouvellet et al.
# Piecewise-constant Rt via 7-day time windows as categorical variables.
# Fit up to a single origin date and compare Rt estimates to naive incidence/lambda_t.

library(EpiEstim) # discr_si()
library(scoringutils)
library(MASS)
library(tidyverse)
library(patchwork)

source("R/compute_lambda.R")

dir.create("outputs/rtglm", recursive = TRUE, showWarnings = FALSE)

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
ggsave("outputs/rtglm/rt_glm.png", p_rt_glm, width = 10, height = 4, dpi = 300, bg = "white")

## Forecast from single origin date -------

horizon <- 28L
test <- dat |> filter(date > origin_date) |> head(horizon)

n_sim <- 100

last_window  <- max(as.integer(as.character(fit_df$time_window)))
newdata_last <- data.frame(
  time_window = factor(last_window, levels = levels(fit_df$time_window)),
  log_Lambda  = 0
)
X_last <- model.matrix(~ 0 + time_window, data = newdata_last)

coef_mean   <- coef(glm_fit)
vcov_matrix <- vcov(glm_fit)

inc_history   <- train$incidence
forecast_list <- vector("list", n_sim)

for (sim in seq_len(n_sim)) {
  coef_sim     <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_matrix)
  log_Rt_const <- as.numeric(X_last %*% coef_sim)
  inc_path     <- inc_history

  for (day in seq_len(horizon)) {
    lambda          <- compute_lambda_last(inc_path, gi_weights)
    expected_cases  <- exp(log_Rt_const) * pmax(lambda, 1e-8)
    simulated_cases <- rpois(1, lambda = expected_cases)
    inc_path        <- c(inc_path, simulated_cases)
  }

  forecast_list[[sim]] <- data.frame(
    sample_id   = sim,
    horizon     = seq_len(horizon),
    predicted   = inc_path[seq(length(inc_history) + 1, length(inc_path))],
    target_date = test$date,
    observed    = test$incidence
  )
}

forecasts <- bind_rows(forecast_list) |>
  mutate(
    model         = "GLM-renewal-baseline",
    forecast_date = origin_date
  )

fc <- as_forecast_sample(
  forecasts,
  forecast_unit = c("model", "forecast_date", "target_date", "horizon")
)

scores_daily <- score(fc)
summarise_scores(scores_daily, by = "model")
summarise_scores(scores_daily, by = c("model", "horizon"))

forecast_quantiles <- as_forecast_quantile(
  fc,
  quantile_level = c(0.05, 0.25, 0.50, 0.75, 0.95)
) |>
  as_tibble() |>
  pivot_wider(names_from = quantile_level, values_from = predicted) |>
  rename(q05 = `0.05`, q25 = `0.25`, q50 = `0.5`, q75 = `0.75`, q95 = `0.95`)

p_forecast_glm <- ggplot() +
  geom_line(data = rt_glm, aes(x = date, y = incidence),
            colour = "grey60", linewidth = 0.4) +
  geom_ribbon(data = forecast_quantiles,
              aes(x = target_date, ymin = q05, ymax = q95),
              fill = "firebrick", alpha = 0.15) +
  geom_ribbon(data = forecast_quantiles,
              aes(x = target_date, ymin = q25, ymax = q75),
              fill = "firebrick", alpha = 0.25) +
  geom_line(data = forecast_quantiles,
            aes(x = target_date, y = q50),
            colour = "firebrick", linewidth = 0.8) +
  geom_point(data = forecast_quantiles,
             aes(x = target_date, y = observed),
             colour = "grey20", size = 1) +
  geom_vline(xintercept = origin_date, linetype = "dashed", colour = "grey40") +
  labs(title    = "GLM renewal model: fit and forecast",
       subtitle = "Grey = observed (training) | Red = forecast (50% and 90% PI) | Points = observed",
       x = "Date", y = "Infections / day") +
  theme_minimal()

p_forecast_glm
ggsave("outputs/rtglm/forecast_single.png", p_forecast_glm, width = 10, height = 4, dpi = 300, bg = "white")

## Rolling-origin forecast function ---------------

fit_and_forecast_glm <- function(dat, origin_date, horizon = 28L, gi_weights,
                                 n_sim = 100L) {

  train <- dat |> filter(date <= origin_date)
  test  <- dat |> filter(date > origin_date) |> head(horizon)

  if (nrow(test) < horizon) return(NULL)

  fit_df <- train |>
    filter(!is.na(Lambda_t)) |>
    mutate(
      days_since_start = as.integer(date - min(date)),
      time_window      = factor(days_since_start %/% 7)
    )

  m <- glm(incidence ~ 0 + time_window + offset(log_Lambda),
            family = poisson(link = "log"),
            data   = fit_df)

  last_window  <- max(as.integer(as.character(fit_df$time_window)))
  newdata_last <- data.frame(
    time_window = factor(last_window, levels = levels(fit_df$time_window)),
    log_Lambda  = 0
  )
  X_last <- model.matrix(~ 0 + time_window, data = newdata_last)

  coef_mean   <- coef(m)
  vcov_matrix <- vcov(m)

  inc_history   <- train$incidence
  forecast_list <- vector("list", n_sim)

  for (sim in seq_len(n_sim)) {
    coef_sim     <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_matrix)
    log_Rt_const <- as.numeric(X_last %*% coef_sim)
    inc_path     <- inc_history

    for (day in seq_len(horizon)) {
      lambda          <- compute_lambda_last(inc_path, gi_weights)
      expected_cases  <- exp(log_Rt_const) * pmax(lambda, 1e-8)
      simulated_cases <- rpois(1, lambda = expected_cases)
      inc_path        <- c(inc_path, simulated_cases)
    }

    forecast_list[[sim]] <- data.frame(
      sample_id   = sim,
      horizon     = seq_len(horizon),
      predicted   = inc_path[seq(length(inc_history) + 1, length(inc_path))],
      target_date = test$date,
      observed    = test$incidence
    )
  }

  bind_rows(forecast_list) |>
    mutate(
      model         = "GLM-renewal-baseline",
      forecast_date = origin_date
    )
}


## Rolling-origin forecasts -----------

origins <- seq(min(dat$date) + 90, as.Date("2021-10-11"), by = "1 month")

all_forecasts <- bind_rows(lapply(origins, function(o)
  fit_and_forecast_glm(dat, origin_date = o, horizon = 28L,
                       gi_weights = gi_weights, n_sim = 100L)))

fc_all <- as_forecast_sample(
  all_forecasts,
  forecast_unit = c("model", "forecast_date", "target_date", "horizon")
)

scores_all <- score(fc_all)
summarise_scores(scores_all, by = "model")
summarise_scores(scores_all, by = c("model", "horizon"))
summarise_scores(scores_all, by = c("model", "forecast_date"))

all_quantiles <- as_forecast_quantile(
  fc_all,
  quantile_level = c(0.05, 0.25, 0.50, 0.75, 0.95)
) |>
  as_tibble() |>
  pivot_wider(names_from = quantile_level, values_from = predicted) |>
  rename(q05 = `0.05`, q25 = `0.25`, q50 = `0.5`, q75 = `0.75`, q95 = `0.95`)

p_rolling_glm <- ggplot() +
  geom_ribbon(data = all_quantiles,
              aes(x = target_date, ymin = q05, ymax = q95,
                  group = forecast_date),
              fill = "firebrick", alpha = 0.05) +
  geom_ribbon(data = all_quantiles,
              aes(x = target_date, ymin = q25, ymax = q75,
                  group = forecast_date),
              fill = "firebrick", alpha = 0.10) +
  geom_line(data = all_quantiles,
            aes(x = target_date, y = q50, group = forecast_date),
            colour = "firebrick", linewidth = 0.3, alpha = 0.5) +
  geom_line(data = dat |> filter(date <= as.Date("2021-10-11") + 28),
            aes(x = date, y = incidence),
            colour = "grey30", linewidth = 0.5) +
  coord_cartesian(xlim = c(min(dat$date), as.Date("2021-10-11") + 28)) +
  labs(title    = "GLM renewal model: rolling-origin forecasts",
       subtitle = "Red = forecast fans (50% and 90% PI) | Grey = observed incidence",
       x = "Date", y = "Infections / day") +
  theme_minimal()

p_rolling_glm
ggsave("outputs/rtglm/rolling_forecasts.png", p_rolling_glm, width = 10, height = 4, dpi = 300, bg = "white")
