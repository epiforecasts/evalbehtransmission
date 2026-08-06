# Renewal equation embedded in a Poisson GAM following Nouvellet et al.
# Autoregressive model only with no behavioural covariates, comparing naive with GLM and GAM Rt

library(mgcv) # gam() and nb()
library(EpiEstim) # discr_si() and overall_infectivity()
library(scoringutils) # as_forecast_sample(), score(), summarise_scores()
library(MASS) # mvrnorm() for parameter uncertainty draws
library(tidyverse)
library(patchwork)

source("R/compute_lambda.R")

dir.create("outputs/rtgam", recursive = TRUE, showWarnings = FALSE)

# Create generation interval, normalised - CHANGE name (serial or generation?)
si_distr   <- discr_si(k = 0:21, mu = 5.5, sigma = 2.1) # Depends on variant
si_distr   <- si_distr / sum(si_distr)
gi_weights <- si_distr[-1] # Drop lag-0 weight, giving weights for lags 1-21

# Prepare incidence data with log scale and total infectivity

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
    Lambda_t = compute_lambda(incidence, gi_weights),
    Lambda_t = if_else(Lambda_t == 0, NA_real_, Lambda_t), # Replace 0s with NA
    log_Lambda = log(Lambda_t),
    t = seq_len(n())
  )


## Fit for a SINGLE origin date - split into train and test, with 28 day horizon -------

origin_date <- as.Date("2021-10-11")
horizon <- 28L

train <- dat |> filter(date <= origin_date)
test <- dat |> filter(date > origin_date) |> head(horizon) # Keep truth up to horizon

fit_df <- train |> filter(!is.na(Lambda_t))

# Fit GAM for this date, using formula log(E[I(t)]) = s(t) + log(lambda_t)
# where lambda_t = sum_i{i=1}^s_max I_{t-s}w_s

first_model <- gam(incidence ~ 0 + s(t, k = 40) + offset(log_Lambda), # 0 intercept, s(t) captures this, easier to extract R(t)
                   family = poisson(link = "log"), # Assume Poisson. Check NegBin for sensitivity or if variance poorly captured, but introduces dispersion param to handle
                   data = fit_df,
                   method = "REML") # Handles smoothing and complexity automatically

mgcv::k.check(first_model) # Validate basis dimensions
summary(first_model)

# Recover Rt estimates on training data

# type = "link" returns full linear predictor - fine when Rt is just one term
pred <- predict(first_model, newdata = fit_df, type = "link", se.fit = TRUE)
# Extract s(t) directly. Note this is centered so R(t) = relative to geom. mean R over fitting period, not absolute value
pred_terms <- predict(first_model, newdata = fit_df, type = "terms", se.fit = TRUE)
centre_const <- attr(predict(first_model, type = "terms"), "constant")

rt_train <- fit_df |>
  dplyr::select(date, incidence, Lambda_t) |>
  mutate(
    Rt       = exp(pred$fit) / Lambda_t, # pred$fit gives full linear predictor, including log(Lambda_t), hence divide to isolate Rt
    Rt_lower = exp(pred$fit - 1.96 * pred$se.fit) / Lambda_t,
    Rt_upper = exp(pred$fit + 1.96 * pred$se.fit) / Lambda_t,
    fitted   = exp(pred$fit) # Expected incidence based on model fit
  ) |>
  mutate(Rt_naive = incidence / Lambda_t)


# Visualise fit within training data (before forecast date)

p_incidence <- rt_train |> filter(date >= min(date) + lubridate::days(6)) |>
  ggplot(aes(x = date)) +
  geom_line(aes(y = incidence), colour = "grey30", linewidth = 0.4) +
  geom_line(aes(y = fitted), colour = "steelblue", linewidth = 0.8) +
  labs(title = "Observed vs fitted incidence (training window)",
       x = NULL, y = "Infections / day") +
  theme_minimal()

p_rt <- rt_train |> filter(date >= min(date) + lubridate::days(6)) |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = Rt), colour = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  labs(title = "Estimated Rt (training window)",
       x = "Date", y = expression(R[t])) +
  theme_minimal()   

p_incidence / p_rt
ggsave("outputs/rtgam/training_fit.png", p_incidence / p_rt, width = 10, height = 6, dpi = 300, bg = "white")

# Compare to naive Rt estimate formula using deterministic renewal formula

p_rt_naive <- rt_train |> filter(date >= min(date) + lubridate::days(6)) |>
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = Rt), colour = "steelblue", linewidth = 0.8) +
  geom_line(aes(y = Rt_naive), colour = "firebrick", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  labs(title = "Estimated in-sample Rt",
       subtitle = expression("Blue = GAM estimate. Red dashed = incidence / " * lambda[t]),
       x = "Date", y = expression(R[t])) +
  theme_minimal()

p_rt_naive
ggsave("outputs/rtgam/rt_naive.png", p_rt_naive, width = 10, height = 4, dpi = 300, bg = "white")

# Create forecast, sampling from vcov matrix

n_sim <- 100 # Number of simulated trajectories
last_train_t <- max(fit_df$t) # Final date in training

# Linear predictor matrix (lpmatrix) at last training t. Rt held constant over horizon
basis_last_t <- predict(first_model,
                        newdata = data.frame(t = last_train_t, log_Lambda = 0), # Add dummy value as log_Lambda expected, but coefficient fixed at 1
                        type = "lpmatrix")

coef_mean <- coef(first_model)
vcov_matrix <- vcov(first_model) # Need to understand this


# Create n_sim forecasts and bind afterwards

inc_history <- train$incidence

forecast_list <- vector("list", n_sim)

for (sim in seq_len(n_sim)) {
  coef_sim <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_matrix)
  log_Rt_const <- as.numeric(basis_last_t %*% coef_sim)
  inc_path <- inc_history # Gives history for renewal calculations
  
  for (day in seq_len(horizon)) {
    lambda          <- compute_lambda_last(inc_path, gi_weights)
    expected_cases  <- exp(log_Rt_const) * pmax(lambda, 1e-8) # Avoids negative values
    simulated_cases <- rpois(1, lambda = expected_cases) # Sample Poisson - could REMOVE (check climateR0 reasoning)
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
    model         = "GAM-renewal-baseline",
    forecast_date = origin_date
  )

# Score using {scoringutils}

fc <- as_forecast_sample(
  forecasts,
  forecast_unit = c("model", "forecast_date", "target_date", "horizon")
)

scores_daily <- score(fc) # Predictions are integer valued. Check scoringRules for discrete prob distributions
summarise_scores(scores_daily, by = "model")
summarise_scores(scores_daily, by = c("model", "horizon"))



# Compute forecast quantiles across samples for each horizon day
forecast_quantiles <- as_forecast_quantile(
  fc,
  quantile_level = c(0.05, 0.25, 0.50, 0.75, 0.95)
) |>
  as_tibble() |>
  pivot_wider(names_from = quantile_level, values_from = predicted) |>
  rename(q05 = `0.05`, q25 = `0.25`, q50 = `0.5`, q75 = `0.75`, q95 = `0.95`)

p_forecast <- ggplot() +
  geom_line(data = rt_train, aes(x = date, y = incidence),
            colour = "grey60", linewidth = 0.4) +
  geom_line(data = rt_train, aes(x = date, y = fitted),
            colour = "steelblue", linewidth = 0.6) +
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
  labs(title    = "GAM renewal model: fit and forecast",
       subtitle = "Blue = fitted training | Red = forecast (50% and 90% PI) | Points = observed",
       x = "Date", y = "Infections / day") +
  theme_minimal()

p_forecast # Low uncertainty given huge training data, and infection counts are large compared to standard deviation
ggsave("outputs/rtgam/forecast_single.png", p_forecast, width = 10, height = 4, dpi = 300, bg = "white")



## Rolling-origin forecast function ---------------

fit_and_forecast <- function(dat, origin_date, horizon = 28L, gi_weights,
                             n_sim = 100L, k = 40L) {

  # Train before forecast date, and test on horizon dates
  train <- dat |> filter(date <= origin_date)
  test  <- dat |> filter(date > origin_date) |> head(horizon)

  if (nrow(test) < horizon) return(NULL)  # Not enough true dates to score

  fit_df <- train |> filter(!is.na(Lambda_t))

  # Fit GAM model, with k knots
  m <- gam(incidence ~ 0 + s(t, k = k) + offset(log_Lambda),
            family = poisson(link = "log"),
            data   = fit_df,
            method = "REML")

  # Extract latest time in training data
  last_train_t <- max(fit_df$t)
  basis_last_t <- predict(m,
                          newdata = data.frame(t = last_train_t, log_Lambda = 0),
                          type    = "lpmatrix")

  coef_mean   <- coef(m)
  vcov_matrix <- vcov(m)

  inc_history   <- train$incidence
  forecast_list <- vector("list", n_sim)

  # Simulate n_sim trajectories, using renewal equation and sampled coefficients
  for (sim in seq_len(n_sim)) {
    coef_sim     <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_matrix)
    log_Rt_const <- as.numeric(basis_last_t %*% coef_sim)
    inc_path     <- inc_history

    for (day in seq_len(horizon)) {
      lambda          <- compute_lambda_last(inc_path, gi_weights)
      expected_cases  <- exp(log_Rt_const) * pmax(lambda, 1e-8)
      simulated_cases <- rpois(1, lambda = expected_cases) # Sample Poisson - could REMOVE (check climateR0 reasoning)
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

  # Return data frame with information required for {scoringutils}
  bind_rows(forecast_list) |>
    mutate(
      model         = "GAM-renewal-baseline",
      forecast_date = origin_date
    )
}


## Rolling-origin forecasts -----------


origins <- seq(min(dat$date) + 90, as.Date("2021-10-11"), by = "1 month") # Fit each month from 3 months in

all_forecasts <- bind_rows(lapply(origins, function(o)
  fit_and_forecast(dat, origin_date = o, horizon = 28L,
                   gi_weights = gi_weights, n_sim = 100L)))

# Score all origins
fc_all <- as_forecast_sample(
  all_forecasts,
  forecast_unit = c("model", "forecast_date", "target_date", "horizon")
)

scores_all <- score(fc_all)
summarise_scores(scores_all, by = "model")
summarise_scores(scores_all, by = c("model", "horizon"))
summarise_scores(scores_all, by = c("model", "forecast_date"))


# Compute quantiles per forecast date and target date
all_quantiles <- as_forecast_quantile(
  fc_all,
  quantile_level = c(0.05, 0.25, 0.50, 0.75, 0.95)
) |>
  as_tibble() |>
  pivot_wider(names_from = quantile_level, values_from = predicted) |>
  rename(q05 = `0.05`, q25 = `0.25`, q50 = `0.5`, q75 = `0.75`, q95 = `0.95`)

p_rolling <- ggplot() +
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
  labs(title    = "GAM renewal model: rolling-origin forecasts",
       subtitle = "Red = forecast fans (50% and 90% PI) | Grey = observed incidence",
       x = "Date", y = "Infections / day") +
  theme_minimal()

p_rolling # Very little spread, given Rt held constant and large infection counts
ggsave("outputs/rtgam/rolling_forecasts.png", p_rolling, width = 10, height = 4, dpi = 300, bg = "white")
