# forecast_comparison.R
# 3-model rolling-origin forecast figure for poster (Fig 3)
# Models: Baseline GAM | +CoMix | +Mobility
# Origins: 3 illustrative dates including pre-Delta (overconfidence story)

library(mgcv)
library(MASS)
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

source("R/compute_lambda.R")

# Generation interval weights (gamma, matching temp_mvp.R)
gi_mean  <- 5.5; gi_sd <- 2.1; max_lag <- 14L
gi_shape <- (gi_mean / gi_sd)^2
gi_rate  <- gi_mean / gi_sd^2
w_raw    <- diff(pgamma(0:max_lag, shape = gi_shape, rate = gi_rate))
gi_weights <- w_raw / sum(w_raw)

# Load harmonised data and compute Lambda_t
dat <- readRDS("stashed/stash_models/harmonised_daily.rds") |>
  arrange(date) |>
  mutate(
    Lambda_t = {
      n <- n()
      lag_mat <- embed(incidence, max_lag + 1L)[, -1, drop = FALSE]
      c(rep(NA_real_, max_lag), as.vector(lag_mat %*% gi_weights))
    }
  )

# ---- Forecast function -------------------------------------------------------
# Refit model at origin_date, hold Rt (+ covariates) constant over horizon,
# propagate n_sim trajectories via renewal equation.

fit_and_forecast <- function(dat, origin_date, model_name,
                             horizon = 28L, n_sim = 300L, k = 40L) {

  # Full training set (used for incidence history in renewal)
  train_full <- dat |>
    filter(date <= origin_date, !is.na(Lambda_t), Lambda_t > 0) |>
    mutate(t = seq_len(n()))

  test <- dat |> filter(date > origin_date) |> head(horizon)
  if (nrow(test) < horizon) return(NULL)

  # Model-specific covariate filtering and formula
  if (model_name == "Baseline") {
    train_fit    <- train_full
    newdata_last <- data.frame(t = max(train_fit$t), Lambda_t = 1)
    formula_gam  <- incidence ~ s(t, k = k) + offset(log(Lambda_t))
    last_cov     <- list()
  } else if (model_name == "CoMix") {
    train_fit <- train_full |> filter(!is.na(lambda1_std))
    if (nrow(train_fit) < 30) return(NULL)
    last_cov     <- list(lambda1_std = tail(train_fit$lambda1_std, 1))
    newdata_last <- data.frame(t = max(train_fit$t), lambda1_std = last_cov$lambda1_std, Lambda_t = 1)
    formula_gam  <- incidence ~ s(t, k = k) + lambda1_std + offset(log(Lambda_t))
  } else if (model_name == "Mobility") {
    train_fit <- train_full |> filter(!is.na(mobility_retail_std))
    if (nrow(train_fit) < 30) return(NULL)
    last_cov     <- list(mobility_retail_std = tail(train_fit$mobility_retail_std, 1))
    newdata_last <- data.frame(t = max(train_fit$t), mobility_retail_std = last_cov$mobility_retail_std, Lambda_t = 1)
    formula_gam  <- incidence ~ s(t, k = k) + mobility_retail_std + offset(log(Lambda_t))
  }

  m <- tryCatch(
    gam(formula_gam, family = poisson(), data = train_fit, method = "REML"),
    error = function(e) { message("GAM failed: ", e$message); return(NULL) }
  )
  if (is.null(m)) return(NULL)

  # Linear predictor (log Rt) at last training point — held constant over horizon
  basis_last   <- predict(m, newdata = newdata_last, type = "lpmatrix")
  coef_mean    <- coef(m)
  vcov_mat     <- vcov(m)
  inc_hist     <- train_full$incidence   # full history for renewal equation

  forecast_list <- vector("list", n_sim)
  for (sim in seq_len(n_sim)) {
    coef_sim     <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_mat)
    log_Rt_const <- as.numeric(basis_last %*% coef_sim)
    inc_path     <- inc_hist
    for (day in seq_len(horizon)) {
      lam      <- compute_lambda_last(inc_path, gi_weights)
      expected <- exp(log_Rt_const) * pmax(lam, 1e-8)
      inc_path <- c(inc_path, rpois(1, lambda = expected))
    }
    forecast_list[[sim]] <- data.frame(
      sample_id   = sim,
      target_date = test$date,
      predicted   = tail(inc_path, horizon),
      observed    = test$incidence
    )
  }

  bind_rows(forecast_list) |>
    mutate(model = model_name, origin = as.character(origin_date))
}

# ---- Run forecasts -----------------------------------------------------------
# Origins chosen to show:
#   2020-09-15 — autumn wave onset, model working reasonably
#   2021-06-01 — pre-Delta peak: overconfidence story (PI fails to cover rapid rise)
#   2021-09-01 — post-Delta, heading into autumn

origins <- as.Date(c("2020-09-15", "2021-06-01", "2021-09-01"))
models  <- c("Baseline", "CoMix", "Mobility")

message("Fitting 9 models (3 origins × 3 models)...")
all_sims <- bind_rows(lapply(origins, function(o) {
  bind_rows(lapply(models, function(mod) {
    message("  ", mod, " @ ", o)
    fit_and_forecast(dat, o, mod)
  }))
}))

# Quantile ribbons
quantiles <- all_sims |>
  group_by(model, origin, target_date, observed) |>
  summarise(
    q05 = quantile(predicted, 0.05),
    q25 = quantile(predicted, 0.25),
    q50 = quantile(predicted, 0.50),
    q75 = quantile(predicted, 0.75),
    q95 = quantile(predicted, 0.95),
    .groups = "drop"
  ) |>
  mutate(
    model  = factor(model,  levels = c("Baseline", "CoMix", "Mobility")),
    origin = as.Date(origin)
  )

# ---- Plot -------------------------------------------------------------------

lockdowns <- data.frame(
  xmin = as_date(c("2020-03-23", "2020-11-05", "2021-01-05")),
  xmax = as_date(c("2020-06-01", "2020-12-02", "2021-03-08"))
)

colours <- c(
  "Baseline" = "#0072B2",
  "CoMix"    = "#009E73",
  "Mobility" = "#E69F00"
)

x_lim    <- as_date(c("2020-03-23", "2021-10-29"))
x_breaks <- seq(as_date("2020-04-01"), as_date("2021-10-01"), by = "3 months")
date_labels_fn <- function(x) ifelse(month(x) == 1, format(x, "%b\n%Y"), format(x, "%b"))

obs_line <- dat |> filter(date >= x_lim[1], date <= x_lim[2])

# Dummy data for observed line in legend
legend_obs <- data.frame(x = as_date(NA), y = NA_real_, label = "Observed")

p <- ggplot() +
  geom_rect(data = lockdowns,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  # 90% PI
  geom_ribbon(data = quantiles,
              aes(x = target_date, ymin = q05, ymax = q95,
                  fill = model, group = interaction(model, origin)),
              alpha = 0.12) +
  # 50% PI
  geom_ribbon(data = quantiles,
              aes(x = target_date, ymin = q25, ymax = q75,
                  fill = model, group = interaction(model, origin)),
              alpha = 0.28) +
  # Forecast medians
  geom_line(data = quantiles,
            aes(x = target_date, y = q50, colour = model,
                group = interaction(model, origin)),
            linewidth = 0.9) +
  # Origin date markers
  geom_vline(xintercept = origins, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  # Observed incidence (on top)
  geom_line(data = obs_line,
            aes(x = date, y = incidence),
            colour = "grey25", linewidth = 0.65) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  scale_x_date(limits = x_lim, breaks = x_breaks, labels = date_labels_fn) +
  scale_y_log10(
    breaks = c(1000, 5000, 10000, 50000, 100000, 200000),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE)
  ) +
  labs(x = NULL, y = "Infections / day") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.4),
    axis.line.y.left   = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    plot.margin        = margin(t = 4, r = 8, b = 4, l = 4)
  )

out_path <- "stashed/forecast_comparison.png"
ggsave(out_path, p, width = 12, height = 4.5, dpi = 300, bg = "white")
message("Saved: ", out_path)
