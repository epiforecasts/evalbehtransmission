# rolling_forecast_comparison.R
# Rolling-origin forecasts for all three models (Baseline / +CoMix / +Mobility)
# Produces:
#   outputs/rolling_forecast_plot.png   — styled poster figure
#   outputs/forecast_metrics.csv        — CRPS, bias, coverage by model (+by week)

library(mgcv)
library(MASS)
library(scoringutils)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)

source("R/compute_lambda.R")

dir.create("stashed", showWarnings = FALSE)

# Generation interval (gamma, matching temp_mvp.R)
gi_mean  <- 5.5; gi_sd <- 2.1; max_lag <- 14L
w_raw    <- diff(pgamma(0:max_lag,
                        shape = (gi_mean/gi_sd)^2,
                        rate  =  gi_mean/gi_sd^2))
gi_weights <- w_raw / sum(w_raw)

# Load harmonised data and compute Lambda_t
dat <- readRDS("stashed/stash_models/harmonised_daily.rds") |>
  arrange(date) |>
  mutate(
    Lambda_t = {
      lag_mat <- embed(incidence, max_lag + 1L)[, -1, drop = FALSE]
      c(rep(NA_real_, max_lag), as.vector(lag_mat %*% gi_weights))
    }
  )

# ---- Forecast function -------------------------------------------------------

fit_and_forecast <- function(dat, origin_date, model_name,
                             horizon = 28L, n_sim = 200L, k = 40L) {

  train_full <- dat |>
    filter(date <= origin_date, !is.na(Lambda_t), Lambda_t > 0) |>
    mutate(t = seq_len(n()))

  test <- dat |> filter(date > origin_date) |> head(horizon)
  if (nrow(test) < horizon) return(NULL)

  if (model_name == "Baseline") {
    train_fit    <- train_full
    formula_gam  <- incidence ~ s(t, k = k) + offset(log(Lambda_t))
    newdata_last <- data.frame(t = max(train_fit$t), Lambda_t = 1)
  } else if (model_name == "CoMix") {
    train_fit <- train_full |> filter(!is.na(lambda1_std))
    if (nrow(train_fit) < 30) return(NULL)
    formula_gam  <- incidence ~ s(t, k = k) + lambda1_std + offset(log(Lambda_t))
    newdata_last <- data.frame(t          = max(train_fit$t),
                               lambda1_std = tail(train_fit$lambda1_std, 1),
                               Lambda_t   = 1)
  } else if (model_name == "Mobility") {
    train_fit <- train_full |> filter(!is.na(mobility_retail_std))
    if (nrow(train_fit) < 30) return(NULL)
    formula_gam  <- incidence ~ s(t, k = k) + mobility_retail_std + offset(log(Lambda_t))
    newdata_last <- data.frame(t                   = max(train_fit$t),
                               mobility_retail_std  = tail(train_fit$mobility_retail_std, 1),
                               Lambda_t             = 1)
  }

  m <- tryCatch(
    gam(formula_gam, family = poisson(), data = train_fit, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  basis_last <- predict(m, newdata = newdata_last, type = "lpmatrix")
  coef_mean  <- coef(m)
  vcov_mat   <- vcov(m)
  inc_hist   <- train_full$incidence

  sims <- vector("list", n_sim)
  for (sim in seq_len(n_sim)) {
    coef_sim     <- MASS::mvrnorm(1, mu = coef_mean, Sigma = vcov_mat)
    log_Rt_const <- as.numeric(basis_last %*% coef_sim)
    inc_path     <- inc_hist
    for (day in seq_len(horizon)) {
      lam      <- compute_lambda_last(inc_path, gi_weights)
      inc_path <- c(inc_path, rpois(1, exp(log_Rt_const) * pmax(lam, 1e-8)))
    }
    sims[[sim]] <- data.frame(
      sample_id   = sim,
      horizon     = seq_len(horizon),
      target_date = test$date,
      predicted   = tail(inc_path, horizon),
      observed    = test$incidence
    )
  }

  bind_rows(sims) |>
    mutate(model = model_name, forecast_date = origin_date)
}

# ---- Rolling-origin evaluation ----------------------------------------------

origins <- seq(min(dat$date, na.rm = TRUE) + 90,
               as.Date("2021-10-11"),
               by = "1 month")
models  <- c("Baseline", "CoMix", "Mobility")

message("Running rolling-origin forecasts (", length(origins), " origins × 3 models)...")
all_forecasts <- bind_rows(lapply(origins, function(o) {
  bind_rows(lapply(models, function(mod) {
    message("  ", mod, " @ ", o)
    fit_and_forecast(dat, o, mod)
  }))
}))

write_csv(all_forecasts, "stashed/rolling_forecasts_raw.csv")
message("Raw forecasts saved.")

# ---- Score with scoringutils ------------------------------------------------

fc <- as_forecast_sample(
  all_forecasts,
  forecast_unit = c("model", "forecast_date", "target_date", "horizon")
)

scores <- score(fc)

# Overall metrics by model
metrics_model <- summarise_scores(scores, by = "model") |>
  select(model, crps, bias, mad) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

# CRPS by model and forecast week (1–4)
metrics_week <- all_forecasts |>
  mutate(week = paste0("Week ", ceiling(horizon / 7))) |>
  as_forecast_sample(
    forecast_unit = c("model", "forecast_date", "target_date", "horizon", "week")
  ) |>
  score() |>
  summarise_scores(by = c("model", "week")) |>
  select(model, week, crps) |>
  mutate(crps = round(crps, 1)) |>
  pivot_wider(names_from = week, values_from = crps)

# 50% and 90% empirical coverage
coverage <- all_forecasts |>
  group_by(model, forecast_date, target_date) |>
  summarise(
    q05 = quantile(predicted, 0.05),
    q25 = quantile(predicted, 0.25),
    q75 = quantile(predicted, 0.75),
    q95 = quantile(predicted, 0.95),
    obs = first(observed),
    .groups = "drop"
  ) |>
  mutate(
    cov50 = obs >= q25 & obs <= q75,
    cov90 = obs >= q05 & obs <= q95
  ) |>
  group_by(model) |>
  summarise(
    `50% coverage` = round(mean(cov50, na.rm = TRUE), 2),
    `90% coverage` = round(mean(cov90, na.rm = TRUE), 2),
    .groups = "drop"
  )

metrics_full <- metrics_model |>
  left_join(metrics_week,  by = "model") |>
  left_join(coverage,      by = "model") |>
  mutate(model = factor(model, levels = c("Baseline", "CoMix", "Mobility"))) |>
  arrange(model)

write_csv(metrics_full, "stashed/forecast_metrics.csv")
message("Metrics saved to stashed/forecast_metrics.csv")
print(metrics_full)

# ---- Rolling forecast plot --------------------------------------------------

# Only LD2 and LD3 are in view (LD1 is before Sep 2020)
lockdowns <- data.frame(
  xmin    = as_date(c("2020-11-05", "2021-01-05")),
  xmax    = as_date(c("2020-12-02", "2021-03-08")),
  label   = c("Lockdown 2", "Lockdown 3"),
  label_x = as_date(c("2020-11-18", "2021-02-04"))
)

colours <- c(
  "Observed" = "grey45",
  "Baseline" = "grey55",
  "CoMix"    = "#009E73",
  "Mobility" = "#E69F00"
)

linetypes <- c(
  "Observed" = "solid",
  "Baseline" = "dashed",
  "CoMix"    = "solid",
  "Mobility" = "solid"
)

colour_labels <- c(
  "Observed" = "Observed",
  "Baseline" = "Baseline GAM",
  "CoMix"    = expression("+ CoMix " * lambda[1]),
  "Mobility" = "+ Mobility"
)

x_lim    <- as_date(c("2020-09-01", "2021-09-01"))
x_breaks <- seq(as_date("2020-09-01"), as_date("2021-09-01"), by = "2 months")
date_labels_fn <- function(x) ifelse(month(x) == 1, format(x, "%b\n%Y"), format(x, "%b"))

# Per-origin quantiles including median
# group = interaction(model, forecast_date) keeps each 28-day segment separate
rolling_quantiles <- all_forecasts |>
  group_by(model, forecast_date, target_date) |>
  summarise(
    q05 = quantile(predicted, 0.05),
    q25 = quantile(predicted, 0.25),
    q50 = quantile(predicted, 0.50),
    q75 = quantile(predicted, 0.75),
    q95 = quantile(predicted, 0.95),
    .groups = "drop"
  ) |>
  mutate(model = factor(model, levels = c("Baseline", "CoMix", "Mobility")))

obs_line <- dat |> filter(date >= x_lim[1], date <= x_lim[2])

p <- ggplot() +
  # Lockdown shading
  geom_rect(data = lockdowns,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  # Lockdown text labels — hjust=1 anchors top of rotated text at y (near top of box)
  geom_text(data = lockdowns,
            aes(x = label_x, y = 220000, label = label),
            angle = 90, colour = "grey55", size = 3.8,
            hjust = 1, inherit.aes = FALSE) +
  # 90% PI per origin
  geom_ribbon(data = rolling_quantiles,
              aes(x = target_date, ymin = q05, ymax = q95,
                  fill = model, group = interaction(model, forecast_date)),
              alpha = 0.20) +
  # 50% PI per origin
  geom_ribbon(data = rolling_quantiles,
              aes(x = target_date, ymin = q25, ymax = q75,
                  fill = model, group = interaction(model, forecast_date)),
              alpha = 0.45) +
  # Median lines: Baseline dashed grey, behavioural models solid coloured
  geom_line(data = rolling_quantiles,
            aes(x = target_date, y = q50, colour = model, linetype = model,
                group = interaction(model, forecast_date)),
            linewidth = 0.65, alpha = 0.85) +
  # Observed — solid, slightly lighter
  geom_line(data = obs_line,
            aes(x = date, y = incidence, colour = "Observed", linetype = "Observed"),
            linewidth = 0.75) +
  scale_colour_manual(
    values  = colours,
    labels  = colour_labels,
    breaks  = c("Observed", "Baseline", "CoMix", "Mobility"),
    name    = NULL
  ) +
  scale_linetype_manual(
    values  = linetypes,
    labels  = colour_labels,
    breaks  = c("Observed", "Baseline", "CoMix", "Mobility"),
    name    = NULL
  ) +
  scale_fill_manual(
    values  = c("Baseline" = "grey70", "CoMix" = "#009E73", "Mobility" = "#E69F00"),
    guide   = "none"
  ) +
  scale_x_date(limits = x_lim, breaks = x_breaks, labels = date_labels_fn) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL, y = "Daily infections") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid         = element_blank(),
    panel.background   = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.4),
    axis.line.y.left   = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    legend.position    = c(0.83, 0.85),
    legend.background  = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.margin      = margin(4, 8, 4, 8),
    legend.key         = element_blank(),
    legend.text        = element_text(size = 13),
    plot.margin        = margin(4, 8, 2, 4)
  )

# ---- Table panel ------------------------------------------------------------

# Dev.expl: fit each model on training data up to end of plot range
train_base <- dat |>
  filter(date <= x_lim[2], !is.na(Lambda_t), Lambda_t > 0) |>
  mutate(t = seq_len(n()))

m_base  <- gam(incidence ~ s(t, k = 40L) + offset(log(Lambda_t)),
               family = poisson(), data = train_base, method = "REML")
dev_base  <- paste0(round(summary(m_base)$dev.expl * 100, 1), "%")

train_comix <- train_base |> filter(!is.na(lambda1_std))
m_comix <- gam(incidence ~ s(t, k = 40L) + lambda1_std + offset(log(Lambda_t)),
               family = poisson(), data = train_comix, method = "REML")
dev_comix <- paste0(round(summary(m_comix)$dev.expl * 100, 1), "%")

train_mob <- train_base |> filter(!is.na(mobility_retail_std))
m_mob   <- gam(incidence ~ s(t, k = 40L) + mobility_retail_std + offset(log(Lambda_t)),
               family = poisson(), data = train_mob, method = "REML")
dev_mob   <- paste0(round(summary(m_mob)$dev.expl * 100, 1), "%")

# CRPS: target dates within the plot x range
crps_range <- all_forecasts |>
  filter(target_date >= x_lim[1], target_date <= x_lim[2]) |>
  as_forecast_sample(forecast_unit = c("model", "forecast_date", "target_date", "horizon")) |>
  score() |>
  summarise_scores(by = "model") |>
  select(model, crps)
crps_vals <- setNames(format(round(crps_range$crps, 0), big.mark = ",", trim = TRUE),
                      crps_range$model)

tab_data <- data.frame(
  model_parse = c('"Baseline GAM"', '"+ CoMix "*lambda[1]', '"+ Mobility"'),
  colour      = c("grey55", "#009E73", "#E69F00"),
  dev_expl    = c(dev_base, dev_comix, dev_mob),
  crps        = c(crps_vals["Baseline"], crps_vals["CoMix"], crps_vals["Mobility"]),
  y           = c(3, 2, 1),
  stringsAsFactors = FALSE
)

p_table <- ggplot() +
  annotate("segment", x = 0, xend = 1, y = 4.55, yend = 4.55, colour = "grey30", linewidth = 0.4) +
  annotate("segment", x = 0, xend = 1, y = 3.45, yend = 3.45, colour = "grey30", linewidth = 0.4) +
  annotate("segment", x = 0, xend = 1, y = 0.55, yend = 0.55, colour = "grey30", linewidth = 0.4) +
  annotate("text", x = 0.02, y = 4, label = "Model",              fontface = "bold", hjust = 0, size = 4.8, colour = "grey20") +
  annotate("text", x = 0.42, y = 4, label = "Explained deviance", fontface = "bold", hjust = 0, size = 4.8, colour = "grey20") +
  annotate("text", x = 0.74, y = 4, label = "Mean CRPS",          fontface = "bold", hjust = 0, size = 4.8, colour = "grey20") +
  geom_text(data = tab_data, aes(x = 0.02, y = y, label = model_parse, colour = colour),
            hjust = 0, size = 4.4, parse = TRUE) +
  geom_text(data = tab_data, aes(x = 0.42, y = y, label = dev_expl),
            hjust = 0, size = 4.4, colour = "grey20") +
  geom_text(data = tab_data, aes(x = 0.74, y = y, label = crps),
            hjust = 0, size = 4.4, colour = "grey20") +
  annotate("text", x = 0.02, y = 0.15,
           label = "Explained deviance from model fit on training data to Sep 2021. CRPS averaged over target dates within plot range.",
           hjust = 0, vjust = 1, size = 3.3, colour = "grey50", fontface = "italic") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.6, 5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(4, 8, 6, 4))

combined <- p / p_table + plot_layout(heights = c(3.5, 1))
ggsave("stashed/rolling_forecast_plot.png", combined,
       width = 12, height = 6.5, dpi = 300, bg = "white")
message("Plot saved: stashed/rolling_forecast_plot.png")
