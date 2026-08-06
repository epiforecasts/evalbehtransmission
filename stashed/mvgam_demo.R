# =============================================================================
# mvgam_demo.R  (stashed — temporary exploration only)
#
# Demonstrates the MVGAM approach as a principled forecasting replacement for
# the renewal GAM s(t) spline used in temp_mvp.R.
#
# Key idea: replace s(t) with a latent random walk (RW) on log(Rt). This gives:
#   - A proper time-series process that forecasts forward naturally
#   - Widening uncertainty at longer horizons (vs fixed spline extrapolation)
#   - Full posterior distributions over Rt at every time point
#
# Models fitted:
#   (A) Baseline MVGAM  — RW trend only, no behavioural covariates
#   (B) Mobility MVGAM  — RW trend + composite mobility (retail + transit + work)
#
# Mobility choice: composite of retail & recreation, transit stations, and
# workplaces — the three indoor/away-from-home streams most plausibly linked
# to transmission-relevant contacts. Already on % change from baseline scale
# (approximately log-relative), so no further log transformation applied.
#
# Data resolution: aggregated to weekly. Daily resolution makes Stan
# sampling prohibitively slow; weekly (~130 obs) is feasible and natural
# for ONS CIS incidence.
#
# Outputs (stashed/stash_plots/):
#   mvgam_latent_rt.png        — latent Rt posterior vs observed ONS Rt
#   mvgam_insample_fit.png     — posterior mean fit vs observed, both models
#   mvgam_forecasts.png        — probabilistic forecasts with 50/95% intervals
#   mvgam_crps_horizon.png     — CRPS by forecast horizon (1-8 weeks)
#   mvgam_coef_mobility.png    — posterior distribution of mobility coefficient
#
# NOTE: mvgam uses Stan. Each model takes ~3-8 min on a modern laptop.
# Install if needed: install.packages("mvgam")
# =============================================================================

library(mvgam)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scoringutils)  # v2.2.0; crps_sample() confirmed available


# =============================================================================
# SECTION 1 — Load data from data-processed/ and compute daily Lambda_t
# =============================================================================
# Read directly from the canonical processed layer, not from stashed/ intermediates.

incidence_raw <- read_csv("data-processed/incidence_national.csv",
                           show_col_types = FALSE) |>
  rename(incidence = median) |>
  select(date, incidence)

rt_raw <- read_csv("data-processed/rt_national.csv",
                   show_col_types = FALSE) |>
  rename(Rt = median) |>
  select(date, Rt)

mob_raw <- read_csv("data-processed/google_mobility_UK.csv",
                    show_col_types = FALSE) |>
  select(
    date,
    retail    = retail_and_recreation_percent_change_from_baseline,
    transit   = transit_stations_percent_change_from_baseline,
    workplace = workplaces_percent_change_from_baseline
  )

dat <- incidence_raw |>
  left_join(rt_raw,  by = "date") |>
  left_join(mob_raw, by = "date") |>
  mutate(
    # Composite: mean of three indoor/non-household streams.
    # Parks excluded (weather-driven); residential excluded (rises during
    # lockdown but inversely related to transmission-relevant contacts).
    mob_composite = (retail + transit + workplace) / 3
  )

# Discretise gamma generation interval (mean = 5.5 days, SD = 2.1 days)
gi_mean  <- 5.5; gi_sd <- 2.1; max_lag <- 14L
gi_shape <- (gi_mean / gi_sd)^2
gi_rate  <- gi_mean / gi_sd^2
w        <- diff(pgamma(0:max_lag, shape = gi_shape, rate = gi_rate))
w        <- w / sum(w)

lag_mat      <- embed(dat$incidence, max_lag + 1L)[, -1, drop = FALSE]
dat$Lambda_t <- c(rep(NA_real_, max_lag), as.vector(lag_mat %*% w))

cat(sprintf("Daily data: %d rows | %s to %s\n",
    nrow(dat), as.character(min(dat$date)), as.character(max(dat$date))))


# =============================================================================
# SECTION 2 — Aggregate to weekly
# =============================================================================
# Weekly incidence = sum of daily infections.
# Weekly Lambda_t  = sum of daily Lambda_t values (correct renewal offset:
#   E[I_week] = Rt * sum_day(Lambda_t_day) for approximately constant weekly Rt).
# Weekly Rt        = mean of daily ONS Rt (for visualisation only, not modelled).
# Weekly mobility  = mean of daily values.

dat_weekly <- dat |>
  mutate(year_iso = isoyear(date), week_iso = isoweek(date)) |>
  group_by(year_iso, week_iso) |>
  summarise(
    date          = min(date),
    incidence     = round(sum(incidence, na.rm = TRUE)),
    Lambda_t      = sum(Lambda_t, na.rm = TRUE),
    mob_composite = mean(mob_composite, na.rm = TRUE),
    Rt            = mean(Rt, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(date) |>
  filter(
    !is.na(incidence),     incidence > 0,
    !is.na(Lambda_t),      Lambda_t > 0,
    !is.na(mob_composite), is.finite(mob_composite)
  ) |>
  mutate(
    mob_std = (mob_composite - mean(mob_composite)) / sd(mob_composite),
    time    = seq_len(n()),
    series  = factor("national")
  )

cat(sprintf("Weekly data: %d weeks | %s to %s\n",
    nrow(dat_weekly), as.character(min(dat_weekly$date)), as.character(max(dat_weekly$date))))
cat(sprintf("Rt available for %d / %d weeks\n",
    sum(!is.na(dat_weekly$Rt)), nrow(dat_weekly)))


# =============================================================================
# SECTION 3 — Train / test split
# =============================================================================
# Single illustrative split: last 8 weeks as test set.
# Full rolling-origin evaluation is the correct approach for the main analysis
# but is too slow for this demo (requires refitting once per forecast origin).

n_test  <- 8L
n_train <- nrow(dat_weekly) - n_test

train_data <- dat_weekly[seq_len(n_train), ]
test_data  <- dat_weekly[(n_train + 1L):nrow(dat_weekly), ]

cat(sprintf("Train: %d weeks (%s to %s)\n",
    nrow(train_data), as.character(min(train_data$date)), as.character(max(train_data$date))))
cat(sprintf("Test:  %d weeks (%s to %s)\n",
    nrow(test_data), as.character(min(test_data$date)), as.character(max(test_data$date))))


# =============================================================================
# SECTION 4 — Fit baseline MVGAM (RW trend, no covariates)
# =============================================================================
# log E[I_t] = log(Lambda_t)  [offset]  +  eta_t
# eta_t = eta_{t-1} + epsilon_t,   epsilon ~ N(0, sigma^2)
#
# eta_t is the latent log(Rt). At forecast time it follows a random walk from
# its last estimated state, with uncertainty growing as sqrt(h) * sigma.
# This is the natural Bayesian baseline: no behavioural input, transmission
# evolves as a random walk.
#
# newdata = test_data is passed at fit time so MVGAM pre-computes forecasts
# during sampling (more efficient than calling forecast() separately).

cat("\nFitting baseline MVGAM (RW trend, no covariates)...\n")
cat("Stan sampling: expect 3-8 minutes.\n\n")

fit_baseline <- mvgam(
  formula     = incidence ~ 0 + offset(log(Lambda_t)),
  trend_model = RW(),
  family      = poisson(),
  data        = train_data,
  chains      = 1,
  iter        = 3000,
  warmup      = 1000,
  control     = list(max_treedepth = 14)
)

cat("\nBaseline MVGAM summary:\n")
print(summary(fit_baseline))
saveRDS(fit_baseline, "stashed/stash_models/fit_baseline_mvgam.rds")
cat("Saved: stashed/stash_models/fit_baseline_mvgam.rds\n")


# =============================================================================
# SECTION 5 — Fit mobility-augmented MVGAM
# =============================================================================
# log E[I_t] = beta * mob_std_t  +  log(Lambda_t)  +  eta_t
#
# beta: log-multiplicative effect of 1-SD increase in composite mobility on Rt.
# Expected direction: beta > 0 (more activity -> more contacts -> higher Rt).
# The RW trend eta_t absorbs residual Rt variation not explained by mobility
# (e.g., immunity changes, variant emergence). If mobility explains Rt well,
# the RW variance sigma^2 should be smaller than in the baseline.

cat("\nFitting mobility MVGAM (RW trend + composite mobility)...\n")
cat("Stan sampling: expect 3-8 minutes.\n\n")

fit_mobility <- mvgam(
  formula     = incidence ~ 0 + mob_std + offset(log(Lambda_t)),
  trend_model = RW(),
  family      = poisson(),
  data        = train_data,
  chains      = 1,
  iter        = 3000,
  warmup      = 1000,
  control     = list(max_treedepth = 14)
)

cat("\nMobility MVGAM summary:\n")
print(summary(fit_mobility))
saveRDS(fit_mobility, "stashed/stash_models/fit_mobility_mvgam.rds")
cat("Saved: stashed/stash_models/fit_mobility_mvgam.rds\n")

# To reload without refitting:
# fit_baseline <- readRDS("stashed/stash_models/fit_baseline_mvgam.rds")
# fit_mobility <- readRDS("stashed/stash_models/fit_mobility_mvgam.rds")


# =============================================================================
# SECTION 6 — Plot A: Latent Rt trajectory
# =============================================================================
# hindcast() returns an mvgam_forecast object. The hindcasts slot is a named
# list keyed by series name: hindcasts$national is an (n_draws x n_train) matrix.
# Apply over columns (dim 2) to get per-time-point posterior summaries.
# exp() converts from log(Rt) scale to Rt scale.

# Built-in trend plot — always reliable
cat("\nSaving built-in trend plots...\n")
png("stashed/stash_plots/mvgam_trend_baseline_builtin.png", width = 1400, height = 500, res = 150)
plot(fit_baseline, type = "trend", series = 1)
dev.off()

png("stashed/stash_plots/mvgam_trend_mobility_builtin.png", width = 1400, height = 500, res = 150)
plot(fit_mobility, type = "trend", series = 1)
dev.off()
cat("Saved: mvgam_trend_baseline_builtin.png, mvgam_trend_mobility_builtin.png\n")

# Custom ggplot overlaying both models with observed ONS Rt
hc_base <- hindcast(fit_baseline, type = "trend")  # log(Rt) scale
hc_mob  <- hindcast(fit_mobility, type = "trend")

summarise_trend <- function(hc, label, dates) {
  # hindcasts$national: (n_draws x n_train) matrix
  mat <- hc$hindcasts$national
  tibble(
    date    = dates,
    median  = exp(apply(mat, 2, median)),
    lo95    = exp(apply(mat, 2, quantile, 0.025)),
    hi95    = exp(apply(mat, 2, quantile, 0.975)),
    lo50    = exp(apply(mat, 2, quantile, 0.25)),
    hi50    = exp(apply(mat, 2, quantile, 0.75)),
    model   = label
  )
}

trend_df <- bind_rows(
  summarise_trend(hc_base, "Baseline (RW only)",  train_data$date),
  summarise_trend(hc_mob,  "Mobility-augmented",  train_data$date)
)

p_rt <- ggplot(trend_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = model), alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = model), alpha = 0.25) +
  geom_line(aes(y = median, colour = model), linewidth = 0.8) +
  geom_line(
    data    = train_data |> filter(!is.na(Rt)),
    mapping = aes(y = Rt),
    colour = "black", linewidth = 0.4, linetype = "dashed"
  ) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = c("Baseline (RW only)" = "steelblue",
                                  "Mobility-augmented"  = "firebrick")) +
  scale_fill_manual(  values = c("Baseline (RW only)" = "steelblue",
                                  "Mobility-augmented"  = "firebrick")) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    title    = "Latent Rt: MVGAM posterior (shaded) vs observed ONS Rt (dashed)",
    subtitle = "Bands = 50% and 95% posterior credible intervals",
    x = NULL, y = expression(R[t])
  ) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave("stashed/stash_plots/mvgam_latent_rt.png", p_rt, width = 12, height = 6, dpi = 150)
cat("Saved: mvgam_latent_rt.png\n")


# =============================================================================
# SECTION 7 — Plot B: In-sample fit
# =============================================================================
# hindcast() with type = "response" gives posterior predictive samples on the
# incidence scale. hindcasts$national is (n_draws x n_train).

hc_base_resp <- hindcast(fit_baseline, type = "response")
hc_mob_resp  <- hindcast(fit_mobility, type = "response")

summarise_fit <- function(hc, label, train_df) {
  mat <- hc$hindcasts$national
  tibble(
    date     = train_df$date,
    observed = train_df$incidence,
    fitted   = apply(mat, 2, median),
    lo95     = apply(mat, 2, quantile, 0.025),
    hi95     = apply(mat, 2, quantile, 0.975),
    model    = label
  )
}

insample_df <- bind_rows(
  summarise_fit(hc_base_resp, "Baseline (RW only)", train_data),
  summarise_fit(hc_mob_resp,  "Mobility-augmented", train_data)
)

p_insample <- ggplot(insample_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.15) +
  geom_line(aes(y = observed), colour = "grey50", linewidth = 0.4) +
  geom_line(aes(y = fitted),   colour = "firebrick", linewidth = 0.7) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    title    = "MVGAM: in-sample fit (posterior median + 95% CrI)",
    subtitle = "Grey = observed  |  Red = posterior median  |  Blue = 95% credible interval",
    x = NULL, y = "Weekly infections"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave("stashed/stash_plots/mvgam_insample_fit.png", p_insample, width = 12, height = 7, dpi = 150)
cat("Saved: mvgam_insample_fit.png\n")


# =============================================================================
# SECTION 8 — Plot C: Probabilistic forecasts
# =============================================================================
# forecast() was called at fit time (newdata = test_data passed to mvgam()).
# Retrieve forecast samples from the fitted object directly.
# forecasts$national: (n_draws x n_test) matrix.

# Built-in forecast plot requires test data passed at fit time (not available here).
# Forecasts are plotted via ggplot below using fc_base and fc_mob.

# Generate forecasts by passing test_data to forecast() after fitting
fc_base <- forecast(fit_baseline, newdata = test_data)
fc_mob  <- forecast(fit_mobility, newdata = test_data)

summarise_forecast <- function(fc_obj, label, test_dates) {
  mat <- fc_obj$forecasts$national  # (n_draws x n_test)
  tibble(
    date    = test_dates,
    horizon = seq_along(test_dates),
    median  = apply(mat, 2, median),
    lo50    = apply(mat, 2, quantile, 0.25),
    hi50    = apply(mat, 2, quantile, 0.75),
    lo95    = apply(mat, 2, quantile, 0.025),
    hi95    = apply(mat, 2, quantile, 0.975),
    model   = label
  )
}

fc_df   <- bind_rows(
  summarise_forecast(fc_base, "Baseline (RW only)", test_data$date),
  summarise_forecast(fc_mob,  "Mobility-augmented", test_data$date)
)
context <- train_data |> slice_tail(n = 12)

p_fc <- ggplot(fc_df, aes(x = date)) +
  geom_line(data = context, aes(y = incidence), colour = "grey50", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = model), alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = model), alpha = 0.25) +
  geom_line(aes(y = median, colour = model), linewidth = 0.9) +
  geom_point(data = test_data, aes(y = incidence), colour = "black", size = 2) +
  geom_vline(xintercept = min(test_data$date), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Baseline (RW only)" = "steelblue",
                                  "Mobility-augmented"  = "firebrick")) +
  scale_fill_manual(  values = c("Baseline (RW only)" = "steelblue",
                                  "Mobility-augmented"  = "firebrick")) +
  labs(
    title    = "8-week probabilistic forecast: baseline vs mobility-augmented MVGAM",
    subtitle = "Dashed = forecast origin  |  Dots = observed  |  Bands = 50% and 95% predictive intervals",
    x = NULL, y = "Weekly infections",
    colour = NULL, fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("stashed/stash_plots/mvgam_forecasts.png", p_fc, width = 12, height = 5, dpi = 150)
cat("Saved: mvgam_forecasts.png\n")


# =============================================================================
# SECTION 9 — Plot D: CRPS by forecast horizon
# =============================================================================
# crps_sample(observed, predicted) from scoringutils v2.
# predicted must be a numeric vector of posterior draws.

samps_base <- fc_base$forecasts$national  # (n_draws x n_test)
samps_mob  <- fc_mob$forecasts$national

crps_df <- map_dfr(seq_len(n_test), function(h) {
  bind_rows(
    tibble(horizon = h, model = "Baseline (RW only)",
           crps = crps_sample(observed  = test_data$incidence[h],
                               predicted = samps_base[, h])),
    tibble(horizon = h, model = "Mobility-augmented",
           crps = crps_sample(observed  = test_data$incidence[h],
                               predicted = samps_mob[, h]))
  )
})

cat("\nCRPS by horizon:\n")
print(pivot_wider(crps_df, names_from = model, values_from = crps))

p_crps <- ggplot(crps_df, aes(x = horizon, y = crps, colour = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Baseline (RW only)" = "steelblue",
                                  "Mobility-augmented"  = "firebrick")) +
  scale_x_continuous(breaks = seq_len(n_test),
                     labels = paste0(seq_len(n_test), "w")) +
  labs(
    title    = "CRPS by forecast horizon: baseline vs mobility-augmented",
    subtitle = "Lower = better  |  Red below blue = mobility adds forecast accuracy at that horizon",
    x = "Forecast horizon (weeks ahead)", y = "CRPS",
    colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("stashed/stash_plots/mvgam_crps_horizon.png", p_crps, width = 8, height = 5, dpi = 150)
cat("Saved: mvgam_crps_horizon.png\n")


# =============================================================================
# SECTION 10 — Plot E: Posterior of mobility coefficient
# =============================================================================
# as_draws_df.mvgam has a dispatch bug when there is no intercept.
# Bypass it: pull draws directly from the CmdStan object.
# Fixed effects are stored as parameter 'b'; with one covariate, b[1] = mob_std.

b_draws <- posterior::as_draws_df(
  fit_mobility$model_output$draws(variables = "b")
)
cat("\nFixed-effect draw columns:\n")
print(names(b_draws))

# b[1] = mob_std (only fixed effect in this model)
beta_samps <- b_draws[["b[1]"]]

cat(sprintf("\nMobility coefficient (beta):\n  Median exp(β) = %.3f  |  95%% CrI: [%.3f, %.3f]\n",
    exp(median(beta_samps)),
    exp(quantile(beta_samps, 0.025)),
    exp(quantile(beta_samps, 0.975))))

p_coef <- ggplot(data.frame(beta = beta_samps), aes(x = exp(beta))) +
  geom_histogram(bins = 60, fill = "steelblue", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 1,
             linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = exp(median(beta_samps)),
             colour = "firebrick", linewidth = 1) +
  labs(
    title    = "Posterior: multiplicative effect of composite mobility on Rt",
    subtitle = sprintf("Median exp(β) = %.3f  |  95%% CrI: [%.3f, %.3f]",
                       exp(median(beta_samps)),
                       exp(quantile(beta_samps, 0.025)),
                       exp(quantile(beta_samps, 0.975))),
    x = expression("exp(" * beta[mobility] * ")  — multiplicative effect on R"[t]),
    y = "Posterior samples"
  ) +
  theme_minimal()

ggsave("stashed/stash_plots/mvgam_coef_mobility.png", p_coef,
       width = 8, height = 4, dpi = 150)
cat("Saved: mvgam_coef_mobility.png\n")

cat("\nDone. All outputs in stashed/stash_plots/\n")