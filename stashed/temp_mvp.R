# =============================================================================
# temp_mvp.R
#
# TEMPORARY ALL-IN-ONE MVP SCRIPT
#
# This script is an integrated prototype that runs the full analysis from raw
# data to fitted models in a single file for development and verification
# purposes.  It is NOT intended to be the final form of the code.
#
# Once confirmed correct, each section will be migrated to its proper home:
#   - Data loading/processing  -> data-raw/ or analysis/process_*.R
#   - Harmonisation            -> analysis/harmonise.R
#   - Model fitting            -> analysis/model_rtglm.R
#
# Do not run this script as part of the main pipeline (analysis/run_pipeline.R).
# =============================================================================

library(socialmixr)
library(data.table)
library(tidyverse)


# -----------------------------------------------------------------------------
# SECTION 1 — CoMix: load contact survey data and compute dominant eigenvalue
# -----------------------------------------------------------------------------
# For each survey round, build the NGM: diag(s_vec) %*% C(t) %*% diag(i_vec),
# where C(t) is the age-stratified contact matrix, s_vec are age-specific
# susceptibilities, and i_vec are age-specific infectiousness.
# MVP placeholder: s_vec = i_vec = 1 (identity scaling), so NGM = C(t).
# Next step: replace with antibody-informed susceptibilities (ONS CIS) and
# age-specific infectiousness following Munday et al. 2023.
# Output: data frame lambda1_waves with columns wave, date, lambda1.

comix_raw <- tryCatch(
  socialmixr::get_survey("https://doi.org/10.5281/zenodo.13684044"),
  error = function(e) {
    message("Remote fetch failed, falling back to local CSVs: ", e$message)
    socialmixr::survey(
      participants = data.table::fread("data-raw/CoMix/participants.csv"),
      contacts     = data.table::fread("data-raw/CoMix/contacts.csv")
    )
  }
)

# Inspect column names before assuming anything
cat("Participant columns:\n"); print(names(comix_raw$participants))
cat("Contact columns:\n");     print(names(comix_raw$contacts))

# Identify the participant ID and wave columns from what is actually present
part_id_col <- intersect(c("part_id", "participant_id"), names(comix_raw$participants))[1]
# Use survey_round (not wave): wave is a recruitment cohort label spanning months,
# survey_round is the actual per-round time window (~weekly)
wave_col    <- intersect(c("survey_round", "wave"), names(comix_raw$participants))[1]
if (is.na(wave_col)) stop("Cannot identify survey round column: expected 'survey_round' or 'wave' in participants")

# Standardise participant ID column to 'part_id' in both data frames so
# socialmixr::survey() can find it by its expected name
names(comix_raw$participants)[names(comix_raw$participants) == part_id_col] <- "part_id"
names(comix_raw$contacts)[names(comix_raw$contacts)    == part_id_col] <- "part_id"

# Cap contacts at 50 per participant to reduce leverage from outliers
# (wave is only in participants, not contacts — group by participant ID only)
# Base R only: avoids data.table j-expression scoping errors entirely
contacts_df <- as.data.frame(comix_raw$contacts)
contacts_df <- do.call(rbind, lapply(
  split(contacts_df, contacts_df[["part_id"]]),
  function(x) x[seq_len(min(nrow(x), 50)), ]
))
comix_raw$contacts <- contacts_df

uk_pop <- data.frame(
  lower.age.limit = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75),
  population      = c(3950000, 4100000, 3900000, 3750000, 4200000, 4300000,
                      4500000, 4700000, 4600000, 4400000, 4500000, 4200000,
                      3800000, 3200000, 2800000, 3100000)
)

waves <- sort(unique(comix_raw$participants[[wave_col]]))

lambda1_waves <- lapply(waves, function(w) {

  # Subset participants to this wave, then contacts by matching participant IDs
  p <- comix_raw$participants[comix_raw$participants[[wave_col]] == w, ]
  wave_ids <- p[["part_id"]]
  c <- comix_raw$contacts[comix_raw$contacts[["part_id"]] %in% wave_ids, ]
  wave_survey <- socialmixr::survey(participants = p, contacts = c)

  # Representative date: median survey date minus 1 day (contacts = day prior)
  # lubridate::as_date() handles non-standard formats more robustly than as.Date()
  date_col <- intersect(c("date", "survey_date", "sday_id"), names(p))[1]
  wave_date <- median(lubridate::as_date(p[[date_col]]), na.rm = TRUE) - 1

  lambda1 <- tryCatch({
    cm <- socialmixr::contact_matrix(
      wave_survey,
      age_limits      = c(0, 18, 30, 45, 60, 75),
      symmetric       = TRUE,
      weigh_dayofweek = TRUE,
      survey_pop      = uk_pop
    )$matrix

    if (is.null(cm) || all(is.na(cm))) stop("all-NA matrix")

    # MVP: uniform s_vec and i_vec (= 1); NGM reduces to C itself
    # Next step: scale by ONS seroprevalence following Munday et al. 2023
    NGM <- diag(rep(1, nrow(cm))) %*% cm %*% diag(rep(1, nrow(cm)))
    max(Re(eigen(NGM, only.values = TRUE)$values))
  }, error = function(e) {
    message("Wave ", w, " skipped: ", e$message)
    NA_real_
  })

  data.frame(wave = w, date = wave_date, lambda1 = lambda1)
})

lambda1_waves <- do.call(rbind, lambda1_waves)
print(lambda1_waves)


# -----------------------------------------------------------------------------
# SECTION 2 — Extend CoMix eigenvalue to daily series; compare with ONS Rt
# -----------------------------------------------------------------------------
# Extend wave-level lambda1 to a daily step function (value held constant from
# each survey date until the next).  Load rt_national.csv (from process_ons.R)
# and plot both series on the same time axis.

library(zoo)
library(patchwork)

# 1. Daily date sequence over the full span of lambda1_waves
daily_dates <- data.frame(
  date = seq(min(lambda1_waves$date, na.rm = TRUE),
             max(lambda1_waves$date, na.rm = TRUE),
             by = "day")
)

# 2. Step-function extension: left-join brings lambda1 onto survey dates only;
#    na.locf carries each value forward to the next survey date.
#    na.rm = FALSE preserves leading NAs (early sparse rounds — do not impute).
lambda1_daily <- daily_dates |>
  left_join(lambda1_waves |> select(date, lambda1), by = "date") |>
  mutate(lambda1 = zoo::na.locf(lambda1, na.rm = FALSE))

cat("\nlambda1_daily (first 10):\n"); print(head(lambda1_daily, 10))
cat("\nlambda1_daily (last 10):\n");  print(tail(lambda1_daily, 10))

# 3. Load Rt following the same extraction used in process_ons.R:
#    filter ons_estimates to name == "R", take q50 as the daily median.
#    Source process_ons.R to populate ons_estimates if not already present.
if (!exists("ons_estimates")) {
  source(here::here("analysis", "process_ons.R"))
}
ons_estimates <- read_csv("data-raw/ONS/estimates_national.csv") |>
  filter(variable == "England")
rt_national <- ons_estimates |>
  filter(name == "R") |>
  select(date, Rt = q50) |>
  arrange(date)

# 4. Plot — restrict both series to their overlapping date range
plot_start <- max(min(lambda1_daily$date, na.rm = TRUE),
                  min(rt_national$date,   na.rm = TRUE))
plot_end   <- min(max(lambda1_daily$date, na.rm = TRUE),
                  max(rt_national$date,   na.rm = TRUE))

# Dates of skipped rounds (NA lambda1) — shown as vertical dashed lines
na_dates <- lambda1_waves$date[is.na(lambda1_waves$lambda1)]

p_lambda1 <- ggplot(
    lambda1_daily |> filter(date >= plot_start, date <= plot_end),
    aes(x = date, y = lambda1)
  ) +
  geom_step(na.rm = TRUE) +
  geom_vline(xintercept = na_dates, linetype = "dashed", colour = "grey60") +
  labs(y = expression(lambda[1]), x = NULL,
       title = "CoMix dominant eigenvalue (daily step function, MVP: NGM = C(t))") +
  theme_minimal()

p_rt <- ggplot(
    rt_national |> filter(date >= plot_start, date <= plot_end),
    aes(x = date, y = Rt)
  ) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  labs(y = expression(R[t]), x = "Date", title = "ONS Rt (daily median, inc2prev)") +
  theme_minimal()

p_combined <- p_lambda1 / p_rt
ggsave("stashed/stash_plots/lambda1_vs_Rt.png", p_combined, width = 10, height = 6, dpi = 150)
cat("Plot saved to stashed/stash_plots/lambda1_vs_Rt.png\n")


# -----------------------------------------------------------------------------
# SECTION 3 — Google Mobility: load retail & recreation stream
# -----------------------------------------------------------------------------
# Load data-processed/google_mobility_UK.csv (from process_mobility.R) and
# retain date and retail_and_recreation_percent_change_from_baseline.
# No aggregation needed — data are already daily.

mobility_raw <- read_csv("data-processed/google_mobility_UK.csv",
                         show_col_types = FALSE)

cat("Mobility columns:\n"); print(names(mobility_raw))

# process_mobility.R already filtered to national level (sub_region_1 == NA),
# so no further geographic filtering is needed here.
mobility_daily <- mobility_raw |>
  select(date,
         mobility_retail = retail_and_recreation_percent_change_from_baseline)

cat("\nmobility_daily (first 5):\n"); print(head(mobility_daily, 5))
cat("\nmobility_daily (last 5):\n");  print(tail(mobility_daily, 5))
cat("\nDate range:", as.character(range(mobility_daily$date, na.rm = TRUE)), "\n")


# -----------------------------------------------------------------------------
# SECTION 4 — ONS: load daily incidence and Rt
# -----------------------------------------------------------------------------
# Load data-raw/ONS/estimates_national.csv, filter to variable == "England".
# Following process_ons.R: extract incidence (name == "infections", q50) and
# Rt (name == "R", q50).

# ons_estimates is already in the environment if section (2) ran first;
# source process_ons.R only if it is missing (e.g. running section (4) alone)
if (!exists("ons_estimates")) {
  source(here::here("analysis", "process_ons.R"))
}

incidence_national <- ons_estimates |>
  filter(name == "infections") |>
  select(date, incidence = q50) |>
  arrange(date)

rt_national <- ons_estimates |>
  filter(name == "R") |>
  select(date, Rt = q50) |>
  arrange(date)

cat("incidence_national (first 5):\n"); print(head(incidence_national, 5))
cat("\nincidence_national (last 5):\n");  print(tail(incidence_national, 5))
cat("\nIncidence date range:", as.character(range(incidence_national$date)), "\n")

cat("\nrt_national (first 5):\n"); print(head(rt_national, 5))
cat("\nrt_national (last 5):\n");  print(tail(rt_national, 5))
cat("\nRt date range:", as.character(range(rt_national$date)), "\n")


# -----------------------------------------------------------------------------
# SECTION 5 — Harmonise into a single daily data frame
# -----------------------------------------------------------------------------
# Left-join mobility_daily, lambda1_daily, and rt_national onto incidence_national
# (ONS date spine).  Z-score behavioural covariates using non-NA rows only.
# Final columns: date, incidence, Rt, lambda1, lambda1_std,
#                mobility_retail, mobility_retail_std.

# Use incidence_national as the date spine so the range is driven by ONS data.
# Left joins mean rows outside the CoMix or mobility windows get NAs naturally.
dat <- incidence_national |>
  left_join(rt_national,       by = "date") |>
  left_join(lambda1_daily,     by = "date") |>
  left_join(mobility_daily,    by = "date")

# Z-score behavioural covariates using only non-NA rows for mean/SD.
# Leading lambda1 NAs (rounds 1-6, strict lockdown) are left as-is — do not
# backfill, as those contacts are genuinely unobserved and structurally
# different from the first measured wave.
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

dat <- dat |>
  mutate(
    lambda1_std          = z_score(lambda1),
    mobility_retail_std  = z_score(mobility_retail)
  ) |>
  select(date, incidence, Rt, lambda1, lambda1_std,
         mobility_retail, mobility_retail_std)

# Summary diagnostics
cat("Date range of dat:", as.character(range(dat$date)), "\n")
cat("Rows with non-NA incidence:           ", sum(!is.na(dat$incidence)),          "\n")
cat("Rows with non-NA lambda1_std:         ", sum(!is.na(dat$lambda1_std)),         "\n")
cat("Rows with non-NA mobility_retail_std: ", sum(!is.na(dat$mobility_retail_std)), "\n")
cat("Complete cases (incidence + lambda1_std + mobility_retail_std):",
    sum(complete.cases(dat[, c("incidence", "lambda1_std", "mobility_retail_std")])), "\n")

saveRDS(dat, "stashed/stash_models/harmonised_daily.rds")
write_csv(dat, "stashed/harmonised_daily.csv")
cat("Saved to stashed/stash_models/harmonised_daily.rds and stashed/harmonised_daily.csv\n")


# -----------------------------------------------------------------------------
# SECTION 6 — Plot each harmonised variable against time
# -----------------------------------------------------------------------------
# Four separate panels (one each for incidence, Rt, lambda1, mobility_retail)
# plotted against time — no dual axes, unstandardised scales.

if (!exists("dat")) dat <- readRDS("stashed/stash_models/harmonised_daily.rds")

# Date range where lambda1 is NA (early rounds 1-6) — shaded as a grey rect
lambda1_na_start <- min(dat$date)
lambda1_na_end   <- dat$date[which(!is.na(dat$lambda1))[1] - 1L]

# Helper: grey shading + clean theme applied to every panel
na_rect <- annotate("rect",
  xmin = lambda1_na_start, xmax = lambda1_na_end,
  ymin = -Inf, ymax = Inf,
  fill = "grey85", alpha = 0.6
)

p_inc <- ggplot(dat, aes(x = date, y = incidence)) +
  na_rect +
  geom_line(na.rm = TRUE) +
  labs(title = "Daily incidence (ONS median)", x = NULL, y = "Infections / day") +
  theme_minimal()

p_rt <- ggplot(dat, aes(x = date, y = Rt)) +
  na_rect +
  geom_line(na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(title = "Rt (ONS median)", x = NULL, y = expression(R[t])) +
  theme_minimal()

p_l1 <- ggplot(dat, aes(x = date, y = lambda1)) +
  na_rect +
  geom_step(na.rm = TRUE) +
  labs(title = expression("CoMix dominant eigenvalue "*(lambda[1])*", unstandardised"),
       x = NULL, y = expression(lambda[1])) +
  theme_minimal()

p_mob <- ggplot(dat, aes(x = date, y = mobility_retail)) +
  na_rect +
  geom_line(na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "Google Mobility — retail & recreation (% change from baseline)",
       x = "Date", y = "% change") +
  theme_minimal()

p_overview <- p_inc / p_rt / p_l1 / p_mob
ggsave("stashed/stash_plots/harmonised_daily_overview.png", p_overview,
       width = 10, height = 8, dpi = 150)
cat("Plot saved to stashed/stash_plots/harmonised_daily_overview.png\n")


# -----------------------------------------------------------------------------
# SECTION 7 — Baseline renewal GLM (no behavioural covariates)
# -----------------------------------------------------------------------------
# Poisson GLM following Nouvellet et al.:
#   log(E[I_t]) = alpha_0 + log(Lambda_t)
# where Lambda_t = sum_s w_s * I_{t-s} (incidence convolved with generation
# interval) enters as an offset.  Serves as the benchmark model.

if (!exists("dat")) dat <- readRDS("stashed/stash_models/harmonised_daily.rds")

# 2. Discretise gamma generation interval (mean = 5.5 days, SD = 2.1 days)
#    to get daily weights w_s for lags s = 1 ... 14
gi_mean <- 5.5; gi_sd <- 2.1
gi_shape <- (gi_mean / gi_sd)^2
gi_rate  <- gi_mean / gi_sd^2
max_lag  <- 14L
# pgamma differences give probability mass in each 1-day bin
w_raw <- diff(pgamma(0:max_lag, shape = gi_shape, rate = gi_rate))
w     <- w_raw / sum(w_raw)          # normalise to sum = 1

# 3. Compute Lambda_t as weighted lag-convolution of incidence.
#    embed() creates a matrix where row i holds dat$incidence[i], [i-1], ..., [i-max_lag].
#    Column 1 is the current value; columns 2:(max_lag+1) are lags 1:max_lag.
#    Drop col 1, multiply each row by w, then prepend max_lag NAs for the
#    rows that lack a full lag window.  Fully vectorised — no loop.
n          <- nrow(dat)
lag_matrix <- embed(dat$incidence, max_lag + 1L)[, -1, drop = FALSE]
dat$Lambda_t <- c(rep(NA_real_, max_lag), as.vector(lag_matrix %*% w))
cat("Lambda_t: NAs =", sum(is.na(dat$Lambda_t)),
    "| non-NAs =", sum(!is.na(dat$Lambda_t)), "\n")

# 4. Prepare model data frame
dat_model <- dat |>
  slice(-(1:max_lag)) |>                       # drop rows with no lag history
  filter(!is.na(incidence), !is.na(Lambda_t),  # require both outcome and offset
         Lambda_t > 0) |>                       # log offset undefined at 0
  mutate(t = seq_len(n()))                      # numeric time index for GAM smooth

cat("dat_model rows:", nrow(dat_model), "\n")
cat("dat_model date range:", as.character(range(dat_model$date)), "\n")

# 5a. Model A — GLM: single constant log-transmission rate
#     Fitted Rt is flat across 2020-2022 — unrealistic, but illustrates
#     why time-varying structure (Model B) is needed
fit_baseline_glm <- glm(
  incidence ~ offset(log(Lambda_t)),
  family = poisson(),
  data   = dat_model
)
cat("\nModel A (GLM) AIC:", AIC(fit_baseline_glm), "\n")
cat("Model A intercept (= log implied constant Rt):",
    coef(fit_baseline_glm)[["(Intercept)"]], "\n")
cat("Implied constant Rt:", exp(coef(fit_baseline_glm)[["(Intercept)"]]), "\n")

# 5b. Model B — GAM: s(t) captures residual time-varying transmission.
#     k = 40 allows flexibility across ~2 years without overfitting;
#     this is a tuning choice worth revisiting in the final analysis.
fit_baseline_gam <- mgcv::gam(
  incidence ~ s(t, k = 40) + offset(log(Lambda_t)),
  family = poisson(),
  data   = dat_model
)
cat("\nModel B (GAM) AIC:", AIC(fit_baseline_gam), "\n")
print(summary(fit_baseline_gam))

# 6. Plot observed vs fitted from both models
dat_model$fitted_glm <- fitted(fit_baseline_glm)
dat_model$fitted_gam <- fitted(fit_baseline_gam)

p_fit <- ggplot(dat_model, aes(x = date)) +
  geom_line(aes(y = incidence),   colour = "grey60",  linewidth = 0.4) +
  geom_line(aes(y = fitted_glm),  colour = "steelblue", linewidth = 0.8,
            linetype = "dashed") +
  geom_line(aes(y = fitted_gam),  colour = "firebrick", linewidth = 0.8) +
  labs(title = "Baseline models: observed vs fitted incidence",
       subtitle = "Grey = observed  |  Blue dashed = GLM (constant Rt)  |  Red = GAM (time-smooth)",
       x = "Date", y = "Infections / day") +
  theme_minimal()

ggsave("stashed/stash_plots/fit_baseline_glm_vs_gam.png", p_fit,
       width = 10, height = 4, dpi = 150)
cat("Plot saved to stashed/stash_plots/fit_baseline_glm_vs_gam.png\n")

# 7. Save model objects
saveRDS(fit_baseline_glm, "stashed/stash_models/fit_baseline_glm.rds")
saveRDS(fit_baseline_gam, "stashed/stash_models/fit_baseline_gam.rds")
cat("Models saved.\n")


# -----------------------------------------------------------------------------
# SECTION 8 — CoMix model: add lambda1_std as covariate
# -----------------------------------------------------------------------------
# Extend baseline with:
#   log(E[I_t]) = alpha_0 + beta_comix * lambda1_std_{t-l} + log(Lambda_t)

# 1. Load harmonised data and ensure Lambda_t is present
if (!exists("dat")) dat <- readRDS("stashed/stash_models/harmonised_daily.rds")

# Recompute Lambda_t if running section 8 standalone (section 7 not yet run)
if (!"Lambda_t" %in% names(dat)) {
  gi_mean  <- 5.5; gi_sd <- 2.1; max_lag <- 14L
  gi_shape <- (gi_mean / gi_sd)^2
  gi_rate  <- gi_mean / gi_sd^2
  w        <- diff(pgamma(0:max_lag, shape = gi_shape, rate = gi_rate))
  w        <- w / sum(w)
  lag_matrix   <- embed(dat$incidence, max_lag + 1L)[, -1, drop = FALSE]
  dat$Lambda_t <- c(rep(NA_real_, max_lag), as.vector(lag_matrix %*% w))
}

# Rebuild dat_model if section 7 did not run in this session
if (!exists("dat_model")) {
  dat_model <- dat |>
    slice(-(1:14L)) |>
    filter(!is.na(incidence), !is.na(Lambda_t), Lambda_t > 0) |>
    mutate(t = seq_len(n()))
}

# 2. Restrict to rows where lambda1_std is observed (drops early lockdown rounds)
dat_model_comix <- dat_model |>
  filter(!is.na(lambda1_std))
cat("dat_model_comix rows retained:", nrow(dat_model_comix), "\n")

# 3. Fit CoMix GAM
#    lambda1_std enters linearly: log-linear effect on Rt is the natural
#    functional form (eigenvalue scales multiplicatively with transmission).
#    s(t) absorbs residual time-varying transmission not explained by lambda1_std.
fit_comix <- mgcv::gam(
  incidence ~ s(t, k = 40) + lambda1_std + offset(log(Lambda_t)),
  family = poisson(),
  data   = dat_model_comix
)

# 4. Diagnostics
cat("\nCoMix GAM AIC:", AIC(fit_comix), "\n")

# exp(beta) = multiplicative effect on Rt per 1-SD increase in lambda1
beta_comix <- coef(fit_comix)["lambda1_std"]
se_comix   <- sqrt(diag(vcov(fit_comix)))["lambda1_std"]
cat(sprintf("lambda1_std: exp(coef) = %.3f  (95%% CI: %.3f -- %.3f)\n",
            exp(beta_comix),
            exp(beta_comix - 1.96 * se_comix),
            exp(beta_comix + 1.96 * se_comix)))

# Delta-AIC vs baseline GAM — approximate: the two models differ in row set
# (CoMix model drops early rounds where lambda1 is NA)
if (exists("fit_baseline_gam")) {
  cat(sprintf("Delta-AIC vs baseline GAM: %.1f  (NOTE: approximate — different row sets)\n",
              AIC(fit_comix) - AIC(fit_baseline_gam)))
} else {
  cat("fit_baseline_gam not found in environment — load stashed/stash_models/fit_baseline_gam.rds to compare\n")
}

# 5. Plot observed vs fitted
dat_model_comix$fitted_comix <- fitted(fit_comix)

p_fit_comix <- ggplot(dat_model_comix, aes(x = date)) +
  geom_line(aes(y = incidence),    colour = "grey60",    linewidth = 0.4) +
  geom_line(aes(y = fitted_comix), colour = "darkorange", linewidth = 0.8) +
  labs(title = "CoMix GAM: observed vs fitted incidence",
       subtitle = "Grey = observed  |  Orange = fitted (s(t) + lambda1_std)",
       x = "Date", y = "Infections / day") +
  theme_minimal()

ggsave("stashed/stash_plots/fit_comix.png", p_fit_comix, width = 10, height = 4, dpi = 150)
cat("Plot saved to stashed/stash_plots/fit_comix.png\n")

# 6. Save model object
saveRDS(fit_comix, "stashed/stash_models/fit_comix.rds")
cat("Model saved to stashed/stash_models/fit_comix.rds\n")


# -----------------------------------------------------------------------------
# SECTION 9 — Mobility model: add mobility_retail_std as covariate
# -----------------------------------------------------------------------------
# Extend baseline with:
#   log(E[I_t]) = alpha_0 + beta_mob * mobility_retail_std_{t-l} + log(Lambda_t)

# 1. Load harmonised data and ensure Lambda_t is present
if (!exists("dat")) dat <- readRDS("stashed/stash_models/harmonised_daily.rds")

# Recompute Lambda_t if running section 9 standalone (section 7 not yet run)
if (!"Lambda_t" %in% names(dat)) {
  gi_mean  <- 5.5; gi_sd <- 2.1; max_lag <- 14L
  gi_shape <- (gi_mean / gi_sd)^2
  gi_rate  <- gi_mean / gi_sd^2
  w        <- diff(pgamma(0:max_lag, shape = gi_shape, rate = gi_rate))
  w        <- w / sum(w)
  lag_matrix   <- embed(dat$incidence, max_lag + 1L)[, -1, drop = FALSE]
  dat$Lambda_t <- c(rep(NA_real_, max_lag), as.vector(lag_matrix %*% w))
}

# Rebuild dat_model if section 7 did not run in this session
if (!exists("dat_model")) {
  dat_model <- dat |>
    slice(-(1:14L)) |>
    filter(!is.na(incidence), !is.na(Lambda_t), Lambda_t > 0) |>
    mutate(t = seq_len(n()))
}

# 2. Restrict to rows where mobility_retail_std is observed
dat_model_mobility <- dat_model |>
  filter(!is.na(mobility_retail_std))
cat("dat_model_mobility rows retained:", nrow(dat_model_mobility), "\n")

# 3. Fit mobility GAM
#    mobility_retail_std enters linearly on log(Rt).
#    NOTE: mobility_retail is negatively coded (percent change from baseline,
#    so stricter lockdown = more negative values). We therefore expect a
#    positive coefficient: higher retail activity -> higher Rt.
#    s(t) absorbs residual time-varying transmission not explained by mobility.
fit_mobility <- mgcv::gam(
  incidence ~ s(t, k = 40) + mobility_retail_std + offset(log(Lambda_t)),
  family = poisson(),
  data   = dat_model_mobility
)

# 4. Diagnostics
cat("\nMobility GAM AIC:", AIC(fit_mobility), "\n")

# exp(beta) = multiplicative effect on Rt per 1-SD increase in mobility_retail_std
beta_mob <- coef(fit_mobility)["mobility_retail_std"]
se_mob   <- sqrt(diag(vcov(fit_mobility)))["mobility_retail_std"]
cat(sprintf("mobility_retail_std: exp(coef) = %.3f  (95%% CI: %.3f -- %.3f)\n",
            exp(beta_mob),
            exp(beta_mob - 1.96 * se_mob),
            exp(beta_mob + 1.96 * se_mob)))

# Delta-AIC vs baseline GAM — approximate: models may differ in row set
# if mobility has leading/trailing NAs not present in the baseline sample
if (exists("fit_baseline_gam")) {
  cat(sprintf("Delta-AIC vs baseline GAM: %.1f  (NOTE: approximate — different row sets)\n",
              AIC(fit_mobility) - AIC(fit_baseline_gam)))
} else {
  cat("fit_baseline_gam not found in environment — load stashed/stash_models/fit_baseline_gam.rds to compare\n")
}

# 5. Plot observed vs fitted
dat_model_mobility$fitted_mobility <- fitted(fit_mobility)

p_fit_mobility <- ggplot(dat_model_mobility, aes(x = date)) +
  geom_line(aes(y = incidence),       colour = "grey60",   linewidth = 0.4) +
  geom_line(aes(y = fitted_mobility), colour = "steelblue", linewidth = 0.8) +
  labs(title = "Mobility GAM: observed vs fitted incidence",
       subtitle = "Grey = observed  |  Blue = fitted (s(t) + mobility_retail_std)",
       x = "Date", y = "Infections / day") +
  theme_minimal()

ggsave("stashed/stash_plots/fit_mobility.png", p_fit_mobility, width = 10, height = 4, dpi = 150)
cat("Plot saved to stashed/stash_plots/fit_mobility.png\n")

# 6. Save model object
saveRDS(fit_mobility, "stashed/stash_models/fit_mobility.rds")
cat("Model saved to stashed/stash_models/fit_mobility.rds\n")


# -----------------------------------------------------------------------------
# SECTION 10 — Model comparison: baseline GAM, CoMix GAM, mobility GAM
# -----------------------------------------------------------------------------
# Compare all three fitted renewal GAMs on AIC and fitted-vs-observed plots.
# NOTE: delta-AIC comparisons across models fitted on different row sets
# (CoMix drops early NAs; mobility may have trailing NAs) are approximate
# and should be interpreted cautiously — they are informative only within
# the overlapping period.

# 1. Load model objects if not already in the environment
if (!exists("fit_baseline_gam")) fit_baseline_gam <- readRDS("stashed/stash_models/fit_baseline_gam.rds")
if (!exists("fit_comix"))        fit_comix        <- readRDS("stashed/stash_models/fit_comix.rds")
if (!exists("fit_mobility"))     fit_mobility     <- readRDS("stashed/stash_models/fit_mobility.rds")

# nobs() gives the number of rows each model was actually fitted on
aic_table <- data.frame(
  model     = c("baseline GAM", "CoMix GAM", "mobility GAM"),
  n_obs     = c(nobs(fit_baseline_gam), nobs(fit_comix), nobs(fit_mobility)),
  AIC       = c(AIC(fit_baseline_gam),  AIC(fit_comix),  AIC(fit_mobility))
)
aic_table$delta_AIC <- aic_table$AIC - AIC(fit_baseline_gam)

# 2. Print comparison table
cat("\nModel comparison (delta_AIC relative to baseline GAM):\n")
print(aic_table, row.names = FALSE, digits = 4)
cat("(!) delta_AIC comparisons across models with different row sets are approximate\n")

# 3. Fitted-vs-observed plot — three panels, shared y-axis scale
#    Ensure dat_model, dat_model_comix, dat_model_mobility are present with
#    a date column.  gam()$model stores only the model matrix (no date), so
#    we reconstruct from harmonised_daily.rds if sections 7-9 did not run.
if (!exists("dat_model") || !"date" %in% names(dat_model)) {
  if (!exists("dat")) dat <- readRDS("stashed/stash_models/harmonised_daily.rds")
  if (!"Lambda_t" %in% names(dat)) {
    gi_mean <- 5.5; gi_sd <- 2.1; max_lag <- 14L
    w <- diff(pgamma(0:max_lag, shape = (gi_mean/gi_sd)^2, rate = gi_mean/gi_sd^2))
    w <- w / sum(w)
    dat$Lambda_t <- c(rep(NA_real_, max_lag),
                      as.vector(embed(dat$incidence, max_lag + 1L)[, -1, drop = FALSE] %*% w))
  }
  dat_model <- dat |>
    slice(-(1:14L)) |>
    filter(!is.na(incidence), !is.na(Lambda_t), Lambda_t > 0) |>
    mutate(t = seq_len(n()))
}
if (!exists("dat_model_comix") || !"date" %in% names(dat_model_comix)) {
  dat_model_comix <- dat_model |> filter(!is.na(lambda1_std))
}
if (!exists("dat_model_mobility") || !"date" %in% names(dat_model_mobility)) {
  dat_model_mobility <- dat_model |> filter(!is.na(mobility_retail_std))
}

fitted_baseline <- dat_model |>
  select(date, incidence) |>
  mutate(fitted = fitted(fit_baseline_gam), model = "Baseline GAM")

fitted_comix <- dat_model_comix |>
  select(date, incidence) |>
  mutate(fitted = fitted(fit_comix), model = "CoMix GAM")

fitted_mobility <- dat_model_mobility |>
  select(date, incidence) |>
  mutate(fitted = fitted(fit_mobility), model = "Mobility GAM")

all_fitted <- bind_rows(fitted_baseline, fitted_comix, fitted_mobility) |>
  mutate(model = factor(model, levels = c("Baseline GAM", "CoMix GAM", "Mobility GAM")))

# Shared y-axis limits across panels so fits are visually comparable
y_max <- max(all_fitted$incidence, all_fitted$fitted, na.rm = TRUE)

p_comparison <- ggplot(all_fitted, aes(x = date)) +
  geom_line(aes(y = incidence), colour = "grey60",   linewidth = 0.4) +
  geom_line(aes(y = fitted),    colour = "firebrick", linewidth = 0.7) +
  facet_wrap(~ model, ncol = 1) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(title = "Renewal GAMs: observed vs fitted incidence",
       subtitle = "Grey = observed  |  Red = fitted",
       x = "Date", y = "Infections / day") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave("stashed/stash_plots/fit_comparison.png", p_comparison,
       width = 12, height = 8, dpi = 150)
cat("Plot saved to stashed/stash_plots/fit_comparison.png\n")

# 4. Pearson residuals plot — three panels, reference line at zero
#    Pearson residual = (observed - fitted) / sqrt(fitted) under Poisson
all_fitted <- all_fitted |>
  mutate(pearson_resid = (incidence - fitted) / sqrt(fitted))

p_residuals <- ggplot(all_fitted, aes(x = date, y = pearson_resid)) +
  geom_line(colour = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ model, ncol = 1) +
  labs(title = "Renewal GAMs: Pearson residuals over time",
       subtitle = "Systematic departures flag periods where the model fits poorly",
       x = "Date", y = "Pearson residual") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave("stashed/stash_plots/residuals_comparison.png", p_residuals,
       width = 12, height = 8, dpi = 150)
cat("Plot saved to stashed/stash_plots/residuals_comparison.png\n")
# =============================================================================
