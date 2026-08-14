# Chapter 1: explanatory outputs from the rolling windows, for the supplement
# Coefficient trajectories show how each behavioural slope moves across origins
# Deviance explained shows how much s(t) absorbs relative to the covariates

source("analysis/ch1_gam.R")
source("R/ch1_period_scale.R") # Shared period colours, so shading matches the other figures

library(ggplot2)

## Config ----------------------------------------------------------------------

ch1_window_config <- list(
  coef_path     = "data-processed/ch1_window_coefficients.csv",
  deviance_path = "data-processed/ch1_window_deviance.csv",
  periods_path  = "data-processed/ch1_periods.csv",
  output_dir    = "outputs/ch1",

  # Single-covariate models only, as the combined fit infers two slopes jointly
  single_covariate_models = c("contacts", "mobility")
)

## Load ------------------------------------------------------------------------

periods <- read_csv(ch1_window_config$periods_path, show_col_types = FALSE)

window_coefficients <- read_csv(ch1_window_config$coef_path, show_col_types = FALSE)
window_deviance     <- read_csv(ch1_window_config$deviance_path, show_col_types = FALSE)

# Origins run weekly, so this is the date range every window figure shares
# Padded so the points at the first and last origin do not sit against the panel edge
origin_range <- range(window_coefficients$origin) + c(-7, 7)

# Bands are clamped to that range, as ggplot drops any rect reaching outside it
period_bands <- periods |>
  group_by(period) |>
  summarise(start = min(date), end = max(date), .groups = "drop") |>
  arrange(start) |>
  mutate(period = factor(period, levels = ch1_period_levels)) |>
  filter(end >= origin_range[1], start <= origin_range[2]) |>
  mutate(start = pmax(start, origin_range[1]),
         end   = pmin(end,   origin_range[2]))

dir.create(ch1_window_config$output_dir, recursive = TRUE, showWarnings = FALSE)

## Coefficient trajectories ----------------------------------------------------
# One panel per covariate, taken from the model fitting that covariate alone
# Intervals are Wald (normal approximation, estimate +/- 1.96 se) from the mgcv parametric table

coefficient_trajectories <- window_coefficients |>
  filter(!used_smooth,
         model %in% ch1_window_config$single_covariate_models,
         term == model) |>
  mutate(model = factor(model, levels = ch1_window_config$single_covariate_models))

fig_coefficients_by_origin <- ggplot(coefficient_trajectories, aes(x = origin, y = estimate)) +
  geom_rect(data = period_bands, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = period),
            alpha = 0.08) +
  scale_fill_period() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower, ymax = upper), colour = "grey30") +
  geom_point(size = 1.2, colour = "grey10") +
  # Shared y-axis, as both covariates are z-scored and their slopes are directly comparable
  facet_wrap(~model, ncol = 1) +
  scale_x_date(limits = origin_range, expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.6))) + # Period key on one row, drawn more opaque than the bands
  labs(title = expression("Behavioural coefficients"~(beta[k])~"across forecast origins"),
       subtitle = "One 8-week training window per point, fitted without s(t), with 95% Wald intervals",
       x = "Forecast origin", y = expression("Coefficient on log"~R[t]~"("*beta[k]*")"), fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(file.path(ch1_window_config$output_dir, "fig_window_coefficients.png"),
       fig_coefficients_by_origin, width = 10, height = 6, dpi = 300, bg = "white")

## Deviance explained ----------------------------------------------------------
# Both smooth settings are shown, as s(t) alone absorbs most of what the covariates explain
# Baseline without s(t) is exactly zero, holding only an intercept and the offset

fig_deviance_explained <- window_deviance |>
  mutate(model       = factor(model, levels = names(ch1_models)),
         dev_expl    = 100 * dev_expl,
         smooth_term = factor(used_smooth, levels = c(FALSE, TRUE),
                              labels = c("without s(t)", "with s(t)"))) |>
  ggplot(aes(x = model, y = dev_expl, fill = smooth_term)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  scale_fill_manual(values = c("without s(t)" = "grey75", "with s(t)" = "steelblue")) +
  labs(title = "Deviance explained across the 8-week training windows",
       subtitle = "One box per model, with and without the time smooth",
       x = NULL, y = "Deviance explained (%)", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(ch1_window_config$output_dir, "fig_window_deviance.png"),
       fig_deviance_explained, width = 8, height = 5, dpi = 300, bg = "white")

## Checks ----------------------------------------------------------------------

coefficient_summary <- coefficient_trajectories |>
  group_by(model) |>
  summarise(origins       = n(),
            median_beta   = round(median(estimate), 3),
            min_beta      = round(min(estimate), 3),
            max_beta      = round(max(estimate), 3),
            pct_excl_zero = round(100 * mean(lower > 0 | upper < 0)),
            .groups = "drop")

cat("\n--- Coefficient summary, single-covariate models without s(t) ---\n")
print(as.data.frame(coefficient_summary), row.names = FALSE)

write_csv(coefficient_summary,
          file.path(ch1_window_config$output_dir, "table_window_coefficients.csv"))

message("Saved window plots to ", ch1_window_config$output_dir)
