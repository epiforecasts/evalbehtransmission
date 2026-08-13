# Run the Chapter 1 pipeline in order.
# Each script reads what an earlier one wrote, so the order matters and
# data-processed/ is empty on a fresh clone.
# Every source is fetched remotely, so a fresh clone needs no local data.

ch1_pipeline <- c(
  "analysis/process_mobility.R",   # Google Mobility, per-country files
  "analysis/process_comix.R",      # CoMix contact matrices, slow to download
  "analysis/ch1_data.R",
  "analysis/ch1_covariates.R",
  "analysis/ch1_periods.R",
  "analysis/ch1_diagnostics.R",
  "analysis/ch1_rolling.R",        # Slowest, refitting every window
  "analysis/ch1_scoring.R",
  "analysis/ch1_window_plots.R",
  "analysis/ch1_forecast_plots.R", # Needs the scores written above
  "analysis/ch1_descriptive.R"
)

# Each in its own environment, so no script relies on another's objects
for (script in ch1_pipeline) {
  message("\n=== ", script, " ===")
  sys.source(script, envir = new.env(parent = globalenv()))
}

message("\nPipeline complete")
