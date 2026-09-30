# Runs the Appendix A.1 scripts
# Refits GAMs on the Chapter 1 input, and A.1.5 also reads its stored forecasts and scores,
# so run analysis/run_pipeline.R first
# Writes only to analysis/appendix/outputs/

appendix_scripts <- c(
  "analysis/appendix/a1_1_overdispersion.R",  # Mean-variance check on the baseline
  "analysis/appendix/a1_4_smooth_acf.R",      # In-sample fit with and without s(t), and residual autocorrelation
  "analysis/appendix/a1_5_single_origin.R"    # In-sample fit and forecast at two origins, compared across models
)

# Each in its own environment, so no script relies on another's objects
for (script in appendix_scripts) {
  message("\n=== ", script, " ===")
  sys.source(script, envir = new.env(parent = globalenv()))
}

message("\nAppendix outputs written to analysis/appendix/outputs")
