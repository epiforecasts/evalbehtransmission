# Runs the Appendix A.1 scripts, in the same style as analysis/run_pipeline.R
# Reads what the Chapter 1 pipeline has already written, so run that first
# Writes only to stashed/appendix/outputs/, never to outputs/ch1/ or data-processed/

appendix_scripts <- c(
  "stashed/appendix/a1_1_overdispersion.R",  # Mean-variance check on the baseline
  "stashed/appendix/a1_4_smooth_acf.R",      # What s(t) absorbs, and residual autocorrelation
  "stashed/appendix/a1_5_single_origin.R"    # Fit and forecast at one origin
)

# A.1.3 is documentation only, in stashed/appendix/a1_3_eigenvalue_construction.md

for (script in appendix_scripts) {
  message("\n=== ", script, " ===")
  sys.source(script, envir = new.env(parent = globalenv()))
}

message("\nAppendix outputs written to stashed/appendix/outputs")
