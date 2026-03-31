# Evaluate forecast accuracy and explanatory power across model variants.
#
# Metrics:
#   - Incremental partial R² (variance explained beyond autoregressive baseline)
#   - CRPS via {scoringutils} at 1–4 week horizons
#   - Bias and coverage
#   Stratified by pandemic period and age group.