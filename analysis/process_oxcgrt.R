# Reads raw OxCGRT policy stringency data from data-raw/OxCGRT/,
# filters to UK national level, and saves to data-processed/.
#
# OxCGRT variables of interest:
#   - StringencyIndex_Average: composite policy stringency score (0–100)
#   - Individual policy indicators (school closures, workplace closures, etc.)
#
# Downstream use: policy stringency covariate in transmission models
