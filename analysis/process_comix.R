# Reads raw CoMix social contact survey data from data-raw/CoMix/,
# computes age-stratified contact matrices, and saves to data-processed/.
#
# CoMix variables of interest:
#   - Age-stratified contact matrices (home, work, other, all settings)
#   - Weighted by participant weights and age-group population denominators
#
# Downstream use: next-generation matrix construction