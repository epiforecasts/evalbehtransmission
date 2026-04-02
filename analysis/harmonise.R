# Reads processed data streams from data-processed/ and aligns them
# to a common temporal grid. Standardises (z-scores) behavioural covariates
# prior to model fitting so that coefficients are directly comparable.
#
# TODO: decide on data architecture before implementing:
#   Option A — one consolidated master dataset with all covariates joined
#               (simpler pipeline, but forces all models to share the same
#               observation window where all streams overlap)
#   Option B — separate datasets per model, each containing only the
#               covariates relevant to that model variant
#               (more flexible, handles differing temporal coverage cleanly)