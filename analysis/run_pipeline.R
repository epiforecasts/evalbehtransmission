# Run full analysis pipeline in order.
# Each script reads from data-raw/ or data-processed/ and writes its outputs
# before the next script runs.

source("analysis/process_ons.R")
source("analysis/process_mobility.R")
source("analysis/process_comix.R")
source("analysis/harmonise.R")
source("analysis/estimate_rt.R")
source("analysis/model_rtglm.R")
source("analysis/forecast.R")
source("analysis/evaluate.R")
