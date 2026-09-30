# Chapter 1, issue #51: compare forecast scores across generation intervals
# Reads the records written by ch1_scoring.R, one per generation interval

library(dplyr)
library(readr)

gi_records <- list.files("data-processed", pattern = "^ch1_gi_.*\\.csv$", full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |>
  bind_rows()

# One table per generation interval, printed for comparison
for (gi_name in unique(gi_records$gi)) {
  cat("\n---", gi_name, "---\n")
  print(as.data.frame(filter(gi_records, gi == gi_name)), row.names = FALSE)
}
