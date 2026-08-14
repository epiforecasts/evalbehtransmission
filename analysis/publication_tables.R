# Renders the score tables from ch1_scoring.R as Word tables for the report.
# Presentation only: values are read from outputs/ch1/ and nothing is recomputed, so review
# need go no further than checking the rendered values against the source CSVs.

library(dplyr)
library(readr)
library(flextable)

## Config ----------------------------------------------------------------------

ch1_table_config <- list(
  overall_path   = "outputs/ch1/table_crps_overall.csv",
  horizon_path   = "outputs/ch1/table_crps_horizon.csv",
  period_path    = "outputs/ch1/table_crps_period.csv",
  horizon_output = "outputs/ch1/table_crps_main.docx",
  period_output  = "outputs/ch1/table_crps_period.docx",

  model_levels   = c("baseline", "contacts", "mobility", "combined"),
  horizon_levels = c("Overall", "1 week ahead", "2 weeks ahead",
                     "3 weeks ahead", "4 weeks ahead"),

  # Components tie at three decimals and coverage sits far below nominal, so neither is ranked
  bold_columns = c("crps", "rel_crps"),

  # Signed errors and proportions do not carry three decimals
  two_decimals = c("bias", "interval_coverage_50", "interval_coverage_90")
)

header_labels <- list(
  model = "Model", crps = "log-CRPS", rel_crps = "rCRPS",
  overprediction = "Over-prediction", underprediction = "Under-prediction",
  dispersion = "Dispersion", bias = "Bias",
  interval_coverage_50 = "50% coverage", interval_coverage_90 = "90% coverage"
)

## Build one table -------------------------------------------------------------
# scores must hold a grouping column, model, and the value columns to display

build_score_table <- function(scores, group_column, value_columns, heading, caption) {

  bold_columns <- intersect(ch1_table_config$bold_columns, value_columns)

  scores <- scores |>
    arrange(.data[[group_column]], model) |>
    select(all_of(c(group_column, "model", value_columns)))

  # Flagged before rounding, so ties are decided on full precision
  is_best <- scores |>
    group_by(.data[[group_column]]) |>
    mutate(across(all_of(bold_columns), \(x) x == min(x))) |>
    ungroup()

  # rel_crps is 1 by definition for the baseline, so it is neither shown nor bolded
  baseline <- scores$model == "baseline"
  is_best$rel_crps[baseline] <- FALSE
  scores$rel_crps[baseline]  <- NA_real_

  # as_grouped_data inserts a label row above each block, so model is NA on those rows
  grouped    <- as_grouped_data(scores, groups = group_column)
  label_rows <- which(is.na(grouped$model))
  data_rows  <- which(!is.na(grouped$model))

  stopifnot(length(data_rows) == nrow(scores))

  table <- as_flextable(grouped, hide_grouplabel = TRUE) |>
    set_header_labels(values = header_labels) |>
    colformat_double(j = setdiff(value_columns, ch1_table_config$two_decimals),
                     digits = 3, na_str = "—") |>
    colformat_double(j = intersect(value_columns, ch1_table_config$two_decimals),
                     digits = 2) |>
    align(j = value_columns, align = "right", part = "all") |>
    add_header_lines(heading, top = TRUE) |>
    border_remove() |>
    hline_top(part = "header", border = fp_border_default(width = 1.5)) |>
    hline(i = 1, part = "header", border = fp_border_default(width = 1)) |>
    hline_bottom(part = "header", border = fp_border_default(width = 1)) |>
    bold(part = "header") |>
    align(i = 1, align = "left", part = "header") |>
    bold(i = label_rows) |>
    add_footer_lines(caption) |>
    hline_top(part = "footer", border = fp_border_default(width = 1.5)) |>
    font(fontname = "Calibri", part = "all") |>
    fontsize(size = 8, part = "all") |>
    fontsize(size = 7, part = "footer") |>
    padding(padding.top = 1, padding.bottom = 1,
            padding.left = 3, padding.right = 3, part = "all") |>
    # Fills the text width of whatever page it is pasted into, wrapping headers as needed
    set_table_properties(layout = "autofit", width = 1)

  for (column in bold_columns) {
    table <- bold(table, i = data_rows[is_best[[column]]], j = column)
  }

  table
}

## By horizon ------------------------------------------------------------------

horizon_scores <- bind_rows(
  read_csv(ch1_table_config$overall_path, show_col_types = FALSE) |>
    mutate(horizon = "Overall"),
  read_csv(ch1_table_config$horizon_path, show_col_types = FALSE) |>
    mutate(horizon = paste0(horizon, if_else(horizon == 1, " week", " weeks"), " ahead"))
) |>
  mutate(model   = factor(model, levels = ch1_table_config$model_levels),
         horizon = factor(horizon, levels = ch1_table_config$horizon_levels))

table_by_horizon <- build_score_table(
  horizon_scores, "horizon",
  value_columns = c("crps", "rel_crps", "overprediction", "underprediction",
                    "dispersion", "interval_coverage_50", "interval_coverage_90"),
  heading = paste("Evaluation scores for incidence-only, contacts, mobility and",
                  "combined models, England"),
  caption = paste("Forecast performance by model and horizon, on the log scale, averaged",
                  "over 27 weekly forecast origins from 01 Jul 2020 to 30 Dec 2020.",
                  "rCRPS is the ratio to the incidence-only baseline, so values below one",
                  "indicate improvement. Bold indicates the best performing model at each",
                  "horizon.")
)

## By pandemic period ----------------------------------------------------------
# ch1_scoring.R writes periods in chronological order, and omits any holding no origins

period_scores <- read_csv(ch1_table_config$period_path, show_col_types = FALSE) |>
  mutate(model  = factor(model, levels = ch1_table_config$model_levels),
         # The origin count belongs with the period rather than repeated on every row
         period = sprintf("%s (%d origins)", period, n_origins),
         period = factor(period, levels = unique(period)))

table_by_period <- build_score_table(
  period_scores, "period",
  value_columns = c("crps", "rel_crps", "overprediction", "underprediction", "bias"),
  heading = "Evaluation scores by pandemic period, England",
  caption = paste("Forecast performance by model and pandemic period, on the log scale,",
                  "averaged over horizons 1 to 4 and over the origins falling in each",
                  "period. Periods containing no forecast origins are omitted. rCRPS is the",
                  "ratio to the incidence-only baseline, so values below one indicate",
                  "improvement. Bias is the mean signed error, where negative indicates",
                  "under-prediction. Bold indicates the best performing model in each period.")
)

## Save ------------------------------------------------------------------------

save_as_docx(table_by_horizon, path = ch1_table_config$horizon_output)
save_as_docx(table_by_period,  path = ch1_table_config$period_output)
message("Saved ", ch1_table_config$horizon_output, " and ", ch1_table_config$period_output)
