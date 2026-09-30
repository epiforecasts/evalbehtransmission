# Reads the CoMix social contact survey from Zenodo, or from data-raw/CoMix/, and writes
# two contact series to data-processed/, both used as covariates in Chapter 1:
#   - the dominant eigenvalue of the age-stratified contact matrix, per survey round
#   - mean contacts per participant per day, as a daily trailing mean
#
# Contacts are aggregated across all settings, and matrices are symmetrised against
# UK population denominators. No susceptibility scaling is applied, so the
# eigenvalue is rho(C) rather than rho(K)
#
# Public CoMix has no region, so the sample is UK-wide and is weighted to the UK population
# The outcome (inc2prev) is England only, so covariates and outcome currently differ in geography
# Private CoMix adds region, allowing England-only or per-nation weighting (see #54)

library(socialmixr)
library(data.table)
library(tidyverse)

source("R/inc2prev_path.R")

## Define paths and parameters -------------------------------------------------

data_dir <- "data-raw/CoMix"
eigenvalues_path <- "data-processed/comix_eigenvalues.csv"
mean_contacts_path <- "data-processed/comix_mean_contacts.csv"

# A pinned version, so a later Zenodo release cannot move the series
comix_record_url <- "https://zenodo.org/records/13684044/files"

# FALSE reads CoMix from data_dir and inc2prev from a local clone, matching ch1_data.R
use_remote <- TRUE

# Read file by file rather than through socialmixr::get_survey(), which joins each extra
# file on every shared column name, so sday and hh_common collide on their row number X
comix_path <- function(file, use_remote) {
  if (use_remote) paste(comix_record_url, file, sep = "/") else file.path(data_dir, file)
}

# Ages are sampled within their reported band, so fix the seed to make a run reproducible
seed <- 42

# Contacts truncated per participant per day, following the CMMID CoMix weekly reports
# A few participants report into the thousands, so the untruncated mean is unusable
contact_cap <- 50

# Two-week pooling of the alternating panels, trailing rather than centred
contact_window_days <- 14

# CoMix participant age bands, as recorded, and their lower limits
# Children join from May 2020 and their share of respondents then swings, so the mean
# is weighted to the population rather than taken over whoever responded
participant_age_bands <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49",
                           "50-59", "60-69", "70-120")
participant_age_lower <- c(0, 5, 12, 18, 30, 40, 50, 60, 70)
child_age_bands       <- c("0-4", "5-11", "12-17")

# UK population by single year of age and sex, ONS mid-2020 estimates
# The file inc2prev uses, read as in its data-raw/update_cis.R but from the sex-specific sheets
load_uk_population <- function(use_remote) {
  pop_file <- inc2prev_path("data-raw/uk_pop.xls", use_remote)
  if (use_remote) { # readxl cannot read from a URL
    local_file <- tempfile(fileext = ".xls")
    download.file(pop_file, local_file, mode = "wb", quiet = TRUE)
    pop_file <- local_file
  }
  rbindlist(lapply(c(M = "Males", F = "Females"), function(sex) {
    sheet <- readxl::read_excel(pop_file, sheet = paste("MYE2 -", sex), skip = 7)
    uk    <- unlist(sheet[sheet$Name == "UNITED KINGDOM", c(as.character(0:89), "90+")])
    data.table(age = 0:90, population = uk) # 90 stands for 90+
  }), idcol = "sex")
}

uk_population <- load_uk_population(use_remote)

# Population share of each participant age band x sex cell, the post-stratification weights
cell_population <- uk_population[
  , .(population = sum(population)),
  by = .(part_band = participant_age_bands[findInterval(age, participant_age_lower)], sex)
][, share := population / sum(population)]

# Age limits: align with ONS CIS antibody data bins and CoMix ranges
age_limits <- c(2, 11, 16, 25, 35, 50, 70)
# Survey population on age_limits from the same UK table, dropping ages below the first limit
# Do not use survey_pop = "United Kingdom" — socialmixr's bundled WPP data only goes to 2015
# and uses 5-year bands that don't align to age_limits, causing interpolation artefacts.
survey_pop <- uk_population[
  , .(population = sum(population)),
  by = .(lower.age.limit = c(NA, age_limits)[findInterval(age, age_limits) + 1])
][!is.na(lower.age.limit)][order(lower.age.limit)]


## Load and merge data ---------------------------------------------------------

load_merge_comix <- function(use_remote) {

  message("Loading CoMix data...")

  # The contact file is about 100 MB, too slow for R's 60 second download timeout
  if (use_remote) options(timeout = max(600, getOption("timeout")))

  participants_raw <- fread(comix_path("CoMix_uk_participant_common.csv", use_remote))
  contacts_raw     <- fread(comix_path("CoMix_uk_contact_common.csv", use_remote))
  sday_raw         <- fread(comix_path("CoMix_uk_sday.csv", use_remote))
  extra_raw        <- fread(comix_path("CoMix_uk_participant_extra.csv", use_remote))

  # sday_id dates each response, being the day the survey was filled in
  # survey_round groups responses for the matrices below, taking its date from sday_id
  # wave is the panel's nth round, carrying no date, so it is never used
  participants_raw <- participants_raw |>
    merge(sday_raw[, .(part_id, sday_id, dayofweek)], by = "part_id", all.x = TRUE) |>
    merge(extra_raw[, .(part_id, survey_round)], by = "part_id", all.x = TRUE)

  participants_raw <- participants_raw[!is.na(survey_round)]

  list(participants = participants_raw, contacts = contacts_raw)
}

## Contact matrices by survey round ---------------------------------------------

compute_round_matrices <- function(data_list, age_limits, survey_pop) {

  message("Computing contact matrices by survey round...")

  participants <- data_list$participants
  contacts <- data_list$contacts

  rounds <- sort(unique(participants$survey_round))

  lapply(rounds, function(this_round) {
    round_participants <- participants[survey_round == this_round]
    round_contacts     <- contacts[part_id %in% round_participants$part_id]

    # Same cap as the mean contacts below, so the two series truncate identically
    round_contacts <- round_contacts[, head(.SD, contact_cap), by = part_id]

    # Date the round by the median day its participants responded, less one day as
    # contacts are reported for the previous day. Rounds are unevenly spaced
    response_dates <- as.Date(round_participants$sday_id, format = "%Y.%m.%d")
    round_date     <- median(response_dates, na.rm = TRUE) - 1

    # dayofweek records the day the survey was filled in, but contacts happened the day before
    # Shifting it means the weekday weighting applies to the day the contacts occurred
    round_participants <- copy(round_participants)[
      , dayofweek := as.POSIXlt(response_dates - 1)$wday]

    # Build local survey object and clean (parses part_age string bands into numeric)
    round_survey <- socialmixr::clean(socialmixr::survey(
      participants = as.data.frame(round_participants),
      contacts     = as.data.frame(round_contacts)
    ))

    # Use tryCatch to handle errors in occasional rounds e.g. sparse data early on
    round_matrix <- tryCatch({
      socialmixr::contact_matrix(
        round_survey,
        age_limits = age_limits,
        symmetric = TRUE,
        weigh_dayofweek = TRUE,
        survey_pop = survey_pop,
        # Almost every age is reported as a band, so socialmixr draws a year uniformly
        # between est_min and est_max, rather than taking the band midpoint
        estimated_participant_age = "sample",
        estimated_contact_age     = "sample",
        # A contact with no age is dropped on its own, and the participant's other contacts stay
        missing_contact_age       = "ignore",
        missing_participant_age   = "remove"
      )$matrix
    }, error = function(e) {
      warning(sprintf("Matrix computation failed for round %s: %s", this_round, e$message))
      NULL
    })

    list(
      survey_round = this_round,
      date         = round_date,
      matrix       = round_matrix
    )
  })
}


## Mean contacts per participant per day ---------------------------------------
# Follows the CMMID CoMix weekly reports: truncate per participant per day, then take the
# mean over participants, pooled across two weeks so both panels contribute
# Their window is centred on a survey round, this one is trailing to avoid future leakage
# The mean is post-stratified by age band x sex to the UK population, as Gimma et al. (2022)
# weight by age and gender

# Population-weighted mean over the cells observed in a window
weighted_cell_mean <- function(mean_contacts, share) {
  observed <- !is.na(mean_contacts)
  if (!any(observed)) return(NA_real_)
  sum(mean_contacts[observed] * share[observed]) / sum(share[observed])
}

compute_mean_contacts <- function(data_list, window_length = contact_window_days,
                                  cap = contact_cap) {

  message("Computing mean contacts per day...")

  # Contacts are reported for the previous day, as in the wave dates above
  responses <- data.table(
    part_id      = data_list$participants$part_id,
    part_age     = data_list$participants$part_age,
    sex          = data_list$participants$part_gender,
    contact_date = as.Date(data_list$participants$sday_id, format = "%Y.%m.%d") - 1
  )[!is.na(contact_date)]

  # Under 1 joins 0-4 and 18-19 joins 18-29, while other bands crossing these are dropped
  responses[, part_band := fcase(part_age == "Under 1", "0-4",
                                 part_age == "18-19",   "18-29",
                                 part_age %in% participant_age_bands, part_age,
                                 default = NA_character_)]

  # Missing or unmappable ages, and other or missing gender, have no population to weight against
  n_missing_age  <- responses[is.na(part_age), .N]
  n_unmapped_age <- responses[!is.na(part_age) & is.na(part_band), .N]
  n_no_sex       <- responses[!is.na(part_band) & !sex %in% c("M", "F"), .N]
  responses <- responses[!is.na(part_band) & sex %in% c("M", "F")]

  # Participants reporting nothing are absent from the contact table and count as zero
  counts <- data_list$contacts[, .(n_contacts = pmin(.N, cap)), by = part_id]
  responses <- merge(responses, counts, by = "part_id", all.x = TRUE)
  responses[is.na(n_contacts), n_contacts := 0]

  days <- seq(min(responses$contact_date), max(responses$contact_date), by = "day")

  # Mean per age band x sex cell over each trailing window, one row per day and cell
  cell_means <- rbindlist(lapply(days, function(day) {
    responses[contact_date > day - window_length & contact_date <= day,
              .(mean_contacts = mean(n_contacts), n = .N), by = .(part_band, sex)][, date := day]
  }))
  cell_means <- cell_means[CJ(date = days, part_band = participant_age_bands, sex = c("M", "F")),
                           on = .(date, part_band, sex)]
  cell_means <- merge(cell_means, cell_population, by = c("part_band", "sex"))
  cell_means[, is_child := part_band %in% child_age_bands]
  setorder(cell_means, part_band, sex, date)

  # Children enter CoMix on 7 May 2020, so earlier days take each child cell's first mean
  # No filled day enters a forecast window, and schools were closed throughout (see #53)
  # Days between child panels in late July carry each child cell's mean forward
  cell_means[, mean_filled := if (is_child[1]) nafill(nafill(mean_contacts, "locf"), "nocb")
                              else mean_contacts,
             by = .(part_band, sex)]

  # One population-weighted mean over all cells, so the series tracks behaviour rather than
  # recruitment. Adult, child (before filling) and unweighted sample means are kept for comparison
  cell_means[, .(
    mean_contacts_standardised = weighted_cell_mean(mean_filled, share),
    mean_contacts_adult        = weighted_cell_mean(mean_contacts[!is_child], share[!is_child]),
    mean_contacts_child        = weighted_cell_mean(mean_contacts[is_child],  share[is_child]),
    mean_contacts_sample       = sum(mean_contacts * n, na.rm = TRUE) / sum(n, na.rm = TRUE),
    n_responses                = sum(n, na.rm = TRUE)
  ), keyby = date]
}


## Convert to next-generation matrix using antibody data -----------------------

matrix_to_nextgen <- function(cm) {
  # TODO: implement NGM scaling
  cm
}


## Extract dominant eigenvalue from matrices -----------------------------------

calculate_dom_eigenvalues <- function(round_data) {

  results <- lapply(round_data, function(round) {

    lambda1 <- NA_real_

    # Skip rounds where the matrix contains NA/Inf — sparse early rounds; lambda1 stays NA
    # matrix_to_nextgen() is a pass-through, so this is rho(C), not rho(K)
    if (!is.null(round$matrix) && all(is.finite(round$matrix))) {
      lambda1 <- max(Re(eigen(round$matrix, only.values = TRUE)$values))
    }

    data.frame(
      survey_round = round$survey_round,
      date         = round$date,
      lambda1      = lambda1
    )
  })

  bind_rows(results)
}


## Run full execution ----------------------------------------------------------

main <- function() {

  set.seed(seed)

  # Load and merge data
  comix_data <- load_merge_comix(use_remote)
  
  # One contact matrix per survey round
  round_matrices <- compute_round_matrices(
    data_list  = comix_data,
    age_limits = age_limits,
    survey_pop = survey_pop
  )

  # Currently a pass-through, so no susceptibility scaling is applied
  round_matrices <- lapply(round_matrices, function(round) {
    round$matrix <- matrix_to_nextgen(round$matrix)
    round
  })

  final_eigenvalues <- calculate_dom_eigenvalues(round_matrices)

  # Already daily, so it needs no anchoring to a round
  mean_contacts <- compute_mean_contacts(comix_data)

  if (!dir.exists(dirname(eigenvalues_path))) {
    dir.create(dirname(eigenvalues_path), recursive = TRUE)
  }

  write.csv(final_eigenvalues, eigenvalues_path, row.names = FALSE)
  write.csv(mean_contacts, mean_contacts_path, row.names = FALSE)
  message(sprintf("Survey round eigenvalues saved to %s", eigenvalues_path))
  message(sprintf("Daily mean contacts saved to %s", mean_contacts_path))

  # Days where the child mean was filled rather than observed
  gaps <- mean_contacts[is.na(mean_contacts_child), .(date)]
  cat("\n--- Days with no children in the trailing window ---\n")
  cat(nrow(gaps), "of", nrow(mean_contacts), "days, first child response",
      format(min(mean_contacts$date[!is.na(mean_contacts$mean_contacts_child)])), "\n")
  print(as.data.frame(gaps[date > min(mean_contacts$date[!is.na(mean_contacts$mean_contacts_child)])]),
        row.names = FALSE)
}

# Run pipeline
main()
