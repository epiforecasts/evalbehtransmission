# Reads raw CoMix social contact survey data from data-raw/CoMix/,
# computes age-stratified contact matrices, and saves to data-processed/.
#
# CoMix variables of interest:
#   - Age-stratified contact matrices (home, work, other, all settings)
#   - Weighted by participant weights and age-group population denominators
#
# Downstream use: next-generation matrix construction

install.packages("socialmixr")
library(socialmixr)
library(data.table)
library(tidyverse)

## Define paths and parameters -------------------------------------------------

data_dir <- "data-raw/CoMix"
output_dir <- "data-processed/comix_eigenvalues.csv"

# Age limits: align with ONS CIS antibody data bins and CoMix ranges
age_limits <- c(2, 11, 16, 25, 35, 50, 70)
# Survey population: England 2020, from inc2prev populations.csv, aligned to age_limits above.
# Do NOT use survey_pop = "United Kingdom" — socialmixr's bundled WPP data only goes to 2015
# and uses 5-year bands that don't align to age_limits, causing interpolation artefacts.
survey_pop <- read.csv("data-raw/inc2prev-main/data-processed/populations.csv") |>
  dplyr::filter(level == "age_school", geography == "England") |>
  dplyr::select(lower.age.limit = lower_age_limit, population)


## Load and merge data ---------------------------------------------------------

load_merge_comix <- function(data_dir) {
  
  message("Loading raw CoMix data...")
  
  participants_raw <- fread(file.path(data_dir, "CoMix_uk_participant_common.csv"))
  contacts_raw     <- fread(file.path(data_dir, "CoMix_uk_contact_common.csv"))
  sday_raw         <- fread(file.path(data_dir, "CoMix_uk_sday.csv"))
  extra_raw        <- fread(file.path(data_dir, "CoMix_uk_participant_extra.csv"))

  # Merge sday onto participants: provides sday_id and dayofweek
  # NOTE: wave in sday is panel-specific (participant's nth survey), NOT a
  # global time index. Do not group by wave.
  participants_raw <- merge(
    participants_raw,
    sday_raw[, .(part_id, sday_id, dayofweek)],
    by = "part_id",
    all.x = TRUE
  )

  # Merge participant_extra to get survey_round — the global weekly time unit
  # (survey_round 1-101 each span ~1 calendar week across all panels)
  participants_raw <- merge(
    participants_raw,
    extra_raw[, .(part_id, survey_round)],
    by = "part_id",
    all.x = TRUE
  )

  participants_raw <- participants_raw[!is.na(survey_round)]

  list(participants = participants_raw, contacts = contacts_raw)
}

## Compute weekly contact matrices based on wave -------------------------------

compute_wave_matrices <- function(data_list, age_limits, survey_pop) {
  
  message("Computing contact matrices by wave...")
  
  participants <- data_list$participants
  contacts <- data_list$contacts
  
  waves <- sort(unique(participants$survey_round))

  wave_data <- lapply(waves, function(w) {
    # Subset to current survey_round (global weekly time unit)
    p_wave <- participants[survey_round == w]
    c_wave <- contacts[part_id %in% p_wave$part_id]
    
    # Calculate median contact date for this wave
    dates <- as.Date(p_wave$sday_id, format = "%Y.%m.%d")
    wave_date <- median(dates, na.rm = TRUE) - 1
    
    # Build local survey object and clean (parses part_age string bands into numeric)
    wave_survey <- socialmixr::survey(
      participants = as.data.frame(p_wave),
      contacts = as.data.frame(c_wave)
    )
    wave_survey <- socialmixr::clean(wave_survey)

    # Compute the matrix using socialmixr
    # Use tryCatch to handle errors in occasional surveys e.g. sparse data early on
    matrix_all <- tryCatch({
      socialmixr::contact_matrix(
        wave_survey,
        age_limits = age_limits,
        symmetric = TRUE,
        weigh_dayofweek = TRUE,
        survey_pop = survey_pop
      )$matrix
    }, error = function(e) {
      warning(sprintf("Matrix computation failed for wave %s: %s", w, e$message))
      NULL
    })
    
    list(
      wave = w,
      date = wave_date,
      matrix = matrix_all
    )
  })
  
  return(wave_data)
}


## Convert to next-generation matrix using antibody data -----------------------

matrix_to_nextgen <- function(cm) {
  # TODO: implement NGM scaling
  cm
}


## Extract dominant eigenvalue from matrices -----------------------------------

calculate_dom_eigenvalues <- function(wave_data) {
  
  results <- lapply(wave_data, function(wd) {
    
    lambda1 <- NA_real_
    
    # Skip waves where matrix contains NA/Inf — sparse early waves; lambda1 stays NA
    if (!is.null(wd$matrix) && all(is.finite(wd$matrix))) {
      
      NGM <- wd$matrix
      
      lambda1 <- max(Re(eigen(NGM, only.values = TRUE)$values))
    }
    
    data.frame(
      wave = wd$wave,
      date = wd$date, 
      lambda1 = lambda1
    )
  })
  
  bind_rows(results)
}


## Run full execution ----------------------------------------------------------

main <- function() {
  
  # Load and merge data
  comix_data <- load_merge_comix(data_dir)
  
  # Compute matrices
  wave_matrices <- compute_wave_matrices(
    data_list = comix_data,
    age_limits = age_limits,
    survey_pop = survey_pop
  )

  # Apply NGM scaling
  wave_matrices <- lapply(wave_matrices, function(wd) {
    wd$matrix <- matrix_to_nextgen(wd$matrix)
    wd
  })

  # Extract eigenvalues
  final_eigenvalues <- calculate_dom_eigenvalues(wave_matrices)
  
  # Save outputs
  if (!dir.exists(dirname(output_dir))) {
    dir.create(dirname(output_dir), recursive = TRUE)
  }
  
  write.csv(final_eigenvalues, output_dir, row.names = FALSE)
  message(sprintf("Wave eigenvalues saved to %s", output_dir))
}

# Run pipeline
main()














## Single-wave diagnostic — run line by line -------

# Step 1: Load raw files exactly as load_merge_comix() does
d_participants <- fread(file.path(data_dir, "CoMix_uk_participant_common.csv"))
d_contacts     <- fread(file.path(data_dir, "CoMix_uk_contact_common.csv"))
d_sday         <- fread(file.path(data_dir, "CoMix_uk_sday.csv"))
d_extra        <- fread(file.path(data_dir, "CoMix_uk_participant_extra.csv"))

# Step 2: Check for duplicates in sday — if > 1 row per part_id the merge will expand participants
cat("Rows in sday:", nrow(d_sday), "\n")
cat("Unique part_ids in sday:", length(unique(d_sday$part_id)), "\n")

# Step 3: Check for duplicates in extra
cat("Rows in extra:", nrow(d_extra), "\n")
cat("Unique part_ids in extra:", length(unique(d_extra$part_id)), "\n")

# Step 4: Merge sday and check row count — should equal nrow(d_participants)
p_merged <- merge(d_participants, d_sday[, .(part_id, sday_id, dayofweek)],
                  by = "part_id", all.x = TRUE)
cat("Participants before sday merge:", nrow(d_participants), "\n")
cat("Participants after sday merge: ", nrow(p_merged), "\n")

# Step 5: Merge extra and check again
p_merged <- merge(p_merged, d_extra[, .(part_id, survey_round)],
                  by = "part_id", all.x = TRUE)
p_merged <- p_merged[!is.na(survey_round)]
cat("Participants after extra merge and NA drop:", nrow(p_merged), "\n")
cat("Unique part_ids after merges:", length(unique(p_merged$part_id)), "\n")

# Step 6: Subset to survey_round 1 and inspect
p1 <- p_merged[survey_round == 10]
c1 <- d_contacts[part_id %in% p1$part_id]
cat("survey_round 1 — participants:", nrow(p1), " contacts:", nrow(c1), "\n")
cat("dayofweek present:", "dayofweek" %in% names(p1), "\n")
cat("sday_id sample:", head(p1$sday_id), "\n")

# Step 7: Build socialmixr survey object and clean (parses part_age string bands into numeric)
s1 <- socialmixr::survey(participants = as.data.frame(p1), contacts = as.data.frame(c1))
s1 <- socialmixr::clean(s1)

# Step 8: Compute contact matrix — no tryCatch so real errors surface
cm1 <- contact_matrix(
  s1,
  age_limits      = age_limits,
  symmetric       = TRUE,
  weigh_dayofweek = TRUE,
  survey_pop      = survey_pop
)
cat("Matrix result:\n"); print(cm1$matrix)

# Step 9: Check matrix is finite before computing eigenvalue
NGM1 <- cm1$matrix
cat("Any NA:", anyNA(NGM1), "\n")
cat("Any Inf/NaN:", any(!is.finite(NGM1[!is.na(NGM1)])), "\n")

# Step 10: Compute dominant eigenvalue with guard
if (is.null(NGM1) || !all(is.finite(NGM1))) {
  cat("Cannot compute eigenvalue — matrix contains non-finite values\n")
} else {
  lambda1 <- max(Re(eigen(NGM1, only.values = TRUE)$values))
  cat("lambda1:", lambda1, "\n")
}

## End single-wave diagnostic ---------------------------------------------------