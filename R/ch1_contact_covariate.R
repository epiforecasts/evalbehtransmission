# Which CoMix series is used as the contact covariate, and how it is named in figures
# Sourced by ch1_covariates.R and ch1_descriptive.R so both stay in step

contact_covariate <- "mean_contacts"   # or "eigenvalue"

contact_covariate_labels <- c(
  eigenvalue    = "CoMix contact matrix dominant eigenvalue",
  mean_contacts = "CoMix mean contacts per participant per day"
)

contact_covariate_label <- function(which = contact_covariate) {
  unname(contact_covariate_labels[which])
}
