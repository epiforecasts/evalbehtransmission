#' Compute total infectivity (lambda_t) from an incidence vector
#'
#' @param incidence Numeric vector of daily incidence values
#' @param gi_weights Numeric vector of generation interval weights (lag 1 onwards, normalised)
#'
#' @returns Numeric vector of lambda_t values, same length as incidence
compute_lambda <- function(incidence, gi_weights) {
  max_lag <- length(gi_weights)
  vapply(seq_along(incidence), function(t) {
    lags <- seq_len(min(t - 1, max_lag))
    if (length(lags) == 0) return(0)
    sum(incidence[t - lags] * gi_weights[lags])
  }, numeric(1))
}

#' Compute total infectivity at the final time point only
#'
#' @param incidence Numeric vector of daily incidence values
#' @param gi_weights Numeric vector of generation interval weights (lag 1 onwards, normalised)
#'
#' @returns Single numeric value: lambda at time length(incidence)
compute_lambda_last <- function(incidence, gi_weights) {
  t <- length(incidence)
  max_lag <- length(gi_weights)
  lags <- seq_len(min(t - 1, max_lag))
  if (length(lags) == 0) return(0)
  sum(incidence[t - lags] * gi_weights[lags])
}
