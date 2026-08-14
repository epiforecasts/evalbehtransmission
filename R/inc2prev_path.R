#' Path to a file in the inc2prev repository
#'
#' @param path Path to the file within the inc2prev repository
#' @param use_remote If TRUE read from GitHub, otherwise from a local clone in
#'   data-raw/inc2prev-main/
#'
#' @returns A URL or local file path
inc2prev_path <- function(path, use_remote = TRUE) {
  if (use_remote) {
    paste0("https://raw.githubusercontent.com/epiforecasts/inc2prev/refs/heads/main/", path)
  } else {
    file.path("data-raw/inc2prev-main", path)
  }
}
