#' Create the location of a new case
#'
#' This function creates the location of a new case, with or without a known
#' infector. If it is missing (NULL, default), the case is imported according to
#' parameters in \code{config}.
#'
#' @export
#'
#' @param location An optional location of the infector; if provided, it must be
#'     a list with the following slots: 'location', 'dna'
#'
#' @param ... further arguments passed to \code{\link{new_config}}
#'
#' @param config optionally, a config returned by \code{\link{new_config}}
#'
#' @details This function uses a Normal spatial kernel for all spatial
#'     dimensions, with a single parameter 'sd_spatial'.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#'
#' ## without initial location
#' x_loc <- new_location()
#'
#' ## dispersal from this location, different sd
#' new_locs_sd1 <- t(replicate(100, new_location(x_loc)))
#' new_locs_sd5 <- t(replicate(100, new_location(x_loc, 5)))
#'
#' plot(matrix(x_loc, ncol = 2), xlim = c(-20, 20), ylim = c(-20, 20), pch = "x")
#' points(new_locs_sd1, col = "blue")
#' points(new_locs_sd5, col = "red")
#'
#' ## using a custom spatial kernel
#' custom_kernel <- function(x) {
#'   out <- x + runif(length(x), min = -2, max = 20)
#'   out
#' }
#'
#' new_location(c(0,0)) # default kernel
#' new_location(c(0,0), spatial_kernel = custom_kernel) # custom kernel
#'


new_location <- function(location = NULL, ..., config = new_config(...)) {
  if (is.null(location)) {
    out <- c(stats::runif(1, min = config$x_min, max = config$x_max),
             stats::runif(1, min = config$y_min, max = config$y_max))
  } else {
    out <- config$spatial_kernel(location)
  }
  return(out)
}
