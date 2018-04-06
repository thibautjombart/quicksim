#' Create a new case
#'
#' This function creates a new case, with or without a known infector. If it is
#' missing (NULL, default), the case is imported according to parameters in
#' \code{config}.
#'
#' It calls upon \code{\link{new_location}} and \code{\link{new_dna}} to
#' determine locations and DNA sequences.
#'
#' @export
#'
#' @param infector An optional infector; if provided, it must be a list with the
#'     following slots: 'location', 'dna'
#'
#' @param date The date at which the new case is infected; if the case is
#'   imported and this is not provided, defaults to 0.
#'
#' @param ... further arguments passed to \code{\link{new_config}}
#'
#' @param config optionally, a config returned by \code{\link{new_config}}
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#'
#' ## imported cases
#' A <- new_case()
#' B <- new_case()
#' A
#' B
#'
#' ## transmission chain A -> C -> D
#' C <- new_case(A, date = 3)
#' D <- new_case(C, date = 10)
#'
#' A
#' B
#' C
#' D

new_case <- function(infector = NULL, date = NULL,
                     ..., config = new_config(...)) {
  out <- list()

  if (is.null(infector)) {
    if (is.null(date)) {
      date <- 0L
    }
    
    out$date <- date
    out$location <- new_location(NULL, config = config)
    out$dna <- new_dna(NULL, config = config)
    
  } else {

    if (is.null(date)) {
      stop("date of infection must be provided")
    }

    out$date <- date
    delay <- date - infector$date

    if (delay < 1L) {
      stop("Generation time < 1")
    }
    
    out$dna <- new_dna(infector$dna, delay, config = config)
    out$location <- new_location(infector$location, config = config)
    
  }

  class(out) <- c("case", "list")
  attr(out, "config") <- config
  return(out)
}
