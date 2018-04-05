#' Create a new case
#'
#' This function creates a new case, with or without a known infector. If it is
#' missing (NULL, default), the case is imported according to parameters in
#' \code{config}.
#'
#' @export
#'
#' @param infector An optional infector; if provided, it must be a list with the
#'     following slots: 'location', 'dna'
#'
#' @param ... further arguments passed to \code{new_config}
#'

new_case <- function(infector = NULL, ..., config = new_config(...)) {
    out <- list()
    
    out$location <- new_location(infector$location, config = config)

    out$dna <- new_dna(infector$dna, config = config)

    return(out)
}
