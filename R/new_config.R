#' Create a simulation configuration
#'
#' This function creates a list storing simulation settings - see details.
#'
#' @export
#'
#' @param ... named simulation settings - see details
#'
#' @details
#'
#' Simulation settings include:
#' \itemize{
#'   \item \code{x_min}: lower bound of spatial coordinates (x-axis)
#'   \item \code{x_max}: upper bound of spatial coordinates (x-axis)
#'   \item \code{y_min}: lower bound of spatial coordinates (y-axis)
#'   \item \code{y_max}: upper bound of spatial coordinates (y-axis)
#'   \item \code{sd_spatial}: standard deviation of the Normal kernel use in
#'     spatial dispersion
#'   \item \code{genome_length}: number of nucleotides in the pathogen genome
#'   \item \code{mutation_rate}: per nucleotide and day
#'   \item \code{separation_lineages}: number of days to ancestral lineage for
#'     newly introduced pathogens 
#' }
#'
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#'
#' new_config()
#' new_config(genome_length = 100, xmax = 10)
#' 
new_config <- function(...) {
    defaults <- list(x_min = 0, # min x, spatial coords
                     x_max = 100, # max x, spatial coords
                     y_min = 0, # min y, spatial coords
                     y_max = 100, # max y, spatial coords
                     sd_spatial = 1, # sd of normal spatial kernel
                     genome_length = 3e4, # genome length
                     mutation_rate = 1e-5, # mutation rate
                     separation_lineages = 365 # nb days to ancestral lineage
                     )
    
    args <- list(...)
    
    modify_defaults(defaults, args)    
}


