#' Create the location of a new case
#'
#' This function creates the DNA sequence for a new case, with or without a known
#' infector. If the infector is missing (NULL, default), the case is imported according to
#' parameters in \code{config}.
#'
#' @export
#'
#' @param dna An optional location of the infector; if provided, it must be
#'     a list with the following slots: 'location', 'dna'
#'
#' @param generation_time The number of days separating the dates of infection
#'   of the infector, and the infectee. Only needed if the infector is provided.
#'
#' @param ... further arguments passed to \code{\link{new_config}}
#'
#' @param config optionally, a config returned by \code{\link{new_config}}
#'
#' @details We consider that all markers are biallelic, and store only
#'   mutations.
#'
#' @seealso \code{\link{dist_dna}} to compute pairwise genetic distances from
#'   vectors of mutations.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#'
#' ## without initial location
#' x_dna <- new_dna()
#' x_dna
#'
#' ## drift from this sequence, different nb of generations
#' new_dna(x_dna, 1L) # 1 generation
#' new_dna(x_dna, 12L) # 12 generations
#'
#' ## 100 generations, other genetic params
#' new_dna(x_dna, 100, mutation_rate = 1e-6, genome_length = 2e6)
#'
#' ## illustrate reverse mutations
#' my_config <- new_config(mutation_rate = 0.1, genome_length = 10)
#' x <- integer(0)
#' set_seed(1)
#' x <- new_dna(x, 5); x
#' x <- new_dna(x, 5); x
#' x <- new_dna(x, 5); x
#' x <- new_dna(x, 5); x
#' x <- new_dna(x, 5); x


new_dna <- function(dna = NULL, generation_time = NULL,
                    ..., config = new_config(...)) {

    ## generate new lineage
    
    if (is.null(dna)) {
        out <- new_dna(integer(0),
                     generation_time = config$separation_lineages,
                     config = config)
    } else {

        if (is.null(generation_time)) {
            stop("generation_time must be specified (currently NULL)")
        }
        
        ## generate new mutations
        lambda <- generation_time * config$mutation_rate * config$genome_length
        n_mutations <- stats::rpois(1, lambda)
      new_mutations <- stats::runif(n_mutations,
                                    min = 1,
                                    max = config$genome_length)
        new_mutations <- as.integer(round(new_mutations))
        out <- c(dna, new_mutations)

        ## cleanup: remove reverse mutations
        reverse_mutations <- as.integer(names(which(table(out) %% 2 == 0L)))
        out <- setdiff(out, reverse_mutations)
    }
    
    return(out)
}

