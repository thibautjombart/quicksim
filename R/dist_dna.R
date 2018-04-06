## COMPUTE GENETIC DISTANCES BETWEEN A PAIR OF CASES

## Using the storage of mutations defined above, i.e. one individual is a vector
## of mutations. See 'dist_dna' below for the version computing all pairwise
## comparisons from a list of 'DNA sequences' (i.e. list of vectors of
## mutations).

## x: a 'DNA sequence', stored as a vector of integers indicating mutations
## y: a 'DNA sequence', stored as a vector of integers indicating mutations

dist_dna_pair <- function(x, y) {
    common_mutations <- intersect(x, y)
    length(setdiff(c(x, y), common_mutations))
}


## examples
dist_dna_pair(NULL, NULL)
dist_dna_pair(NULL, 1)
dist_dna_pair(1, 1)
dist_dna_pair(1, 1:10)
dist_dna_pair(1:3, 1:10)





#' Compute pairwise genetic distances
#'
#' This function Compute pairwise genetic distances from vectors of mutations as
#' generated in the simulations. The computed distance is the Hamming distance,
#' i.e. the number of nucleotide differences between sequences.
#'
#' @export
#'
#' @param x A list of integer vectors, each of which represents the genome of
#'   a pathogen as positions of mutations.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#' dist_dna(list(integer(0), integer(0), 1, 1))
#' dist_dna(list(integer(0), integer(0), 1, 1, 1:2, 1:3, 1:4, 11:20))
#'

dist_dna <- function(x) {
    N <- length(x)
    all_pairs <- combn(seq_len(N), 2L, simplify = FALSE)
    n_pairs <- ncol(all_pairs)
    out <- dist(seq_len(N))
    out[] <- NA_integer_

    for (i in seq_along(all_pairs)) {
	out[i] <- dist_dna_pair(x[[all_pairs[[i]][1]]],
				x[[all_pairs[[i]][2]]])
    }

    return(out)
}
