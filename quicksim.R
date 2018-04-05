## This script contains a series of functions which determine the state of a
## newly infected case given its infector - e.g. geographic location, genetic
## sequences, etc.


## INTERNAL FUNCTIONS

modify_defaults <- function(defaults, x, strict = TRUE) {
    extra <- setdiff(names(x), names(defaults))
    if (strict && (length(extra) > 0L)) {
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    utils::modifyList(defaults, x, keep.null = TRUE) # keep.null is needed here
}


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






## SPATIAL COORDINATES

## This function uses a Normal spatial kernel for all spatial dimensions, with a
## single parameter 'sd'.

## x: infector's location

new_location <- function(x = NULL, ..., config = new_config(...)) {
    if (is.null(x)) {
        out <- c(runif(1, min = config$x_min, max = config$x_max),
                 runif(1, min = config$y_min, max = config$y_max))
    } else {
        out <- rnorm(length(x), mean = x, sd = config$sd_spatial)
    }
    return(out)
}


## example

x_loc <- c(3, 5)
new_locs_sd1 <- t(replicate(100, new_location(x_loc)))
new_locs_sd5 <- t(replicate(100, new_location(x_loc, 5)))

plot(matrix(x_loc, ncol = 2), xlim = c(-20, 20), ylim = c(-20, 20), pch = "x")
points(new_locs_sd1, col = "blue")
points(new_locs_sd5, col = "red")





## GENETIC DATA

## We consider that all markers are biallelic, and store only mutations.

## x: ancestor's 'DNA', stored as a vector of mutations

## generation_time: the number of days between the infection of x and the new
## case

## mu: mutation rate (per site and per day)

## genome_length: total length of the genome

new_dna <- function(x = NULL, generation_time = NULL, ..., config = new_config(...)) {

    ## generate new lineage
    
    if (is.null(x)) {
        out <- new_dna(integer(0),
                     generation_time = config$separation_lineages,
                     config = config)
    } else {

        if (is.null(generation_time)) {
            stop("generation_time must be specified (currently NULL)")
        }
        
        ## generate new mutations
        lambda <- generation_time * config$mutation_rate * config$genome_length
        n_mutations <- rpois(1, lambda)
        new_mutations <- runif(n_mutations, min = 1, max = config$genome_length)
        new_mutations <- as.integer(round(new_mutations))
        out <- c(x, new_mutations)

        ## cleanup: remove reverse mutations
        reverse_mutations <- as.integer(names(which(table(out) %% 2 == 0L)))
        out <- setdiff(out, reverse_mutations)
    }
    
    return(out)
}


## examples

new_dna(1:10, 10, 0, 10) # no mutation

set.seed(1)
new_dna(1:10, 10, 1, 10) # reverse mutations delete old ones
## [1]  3  4  5  6  8 10




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






## x: a list of integer vectors, each of which represents the genome of a
## pathogen as positions of mutations

dist_dna <- function(x) {
    N <- length(x)
    all_pairs <- combn(seq_len(N), 2L, simplify = FALSE)
    n_pairs <- ncol(all_pairs)
    out <- dist(seq_len(N))
    out[] <- NA_integer_

    ## browser()
    for (i in seq_along(all_pairs)) {
	out[i] <- dist_dna_pair(x[[all_pairs[[i]][1]]],
				x[[all_pairs[[i]][2]]])
    }

    return(out)
}


## examples
dist_dna(list(integer(0), integer(0), 1, 1))
dist_dna(list(integer(0), integer(0), 1, 1, 1:2, 1:3, 1:4, 11:20))






## MAKE NEW INDIVIDUALS, EITHER FROM AN INFECTOR, OR FROM SCRATCH

new_case <- function(infector = NULL, ..., config = new_config(...)) {
    out <- list()
    
    out$location <- new_location(infector, config = config)

    out$dna <- new_dna(infector, config = config)

    return(out)
}
