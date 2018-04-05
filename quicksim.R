

## This script contains a series of functions which determine the state of a
## newly infected case given its infector - e.g. geographic location, genetic
## sequences, etc.




## SPATIAL COORDINATES

## This function uses a Normal spatial kernel for all spatial dimensions, with a
## single parameter 'sd'.

## x: infector's location

new_location <- function(x, sd = 1) {
    out <- rnorm(length(x), mean = x, sd = sd)
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

## mu: mutation rate (per genome and per day)

## genome_length: total length of the genome

new_dna <- function(x, generation_time, mu, genome_length) {
    ## generate new mutations
    n_mutations <- rpois(1, mu * generation_time)
    new_mutations <- runif(n_mutations, min = 1, max = genome_length)
    new_mutations <- as.integer(round(new_mutations))
    out <- c(x, new_mutations)

    ## cleanup: remove reverse mutations
    reverse_mutations <- as.integer(names(which(table(out) %% 2 == 0L)))
    out <- setdiff(out, reverse_mutations)
    return(out)
}


## examples

new_dna(1:10, 10, 0, 10) # no mutation

set.seed(1)
new_dna(1:10, 10, 1, 10) # reverse mutations delete old ones
## [1]  3  4  5  6  8 10




## COMPUTE GENETIC DISTANCES BETWEEN A PAIR OF CASES

## Using the storage of mutations defined above, i.e. one individual is a vector
## of mutations.

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
