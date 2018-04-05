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
