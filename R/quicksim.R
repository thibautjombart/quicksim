## This script contains a series of functions which determine the state of a
## newly infected case given its infector - e.g. geographic location, genetic
## sequences, etc.


## MAKE NEW INDIVIDUALS, EITHER FROM AN INFECTOR, OR FROM SCRATCH

new_case <- function(infector = NULL, ..., config = new_config(...)) {
    out <- list()
    
    out$location <- new_location(infector, config = config)

    out$dna <- new_dna(infector, config = config)

    return(out)
}
