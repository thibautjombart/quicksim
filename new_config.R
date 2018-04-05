

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


