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
