## Auxiliary, non-exported functions.


## This function implements a Normal spatial kernel with a standard deviation
## defaulting to 1.

kernel_normal <- function(sd) {
  function(x) stats::rnorm(length(x), mean = x, sd = sd)
}
