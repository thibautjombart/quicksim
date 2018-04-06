##' @export
##' @rdname new_case
##' @param x A 'case' object, as returned by \code{\link{new_case}}.

print.case <- function(x, ...) {
  cat("<case object>\n")

  cat("\n/// date of infection ($date): \n")
  print(x$date)

  cat("\n/// place of infection ($location): \n")
  print(x$location)

  cat("\n/// DNA sequence mutations ($dna): \n")
  print(x$dna)
  
  invisible(x)
}
