##' @export
##' @rdname new_case
##' @param x A 'case' object, as returned by \code{\link{new_case}}.

print.case <- function(x, ...) {
  cat("<case object>\n")

  cat("\ndate of infection ($date): ")
  print(x$date)

  cat("\nplace of infection ($location): ")
  print(x$location)

  cat("\nDNA sequence mutations ($dna): ")
  print(x$dna)
  
  invisible(x)
}
