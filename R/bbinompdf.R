#' The beta binomial pdf for sparsity priors
#'
#' @param x Number of neighbors (scalar)
#' @param nsize Maximal number of elements
#' @param a Scalar prior parameter a
#' @param b Scalar prior parameter b
#' @param min_k Minimum number of elements (defaults to 0)
#' @param max_k Maximum number of elements (defaults to nsize + 1)
#'
#' @export bbinompdf
#'
#' @return Density of neighbors
bbinompdf <- function(x, nsize, a, b, min_k = 0, max_k = nsize) {
  x2 <- base::beta(a + x, b + nsize - x) / base::beta(a, b)
  if (any(x < min_k || x > max_k)) {
    x2[x < min_k | x > max_k] <- 0
  }
  return(x2)
}
