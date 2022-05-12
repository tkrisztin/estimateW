#' Return the density of the four-paramater Beta distribution for a given rho value
#'
#' @param rho The scalar value for \eqn{\rho}
#' @param a The distribution parameter; see LeSage and Pace (2019)
#'
#' @return
#' Returns the beta probability.
beta_prob <- function(rho, a) {
  return(1 / beta(a, a) * ((1 + rho)^(a - 1) * (1 - rho)^(a - 1)) / (2^(2 * a - 1)))
}
