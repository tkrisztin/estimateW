#' Return the density of the four-paramater Beta distribution for a given rho value
#'
#' @param rho The scalar value for \eqn{\rho}
#' @param a The first shape parameter of the Beta distribution
#' @param b The second shape parameter of the Beta distribution
#' @param rmin Minimum value of \eqn{\rho}
#' @param rmax Maximum value of \eqn{\rho}
#'
#' @return
#' Returns the beta probability.
beta_prob <- function(rho, a = 1, b = 1, rmin = 0, rmax = 1) {
  return(1 / beta(a, b) * ((rho - rmin)^(a - 1) * (rmax - rho)^(b - 1)) / ((rmax-rmin)^(a + b - 1)))
}
