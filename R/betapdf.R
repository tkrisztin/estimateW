#' The four-parameter Beta probability density function
#'
#' A four-parameter Beta specification as the prior for the spatial autoregressive parameter \eqn{\rho},
#' as proposed by LeSage and Parent (2007) .
#'
#' The prior density is given by:
#'
#' \deqn{ p(\rho) \sim \frac{1}{Beta(a,b)} \frac{(\rho - \underline{\rho}_{min})^{(a-1)} (\underline{\rho}_{max} - \rho)^{(b-1)} }{2^{a + b - 1}} }
#'
#' where \eqn{Beta(a, b)} (\eqn{a,b > 0}) represents the Beta function,
#' \eqn{Beta(a, b)= \int_{0}^{1} t^{a-1} (1-t)^{b-1} dt}.
#'
#' @param rho The scalar value for \eqn{\rho}
#' @param a The first shape parameter of the Beta distribution
#' @param b The second shape parameter of the Beta distribution
#' @param rmin Scalar \eqn{\underline{\rho}_{min}}: the minimum value of \eqn{\rho}
#' @param rmax Scalar \eqn{\underline{\rho}_{max}}: the maximum value of \eqn{\rho}
#'
#' @return Density value evaluated at \code{rho}.
#'
#' @export
#'
#' @references
#'  LeSage, J. P., and Parent, O. (2007) Bayesian model averaging for spatial econometric models.
#'  \emph{Geographical Analysis}, \bold{39(3)}, 241-267.
betapdf <- function(rho, a = 1, b = 1, rmin = 0, rmax = 1) {
  return(1 / beta(a, b) * ((rho - rmin)^(a - 1) * (rmax - rho)^(b - 1)) / ((rmax-rmin)^(a + b - 1)))
}
