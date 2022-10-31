#' The beta binomial probability density function for sparsity priors
#'
#' The beta binomial distribution can be used in \code{\link{W_priors}} as a prior on the
#' number of expected neighbors. Assuming a \emph{fixed} prior for \eqn{\Omega}  implies
#' that the number of neighbors follows a Binomial distribution with a prior
#' expected number of neighbors of \eqn{(N-1)\underline{p}}. However, this prior structure
#' can promote a relatively large number of neighbors. To put more prior weight on
#' parsimonious neighborhood structures and  promote sparsity in \eqn{\Omega},
#' the beta binomial prior accounts for the number of linkages in each row of \eqn{\Omega}.
#'
#' ### TODO: REWRITE
#' This flexible prior structure on the number of neighbors is a truncated beta-binomial distribution
#' with two prior hyperparameters \code{a} and \code{b}. The resulting prior on \eqn{\omega_{ij}} can be written as:
#'
#'  \deqn{
#'  p(\omega_{ij} = 1 | x)\propto \Gamma\left(a+ x \right)\Gamma\left(b+(N-1)-x\right),
#'  }
#'
#'  where \eqn{\Gamma(\cdot )} is the Gamma function, and \eqn{a} and
#'  \eqn{b} are prior hyperparameters. In the case of \eqn{a = b = 1}, the prior takes the
#'  form of a discrete uniform distribution over the number of neighbors. By fixing \eqn{a = 1}
#'  the prior can be anchored around the expected number of neighbors \eqn{m} through
#'  \eqn{b=[(N-1)-m]/m} (see Ley and Steel, 2009).
#'
#'  The prior is truncated when \code{x < min_k} or when \code{x > max_k}.
#'
#' @param x Number of neighbors (scalar)
#' @param nsize Maximal number of elements
#' @param a Scalar prior parameter a
#' @param b Scalar prior parameter b
#' @param min_k Minimum number of elements (defaults to 0)
#' @param max_k Maximum number of elements (defaults to nsize)
#'
#'
#' @export bbinompdf
#'
#' @return Resulting prior density of evaluated at \code{x}.
#'
#' @references
#'   Ley, E., & Steel, M. F. (2009). On the effect of prior assumptions in Bayesian
#'   model averaging with applications to growth regression. \emph{Journal of Applied Econometrics},
#'   \bold{24(4)}. \doi{10.1002/jae.1057}.
bbinompdf <- function(x, nsize, a, b, min_k = 0, max_k = nsize) {
  x2 <- base::beta(a + x, b + nsize - x) / base::beta(a, b)
  if (any(x < min_k || x > max_k)) {
    x2[x < min_k | x > max_k] <- 0
  }
  return(x2)
}
