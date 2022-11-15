#' Probability density for a hierarchical prior setup for the elements of the adjacency matrix based on the beta binomial distribution
#'
#' A hierarchical prior setup can be used in \code{\link{W_priors}} to anchor the prior
#' number of expected neighbors. Assuming a \emph{fixed} prior inclusion probability \eqn{\underline{p}=1/2}
#' for the off-diagonal entries in the binary \eqn{n} by \eqn{n} adjacency matrix \eqn{\Omega} implies
#' that the number of neighbors (i.e. the row sums of \eqn{\Omega}) follows a Binomial distribution
#' with a prior expected number of neighbors for the \eqn{n} spatial observations of \eqn{(n-1)\underline{p}}.
#' However, such a prior structure has the potential undesirable effect of promoting a relatively large
#' number of neighbors. To put more prior weight on parsimonious neighborhood structures and promote sparsity
#' in \eqn{\Omega}, the beta binomial prior accounts for the number of neighbors in each row of \eqn{\Omega}.
#'
#' The beta-binomial distribution is the result of treating the prior inclusion probability \eqn{\underline{p}}
#' as random (rather than being fixed) by placing a hierarchical beta prior on it.
#' For the number of neighbors \eqn{x}, the resulting prior on the elements of \eqn{\Omega}, \eqn{\omega_{ij}},
#' can be written as:
#'
#'  \deqn{
#'  p(\omega_{ij} = 1 | x)\propto \Gamma\left(a+ x \right)\Gamma\left(b+(n-1)-x\right),
#'  }
#'
#'  where \eqn{\Gamma(\cdot )} is the Gamma function, and \eqn{a} and
#'  \eqn{b} are hyperparameters from the beta prior. In the case of \eqn{a = b = 1}, the prior takes the
#'  form of a discrete uniform distribution over the number of neighbors. By fixing \eqn{a = 1}
#'  the prior can be anchored around the expected number of neighbors \eqn{m} through
#'  \eqn{b=[(n-1)-m]/m} (see Ley and Steel, 2009).
#'
#'  The prior can be truncated by setting a minimum (\code{min_k}) and/or a maximum number of
#'  neighbors (\code{max_k}). Values outside this range have zero prior support.
#'
#' @param x Number of neighbors (scalar)
#' @param nsize Number of potential neighbors: \code{nsize}\eqn{=(n-1)}
#' @param a Scalar prior parameter \eqn{a}
#' @param b Scalar prior parameter \eqn{b}
#' @param min_k Minimum prior number of neighbors (defaults to 0)
#' @param max_k Maximum prior number of neighbors (defaults to \code{nsize})
#'
#'
#' @export bbinompdf
#'
#' @return Prior density evaluated at \code{x}.
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
