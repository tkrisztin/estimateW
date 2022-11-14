#' Pace and Barry's log determinant approximation
#'
#' Bayesian estimates of parameters of SAR and SDM type spatial models require the computation
#' of the log-determinant of positive-definite spatial projection matrices of the form
#' \eqn{(I_n - \rho W)}, where \eqn{W} is a \eqn{n} by \eqn{n} spatial weight matrix. However, direct computation
#' of the log-determinant is computationally expensive.
#'
#' This function wraps the log-determinant approximation by Barry and Pace (1999), which
#' can be used to precompute the log-determinants over a grid of \eqn{\rho} values.
#'
#' @param W numeric \eqn{n} by \eqn{n} non, negative spatial weights matrix,
#'   with zeros on the main diagonal.
#' @param length.out single, integer number, has to be at least 51 (due to order
#'   of approximation). Sets how fine the grid approximation is. Default
#'   value is 200.
#' @param rmin single number between -1 and 1. Sets the minimum range of the
#'   spatial autoregressive parameter \eqn{\rho}. Has to be lower than
#'   \code{rmax}. Default value is -1.
#' @param rmax single number between -1 and 1. Sets the maximum range of the
#'   spatial autoregressive parameter \eqn{\rho}. Has to be higher than
#'   \code{rmin}. Default value is 1.
#'
#' @return numeric \code{length.out} by  \code{2} matrix; the first column
#'   contains the approximated log-determinants the second column the \eqn{\rho} values
#'   ranging between \code{rmin} and \code{rmax}.
#' @export
#'
#' @references Barry, R. P., and Pace, R. K. (1999) Monte Carlo estimates of the
#' log determinant of large sparse matrices. \emph{Linear Algebra and its
#' applications}, \bold{289(1-3)}, 41-54.
logdetPaceBarry <-
  function(W,
           length.out = 200,
           rmin = -1,
           rmax = 1) {
    #rmin=-1 # <---- CHANGE: perhaps rmin=1e-5, to produce results only for 0 < rho < 1
    #rmax=1 # range of rho
    order = 50
    iter = 30 # <--- CHANGE: tuning parameters according to LeSage suggestions

    n = dim(W)[1]

    # Exact mom3ents from 1 to oexact
    td = matrix(c(0, sum(W ^ 2) / 2), length(c(0, sum(W ^ 2) / 2)), 1)

    oexact = length(td)

    # stochastic moments
    mavmomi = matrix(0, order, iter)

    for (j in 1:iter) {
      u = matrix(stats::rnorm(n, 0, 1), n, 1)
      v = u
      utu = t(u) %*% u

      for (i in 1:order) {
        v = W %*% v
        mavmomi[i, j] = as.double(n * ((t(u) %*% v) / (i * utu)))

      }
    }
    mavmomi[1:oexact, ] = td[, matrix(1, iter, 1)]

    # averages across iterations
    avmomi = as.matrix(rowMeans(mavmomi))

    # alpha matrix
    alpha = seq(rmin, rmax, length.out = length.out)
    valpha = matrixcalc::vandermonde.matrix(alpha, length(alpha))
    alomat = -valpha[, (2:(order + 1))]

    # Estimated ln|I-aD| using mixture of exact, stochastic moments
    # exact from 1 to oexact, stochastic from (oexact+1) to order
    lndetmat = alomat %*% avmomi

    return(cbind(lndetmat, rho =alpha))
  }
