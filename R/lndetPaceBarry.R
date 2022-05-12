#' Pace and Barry's log determinant approximation
#'
#' @param W An \eqn{n x n} spatial weights matrix
#' @param length.out Length of grid (defaults to 200)
#' @param rmin Minimum range of \eqn{\rho} (defaults to -1)
#' @param rmax Minimum \eqn{\rho} value (defaults to 1)
#'
#' @return A \eqn{length.out x 2} matrix; the first column contains the approximated log-determinants
#' the second column the rho values
#' @export lndetPaceBarry
lndetPaceBarry <-
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
