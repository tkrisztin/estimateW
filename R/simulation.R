#' Simulating from a data generating process
#'
#' This function can be used to generate a random data generating process for SDM,
#' SAR (if no \eqn{\beta_1} and \eqn{\beta_2} are supplied), SLX (if \eqn{\rho = 0}) type models.
#'
#' The generated spatial panel model takes the form
#'
#' \deqn{
#' Y = \rho W Y + X \beta_1 + W X \beta_2 + Z \beta_3 +  \epsilon,
#' }
#'
#' with \eqn{\epsilon \sim N(0,I_n\sigma^2)}. he function generates the \eqn{N \times 1} vector \eqn{Y}.
#' The elements of the explanatory variable matrices \eqn{X}
#' (\eqn{N \times k_1}) and \eqn{Z} (\eqn{N \times k_2}) are randomly generated from a Gaussian
#' distribution with zero mean and unity variance (\eqn{N(0,1)}).
#'
#' The non-negative, row-stochastic \eqn{n} by \eqn{n} matrix \eqn{W} is constructed using a k-nearest neighbor specification
#' based on a randomly generated spatial location pattern, with coordinates sampled from a standard normal distribution.
#'
#' Values for the parameters \eqn{\beta_1}, \eqn{\beta_2}, and \eqn{\beta_3}, as well as
#' \eqn{\rho} and \eqn{\sigma^2} have to be provided by the user.
#'
#' @param n Number of spatial observations \eqn{n}.
#' @param tt Number of time observations \eqn{T}.
#' @param rho The true \eqn{\rho} parameter
#' @param beta1 Vector of dimensions \eqn{k_1 \times 1}. Provides the values for \eqn{\beta_1} Defaults
#' to \code{c()}. Note: has to be fo same length as \eqn{\beta_2}.
#' @param beta2 Vector of dimensions \eqn{k_1 \times 1}. Provides the values for \eqn{\beta_2} Defaults
#' to \code{c()}. Note: has to be fo same length as \eqn{\beta_1}.
#' @param beta3 Vector of dimensions \eqn{k_2 \times 1}. Provides the values for \eqn{\beta_3} Defaults
#' to \code{c()}.
#' @param sigma2 The true \eqn{\sigma^2} parameter for the DGP. Has to be a scalar larger than zero.
#' @param n_neighbor Number of neighbors for the generated \eqn{n \times n} spatial weight \eqn{W} matrix
#' @param do_symmetric Should the generated spatial weight matrix be symmetric? (default: FALSE)
#' @param intercept Should the first column of \eqn{Z} be an intercept? Defaults to \code{FALSE}.
#' If \code{intercept = TRUE}, \eqn{\beta_3} has to be at least of length \code{1}.
#'
#' @return A list with the generated \eqn{X}, \eqn{Y} and \eqn{W} and a list of parameters.
#' @export sim_dgp
#'
#' @examples
#' # SDM data generating process
#' dgp_dat = sim_dgp(n =20, tt = 10, rho = .5, beta1 = c(1,-1),
#'                   beta2 = c(0,.5),beta3 = c(.2),sigma2 = .5)
sim_dgp= function(n, tt, rho, beta1 = c(), beta2 = c(), beta3 = c(),
                    sigma2 = .5, n_neighbor = 7, do_symmetric = FALSE,
                  intercept = FALSE) {
  smallk = length(beta1)
  if (length(beta2) != smallk) {stop("Beta2 has to be same length as beta1!")}
  k_dum = length(beta3)
  if (smallk == 0 && k_dum == 0) {stop("At least beta1 or beta2 has to be specified.")}


  # Randomly generated spatial patterns
  xy <- cbind(stats::runif(n),stats::runif(n))
  dist = as.matrix(dist(xy))
  diag(dist) = Inf
  dist = 1/dist
  W = t(apply(dist,c(1),function(x) {x[x<sort(x,decreasing = TRUE)[n_neighbor]] = 0; x[x>0] = 1; return(x)} ))
  if (do_symmetric) {
    W = t(W) + W
    W[W>0] = 1
  }
  W = as.matrix(W/rowSums(W))
  diag(W) = 0

  A <- as.matrix(diag(n) - rho*W)
  Ainv <- solve(A)
  #X <- cbind(1,matrix(stats::rnorm( n*(k-1) * tt),n*tt,k-1))
  X <- matrix(stats::rnorm( n*(smallk) * tt),n*tt,smallk)
  Z <- matrix(stats::rnorm( n*(k_dum) * tt),n*tt,k_dum)
  if (k_dum > 0 && intercept) {
    Z[,1] = 1
  }
  MU = stats::rnorm(n*tt, mean = 0,sd = sqrt(sigma2))
  if (smallk > 0) {
    MU = MU + cbind(X, kronecker(Matrix::.sparseDiagonal(tt),W) %*% X) %*% c(beta1,beta2)
  }
  if (k_dum > 0) {MU = MU + Z %*% beta3}
  Y = kronecker(Matrix::.sparseDiagonal(tt),Ainv) %*% MU

  ret <- list(Y = as.matrix(Y),
              X = as.matrix(X),
              Z = as.matrix(Z),
              W = W,
              para = list(rho = rho,
                          beta1 = beta1,
                          beta2 = beta2,
                          beta3 = beta3,
                          sigma2 = sigma2,
                          rho = rho,
                          xy = xy,
                          n_neighbor = n_neighbor,
                          do_symmetric = do_symmetric,
                          intercept = intercept))
  class(ret) = "sim_dgp"
  return(ret)
}
