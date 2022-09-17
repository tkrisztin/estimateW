#' Simulating from an SDM data generating process
#'
#' The model takes the form \eqn{Y = \rho W Y + X \beta_1 + W X \beta_2 + Z \beta_3 +  \epsilon}, with \eqn{\epsilon \sim N(0,\sigma^2)}
#'
#' @param n Number of spatial observations
#' @param tt Number of time observations
#' @param rho The true rho parameter
#' @param beta1 The \eqn{k_1 x 1} vector of true \eqn{\beta_1}
#' @param beta2 The \eqn{k_1 x 1} vector of true \eqn{\beta_2} (with \eqn{k_1 = k_2})
#' @param beta3 The \eqn{k_3 x 1} vector of true \eqn{\beta_3}
#' @param sigma2 The true \eqn{\sigma^2} parameter
#' @param n_neighbor Number of neighbours for the \eqn{W} matrix
#' @param do_symmetric Should the generated \eqn{W} matrix be symmetric? (default: FALSE)
#'
#' @return A list with the generated \eqn{X}, \eqn{Y} and \eqn{W} and a list of parameters.
#' @export sim_dgp
#'
#' @examples
#' # SDM data generating process
#' dgp_dat = sim_dgp(n =20, tt = 10, rho = .5, beta1 = c(1,-1),
#'                   beta2 = c(0,.5),beta3 = c(.2),sigma2 = .5)
sim_dgp= function(n, tt, rho, beta1 = c(), beta2 = c(), beta3 = c(),
                    sigma2 = .5, n_neighbor = 7, do_symmetric = FALSE) {
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
                          do_symmetric = do_symmetric))
  class(ret) = "sim_dgp"
  return(ret)
}
