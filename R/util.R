

#' Specify Prior Distributions for Estimate W Models
#'
#' @param X An \eqn{n x q_x} matrix of independent variables that will spatially lagged (can be NULL)
#' @param Z An \eqn{n x q_z} matrix of independent variables that will not be spatially lagged
#' @param tt Number of time observations
#' @param W_prior An \eqn{n x n} matrix of priors for \eqn{W}
#' @param beta_prior_mean A \eqn{k x 1} matrix of \eqn{\beta} prior means
#' @param beta_prior_var A \eqn{k x k} matrix of \eqn{\beta} prior variances
#' @param sigma_a Sigma rate prior parameter (scalar)
#' @param sigma_b Sigma shape prior parameter (scalar)
#' @param rho_pr Prior for the four-part beta ditribution
#' @param griddy_n Number of griddy gibbs steps
#' @param rmin Minimum range of \eqn{\rho}
#' @param rmax Maximum range of \eqn{\rho}
#' @param bbinom_a Parameter a of sparsity prior
#' @param bbinom_b Parameter a of sparsity prior
#' @param bb_pr Should sparsity priors be used? (default: TRUE)
#' @param min_k Minimum number of a priori neighbours
#' @param max_k Maximum number of a priori neighbours
#' @param rjct_pr Should rejection priors be used? (default: TRUE)
#'
#' This function gives access to a larger set of prior distributions in case the default choice is unsatisfactory.
#'
#' @export
init_priors = function(X,Z=NULL,tt = 1,
                       W_prior = matrix(.5,nrow(X) / tt,nrow(X) / tt),
                       beta_prior_mean = matrix(0,
                                    ifelse(is.null(X),0,ncol(X))*2 + ifelse(is.null(Z),0,ncol(Z)) , 1),
                       beta_prior_var =
                         diag(ifelse(is.null(X),0,ncol(X))*2 + ifelse(is.null(Z),0,ncol(Z))) * 10,
                       sigma_a = .1, sigma_b = .1,
                       rho_pr = 1.1, griddy_n = 105, rmin = 0, rmax = 1,
                       bbinom_a = 1, bbinom_b = 1, bb_pr = TRUE,
                       min_k = 1, max_k = min(10, floor(nrow(X)/tt / 2)), rjct_pr = FALSE) {
  return(list(W_prior = W_prior,
              beta_prior_mean = beta_prior_mean,
              beta_prior_var = beta_prior_var,
              sigma_a = sigma_a, sigma_b = sigma_b,
              rho_pr = rho_pr, griddy_n = griddy_n, rmin = rmin, rmax = rmax,
              bbinom_a = bbinom_a, bbinom_b = bbinom_b, bb_pr = bb_pr,
              min_k = min_k, max_k = max_k, rjct_pr = rjct_pr))
}
