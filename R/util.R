

#' Specify Prior Distributions for the spatial weight matrix
#'
#' @param n The number of observations
#' @param W_prior An \eqn{n x n} matrix of priors for \eqn{W}
#' @param SYMMETRIC Should the estimated \eqn{W} matrix be symmetric (default: TRUE)
#' @param ROW_STANDARDIZED Should the estimated \eqn{W} matrix be row-standardized (default: TRUE)
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
W_priors = function(n,
                       W_prior = matrix(.5,n,n),
                       SYMMETRIC = TRUE,ROW_STANDARDIZED = TRUE,
                       bbinom_a = 1, bbinom_b = 1, bb_pr = TRUE,
                       min_k = 1, max_k = min(10, floor(n / 2)), rjct_pr = FALSE) {
  # Ensure diagonal of W_prior is zero
  diag(W_prior) <- 0
  return(list(W_prior = W_prior,SYMMETRIC = SYMMETRIC,ROW_STANDARDIZED = ROW_STANDARDIZED,
              bbinom_a = bbinom_a, bbinom_b = bbinom_b, bb_pr = bb_pr,
              min_k = min_k, max_k = max_k, rjct_pr = rjct_pr))
}

#' Specify Prior Distributions for the spatial autoregressive coefficient
#'
#' @param rho_pr Prior for the four-part beta ditribution
#' @param griddy_n Number of griddy gibbs steps
#' @param rmin Minimum range of \eqn{\rho}
#' @param rmax Maximum range of \eqn{\rho}
#' @param GRIDDY_GIBBS Should griddy-Gibbs be used for \eqn{\rho} estimation?
#' Does not work if \code{ROW_STANDARDIZED = FALSE} is specified in the \eqn{W} prior specification. Main advantage is that less draws are required for \eqn{\rho}
#'
#' This function gives access to a larger set of prior distributions for \eqn{\rho} in case the default choice is unsatisfactory.
#'
#' @export
rho_priors = function(rho_pr = 1.1, griddy_n = 105, rmin = 0, rmax = 1, GRIDDY_GIBBS = TRUE) {
  return(list(rho_pr = rho_pr, griddy_n = griddy_n, rmin = rmin, rmax = rmax, GRIDDY_GIBBS = GRIDDY_GIBBS))
}

#' Specify Prior Distributions for the slope parameters
#'
#' @param k The total number of coefficients in the model.
#' @param beta_prior_mean A \eqn{k x 1} matrix of \eqn{\beta} prior means (default: vector of zeros)
#' @param beta_prior_var A \eqn{k x k} matrix of \eqn{\beta} prior variances (default: \eqn{10})
#'
#' This function allows the user to specify priors for the slope coefficients.
#'
#' @export
beta_priors = function(k,
                        beta_prior_mean = matrix(0,k , 1),
                        beta_prior_var = diag(k) * 10) {
  beta_prior_var_inv <- solve(beta_prior_var)
  return(list(beta_prior_mean = beta_prior_mean,
              beta_prior_var = beta_prior_var,
              beta_prior_var_inv = beta_prior_var_inv))
}

#' Specify Prior Distributions for the error variance
#'
#' @param sigma_a Sigma rate prior parameter (scalar), default: \eqn{.1}
#' @param sigma_b Sigma shape prior parameter (scalar), default: \eqn{.1}
#'
#' This function allows the user to specify priors for the error variance \eqn{\sigma^2}.
#'
#' @export
sigma_priors = function(sigma_a = .1, sigma_b = .1) {
  return(list(sigma_a = sigma_a, sigma_b = sigma_b))
}
