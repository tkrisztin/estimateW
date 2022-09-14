

#' Specify Prior Distributions for the spatial weight matrix
#'
#' @param n The number of observations
#' @param W_prior An \eqn{n x n} matrix of priors for \eqn{W}
#' @param symmetric_prior Should the estimated \eqn{W} matrix be symmetric (default: TRUE)
#' @param row_standardized_prior Should the estimated \eqn{W} matrix be row-standardized (default: TRUE)
#' @param bbinom_a_prior Parameter a of sparsity prior
#' @param bbinom_b_prior Parameter a of sparsity prior
#' @param use_bbinom_prior Should sparsity priors be used? (default: TRUE)
#' @param min_neighbors Minimum number of a priori neighbours (default: 0)
#' @param max_neighbors Maximum number of a priori neighbours (default: n-1)
#'
#' This function gives access to a larger set of prior distributions in case the default choice is unsatisfactory.
#'
#' @export
W_priors = function(n,
                       W_prior = matrix(.5,n,n),
                       symmetric_prior = TRUE,row_standardized_prior = TRUE,
                       bbinom_a_prior = 1, bbinom_b_prior = 1, use_bbinom_prior = TRUE,
                       min_neighbors = 0, max_neighbors = n-1) {
  # Ensure diagonal of W_prior is zero
  diag(W_prior) <- 0
  if (min_neighbors != 0 || max_neighbors != n-1) {use_reject_prior = TRUE} else {use_reject_prior = FALSE}
  return(list(W_prior = W_prior,symmetric_prior = symmetric_prior,row_standardized_prior = row_standardized_prior,
              bbinom_a_prior = bbinom_a_prior, bbinom_b_prior = bbinom_b_prior,
              use_bbinom_prior = use_bbinom_prior,
              min_neighbors = min_neighbors, max_neighbors = max_neighbors,
              use_reject_prior = use_reject_prior))
}

#' Specify Prior Distributions for the spatial autoregressive coefficient
#'
#' @param rho_a_prior Prior for the four-part beta ditribution
#' @param rho_b_prior Prior for the four-part beta ditribution
#' @param rho_min Minimum range of \eqn{\rho} (default: 0)
#' @param rho_max Maximum range of \eqn{\rho} (default: 1)
#' @param init_rho_scale For Metropolis-Hastings step the initial candidate variance (default: 1)
#' @param griddy_n Number of griddy gibbs steps
#' @param use_griddy_gibbs Should griddy-Gibbs be used for \eqn{\rho} estimation?
#'
#' Does not work if \code{row_standardized_prior = FALSE} is specified in the \eqn{W} prior specification. Main advantage is that less draws are required for \eqn{\rho}
#'
#' This function gives access to a larger set of prior distributions for \eqn{\rho} in case the default choice is unsatisfactory.
#'
#' @export
rho_priors = function(rho_a_prior = 1.1, rho_b_prior = 1.1,
                      rho_min = 0, rho_max = 1, init_rho_scale = 1,
                      griddy_n = 105, use_griddy_gibbs = TRUE) {
  init.rho = stats::runif(1,rho_min,rho_max)
  return(list(rho_a_prior = rho_a_prior, rho_b_prior = rho_b_prior,
              griddy_n = griddy_n, rho_min = rho_min, rho_max = rho_max, init_rho_scale,
              use_griddy_gibbs = use_griddy_gibbs,init.rho=init.rho))
}

#' Specify Prior Distributions for the slope parameters
#'
#' @param k The total number of coefficients in the model.
#' @param beta_mean_prior A \eqn{k x 1} matrix of \eqn{\beta} prior means (default: vector of zeros)
#' @param beta_var_prior A \eqn{k x k} matrix of \eqn{\beta} prior variances (default: \eqn{10})
#'
#' This function allows the user to specify priors for the slope coefficients.
#'
#' @export
beta_priors = function(k,
                        beta_mean_prior = matrix(0,k , 1),
                        beta_var_prior = diag(k) * 10) {
  beta_var_prior_inv <- solve(beta_var_prior)
  return(list(beta_mean_prior = beta_mean_prior,
              beta_var_prior = beta_var_prior,
              beta_var_prior_inv = beta_var_prior_inv))
}

#' Specify Prior Distributions for the error variance
#'
#' @param sigma_rate_prior Sigma rate prior parameter (scalar), default: \eqn{.1}
#' @param sigma_shape_prior Sigma shape prior parameter (scalar), default: \eqn{.1}
#'
#' This function allows the user to specify priors for the error variance \eqn{\sigma^2}.
#'
#' @export
sigma_priors = function(sigma_rate_prior = .1, sigma_shape_prior = .1) {
  return(list(sigma_rate_prior = sigma_rate_prior, sigma_shape_prior = sigma_shape_prior))
}
