

#' Set prior specifications for the spatial weight matrix
#'
#' Set prior specifications for the \eqn{n} by \eqn{n} spatial weight matrix \eqn{W=f(\Omega)},
#' where \eqn{\Omega} is an \eqn{n} by \eqn{n} unknown binary adjacency matrix (with zeros on the
#' main diagonal), and \eqn{f()} denotes the (optional) row-standardization function
#'
#'
#' @param n The number of spatial observations
#' @param W_prior An \eqn{n} by \eqn{n} matrix of prior inclusion probabilities for \eqn{W}
#' @param symmetric_prior Binary value. Should the estimated adjacency matrix \eqn{\Omega} be symmetric (default: FALSE)?
#' if TRUE: \eqn{\Omega} is forced symmetric; if FALSE: \eqn{\Omega} not necessarily symmetric.
#' @param row_standardized_prior Binary value. Should the estimated \eqn{W} matrix be row-standardized (default: TRUE)?
#' if TRUE: row-stochastic \eqn{W}; if FALSE: \eqn{W} not row-standardized.
#' @param nr_neighbors_prior An \eqn{n \times 1} vector of prior inclusion probabilities on the number of neighbors.
#' Defaults to a \code{\link{bbinompdf}} prior, with prior parameters \eqn{a = 1}, \eqn{b = 1} and
#' no minimum or maximum restrictions on the number of neighbors. A flat prior would be an
#' \eqn{n \times 1} vector of ones.
#'
#' @export
W_priors = function(n,
                       W_prior = matrix(.5,n,n),
                       symmetric_prior = FALSE,row_standardized_prior = TRUE,
                       nr_neighbors_prior = bbinompdf(0:(n-1), nsize = n - 1,
                                                      a = 1,b = 1,
                                                      min_k = 0,max_k = n-1)
                    ) {
  # Ensure diagonal of W_prior is zero
  diag(W_prior) <- 0
  return(list(W_prior = W_prior,symmetric_prior = symmetric_prior,
              row_standardized_prior = row_standardized_prior,
              nr_neighbors_prior = nr_neighbors_prior))
}

#' Specify prior for the spatial autoregressive parameter and sampling settings
#'
#' Specify prior for the spatial autoregressive parameter and sampling settings
#'
#' @param rho_a_prior Single number. Prior hyperparameter for the four-parameter beta distribution \code{\link{betapdf}}.
#' Defaults to 1.
#' @param rho_b_prior Single number. Prior hyperparameter for the four-parameter beta distribution \code{\link{betapdf}}.
#' Defaults to 1.
#' @param rho_min Minimum value for \eqn{\rho} (default: 0)
#' @param rho_max Maximum value for \eqn{\rho} (default: 1)
#' @param init_rho_scale For Metropolis-Hastings step the initial candidate variance (default: 1)
#' @param griddy_n single integer number. Sets how fine the grid approximation is. Default
#'   value is 60.
#' @param use_griddy_gibbs Binary value. Should griddy-Gibbs be used for \eqn{\rho} estimation?
#' \code{use_griddy_gibbs=TRUE} does not work if \code{row_standardized_prior = FALSE} is specified in the \eqn{W} prior specification.
#' if TRUE: griddy-Gibbs step for sampling \eqn{\rho}; if FALSE: tuned random-walk Metropolis-Hastings step
#'
#' @param mh_tune_low Lower bound of acceptance rate for Metropolis-Hastings tuning
#' (used if \code{use_griddy_gibbs==FALSE})
#' @param mh_tune_high Upper bound of acceptance rate for Metropolis-Hastings tuning
#' (used if \code{use_griddy_gibbs==FALSE})
#' @param mh_tune_scale Scaling factor for Metropolis-Hastings tuning
#' (used if \code{use_griddy_gibbs==FALSE})
#'
#' @export
rho_priors = function(rho_a_prior = 1, rho_b_prior = 1,
                      rho_min = 0, rho_max = 1, init_rho_scale = 1,
                      griddy_n = 60, use_griddy_gibbs = TRUE,
                      mh_tune_low = .4,mh_tune_high = .6, mh_tune_scale = .1) {
  return(list(rho_a_prior = rho_a_prior, rho_b_prior = rho_b_prior,
              griddy_n = griddy_n, rho_min = rho_min, rho_max = rho_max,
              init_rho_scale = init_rho_scale,
              use_griddy_gibbs = use_griddy_gibbs,
              mh_tune_low = mh_tune_low,mh_tune_high = mh_tune_high,
              mh_tune_scale = mh_tune_scale))
}

#' Set prior specifications for the slope parameters
#'
#' This function allows the user to specify custom values for Gaussian priors on the slope parameters.
#'
#' For the slope parameters \eqn{\beta} the package uses common Normal
#' prior specifications. Specifically,  \eqn{p(\beta)\sim\mathcal{N}(\underline{\mu}_\beta,\underline{V}_\beta)}.
#'
#' This function allows the user to specify custom values for the prior hyperparameters \eqn{\underline{\mu}_\beta}
#' and \eqn{\underline{V}_\beta}. The default values correspond to weakly informative Gaussian priors with mean
#' zero and a diagonal prior variance-covariance matrix with \eqn{100} on the main diagonal.
#'
#' @param k The total number of slope parameters in the model.
#' @param beta_mean_prior numeric \eqn{k} by \eqn{1} matrix of prior means \eqn{\underline{\mu}_\beta}.
#' @param beta_var_prior A \eqn{k} by \eqn{k} matrix of prior variances \eqn{\underline{V}_\beta}. Defaults to a
#' diagonal matrix with \code{100} on the main diagonal.
#'
#' @return A list with the prior mean vector (\code{beta_mean_prior}), the prior variance matrix
#' (\code{beta_var_prior}) and the inverse of the prior variance matrix (\code{beta_var_prior_inv}).
#'
#' @export
beta_priors = function(k,
                        beta_mean_prior = matrix(0,k , 1),
                        beta_var_prior = diag(k) * 100) {
  beta_var_prior_inv <- solve(beta_var_prior)
  return(list(beta_mean_prior = beta_mean_prior,
              beta_var_prior = beta_var_prior,
              beta_var_prior_inv = beta_var_prior_inv))
}

#' Set prior specification for the error variance using an inverse Gamma distribution
#'
#' @param sigma_rate_prior Sigma rate prior parameter (scalar), default: \eqn{0.001}.
#' @param sigma_shape_prior Sigma shape prior parameter (scalar), default: \eqn{0.001}.
#'
#' This function allows the user to specify priors for the error variance \eqn{\sigma^2}.
#'
#' @export
sigma_priors = function(sigma_rate_prior = 0.001, sigma_shape_prior = 0.001) {
  return(list(sigma_rate_prior = sigma_rate_prior, sigma_shape_prior = sigma_shape_prior))
}
