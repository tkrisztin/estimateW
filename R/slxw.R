
#' A Markov Chain Monte Carlo (MCMC) sampler for the panel spatial SLX model with unknown spatial weight matrix
#'
#' The sampler uses independent Normal-inverse-Gamma priors for the slope and variance parameters.
#' It is a wrapper around \code{\link{W_sampler}}.
#'
#' The considered spatial panel SLX model with unknown (\eqn{n} by \eqn{n}) spatial weight
#' matrix \eqn{W} takes the form:
#'
#' \deqn{
#'  Y_t = X_t \beta_1 + W X_t \beta_2 + Z \beta_3 + \varepsilon_t,
#'  }
#'
#' with \eqn{\varepsilon_t \sim N(0,I_n \sigma^2)} and \eqn{W = f(\Omega)}. The \eqn{n} by \eqn{n}
#' matrix \eqn{\Omega} is an unknown binary adjacency matrix with zeros on the main diagonal and
#' \eqn{f(\cdot)} is the (optional) row-standardization function.
#'
#' \eqn{Y_t} (\eqn{n \times 1}) collects the \eqn{n} cross-sectional (spatial) observations for time
#' \eqn{t=1,...,T}. \eqn{X_t} (\eqn{n \times k_1}) and \eqn{Z_t} (\eqn{n \times k_2}) are
#' matrices of explanatory variables, where the former will also be spatially lagged. \eqn{\beta_1}
#' (\eqn{k_1 \times 1}), \eqn{\beta_2} (\eqn{k_1 \times 1}) and \eqn{\beta_3} (\eqn{k_2 \times 1})
#' are unknown slope parameter vectors.
#'
#' After vertically staking the \eqn{T} cross-sections  \eqn{Y=[Y_1',...,Y_T']'} (\eqn{N \times 1}),
#' \eqn{X=[X_1',...,X_T']'} (\eqn{N \times k_1}) and \eqn{Z=[Z_1', ..., Z_T']'} (\eqn{N \times k_2}),
#' with \eqn{N=nT}. The final model can be expressed as:
#'
#' \deqn{
#'  Y = X \beta_1 + \tilde{W} X \beta_2 + Z \beta_3 + \varepsilon,
#' }
#'
#' where \eqn{\tilde{W}=I_T \otimes W} and \eqn{\varepsilon \sim N(0,I_N \sigma^2)}. Note that the input
#' data matrices have to be ordered first by the cross-sectional spatial units and then stacked by time.
#'
#' Estimation usually even works well in cases of \eqn{n >> T}. However, note that for applications with \eqn{n > 200} the
#' estimation process becomes computationally demanding and slow. Consider in this case reducing \code{niter} and
#' \code{nretain} and carefully check whether the posterior chains have converged.
#'
#'
#' @inheritParams sdmw
#'
#' @return List with posterior samples for the slope parameters, \eqn{\sigma^2}, \eqn{W},
#' and average direct, indirect, and total effects.
#' @export slxw
#'
#' @examples
#' set.seed(123)
#' n = 20; tt = 10
#' dgp_dat = sim_dgp(n = 20, tt = 10, rho = 0, beta1 = c(1,-1),
#'                   beta2 = c(0,.5), beta3 = c(.2), sigma2 = .5)
#' res = slxw(Y = dgp_dat$Y, tt = tt, X = dgp_dat$X, Z = dgp_dat$Z,
#'                   niter = 20, nretain = 10)
slxw <- function(Y, tt, X = matrix(0,nrow(Y),0),Z = matrix(1,nrow(Y),1), niter = 100, nretain = 50,
                 W_prior = W_priors(n = nrow(Y)/tt),
                 beta_prior = beta_priors(k = ncol(X)*2 + ncol(Z)),sigma_prior = sigma_priors()) {
  if (ncol(X) == 0 && ncol(Z) == 0) {stop("Error: At least either X or Z matrix have to be provided.")}
  origX = X

  smallk = ncol(X)
  if (smallk>0 && is.null(colnames(X))) {
    colnames(X) = paste0("X",1:smallk)
  }
  k_dum = ncol(Z)
  if (k_dum>0 && is.null(colnames(Z))) {
    colnames(Z) = paste0("Z",1:k_dum)
  }
  if (smallk>0) {varnames = c(colnames(X), paste0("W_",colnames(X)))} else {varnames = c()}
  if (k_dum > 0) varnames = c(varnames,colnames(Z))

  ndiscard <- niter - nretain
  k <- smallk*2 + k_dum
  n <- nrow(X) / tt

  # map variable positions for spatial effects
  ind_baseX = ind_WX = ind_lagFX = c()
  if (smallk > 0) {
    ind_baseX = c(1:smallk)
    # the columns of XX that are spatially lagged
    ind_WX = c(1:smallk) + smallk
    # the spatial FX corresponding to these
    ind_lagFX = 1:smallk
  }
  if (k_dum > 0) {ind_baseX = c(ind_baseX,(2*smallk + 1):k)}

  # save the posterior draws here
  postb <- matrix(0, k, nretain)
  rownames(postb) <- varnames
  posts <- matrix(0, 1, nretain); rownames(posts) = "sigma"
  postw <- array(0, c(n, n, nretain))
  postwprob <- array(0, c(n, n, nretain))

  sampler_W = W_sampler$new(W_prior)
  sampler_beta = beta_sampler$new(beta_prior)
  sampler_sigma = sigma_sampler$new(sigma_prior)

  curr_WX = as.matrix(kronecker(Matrix::.sparseDiagonal(tt),sampler_W$curr_w) %*% X)
  tX = cbind(X,curr_WX,Z)
  tY <- matrix(Y, n, tt)
  curr_mu = curr_mu_lag = matrix(0,n,tt)

  ### Gibbs sampling
  pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  for (iter in 1:niter) {

    # draw beta
    sampler_beta$sample(Y,tX,sampler_sigma$curr_sigma)
    curr_xb <- tX %*% sampler_beta$curr_beta
    curr_txb <- matrix(curr_xb, n, tt)
    if (smallk > 0) {
      curr_mu = matrix(tX[,-ind_WX] %*% sampler_beta$curr_beta[-ind_WX],n,tt)
      curr_mu_lag = matrix(X %*% sampler_beta$curr_beta[ind_WX],n,tt)
    } else {
      curr_mu = matrix(curr_xb,n,tt)
    }

    # draw sigma
    sampler_sigma$sample(Y,curr_xb)

    # Gibbs step for W - element-wise
    sampler_W$sample(Y = tY,curr_sigma = sampler_sigma$curr_sigma,
                     mu = curr_mu,lag_mu = curr_mu_lag)
    curr_WX = as.matrix(kronecker(Matrix::.sparseDiagonal(tt),sampler_W$curr_w) %*% X)
    tX = cbind(X,curr_WX,Z)

    # we are past the burn-in, save the draws
    if (iter > ndiscard) {
      s <- iter - ndiscard
      postb[, s] <- sampler_beta$curr_beta
      posts[s] <- sampler_sigma$curr_sigma
      postw[, , s] <- sampler_W$curr_w
    }
    utils::setTxtProgressBar(pb,iter)
  }
  close(pb)

  ret = list( Y = Y, X = X,Z = Z,
             postb = postb, posts = posts, postw = postw,
             W_prior = W_prior,
             beta_prior = beta_prior,sigma_prior = sigma_prior,
             param = list(niter = niter, nretain = nretain)
  )
  class(ret) = "estimateW"
  return(ret)
}
