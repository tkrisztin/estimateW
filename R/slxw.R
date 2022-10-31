
#' A sampler for estimating the W matrix in an SLX type model
#'
#' The model takes the form \eqn{Y = X \beta_1 + f(\Omega)X \beta_2 + Z \beta_3 +  \epsilon}, with \eqn{\epsilon \sim N(0,I\sigma^2)}
#'
#' @inheritParams sdmw
#'
#' @return List with posterior samples for \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, and \eqn{\sigma^2}.
#' @export slxw
#'
#' @examples
#' set.seed(123)
#' n = 20; tt = 10
#' dgp_dat = sim_dgp(n = 20, tt = 10, rho = 0, beta1 = c(1,-1),
#'                   beta2 = c(0,.5), beta3 = c(.2), sigma2 = .5)
#' res = slxw(Y = dgp_dat$Y, tt = tt, X = dgp_dat$X, Z = dgp_dat$Z,
#'                   niter = 20, nretain = 10)
slxw <- function(Y, tt, X = matrix(0,nrow(Y),0),Z = matrix(1,nrow(Y),1), niter = 1000, nretain = 250,
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
