#' A Markov Chain Monte Carlo (MCMC) sampler for the panel spatial autoregressive model (SAR) with unknown spatial weight matrix
#'
#' The sampler uses independent  Normal-inverse-Gamma priors for the slope and variance parameters, as well as a four-parameter
#' beta prior for the spatial autoregressive parameter \eqn{\rho}.
#' This is a wrapper function calling \code{\link{sdmw}} with no spatially lagged exogenous variables.
#'
#' The considered panel spatial autoregressive model (SAR) with unknown (\eqn{n} by \eqn{n}) spatial weight
#' matrix \eqn{W} takes the form:
#'
#' \deqn{
#'  Y_t = \rho W Y_t + Z \beta + \varepsilon_t,
#'  }
#'
#' with \eqn{\varepsilon_t \sim N(0,I_n \sigma^2)} and \eqn{W = f(\Omega)}. The \eqn{n} by \eqn{n}
#' matrix \eqn{\Omega} is an unknown binary adjacency matrix with zeros on the main diagonal and
#' \eqn{f(\cdot)} is the (optional) row-standardization function. \eqn{\rho} is a scalar spatial autoregressive parameter.
#'
#' \eqn{Y_t} (\eqn{n \times 1}) collects the \eqn{n} cross-sectional (spatial) observations for time
#' \eqn{t=1,...,T}. \eqn{Z_t} (\eqn{n \times k_3}) is a matrix of explanatory variables.
#' \eqn{\beta} (\eqn{k_3 \times 1}) is an unknown slope parameter vector.
#'
#' After vertically staking the \eqn{T} cross-sections  \eqn{Y=[Y_1',...,Y_T']'} (\eqn{N \times 1}),
#' and \eqn{Z=[Z_1', ..., Z_T']'} (\eqn{N \times k_3}), with \eqn{N=nT}. The final model can be expressed as:
#'
#' \deqn{
#'  Y = \rho \tilde{W}Y + Z \beta + \varepsilon,
#' }
#'
#' where \eqn{\tilde{W}=I_T \otimes W} and \eqn{\varepsilon \sim N(0,I_N \sigma^2)}. Note that the input
#' data matrices have to be ordered first by the cross-sectional spatial units and then stacked by time.
#'
#' Estimation usually even works well in cases of \eqn{n >> T}. However, note that for applications with \eqn{n > 200} the
#' estimation process becomes computationally demanding and slow. Consider in this case reducing \code{niter} and
#' \code{nretain} and carefully check whether the posterior chains have converged.
#'
#' @inheritParams sdmw
#' @param Z numeric \eqn{N \times k_3} design matrix of independent variables.
#' The default value is a \eqn{N \times 1} vector of ones (i.e. an intercept for the model).
#' @param beta_prior list containing priors for the slope coefficients \eqn{\beta},
#' generated by the smart constructor \code{\link{beta_priors}}.
#'
#' @return List with posterior samples for the slope parameters, \eqn{\rho}, \eqn{\sigma^2}, \eqn{W},
#' and average direct, indirect, and total effects.
#' @export sarw
#'
#' @examples
#' n = 20; tt = 10
#' dgp_dat = sim_dgp(n =n, tt = tt, rho = .5, beta3 = c(.5,1),
#'             sigma2 = .05,n_neighbor = 3,intercept = TRUE)
#' res = sarw(Y = dgp_dat$Y,tt = tt,Z = dgp_dat$Z,niter = 20,nretain = 10)
sarw <- function(Y, tt, Z, niter = 100, nretain = 50,
                 W_prior = W_priors(n = nrow(Y)/tt),rho_prior = rho_priors(),
                 beta_prior = beta_priors(k = ncol(Z)),sigma_prior = sigma_priors()) {
  ret = sdmw(Y = Y, tt =tt, X = matrix(0,nrow(Y),0),Z = Z, niter = niter, nretain = nretain,
             W_prior = W_prior,rho_prior = rho_prior,
             beta_prior = beta_prior,sigma_prior = sigma_prior)
  ret$model_type = "SAR"
  return(ret)
}


#' A Markov Chain Monte Carlo (MCMC) sampler for the panel spatial Durbin model (SDM) with unknown spatial weight matrix
#'
#' The sampler uses independent Normal-inverse-Gamma priors for the slope and variance parameters, as well as a four-parameter
#' beta prior for the spatial autoregressive parameter \eqn{\rho}. It is a wrapper around \code{\link{W_sampler}}.
#'
#' The considered panel spatial Durbin model (SDM) with unknown (\eqn{n} by \eqn{n}) spatial weight
#' matrix \eqn{W} takes the form:
#'
#' \deqn{
#'  Y_t = \rho W Y_t + X_t \beta_1 + W X_t \beta_2 + Z \beta_3 + \varepsilon_t,
#'  }
#'
#' with \eqn{\varepsilon_t \sim N(0,I_n \sigma^2)} and \eqn{W = f(\Omega)}. The \eqn{n} by \eqn{n}
#' matrix \eqn{\Omega} is an unknown binary adjacency matrix with zeros on the main diagonal and
#' \eqn{f(\cdot)} is the (optional) row-standardization function. \eqn{\rho} is a scalar spatial autoregressive parameter.
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
#'  Y = \rho \tilde{W}Y + X \beta_1 + \tilde{W} X \beta_2 + Z \beta_3 + \varepsilon,
#' }
#'
#' where \eqn{\tilde{W}=I_T \otimes W} and \eqn{\varepsilon \sim N(0,I_N \sigma^2)}. Note that the input
#' data matrices have to be ordered first by the cross-sectional spatial units and then stacked by time.
#'
#' Estimation usually even works well in cases of \eqn{n >> T}. However, note that for applications with \eqn{n > 200} the
#' estimation process becomes computationally demanding and slow. Consider in this case reducing \code{niter} and
#' \code{nretain} and carefully check whether the posterior chains have converged.
#'
#' @inheritParams sdm
#' @param niter single number greater or equal to 1, indicating the total number of draws.
#' Will be automatically coerced to integer. The default value is 100.
#' @param nretain single number greater or equal to 0, indicating the number of draws
#' kept after the burn-in. Will be automatically coerced to integer. The default value is 50.
#' @param W_prior list containing prior settings for estimating the spatial weight matrix \eqn{W}.
#' Generated by the smart constructor \code{\link{W_priors}}.
#'
#' @return List with posterior samples for the slope parameters, \eqn{\rho}, \eqn{\sigma^2}, \eqn{W},
#' and average direct, indirect, and total effects.
#'
#' @export sdmw
#'
#' @examples
#' n = 20; tt = 10
#' dgp_dat = sim_dgp(n =n, tt = tt, rho = .75, beta1 = c(.5,1),beta2 = c(-1,.5),
#'             beta3 = c(1.5), sigma2 = .05,n_neighbor = 3,intercept = TRUE)
#' res = sdmw(Y = dgp_dat$Y,tt = tt,X = dgp_dat$X,Z = dgp_dat$Z,niter = 20,nretain = 10)
sdmw <- function(Y, tt, X = matrix(0,nrow(Y),0),Z = matrix(1,nrow(Y),1), niter = 100, nretain = 50,
                  W_prior = W_priors(n = nrow(Y)/tt),rho_prior = rho_priors(),
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

  if (!W_prior$row_standardized_prior && rho_prior$use_griddy_gibbs) {
    stop("No support for rho prior Griddy Gibbs without row-standardization prior for W!")
  }

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
  posts <- matrix(0, 1, nretain); rownames(posts) = "sigma2"
  postr <- matrix(0, 1, nretain); rownames(postr) = "rho"
  postw <- array(0, c(n, n, nretain))
  postwprob <- array(0, c(n, n, nretain))

  post.direct <- matrix(0, smallk + k_dum, nretain)
  post.indirect <- matrix(0, smallk + k_dum, nretain)
  post.total <- matrix(0, smallk + k_dum, nretain)
  rownames(post.direct) <- rownames(post.indirect) <- rownames(post.total) <- varnames[ind_baseX]

  # initialize wdraws
  sampler_rho <- rho_sampler$new(rho_prior)
  sampler_W = W_sampler$new(W_prior,sampler_rho$curr_rho)
  sampler_beta = beta_sampler$new(beta_prior)
  sampler_sigma = sigma_sampler$new(sigma_prior)

  curr.WX = as.matrix(kronecker(Matrix::.sparseDiagonal(tt),sampler_W$curr_w) %*% X)
  tX = cbind(X,curr.WX,Z)
  tY <- matrix(Y, n, tt)
  curr_mu = curr_mu_lag = matrix(0,n,tt)

  ### Gibbs sampling
  pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  for (iter in 1:niter) {
    Ay <- matrix(sampler_W$curr_A %*% tY, n * tt, 1)

    # draw beta
    sampler_beta$sample(Ay,tX,sampler_sigma$curr_sigma)
    curr.xb <- tX %*% sampler_beta$curr_beta
    curr.txb <- matrix(curr.xb, n, tt)
    if (smallk > 0) {
      curr_mu = matrix(tX[,-ind_WX] %*% sampler_beta$curr_beta[-ind_WX],n,tt)
      curr_mu_lag = matrix(X %*% sampler_beta$curr_beta[ind_WX],n,tt)
    } else {
      curr_mu = matrix(curr.xb,n,tt)
    }

    # draw sigma
    sampler_sigma$sample(Ay,curr.xb)

    ## Griddy-Gibbs step for rho
    sampler_rho$setW(newW = sampler_W$curr_w,
                       newLogdet = sampler_W$curr_logdet,
                       newA = sampler_W$curr_A, newAI = sampler_W$curr_AI)
    sampler_rho$sample(tY,curr.txb,sampler_sigma$curr_sigma)
    if (iter > (ndiscard / 2)) {
      sampler_rho$stopMHtune()
    }


    # Gibbs step for W - element-wise
    sampler_W$set_rho(new_rho = sampler_rho$curr_rho,
                      newLogdet = sampler_rho$curr_logdet,
                      newA = sampler_rho$curr_A, newAI = sampler_rho$curr_AI)
    # sampler_W$sample_fast(Y = tY,curr_sigma = sampler_sigma$curr_sigma,
    #                  mu = curr_mu,lag_mu = curr_mu_lag)
    sampler_W$sample(Y = tY,curr_sigma = sampler_sigma$curr_sigma,
                          mu = curr_mu,lag_mu = curr_mu_lag)
    curr.WX = as.matrix(kronecker(Matrix::.sparseDiagonal(tt),sampler_W$curr_w) %*% X)
    tX = cbind(X,curr.WX,Z)

    # we are past the burn-in, save the draws
    if (iter > ndiscard) {
      s <- iter - ndiscard
      postb[, s] <- sampler_beta$curr_beta
      posts[s] <- sampler_sigma$curr_sigma
      postr[s] <- sampler_rho$curr_rho
      postw[, , s] <- sampler_W$curr_w

      post.direct[, s] <- sum(diag(sampler_W$curr_AI)) / n * sampler_beta$curr_beta[ind_baseX]
      post.total[, s] <- sum(sampler_W$curr_AI) / n * sampler_beta$curr_beta[ind_baseX]
      # if we have WX
      if (smallk > 0) {
        post.direct[ind_lagFX,s] = post.direct[ind_lagFX,s] +
          sum(diag(sampler_W$curr_AI))/n * sampler_beta$curr_beta[ind_WX]
        post.total[ind_lagFX,s] = post.total[ind_lagFX,s] +
          sum(sampler_W$curr_AI)/n * sampler_beta$curr_beta[ind_WX]
      }
      post.indirect[, s] <- post.total[, s] - post.direct[, s]
    }
    utils::setTxtProgressBar(pb,iter)
  }
  close(pb)

  ret = list(Y = Y, X = X,Z = Z,
    postb = postb, posts = posts, postr = postr, postw = postw,
    post.direct = post.direct, post.indirect = post.indirect, post.total = post.total,
    W_prior = W_prior,rho_prior = rho_prior,
    beta_prior = beta_prior,sigma_prior = sigma_prior,
    param = list(niter = niter, nretain = nretain),
    model_type = "SDM"
  )
  class(ret) = "estimateW"
  return(ret)
}



