#' A sampler for estimating the W matrix in a SAR type model \eqn{Y = \rho WY + X \beta + \epsilon}, with \eqn{\epsilon \sim N(0,\sigma^2)}
#'
#' @param Y An \eqn{n x 1} matrix of dependent variables
#' @param X An \eqn{n x k} matrix of independent variables
#' @param W_prior An \eqn{n x x} matrix of priors for \eqn{W}
#' @param tt Number of time observations
#' @param niter Total number of iterations
#' @param nretain Number of iterations to retain (must be smaller than niter)
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
#' @param VERBOSE Should detailed diagnostic information be provided (default: FALSE)
#' @param SYMMETRIC Should the estimated \eqn{W} matrix be symmetric (default: TRUE)
#' @param GRIDDY_GIBBS Should griddy-Gibbs be used for \eqn{\rho} estimation?
#' Does not work if \code{ROW_STANDARDIZED = FALSE}. Main advantage is that less draws are required for \eqn{\rho}
#' @param ROW_STANDARDIZED Should the estimated \eqn{W} matrix be row-standardized (default: TRUE)
#'
#' @return List with posterior samples for \eqn{\rho}, \eqn{\beta}, \eqn{\sigma^2}, \eqn{w}, and direct, indirect, and total effects.
#' @export SAR_W
SAR_W <- function(Y, X, W_prior, tt, niter = 1000, nretain = 250,
                   beta_prior_mean = matrix(0, ncol(X) + 1, 1),
                   beta_prior_var = diag(ncol(X) + 1) * 10,
                   sigma_a = .1, sigma_b = .1,
                   rho_pr = 1.1, griddy_n = 105, rmin = 0, rmax = 1,
                   bbinom_a = 1, bbinom_b = 1, bb_pr = TRUE,
                   min_k = 1, max_k = min(10, floor(nrow(W_prior) / 2)), rjct_pr = FALSE,
                   VERBOSE = FALSE, SYMMETRIC = TRUE, GRIDDY_GIBBS = TRUE, ROW_STANDARDIZED = TRUE) {
  if (is.null(colnames(X))) {
    varnames <- paste0("X", 1:ncol(X))
  } else {
    varnames <- colnames(X)
  }
  varnames <- c("(Intercept)", varnames)
  if (!ROW_STANDARDIZED && GRIDDY_GIBBS) {
    stop("No support for Griddy Gibbs without row-standardization!")
  }

  ndiscard <- niter - nretain
  k <- ncol(X) + 1
  smallk <- k - 1
  n <- nrow(X) / tt

  diag(W_prior) <- 0

  rho_scale <- 1
  rho_accept <- 0


  # save the posterior draws here
  postb <- matrix(0, k, nretain)
  rownames(postb) <- varnames
  posts <- matrix(0, 1, nretain)
  postr <- matrix(0, 1, nretain)
  postw <- array(0, c(n, n, nretain))
  postwprob <- array(0, c(n, n, nretain))

  post.direct <- matrix(0, smallk, nretain)
  post.indirect <- matrix(0, smallk, nretain)
  post.total <- matrix(0, smallk, nretain)
  rownames(post.direct) <- rownames(post.indirect) <- rownames(post.total) <- varnames[-1]

  # pre-calculate some terms for faster draws
  beta_prior_var_inv <- solve(beta_prior_var)
  curr.W <- matrix(0, n, n) # not standardized W
  ### generate curr.W from the prior distribution
  if (SYMMETRIC) {
    ii_samples <- sample(2:n, n - 1, replace = F)
  } else {
    ii_samples <- sample(1:n, n, replace = F)
  }
  for (i in ii_samples) {
    if (SYMMETRIC) {
      jj_samples <- sample(c(1:(i - 1)), i - 1, replace = F)
    } else {
      jj_samples <- sample(1:n, n, replace = F)
    }
    for (j in jj_samples) {
      curr.Wpr <- W_prior[i, j]
      if (SYMMETRIC) {
        neighb1 <- sum((curr.W + t(curr.W))[i, ])
      } else {
        neighb1 <- sum(curr.W[i, ])
      }
      if (rjct_pr) {
        if (SYMMETRIC) {
          rjct_n <- max(neighb1, sum((curr.W + t(curr.W))[j, ]))
        } else {
          rjct_n <- neighb1
        }
        if (rjct_n < min_k) {
          curr.Wpr <- 1
        } else if (rjct_n == max_k) {
          curr.Wpr <- 0
        }
      }
      if (bb_pr) {
        bbprior1 <- bbinompdf(neighb1, n - 1, bbinom_a, bbinom_b) * curr.Wpr
        bbprior0 <- (1 - bbinompdf(neighb1, n - 1, bbinom_a, bbinom_b)) * (1 - curr.Wpr)
        bbprior_ <- bbprior1 / (bbprior1 + bbprior0)
      } else {
        bbprior_ <- curr.Wpr
      }
      prob.delta <- bbprior_ / (bbprior_ + (1 - bbprior_))
      if (prob.delta == 1) {
        curr.W[i, j] <- 1
      } else if (prob.delta != 0) {
        curr.W[i, j] <- stats::rbinom(1, 1, prob.delta)
      }
    }
  }
  # # set-up random W, so that all have min_k neighbours
  # nr_gen_neighbours = min((min_k + max_k)/2,ceiling(n/10))
  # iter = 1
  # sample_from = 1:n
  # while (any(rowSums(curr.W) < min_k) && iter < n*100) {
  #   nneigh = sample(sample_from,2)
  #   prop.W = curr.W
  #   prop.W[nneigh[1],nneigh[2]] = 1
  #   prop.W[nneigh[2],nneigh[1]] = 1
  #   if ( all(sum(prop.W[nneigh,]) < max_k) ) {
  #     curr.W = prop.W
  #   }
  #   if (any(rowSums(curr.W[sample_from,]) >= min_k) && length(sample_from) > 1) {
  #     sample_from = c(1:n)[-which(rowSums(curr.W)>=min_k)]
  #   }
  #   iter = iter + 1
  # }
  if (SYMMETRIC) {
    curr.W[upper.tri(curr.W, diag = tt)] <- 0
  }
  # curr.w - row-standardized
  curr.w <- matrix(0, n, n)
  if (ROW_STANDARDIZED) {
    if (SYMMETRIC) {
      curr.w <- (curr.W + t(curr.W)) / rowSums((curr.W + t(curr.W)))
    } else {
      curr.w <- curr.W / rowSums(curr.W)
    }
    curr.w[is.na(curr.w)] <- 0
  } else {
    if (SYMMETRIC) {
      curr.w <- curr.W + t(curr.W)
    } else {
      curr.w <- curr.W
    }
  }

  tX <- cbind(1, X)
  tY <- matrix(Y, n, tt)
  tAy <- matrix(0, n, tt)
  curr.beta <- solve(crossprod(tX)) %*% crossprod(tX, Y)
  curr.sigma <- as.double(crossprod(Y - tX %*% curr.beta)) / (tt * n - k)
  curr.gamma <- matrix(0, n, n)
  curr.rho <- 0
  curr.AI <- as.matrix(solve(Matrix::.sparseDiagonal(n) - curr.rho * curr.w))

  ### Gibbs sampling
  curr.A <- Matrix::.sparseDiagonal(n) - curr.rho * curr.w
  curr.logdet <- log(det(curr.A))
  for (iter in 1:niter) {
    if (VERBOSE) {
      cat(
        "iter:", iter,
        "W:", sum(curr.w > 0),
        "min_k", min(rowSums(curr.w > 0)),
        "sigma:", round(curr.sigma, 2),
        "rho:", curr.rho, "\n"
      )
    }
    curr.A <- Matrix::.sparseDiagonal(n) - curr.rho * curr.w
    Ay <- matrix(curr.A %*% tY, n * tt, 1)

    # draw beta
    # e1 = try({
    V <- solve(beta_prior_var_inv + 1 / curr.sigma * crossprod(tX))
    b <- V %*% (beta_prior_var_inv %*% beta_prior_mean + 1 / curr.sigma * crossprod(tX, Ay))
    # curr.beta = mvrnorm(1,b,V)
    curr.beta <- b + t(chol(V)) %*% stats::rnorm(k)
    curr.xb <- tX %*% curr.beta
    curr.txb <- matrix(curr.xb, n, tt)

    # draw sigma
    curr.ESS <- crossprod(Ay - curr.xb)
    curr.sigma <- 1 / stats::rgamma(1, sigma_a + (tt * n) / 2, sigma_b + as.double(curr.ESS) / 2)

    ## Griddy-Gibbs step for rho
    if (GRIDDY_GIBBS) {
      logdets <- lndetPaceBarry(curr.w, length.out = griddy_n, rmin = rmin, rmax = rmax)[-griddy_n, ]
      wY <- curr.w %*% tY
      # ess.grid1 = sapply(logdets[,2], function(x) sum(dnorm(as.matrix(tY - x*wY),curr.txb,sqrt(curr.sigma),log = tt))  )
      ess.grid <- sapply(logdets[, 2], function(x) -sum(((tY - x * wY) - curr.txb)^2) / (2 * curr.sigma))
      den <- tt * logdets[, 1] + ess.grid + log(beta_prob(logdets[, 2], rho_pr))
      log_cond_post_rho <- den
      log_cond_post_rho <- log_cond_post_rho - max(log_cond_post_rho)
      cond_post_rho <- exp(log_cond_post_rho)
      z <- cumsum(cond_post_rho) / sum(cond_post_rho)
      rnd <- stats::runif(1) #* sum(z)
      ind <- min(which(rnd <= z))
      if (is.integer(ind) && ind <= length(logdets[, 2])) {
        curr.rho <- logdets[ind, 2]
        curr.A <- Matrix::.sparseDiagonal(n) - curr.rho * curr.w
        curr.AI <- as.matrix(solve(curr.A))
        curr.logdet <- log(det(curr.A))
      }
    } else {
      # draw p(rho | .) using MH-step
      accept <- 0
      while (accept != 1) {
        prop.rho <- stats::rnorm(1, curr.rho, rho_scale)
        if (prop.rho < rmax && prop.rho > rmin) {
          accept <- 1
        }
      }
      prop.A <- Matrix::.sparseDiagonal(n) - prop.rho * curr.w
      prop.logdet <- suppressWarnings(log(det(prop.A)))

      if (iter == 1) {
        curr.logdet <- log(det(curr.A))
      }
      post_curr <- tt * curr.logdet +
        sum(stats::dnorm(as.matrix(curr.A %*% tY), curr.txb, sqrt(curr.sigma), log = T)) +
        log(beta_prob(curr.rho, rho_pr))
      post_prop <- tt * prop.logdet +
        sum(stats::dnorm(as.matrix(prop.A %*% tY), curr.txb, sqrt(curr.sigma), log = T)) +
        log(beta_prob(prop.rho, rho_pr))

      acc_prob <- post_prop - post_curr
      if (is.nan(acc_prob) == FALSE) {
        if ((acc_prob) > log(stats::runif(1, 0, 1))) {
          curr.rho <- prop.rho
          curr.A <- prop.A
          curr.AI <- as.matrix(solve(prop.A))
          curr.logdet <- prop.logdet
          rho_accept <- rho_accept + 1
        }
      }
      if (iter < (ndiscard / 2)) {
        # rho tuning
        if ((rho_accept / iter) > 0.3) {
          rho_scale <- 1.1 * rho_scale
        }
        if ((rho_accept / iter) < 0.1) {
          rho_scale <- 0.9 * rho_scale
        }
      }
    }

    # Gibbs step for W - element-wise
    # curr.txb = matrix(curr.xb,n,tt)
    if (SYMMETRIC) {
      ii_samples <- sample(2:n, n - 1, replace = F)
    } else {
      ii_samples <- sample(1:n, n, replace = F)
    }
    for (ii in ii_samples) {
      if (SYMMETRIC) {
        jj_samples <- sample(c(1:(ii - 1)), ii - 1, replace = F)
      } else {
        jj_samples <- sample(1:n, n, replace = F)
      }
      for (jj in jj_samples) {
        if (W_prior[ii, jj] == 0) {
          curr.W[ii, jj] <- 0
        } else if (W_prior[ii, jj] == 1) {
          curr.W[ii, jj] <- 1
        } else {
          if (SYMMETRIC) {
            ch_elmnt <- c(ii, jj)
          } else {
            ch_elmnt <- ii
          }
          W0 <- W1 <- curr.W
          was1 <- (curr.W[ii, jj] == 1)
          if (was1) {
            W0[ii, jj] <- 0
            if (SYMMETRIC) {
              WW0 <- (W0 + t(W0))
            } else {
              WW0 <- W0
            }
            w0 <- w1 <- curr.w
            if (ROW_STANDARDIZED) {
              w0[ch_elmnt, ] <- WW0[ch_elmnt, ] / rowSums(WW0[ch_elmnt, , drop = F])
            } else {
              w0[ch_elmnt, ] <- WW0[ch_elmnt, ]
            }
            w0[is.na(w0)] <- 0
            A0 <- diag(n) - curr.rho * w0
            A1 <- curr.A
            diff0 <- A0[ch_elmnt, , drop = F] - curr.A[ch_elmnt, , drop = F]
            res0 <- update_ldetAI(ch_elmnt, diff0, curr.AI, curr.logdet)
            logdet0 <- res0$logdet
            logdet1 <- curr.logdet
          } else {
            W1[ii, jj] <- 1
            if (SYMMETRIC) {
              WW1 <- (W1 + t(W1))
            } else {
              WW1 <- W1
            }
            w0 <- w1 <- curr.w
            if (ROW_STANDARDIZED) {
              w1[ch_elmnt, ] <- WW1[ch_elmnt, ] / rowSums(WW1[ch_elmnt, , drop = F])
            } else {
              w1[ch_elmnt, ] <- WW1[ch_elmnt, ]
            }
            w1[is.na(w1)] <- 0
            A1 <- diag(n) - curr.rho * w1
            A0 <- curr.A
            diff1 <- A1[ch_elmnt, , drop = F] - curr.A[ch_elmnt, , drop = F]
            logdet0 <- curr.logdet
            res1 <- update_ldetAI(ch_elmnt, diff1, curr.AI, curr.logdet)
            logdet1 <- res1$logdet
          }

          curr.W_prior <- W_prior
          # # rejection prior
          if (rjct_pr) {
            if (SYMMETRIC) {
              W_reject1 <- rowSums(w1[c(ii, jj), ] > 0) > max_k
            } else {
              W_reject1 <- sum(W1[ii, ]) > max_k
            }
            if (sum(W_reject1) > 0) {
              curr.W_prior[ii, jj] <- 0
            }
            if (SYMMETRIC) {
              W_reject0 <- rowSums(w0[c(ii, jj), ] > 0) < min_k
            } else {
              W_reject0 <- sum(W0[ii, ]) < min_k
            }
            if (sum(W_reject0) > 0) {
              curr.W_prior[ii, jj] <- 1
            }
          }
          if (bb_pr) {
            if (SYMMETRIC) {
              neighb0 <- sum((W0 + t(W0))[ii, ])
            } else {
              neighb0 <- sum(W0[ii, ])
            }
            # bbprior0 = bbinompdf(neighb0,n-1,bbinom_a,bbinom_b) * (1 - curr.W_prior[ii,jj])
            bbprior0 <- bbinompdf(neighb0, n - 1, bbinom_a, bbinom_b, min_k, max_k) * (1 - curr.W_prior[ii, jj])
            # if (neighb0 == 0) {bbprior0 = 0}
            if (SYMMETRIC) {
              neighb1 <- sum((W1 + t(W1))[ii, ])
            } else {
              neighb1 <- sum(W1[ii, ])
            }
            # bbprior1 = bbinompdf(neighb1,n-1,bbinom_a,bbinom_b) * curr.W_prior[ii,jj]
            bbprior1 <- bbinompdf(neighb1, n - 1, bbinom_a, bbinom_b, min_k, max_k) * curr.W_prior[ii, jj]
            bbprior_ <- bbprior1 / (bbprior1 + bbprior0)
          } else {
            bbprior_ <- curr.W_prior[ii, jj]
          }

          err1 <- sum((A1[ch_elmnt, ] %*% tY - curr.txb[ch_elmnt, ])^2)
          err0 <- sum((A0[ch_elmnt, ] %*% tY - curr.txb[ch_elmnt, ])^2)
          adj <- min(err0, err1)
          err1 <- err1 - adj
          err0 <- err0 - adj

          # change 23.02.2022
          #p1 =   bbprior_ * exp(logdet1*tt) * dnorm(err1,0,sqrt(curr.sigma))
          #p0 = (1- bbprior_) * exp(logdet0*tt) * dnorm(err0,0,sqrt(curr.sigma))
          p1 <- bbprior_ * exp(logdet1 * tt) * exp(-err1 / (curr.sigma))
          p0 <- (1 - bbprior_) * exp(logdet0 * tt) * exp(-err0 / (curr.sigma))

          prob.delta <- p1 / (p1 + p0)
          if (is.na(prob.delta)) {
            prob.delta <- 0
          }
          rnd_draw <- stats::runif(1)
          if (rnd_draw <= prob.delta) {
            curr.W[ii, jj] <- 1
            if (!was1) {
              curr.logdet <- logdet1
              curr.A <- A1
              curr.AI <- res1$AI
              curr.w <- w1
            }
          } else {
            curr.W[ii, jj] <- 0
            if (was1) {
              curr.logdet <- logdet0
              curr.A <- A0
              curr.AI <- res0$AI
              curr.w <- w0
            }
          }
        }
        curr.gamma[ii, jj] <- curr.W[ii, jj]
      }
    }
    if (ROW_STANDARDIZED) {
      if (SYMMETRIC) {
        curr.w <- (curr.W + t(curr.W)) / rowSums(curr.W + t(curr.W))
      } else {
        curr.w <- curr.W / rowSums(curr.W)
      }
      if (any(is.na(curr.w))) {
        curr.w[is.na(curr.w)] <- 0
      }
    } else {
      if (SYMMETRIC) {
        curr.w <- curr.W + t(curr.W)
      } else {
        curr.w <- curr.W
      }
    }

    # we are past the burn-in, save the draws
    if (iter > ndiscard) {
      s <- iter - ndiscard
      postb[, s] <- as.matrix(curr.beta)
      posts[s] <- curr.sigma
      postr[s] <- curr.rho
      postw[, , s] <- curr.w

      post.direct[, s] <- sum(diag(curr.AI)) / n * curr.beta[-1]
      post.total[, s] <- sum(curr.AI) / n * curr.beta[-1]
      post.indirect[, s] <- post.total[, s] - post.direct[, s]
    }
  }

  return(list(
    postb = postb, posts = posts, postr = postr, postw = postw,
    post.direct = post.direct, post.indirect = post.indirect, post.total = post.total
  ))
}

# ## first let's construct our SAR DGP
# n=20
# tt =30
# smallk = 2
# k = 1 + smallk
# # # 7 nearest neighbour W construction from random pattern
# xy <- cbind(runif(n),runif(n))
# #xy <- cbind(1:n,1:n)
# source("knn.R")
# W<-getWknn(xy,3)
# #W = t(W) + W
# #W = as.matrix(W/rowSums(W))
# #W[W>0] = 1
# diag(W) = 0
# # one-forward, one-behind W
# # W = matrix(0,n,n)
# # W[-1,-n] = W[-1,-n] + diag(n-1)
# # W[-n,-1] = W[-n,-1] + diag(n-1)
# # W = W + W%*%W + W%*%W%*%W
# # W[W>0] = 1; diag(W) = 0
# # W = W/rowSums(W)
# #
# #
# # # test1 = matrix(0,30,20)
# # # rrhos = seq(-.99,.99,length.out = 30)
# # # for (rr in 1:30) {
# # #   for (kk in 1:20) {
# # #     W<-getWknn(xy,ceiling(kk))
# # #     W = t(W) + W
# # #     W = as.matrix(W/rowSums(W))
# # #     diag(W) = 0
# # #     test1[rr,kk] = max(solve(.sparseDiagonal(n)-rrhos[rr]*W))
# # #   }
# # # }
# # # plot(rrhos,test1[,1],type="l")
# # # for (j in 2:20) {
# # #   lines(rrhos,test1[,j])
# # # }
# #
# RHO = .5
# A <- as.matrix(diag(n) - RHO*W); Ainv <- solve(A)
# X <- matrix(rnorm(n*smallk * tt),n*tt,smallk)
# ALPHA <- 0
# BETA = c(1,-1)
# SIGMA = 1
# BETA_TRUE <- c(ALPHA, BETA)
# Y = kronecker(.sparseDiagonal(tt),Ainv) %*% (ALPHA + X %*% BETA + rnorm(n*tt, mean = 0,sd = SIGMA))
#
# INDIRECT = (sum(Ainv) - sum(diag(Ainv))) * BETA / n
#
# W_prior = matrix(.5,n,n);
# res1 = SAR_W1(Y,X,W_prior,tt,niter = 100,nretain = 50,VERBOSE = TRUE,rho_pr = 1.01,
#               SYMMETRIC = TRUE,bb_pr = FALSE,rjct_pr = TRUE,min_k = 1,max_k = 4,sigma_a = 1,sigma_b = 1,
#               rmin = 0,rmax = 1,
#               ROW_STANDARDIZED = TRUE,GRIDDY_GIBBS = TRUE)
#
# res2 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = FALSE,rjct_pr = TRUE,min_k = 5,max_k = 25,sigma_a = 1,sigma_b = 1)
#
# W_prior = matrix(.1,n,n);
# res3 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = FALSE,rjct_pr = TRUE,min_k = 1,max_k = n/2,sigma_a = 1,sigma_b = 1)
#
# W_prior = matrix(.5,n,n);
# mbar = n/10; bbinom_b = ((n-1) - mbar) / mbar
# res4 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,bbinom_b = 1,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = TRUE,min_k = 1,max_k = n/2,sigma_a = 1,sigma_b = 1)
#
# res5 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,bbinom_b = 1,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = TRUE,min_k = 5,max_k = 25,sigma_a = 1,sigma_b = 1)
#
# mbar = n/10; bbinom_b = ((n-1) - mbar) / mbar
# res6 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,bbinom_b = mbar,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = TRUE,min_k = 1,max_k = n/2,sigma_a = 1,sigma_b = 1)
#
# res7 = SAR_W1(Y,X,W_prior,tt,niter = 50,nretain = 30,VERBOSE = TRUE,bbinom_b = mbar,rho_pr = 7.5,
#               SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = TRUE,min_k = 5,max_k = 25,sigma_a = 1,sigma_b = 1)
#
#
# sqrt(mean((res1$post.indirect - INDIRECT)^2))
# sqrt(mean((res2$post.indirect - INDIRECT)^2))
# sqrt(mean((res3$post.indirect - INDIRECT)^2))
# sqrt(mean((res4$post.indirect - INDIRECT)^2))
# sqrt(mean((res5$post.indirect - INDIRECT)^2))
# sqrt(mean((res6$post.indirect - INDIRECT)^2))
# sqrt(mean((res7$post.indirect - INDIRECT)^2))

# # res3 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = FALSE,bb_pr = TRUE,bbinom_b = 5,bbinom_a = .5)
# # res4 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = FALSE,bb_pr = TRUE,rjct_pr = TRUE,max_k = 10)
# #
# # res5 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = TRUE,bb_pr = FALSE,rjct_pr = FALSE)
# # res6 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = TRUE,bb_pr = FALSE,rjct_pr = TRUE,max_k = 10)
# # res7 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = FALSE)
# # res8 = SAR_W1(Y,X,W_prior,tt,niter = 30,nretain = 10,VERBOSE = T,SYMMETRIC = TRUE,bb_pr = TRUE,rjct_pr = TRUE,max_k = 10)
# # #
# # plots
# require(plot.matrix)
# plot(as.matrix(W))
# # par(mfrow=c(2,2))
# # plot(apply(res1$postw>0,c(1,2),mean),main = "non-symmetric, no bbinom, no reject")
# # plot(apply(res2$postw>0,c(1,2),mean),main = "non-symmetric, no bbinom, reject")
# plot(apply(res3$postw>0,c(1,2),mean),main = "non-symmetric, bbinom, no reject")
# # plot(apply(res4$postw>0,c(1,2),mean),main = "non-symmetric, bbinom, reject")
# #
# # plot(apply(res5$postw>0,c(1,2),mean),main = "symmetric, no bbinom, no reject")
# # plot(apply(res6$postw>0,c(1,2),mean),main = "symmetric, no bbinom, reject")
# # plot(apply(res7$postw>0,c(1,2),mean),main = "symmetric, bbinom, no reject")
# # plot(apply(res8$postw>0,c(1,2),mean),main = "symmetric, bbinom, reject")
# # par(mfrow=c(1,1))
