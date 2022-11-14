#' An R6 class for sampling slope parameters
#'
#' This class samples slope parameters with a Gaussian prior from the conditional posterior.
#' Use the \link{beta_priors} class for setup.
#'
#' @field beta_prior The current \code{\link{beta_priors}}
#' @field curr_beta The current value of \eqn{\beta}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
beta_sampler = R6::R6Class("beta_sampler",cloneable = FALSE, public = list(
  beta_prior = NULL,
  curr_beta = NULL,
  #' @param beta_prior The list returned by \code{\link{beta_priors}}
  #'
  #' @export
  initialize = function(beta_prior) {
    self$beta_prior = beta_prior
    k = ncol(self$beta_prior$beta_var_prior_inv)
    self$sample(matrix(0,2,1),matrix(0,2,k),1)
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{1} matrix of responses
  #' @param X The \eqn{n} by \eqn{k} design matrix
  #' @param curr_sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample = function(Y,X,curr_sigma) {
    k = ncol(self$beta_prior$beta_var_prior_inv)
    sig_sqr = sqrt(curr_sigma)
    V <- solve(self$beta_prior$beta_var_prior_inv + crossprod(X / sig_sqr))
    b <- V %*% (self$beta_prior$beta_var_prior_inv %*% self$beta_prior$beta_mean_prior +
                  crossprod(X / sig_sqr, Y / sig_sqr))
    self$curr_beta <- b + t(chol(V)) %*% stats::rnorm(k)
    invisible(self)
  }
))

#' An R6 class for sampling for sampling \eqn{\sigma^2}
#'
#' This class samples nuisance parameter which an inverse Gamma prior from the conditional posterior.
#' Use the \link{sigma_priors} class for setup.
#'
#' @field sigma_prior The current \code{\link{sigma_priors}}
#' @field curr_sigma The current value of \eqn{\beta}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
sigma_sampler = R6::R6Class("sigma_sampler", cloneable = FALSE, public = list(
  sigma_prior = NULL,
  curr_sigma = NULL,
  #' @param sigma_prior The list returned by \code{\link{sigma_priors}}
  #'
  #' @export
  initialize = function(sigma_prior) {
    self$sigma_prior = sigma_prior
    self$curr_sigma = 1 / stats::rgamma(1, self$sigma_prior$sigma_rate_prior ,
                      self$sigma_prior$sigma_shape_prior  )
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{1} matrix of responses
  #' @param mu The \eqn{n} by \eqn{1} matrix of means
  #'
  #' @export
  sample = function(Y,mu) {
    n = nrow(Y)
    curr.err = as.double(crossprod(Y - mu))
    self$curr_sigma <- 1 / stats::rgamma(1, self$sigma_prior$sigma_rate_prior + (n) / 2,
                                         self$sigma_prior$sigma_shape_prior + curr.err / 2)
    invisible(self)
  }
))

#' An R6 class for sampling the spatial autoregressive parameter \eqn{\rho}
#'
#' This class samples the spatial autoregressive parameter using either a tuned random-walk
#' Metropolis-Hastings or a griddy Gibbs step. Use the \code{\link{rho_priors}} class for setup.
#'
#' For the Griddy-Gibbs algorithm see Ritter and Tanner (1992).
#'
#' @field rho_prior The current \code{\link{rho_priors}}
#' @field curr_rho The current value of \eqn{\rho}
#' @field curr_W The current value of the spatial weight matrix \eqn{W}; an \eqn{n} by \eqn{n} matrix.
#' @field curr_A The current spatial projection matrix \eqn{I - \rho W}.
#' @field curr_AI The inverse of \code{curr_A}
#' @field curr_logdet The current log-determinant of \code{curr_A}
#' @field curr_logdets A set of log-determinants for various values of \eqn{\rho}. See the
#'  \code{\link{rho_priors}} function for settings of step site and other parameters of the grid.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#'
#' @references
#'  Ritter, C., and Tanner, M. A. (1992). Facilitating the Gibbs sampler: The Gibbs stopper
#'  and the griddy-Gibbs sampler. \emph{Journal of the American Statistical Association},
#'  \bold{87(419)}, 861-868.
rho_sampler = R6::R6Class("rho_sampler", cloneable = FALSE, public = list(
  rho_prior = NULL,
  curr_rho = NULL,
  curr_W = NULL,
  curr_A = NULL,
  curr_AI = NULL,
  curr_logdet = NULL,
  curr_logdets = NULL,
  #' @param rho_prior The list returned by \code{\link{rho_priors}}
  #' @param W An optional starting value for the spatial weight matrix \eqn{W}
  #'
  #' @export
  initialize = function(rho_prior,W = NULL) {
    self$rho_prior = rho_prior
    self$curr_W = W
    self$curr_rho = stats::runif(1,self$rho_prior$rho_min,self$rho_prior$rho_max)
    if (!self$rho_prior$use_griddy_gibbs) {
      private$curr_rho_scale <- self$rho_prior$init_rho_scale
      private$rho_accept = 0
      private$curr_iter = 0
      private$do_MHtune = TRUE
    }
    if (!is.null(self$curr_W)) {
      self$setW(curr_W)
    }
    invisible(self)
  },
  #' @description
  #' Function to stop the tuning of the Metropolis-Hastings step. The tuning of the
  #' Metropolis-Hastings step is usually carried out until half of the burn-in phase.
  #' Call this function to turn it off.
  #'
  #' @export
  stopMHtune = function() {
    if (!self$rho_prior$use_griddy_gibbs) {
      private$do_MHtune = FALSE
    }
  },
  #' @param newW The updated spatial weight matrix \eqn{W}
  #' @param newLogdet An optional value for the log determinant corresponding to \code{newW} and \code{curr_rho}
  #' @param newA An optional value for the spatial projection matrix using \code{newW} and \code{curr_rho}
  #' @param newAI An optional value for the matrix inverse of \code{newA}
  #'
  #' @export
  setW = function(newW, newLogdet = NULL, newA = NULL, newAI = NULL) {
    self$curr_W = newW
    n = nrow(self$curr_W)
    if (is.null(newA)) {
      self$curr_A <- Matrix::.sparseDiagonal(n) - self$curr_rho * self$curr_W
    } else {self$curr_A = newA}
    if (is.null(newAI)) {
      self$curr_AI <- as.matrix(solve(self$curr_A))
    } else {self$curr_AI = newAI}
    if (is.null(newLogdet)) {
      self$curr_logdet <- log(Matrix::det(self$curr_A))
    } else {self$curr_logdet = newLogdet}
    if (self$rho_prior$use_griddy_gibbs) {
      self$curr_logdets <- logdetPaceBarry(self$curr_W, length.out = self$rho_prior$griddy_n,
                                        rmin = self$rho_prior$rho_min,
                                        rmax = self$rho_prior$rho_max)[-self$rho_prior$griddy_n, ]
    }
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} matrix of responses
  #' @param mu The \eqn{n} by \eqn{tt} matrix of means
  #' @param sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample = function(Y,mu, sigma) {
    if (self$rho_prior$use_griddy_gibbs) {
      self$sample_Griddy(Y,mu,sigma)
    } else {
      self$sample_MH(Y,mu,sigma)
    }
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} matrix of responses
  #' @param mu The \eqn{n} by \eqn{tt} matrix of means
  #' @param sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample_Griddy = function(Y,mu, sigma) {
    n = nrow(Y); tt = ncol(Y)
    wY <- self$curr_W %*% Y
    ess.grid <- sapply(self$curr_logdets[, 2],
                       function(x) -sum(((Y - x * wY) - mu)^2) / (2 * sigma))
    den <- tt * self$curr_logdets[, 1] + ess.grid +
      log(betapdf(self$curr_logdets[, 2],
                  self$rho_prior$rho_a_prior, self$rho_prior$rho_b_prior,
                  self$rho_prior$rho_min, self$rho_prior$rho_max))
    log_cond_post_rho <- den
    log_cond_post_rho <- log_cond_post_rho - max(log_cond_post_rho)
    cond_post_rho <- exp(log_cond_post_rho)
    z <- cumsum(cond_post_rho) / sum(cond_post_rho)
    rnd <- stats::runif(1) #* sum(z)
    ind <- min(which(rnd <= z))
    if (is.integer(ind) && ind <= length(self$curr_logdets[, 2])) {
      self$curr_rho <- self$curr_logdets[ind, 2]
      self$curr_A <- Matrix::.sparseDiagonal(n) - self$curr_rho * self$curr_W
      self$curr_AI <- as.matrix(solve(self$curr_A))
      self$curr_logdet <- log(Matrix::det(self$curr_A))
    }
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} matrix of responses
  #' @param mu The \eqn{n} by \eqn{tt} matrix of means
  #' @param sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample_MH = function(Y,mu, sigma) {
    n = nrow(Y); tt = ncol(Y)
    # Proposal for rho
    prop_rho <- stats::rnorm(1, self$curr_rho, private$curr_rho_scale)
    if (prop_rho < self$rho_prior$rho_max && prop_rho > self$rho_prior$rho_min) {
      prop_A <- Matrix::.sparseDiagonal(n) - prop_rho * self$curr_W
      prop_logdet <- suppressWarnings(log(Matrix::det(prop_A)))

      post_curr <- tt * self$curr_logdet +
        sum(stats::dnorm(as.matrix(self$curr_A %*% Y), mu, sqrt(sigma), log = T)) +
        log(betapdf(self$curr_rho,
                    self$rho_prior$rho_a_prior, self$rho_prior$rho_b_prior,
                    self$rho_prior$rho_min, self$rho_prior$rho_max))
      post_prop <- tt * prop_logdet +
        sum(stats::dnorm(as.matrix(prop_A %*% Y), mu, sqrt(sigma), log = T)) +
        log(betapdf(prop_rho,
                    self$rho_prior$rho_a_prior, self$rho_prior$rho_b_prior,
                    self$rho_prior$rho_min, self$rho_prior$rho_max))
      acc_prob <- post_prop - post_curr
      if (is.nan(acc_prob) == FALSE) {
        if ((acc_prob) > log(stats::runif(1, 0, 1))) {
          self$curr_rho <- prop_rho
          self$curr_A <- prop_A
          self$curr_AI <- as.matrix(solve(prop_A))
          self$curr_logdet <- prop_logdet
          private$rho_accept <- private$rho_accept + 1
        }
      }
    }
    if (private$do_MHtune) {
      private$curr_iter = private$curr_iter + 1
      # rho tuning
      if ((private$rho_accept / private$curr_iter) > self$rho_prior$mh_tune_high) {
        private$curr_rho_scale <- (1 + self$rho_prior$mh_tune_scale) * private$curr_rho_scale
      }
      if ((private$rho_accept / private$curr_iter) < self$rho_prior$mh_tune_low) {
        private$curr_rho_scale <- (1- self$rho_prior$mh_tune_scale) * private$curr_rho_scale
      }
    }
    invisible(self)
  }
),private = list(
  curr_rho_scale = NULL,
  rho_accept = NULL,
  curr_iter = NULL,
  do_MHtune = NULL
))

#' An R6 class for sampling the elements of \eqn{W}
#'
#' This class samples the spatial weight matrix. Use the function \link{W_priors} class for setup.
#'
#' The sampling procedure relies on conditional Bernoulli posteriors outlined in
#' Krisztin and Piribauer (2022).
#'
#' @field W_prior The current \code{\link{W_priors}}
#' @field curr_w numeric, non-negative \eqn{n} by \eqn{n} spatial weight matrix with zeros
#' on the main diagonal. Depending on the \code{\link{W_priors}} settings can be symmetric and/or
#' row-standardized.
#' @field curr_W binary \eqn{n} by \eqn{n} spatial connectivity matrix \eqn{\Omega}
#' @field curr_A The current spatial projection matrix \eqn{I - \rho W}.
#' @field curr_AI The inverse of \code{curr_A}
#' @field curr_logdet The current log-determinant of \code{curr_A}
#' @field curr_rho single number between -1 and 1 or NULL, depending on whether the sampler updates
#'         the spatial autoregressive parameter \eqn{\rho}. Set while invoking \code{initialize}
#'         or using the function \code{set_rho}.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#'
#' @references
#'  Krisztin, T., and Piribauer, P. (2022) A Bayesian approach for the estimation
#'  of weight matrices in spatial autoregressive models. \emph{Spatial Economic Analysis}, 1-20.
W_sampler = R6::R6Class("W_sampler", cloneable = FALSE, public =list(
  W_prior = NULL,
  curr_w = NULL,
  curr_W = NULL,
  curr_A = NULL,
  curr_AI = NULL,
  curr_logdet = NULL,
  curr_rho = NULL,
  #' @param W_prior The list returned by \code{\link{W_priors}}
  #' @param curr_rho optional single number between -1 and 1. Value of the spatial autoregressive
  #'  parameter \eqn{\rho}. Defaults to NULL, in which case no updates of the log-determinant, the spatial
  #'   projection matrix, and its inverse are carried out.
  #'
  #' @export
  initialize = function(W_prior,curr_rho = NULL) {
    self$W_prior = W_prior
    n = nrow(self$W_prior$W_prior)
    self$curr_W <- matrix(0, n, n) # not standardized W
    ### generate curr_W from the prior distribution
    if (self$W_prior$symmetric_prior) {
      ii_samples <- sample(2:n, n - 1, replace = F)
    } else {
      ii_samples <- sample(1:n, n, replace = F)
    }
    for (i in ii_samples) {
      if (self$W_prior$symmetric_prior) {
        jj_samples <- sample(c(1:(i - 1)), i - 1, replace = F)
      } else {
        jj_samples <- sample(1:n, n, replace = F)
      }
      for (j in jj_samples) {
        curr.Wpr <- self$W_prior$W_prior[i, j]
        if (self$W_prior$symmetric_prior) {
          neighb1 <- sum((self$curr_W + t(self$curr_W))[i, ])
        } else {
          neighb1 <- sum(self$curr_W[i, ])
        }
        if (self$W_prior$use_reject_prior) {
          if (self$W_prior$symmetric_prior) {
            rjct_n <- max(neighb1, sum((self$curr_W + t(self$curr_W))[j, ]))
          } else {
            rjct_n <- neighb1
          }
          if (rjct_n < self$W_prior$min_neighbors) {
            curr.Wpr <- 1
          } else if (rjct_n == self$W_prior$max_neighbors) {
            curr.Wpr <- 0
          }
        }
        if (self$W_prior$use_bbinom_prior) {
          bbprior1 <- bbinompdf(neighb1, n - 1, self$W_prior$bbinom_a_prior,
                                self$W_prior$bbinom_b_prior) * curr.Wpr
          bbprior0 <- (1 - bbinompdf(neighb1, n - 1, self$W_prior$bbinom_a_prior,
                                self$W_prior$bbinom_b_prior)) * (1 - curr.Wpr)
          bbprior_ <- bbprior1 / (bbprior1 + bbprior0)
        } else {
          bbprior_ <- curr.Wpr
        }
        prob.delta <- bbprior_ / (bbprior_ + (1 - bbprior_))
        if (prob.delta == 1) {
          self$curr_W[i, j] <- 1
        } else if (prob.delta != 0) {
          self$curr_W[i, j] <- stats::rbinom(1, 1, prob.delta)
        }
      }
    }
    if (self$W_prior$symmetric_prior) {
      self$curr_W[upper.tri(self$curr_W, diag = TRUE)] <- 0
    }
    self$curr_w <- matrix(0, n, n)
    if (self$W_prior$row_standardized_prior) {
      if (self$W_prior$symmetric_prior) {
        self$curr_w <- (self$curr_W + t(self$curr_W)) / rowSums((self$curr_W + t(self$curr_W)))
      } else {
        self$curr_w <- self$curr_W / rowSums(self$curr_W)
      }
      self$curr_w[is.na(self$curr_w)] <- 0
    } else {
      if (self$W_prior$symmetric_prior) {
        self$curr_w <- self$curr_W + t(self$curr_W)
      } else {
        self$curr_w <- self$curr_W
      }
    }
    if (!is.null(curr_rho)) {
      self$set_rho(curr_rho)
    }
  },
  #' @description
  #' If the spatial autoregressive parameter \eqn{\rho} is updated during the sampling procedure the log determinant, the
  #' spatial projection matrix \eqn{I - \rho W} and it's inverse must be updated. This function should be
  #' used for a consistent update. At least the new scalar value for \eqn{\rho} must be supplied.
  #'
  #' @param new_rho single, number; must be between -1 and 1.
  #' @param newLogdet An optional value for the log determinant corresponding to \code{newW} and \code{curr_rho}
  #' @param newA An optional value for the spatial projection matrix using \code{newW} and \code{curr_rho}
  #' @param newAI An optional value for the matrix inverse of \code{newA}
  #'
  #' @export
  set_rho = function(new_rho, newLogdet = NULL, newA = NULL, newAI = NULL) {
    self$curr_rho = new_rho
    n = nrow(self$W_prior$W_prior)
    if (is.null(newA)) {
      self$curr_A <- Matrix::.sparseDiagonal(n) - self$curr_rho * self$curr_w
    } else {self$curr_A = newA}
    if (is.null(newAI)) {
      self$curr_AI <- as.matrix(solve(self$curr_A))
    } else {self$curr_AI = newAI}
    if (is.null(newLogdet)) {
      self$curr_logdet <- log(Matrix::det(self$curr_A))
    } else {self$curr_logdet = newLogdet}
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} matrix of responses
  #' @param curr_sigma The variance parameter \eqn{\sigma^2}
  #' @param mu The \eqn{n} by \eqn{tt} matrix of means.
  #' @param lag_mu \eqn{n} by \eqn{tt} matrix of means that will be spatially lagged with
  #' the estimated \eqn{W}. Defaults to a matrix with zero elements.
  #'
  #' @export
  sample = function(Y,curr_sigma,
                    mu,
                    lag_mu  = matrix(0,nrow(tY),ncol(tY))) {
    n = nrow(Y)
    tt = ncol(Y)
    if (self$W_prior$symmetric_prior) {
      ii_samples <- sample(2:n, n - 1, replace = F)
    } else {
      ii_samples <- sample(1:n, n, replace = F)
    }
    for (ii in ii_samples) {
      if (self$W_prior$symmetric_prior) {
        jj_samples <- sample(c(1:(ii - 1)), ii - 1, replace = F)
      } else {
        jj_samples <- sample(1:n, n, replace = F)
      }
      for (jj in jj_samples) {
        if (self$W_prior$W_prior[ii, jj] == 0) {
          self$curr_W[ii, jj] <- 0
        } else if (self$W_prior$W_prior[ii, jj] == 1) {
          self$curr_W[ii, jj] <- 1
        } else {
          if (self$W_prior$symmetric_prior) {
            ch_elmnt <- c(ii, jj)
          } else {
            ch_elmnt <- ii
          }
          W0 <- W1 <- self$curr_W
          was1 <- (self$curr_W[ii, jj] == 1)
          if (was1) {
            W0[ii, jj] <- 0
            if (self$W_prior$symmetric_prior) {
              WW0 <- (W0 + t(W0))
            } else {
              WW0 <- W0
            }
            w0 <- w1 <- self$curr_w
            if (self$W_prior$row_standardized_prior) {
              w0[ch_elmnt, ] <- WW0[ch_elmnt, ] / rowSums(WW0[ch_elmnt, , drop = F])
            } else {
              w0[ch_elmnt, ] <- WW0[ch_elmnt, ]
            }
            w0[is.na(w0)] <- 0
            if (!is.null(self$curr_rho)) {
              A0 <- diag(n) - self$curr_rho * w0
              A1 <- self$curr_A
              diff0 <- A0[ch_elmnt, , drop = F] - self$curr_A[ch_elmnt, , drop = F]
              res0 <- logdetAinvUpdate(ch_elmnt, diff0, self$curr_AI, self$curr_logdet)
              logdet0 <- res0$logdet
              logdet1 <- self$curr_logdet
            }
          } else {
            W1[ii, jj] <- 1
            if (self$W_prior$symmetric_prior) {
              WW1 <- (W1 + t(W1))
            } else {
              WW1 <- W1
            }
            w0 <- w1 <- self$curr_w
            if (self$W_prior$row_standardized_prior) {
              w1[ch_elmnt, ] <- WW1[ch_elmnt, ] / rowSums(WW1[ch_elmnt, , drop = F])
            } else {
              w1[ch_elmnt, ] <- WW1[ch_elmnt, ]
            }
            w1[is.na(w1)] <- 0
            if (!is.null(self$curr_rho)) {
              A1 <- diag(n) - self$curr_rho * w1
              A0 <- self$curr_A
              diff1 <- A1[ch_elmnt, , drop = F] - self$curr_A[ch_elmnt, , drop = F]
              logdet0 <- self$curr_logdet
              res1 <- logdetAinvUpdate(ch_elmnt, diff1, self$curr_AI, self$curr_logdet)
              logdet1 <- res1$logdet
            }
          }

          curr.W_prior <- self$W_prior$W_prior
          # # rejection prior
          if (self$W_prior$use_reject_prior) {
            if (self$W_prior$symmetric_prior) {
              W_reject1 <- rowSums(w1[c(ii, jj), ] > 0) > self$W_prior$max_neighbors
            } else {
              W_reject1 <- sum(W1[ii, ]) > self$W_prior$max_neighbors
            }
            if (sum(W_reject1) > 0) {
              curr.W_prior[ii, jj] <- 0
            }
            if (self$W_prior$symmetric_prior) {
              W_reject0 <- rowSums(w0[c(ii, jj), ] > 0) < self$W_prior$min_neighbors
            } else {
              W_reject0 <- sum(W0[ii, ]) < self$W_prior$min_neighbors
            }
            if (sum(W_reject0) > 0) {
              curr.W_prior[ii, jj] <- 1
            }
          }
          if (self$W_prior$use_bbinom_prior) {
            if (self$W_prior$symmetric_prior) {
              neighb0 <- sum((W0 + t(W0))[ii, ])
            } else {
              neighb0 <- sum(W0[ii, ])
            }
            # bbprior0 = bbinompdf(neighb0,n-1,bbinom_a_prior,bbinom_b_prior) * (1 - curr.W_prior[ii,jj])
            bbprior0 <- bbinompdf(neighb0, n - 1,self$W_prior$bbinom_a_prior,
                                  self$W_prior$bbinom_b_prior, self$W_prior$min_neighbors,
                                  self$W_prior$max_neighbors) *
              (1 - curr.W_prior[ii, jj])
            # if (neighb0 == 0) {bbprior0 = 0}
            if (self$W_prior$symmetric_prior) {
              neighb1 <- sum((W1 + t(W1))[ii, ])
            } else {
              neighb1 <- sum(W1[ii, ])
            }
            # bbprior1 = bbinompdf(neighb1,n-1,bbinom_a_prior,bbinom_b_prior) * curr.W_prior[ii,jj]
            bbprior1 <- bbinompdf(neighb1, n - 1, self$W_prior$bbinom_a_prior,
                                  self$W_prior$bbinom_b_prior, self$W_prior$min_neighbors,
                                  self$W_prior$max_neighbors) *
              curr.W_prior[ii, jj]
            bbprior_ <- bbprior1 / (bbprior1 + bbprior0)
          } else {
            bbprior_ <- curr.W_prior[ii, jj]
          }

          if (!is.null(self$curr_rho)) {
            err1 <- sum((A1[ch_elmnt, ] %*% Y - mu[ch_elmnt, ] - w1[ch_elmnt,] %*% lag_mu)^2)
            err0 <- sum((A0[ch_elmnt, ] %*% Y - mu[ch_elmnt, ] - w0[ch_elmnt,] %*% lag_mu)^2)
          } else {
            err1 = sum( (Y[ch_elmnt,] - w1[ch_elmnt,] %*% lag_mu - mu[ch_elmnt,])^2 )
            err0 = sum( (Y[ch_elmnt,] - w0[ch_elmnt,] %*% lag_mu - mu[ch_elmnt,])^2 )
          }
          adj <- min(err0, err1)
          err1 <- err1 - adj
          err0 <- err0 - adj

          if (!is.null(self$curr_rho)) {
            # change 23.02.2022
            #p1 =   bbprior_ * exp(logdet1*tt) * dnorm(err1,0,sqrt(curr.sigma))
            #p0 = (1- bbprior_) * exp(logdet0*tt) * dnorm(err0,0,sqrt(curr.sigma))
            p1 <- bbprior_ * exp(logdet1 * tt) * exp(-err1 / (curr_sigma))
            p0 <- (1 - bbprior_) * exp(logdet0 * tt) * exp(-err0 / (curr_sigma))
          } else {
            p1 =   bbprior_ * stats::dnorm(err1,0,sqrt(curr_sigma))
            p0 = (1- bbprior_) * stats::dnorm(err0,0,sqrt(curr_sigma))
          }

          prob.delta <- p1 / (p1 + p0)
          if (is.na(prob.delta)) {
            prob.delta <- 0
          }
          rnd_draw <- stats::runif(1)
          if (rnd_draw <= prob.delta) {
            self$curr_W[ii, jj] <- 1
            if (!was1) {
              self$curr_w <- w1
              if (!is.null(self$curr_rho)) {
                self$curr_logdet <- logdet1
                self$curr_A <- A1
                self$curr_AI <- res1$AI
              }
            }
          } else {
            self$curr_W[ii, jj] <- 0
            if (was1) {
              self$curr_w <- w0
              if (!is.null(self$curr_rho)) {
                self$curr_logdet <- logdet0
                self$curr_A <- A0
                self$curr_AI <- res0$AI
              }
            }
          }
        }
      }
    }
    if (self$W_prior$row_standardized_prior) {
      if (self$W_prior$symmetric_prior) {
        self$curr_w <- (self$curr_W + t(self$curr_W)) / rowSums(self$curr_W + t(self$curr_W))
      } else {
        self$curr_w <- self$curr_W / rowSums(self$curr_W)
      }
      if (any(is.na(self$curr_w))) {
        self$curr_w[is.na(self$curr_w)] <- 0
      }
    } else {
      if (self$W_prior$symmetric_prior) {
        self$curr_w <- self$curr_W + t(self$curr_W)
      } else {
        self$curr_w <- self$curr_W
      }
    }
  }
))
