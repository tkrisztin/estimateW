#' An R6 class for sampling for sampling slope coefficients
#'
#' This class samples slope coefficients with a Gaussian prior. Use the \link{beta_priors} class for setup.
#'
#' @field beta_prior The current \code{\link{beta_priors}}
#' @field curr_beta The current value of \eqn{\beta}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
beta_sampler = R6::R6Class("beta_sampler", list(
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
  #' @param Y The \eqn{n} by \eqn{1} vector of responses
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
#' This class samples slope coefficients with an inverse Gamma. Use the \link{sigma_priors} class for setup.
#'
#' @field sigma_prior The current \code{\link{sigma_priors}}
#' @field curr_sigma The current value of \eqn{\beta}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
sigma_sampler = R6::R6Class("sigma_sampler", list(
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
  #' @param Y The \eqn{n} by \eqn{1} vector of responses
  #' @param mu The \eqn{n} by \eqn{1} vector of means (classicaly \eqn{X\beta})
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

#' An R6 class for sampling \eqn{\rho} using a Griddy-Gibbs step
#'
#' This class samples the spatial autoregressive coefficient. Use the \link{rho_priors} class for setup.
#'
#' The sampling relies on the Griddy-Gibbs algorithm outlines in Ritter and Tanner (1992).
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
rho_sampler = R6::R6Class("rho_sampler", list(
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
    if (!is.null(self$curr_W)) {
      self$setW(curr_W)
    }
    invisible(self)
  },
  #' @param newW The updated spatial weight matrix \eqn{W}
  #' @param newLogdet An optional value for the log determinant corresponding to \çode{newW} and \code{curr_rho}
  #' @param newA An optional value for the spatial projection matrix using \çode{newW} and \code{curr_rho}
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
      self$curr_logdets <- lndetPaceBarry(self$curr_W, length.out = self$rho_prior$griddy_n,
                                     rmin = self$rho_prior$rho_min,
                                     rmax = self$rho_prior$rho_max)[-self$rho_prior$griddy_n, ]
    if (is.null(newLogdet)) {
      self$curr_logdet <- log(Matrix::det(self$curr_A))
    } else {self$logdet = newLogdet}
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} vector of responses
  #' @param mu The \eqn{n} by \eqn{tt} vector of means (usually \eqn{X\beta})
  #' @param sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample = function(Y,mu, sigma) {
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
  }
))


#' An R6 class for sampling \eqn{\rho} using a Metropolis-Hastings step
#'
#' This class samples the spatial autoregressive coefficient. Use the \link{rho_priors} class for setup.
#'
#' The sampling relies on a Metropolis-Hastings step, as outlined in LeSage and Pace (2009).
#'
#' @field rho_prior The current \code{\link{rho_priors}}
#' @field curr_rho The current value of \eqn{\rho}
#' @field curr_W The current value of the spatial weight matrix \eqn{W}; an \eqn{n} by \eqn{n} matrix.
#' @field curr_A The current spatial projection matrix \eqn{I - \rho W}.
#' @field curr_AI The inverse of \code{curr_A}
#' @field curr_logdet The current log-determinant of \code{curr_A}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#'
#' @references
#'  LeSage, J. and Pace, R.K., (2009) Introduction to spatial econometrics. Chapman and Hall/CRC.
rho_samplerMH = R6::R6Class("rho_samplerMH", public = list(
  rho_prior = NULL,
  curr_rho = NULL,
  curr_W = NULL,
  curr_A = NULL,
  curr_AI = NULL,
  curr_logdet = NULL,
  #' @param rho_prior The list returned by \code{\link{rho_priors}}
  #' @param W An optional starting value for the spatial weight matrix \eqn{W}
  #'
  #' @export
  initialize = function(rho_prior,W = NULL) {
    self$rho_prior = rho_prior
    self$curr_W = W
    self$curr_rho = stats::runif(1,self$rho_prior$rho_min,self$rho_prior$rho_max)
    private$curr_rho_scale <- self$rho_prior$init_rho_scale
    private$rho_accept = 0
    private$curr_iter = 0
    private$do_MHtune = TRUE
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
    private$do_MHtune = FALSE
  },
  #' @description
  #' If the spatial weight matrix is updated during the sampling procedure the log determinant, the
  #' spatial projection matrix \eqn{I - \rho W} and it's inverse must be updated. This function should be
  #' used for a consistent update. At least the new spatial weight matrix \code{newW} should be supplied.
  #'
  #' @param newW The updated spatial weight matrix \eqn{W}
  #' @param newLogdet An optional value for the log determinant corresponding to \çode{newW} and \code{curr_rho}
  #' @param newA An optional value for the spatial projection matrix using \çode{newW} and \code{curr_rho}
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
    } else {self$logdet = newLogdet}
    invisible(self)
  },
  #' @param Y The \eqn{n} by \eqn{tt} vector of responses
  #' @param mu The \eqn{n} by \eqn{tt} vector of means (usually \eqn{X\beta})
  #' @param sigma The variance parameter \eqn{\sigma^2}
  #'
  #' @export
  sample = function(Y,mu, sigma) {
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
