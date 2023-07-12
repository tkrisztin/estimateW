#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "logdetAinvUpdate_fast.h"

using namespace Rcpp;
using namespace arma;


// Function to calculate probability delta
double calculateProbDelta(double p1, double p0) {
  double prob_delta = p1 / (p1 + p0);
  if (std::isnan(prob_delta)) {
    prob_delta = 0;
  }
  return prob_delta;
}

//' Function to calculate error
double calculateErr(arma::mat A, arma::mat Y,
                    arma::mat w,
                    arma::uvec ch_elmnt,
                    arma::mat mu, arma::mat lag_mu,
                    double curr_rho,
                    bool spatial_error) {
  double err;
  if (!std::isnan(curr_rho)) {
    if (!spatial_error) {
      // SAR
      err = accu(square(A.rows(ch_elmnt) * Y - mu.rows(ch_elmnt) - w.rows(ch_elmnt) * lag_mu));
    } else {
      // SEM
      err = accu(square(A.rows(ch_elmnt) * (Y - mu - w * lag_mu)));
    }
  } else {
    // SLX
    err = accu(square(Y.rows(ch_elmnt) - w.rows(ch_elmnt) * lag_mu - mu.rows(ch_elmnt)));
  }
  return err;
}

//' Function to calculate posterior probality
double calculateProb( arma::mat A0, arma::mat A1,
                      arma::mat Y,
                      arma::mat w0, arma::mat w1,
                      arma::uvec ch_elmnt,
                      arma::mat mu, arma::mat lag_mu,
                      double bbprior_,
                      double logdet0, double logdet1,
                      double curr_rho, double curr_sigma,
                      arma::uword tt,
                      bool spatial_error
) {

  double err1 = calculateErr(A1, Y, w1, ch_elmnt, mu, lag_mu, curr_rho, spatial_error);
  double err0 = calculateErr(A0, Y, w0, ch_elmnt, mu, lag_mu, curr_rho, spatial_error);

  // Calculate adj
  double adj = std::min(err0, err1);
  err1 -= adj;
  err0 -= adj;

  double p0, p1;
  if (!std::isnan(curr_rho)) {
    p1 = bbprior_ * exp(logdet1 * tt) * exp(-err1 / (2 * curr_sigma));
    p0 = (1 - bbprior_) * exp(logdet0 * tt) * exp(-err0 / (2 * curr_sigma));
  } else {
    p1 = bbprior_ * exp(-err1 / (2 * curr_sigma));
    p0 = (1 - bbprior_) * exp(-err0 / (2 * curr_sigma));
  }
  // Calculate probability delta
  return calculateProbDelta(p1, p0);
}


//' A fast sampling step implemented in C++ for updating the spatial weight matrix.
//'
//' This function is intended to be called from the R6 class \code{\link{W_sampler}}.
//'
//' @param Y The \eqn{n} by \eqn{tt} matrix of responses
//' @param curr_sigma The variance parameter \eqn{\sigma^2}
//' @param mu The \eqn{n} by \eqn{tt} matrix of means.
//' @param lag_mu \eqn{n} by \eqn{tt} matrix of means that will be spatially lagged with
//' the estimated \eqn{W}. Defaults to a matrix with zero elements.
//' @param W_prior The current \code{\link{W_priors}}
//' @param curr_W binary \eqn{n} by \eqn{n} spatial connectivity matrix \eqn{\Omega}
//' @param curr_w Row-standardized \eqn{n} by \eqn{n} spatial weight matrix \eqn{W}
//' @param curr_A The current spatial projection matrix \eqn{I - \rho W}.
//' @param curr_AI The inverse of \code{curr_A}
//' @param curr_logdet The current log-determinant of \code{curr_A}
//' @param curr_rho single number between -1 and 1 or NULL, depending on whether the sampler updates
//' the spatial autoregressive parameter \eqn{\rho}. Set while invoking \code{initialize}
//' or using the function \code{set_rho}.
//' @param nr_neighbors_prior An \eqn{n} dimensional vector of prior weights on the number of neighbors
//' (i.e. the row sums of the adjacency matrix \eqn{\Omega}), where the first element denotes the prior probability
//' of zero neighbors and the last those of \eqn{n-1}. A prior using only fixed inclusion probabilities
//' for the entries in \eqn{\Omega} would be an \eqn{n} dimensional vector of \eqn{1/n}. Defaults to
//' a \code{\link{bbinompdf}} prior, with prior parameters \eqn{a = 1}, \eqn{b = 1}.
//' @param symmetric Binary value. Should the estimated adjacency matrix \eqn{\Omega} be symmetric (default: FALSE)?
//' if TRUE: \eqn{\Omega} is forced symmetric; if FALSE: \eqn{\Omega} not necessarily symmetric.
//' @param spatial_error Should a spatial error model be constructed? Defaults to \code{FALSE}.
//' @param row_standardized Binary value. Should the estimated \eqn{W} matrix be row-standardized (default: TRUE)?
//' if TRUE: row-stochastic \eqn{W}; if FALSE: \eqn{W} not row-standardized.
//'
//' @export
// [[Rcpp::export]]
List sampleW_fast(arma::mat Y,
            double curr_sigma,
            arma::mat mu,
            arma::mat lag_mu,
            arma::mat W_prior,
            arma::mat curr_W,
            arma::mat curr_w,
            arma::mat curr_A,
            arma::mat curr_AI,
            double curr_logdet,
            double curr_rho,
            arma::vec nr_neighbors_prior,
            bool symmetric,
            bool spatial_error,
            bool row_standardized) {

  // Declare variables used
  arma::mat w0, w1, A1, A0, W0, W1, WW0, WW1, diff0, diff1;
  double logdet0, logdet1;
  double bbprior0, bbprior1, bbprior_, prob_delta, rnd_draw;
  List res1, res0;
  arma::uword ii,jj, neighb0, neighb1, ch;
  arma::uvec jj_samples, ii_samples;
  bool was1;

  arma::uword n = Y.n_rows;
  arma::uword tt = Y.n_cols;

  if (symmetric) {
    ii_samples = arma::shuffle(arma::regspace<arma::uvec>(1, n - 1));
  } else {
    ii_samples = arma::shuffle(arma::regspace<arma::uvec>(0, n - 1));
  }

  for (int i = 0; i < ii_samples.n_elem; i++) {
    ii = ii_samples(i);

    if (symmetric) {
      jj_samples = arma::shuffle(arma::regspace<arma::uvec>(0, ii - 1));
    } else {
      jj_samples = arma::shuffle(arma::regspace<arma::uvec>(0, n - 1));
    }

    for (int j = 0; j < jj_samples.n_elem; j++) {
      jj = jj_samples(j);

      if (W_prior(ii, jj) == 0) {
        curr_W(ii, jj) = 0;
      } else if (W_prior(ii, jj) == 1) {
        curr_W(ii, jj) = 1;
      } else {
        arma::uvec ch_elmnt;
        if (symmetric) {
          ch_elmnt = arma::uvec{ii, jj};
        } else {
          ch_elmnt = arma::uvec{ii};
        }

        W0 = curr_W;
        W1 = curr_W;
        w0 = curr_w;
        w1 = curr_w;
        was1 = (curr_W(ii, jj) == 1);

        if (was1) {
          W0(ii, jj) = 0;
          if (symmetric) {
            WW0 = W0 + trans(W0);
          } else {
            WW0 = W0;
          }

          for (int k = 0; k < ch_elmnt.n_elem; k++) {
            ch = ch_elmnt(k);
            if (row_standardized) {
              w0.row(ch) = WW0.row(ch) / sum(WW0.row(ch));
            } else {
              w0.row(ch) = WW0.row(ch);
            }
          }

          w0.elem(find_nonfinite(w0)).zeros();

          if (!std::isnan(curr_rho)) {
            A0 = arma::eye(n, n) - curr_rho * w0;
            A1 = curr_A;
            diff0 = A0.rows(ch_elmnt) - curr_A.rows(ch_elmnt);
            if (symmetric) {
              res0 = logdetAinvUpdate_fast(ch_elmnt, diff0, curr_AI, curr_logdet);
              logdet0 = res0["logdet"];
            } else {
              logdet0 = logdetUpdate_fast(ch_elmnt, diff0, curr_AI, curr_logdet);
            }
            logdet1 = curr_logdet;
          }
        } else {
          W1(ii, jj) = 1;
          if (symmetric) {
            WW1 = W1 + trans(W1);
          } else {
            WW1 = W1;
          }

          for (int k = 0; k < ch_elmnt.n_elem; k++) {
            arma::uword ch = ch_elmnt(k);
            if (row_standardized) {
              w1.row(ch) = WW1.row(ch) / sum(WW1.row(ch));
            } else {
              w1.row(ch) = WW1.row(ch);
            }
          }

          w1.elem(find_nonfinite(w1)).zeros();

          if (!std::isnan(curr_rho)) {
            A1 = arma::eye(n, n) - curr_rho * w1;
            A0 = curr_A;
            diff1 = A1.rows(ch_elmnt) - curr_A.rows(ch_elmnt);
            if (symmetric) {
              res1 = logdetAinvUpdate_fast(ch_elmnt, diff1, curr_AI, curr_logdet);
              logdet1 = res1["logdet"];
            } else {
              logdet1 = logdetUpdate_fast(ch_elmnt, diff1, curr_AI, curr_logdet);
            }
            logdet0 = curr_logdet;
          }
        }

        if (symmetric) {
          neighb0 = sum(W0.row(ii)) + sum(W0.col(ii));
          neighb1 = sum(W1.row(ii)) + sum(W1.col(ii));
        } else {
          neighb0 = sum(W0.row(ii));
          neighb1 = sum(W1.row(ii));
        }

        bbprior0 = nr_neighbors_prior(neighb0) * (1 - W_prior(ii, jj));
        bbprior1 = nr_neighbors_prior(neighb1) * W_prior(ii, jj);
        bbprior_ = bbprior1 / (bbprior1 + bbprior0);


        // Calculate probability delta
        prob_delta = calculateProb(A0, A1,
                                          Y,w0,w1,
                                           ch_elmnt, mu, lag_mu,
                                           bbprior_, logdet0, logdet1,
                                           curr_rho, curr_sigma, tt,
                                           spatial_error);

        // Generate random number
        rnd_draw = arma::randu<double>();

        if (rnd_draw <= prob_delta) {
          curr_W(ii, jj) = 1;
          if (!was1) {
            curr_w = w1;
            if (!std::isnan(curr_rho)) {
              curr_logdet = logdet1;
              curr_A = A1;
              if (symmetric) {
                curr_AI = as<arma::mat>(res1["AI"]);
              } else {
                curr_AI = AinvUpdate_fast(ch_elmnt, diff1, curr_AI);
              }
            }
          }
        } else {
          curr_W(ii, jj) = 0;
          if (was1) {
            curr_w = w0;
            if (!std::isnan(curr_rho)) {
              curr_logdet = logdet0;
              curr_A = A0;
              if (symmetric) {
                curr_AI = as<arma::mat>(res0["AI"]);
              } else {
                curr_AI = AinvUpdate_fast(ch_elmnt, diff0, curr_AI);
              }
            }
          }
        }
      }
    }
  }

  // Calculate row sums
  arma::vec row_sums = sum(curr_W + trans(curr_W), 1);

  if (row_standardized) {
    if (symmetric) {
      curr_w = (curr_W + trans(curr_W)) / arma::repmat(row_sums, 1, curr_W.n_cols);
      curr_w.elem(find_nonfinite(curr_w)).zeros();  // Set NaN values to zero
    } else {
      curr_w = curr_W / arma::repmat(sum(curr_W, 1), 1, curr_W.n_cols);
      curr_w.elem(find_nonfinite(curr_w)).zeros();  // Set NaN values to zero
    }
  } else {
    if (symmetric) {
      curr_w = curr_W + trans(curr_W);
    } else {
      curr_w = curr_W;
    }
  }

  return List::create(Named("curr_W") = curr_W, Named("curr_w") = curr_w,
                      Named("curr_A") = curr_A, Named("curr_AI") = curr_AI,
                      Named("curr_logdet") = curr_logdet);
}
