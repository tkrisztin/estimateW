#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



List logdetAinvUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI, double logdet) {
  int n = ch_ind.size();
  arma::rowvec diffAI;

  for (int i = 0; i < n; i++) {
    int ii = ch_ind[i];
    diffAI = diff.row(i) * AI;

    if (diffAI(ii) <= -1) {
      return List::create(Named("AI") = AI, Named("logdet") = NA_REAL);
    } else {
      logdet += log(1 + diffAI(ii));
    }

    AI = AI - AI.col(ii) * diffAI / (1 + diffAI(ii));
  }

  return List::create(Named("AI") = AI, Named("logdet") = logdet);
}

double logdetUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI, double logdet) {
   int ii = ch_ind[0];
   arma::rowvec diffAI = diff.row(0) * AI;

   if (diffAI(ii) <= -1) {
     return NA_REAL;
   } else {
     logdet += log(1 + diffAI(ii));
   }
   return logdet;
}

 // [[Rcpp::export]]
arma::mat AinvUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI) {
  arma::rowvec diffAI;

  int ii = ch_ind[0];
  diffAI = diff.row(0) * AI;

  if (diffAI(ii) <= -1) {
   return AI;
  }
  AI = AI - AI.col(ii) * diffAI / (1 + diffAI(ii));

  return AI;
}
