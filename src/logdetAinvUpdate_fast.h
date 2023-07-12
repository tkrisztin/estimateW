// Defines a header file containing function signatures for functions in src/

// Protect signatures using an inclusion guard.
#ifndef logdetAinvUpdate_fast_H
#define logdetAinvUpdate_fast_H

Rcpp::List logdetAinvUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI, double logdet);

arma::mat AinvUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI);

double logdetUpdate_fast(arma::uvec ch_ind, arma::mat diff, arma::mat AI, double logdet);

#endif
