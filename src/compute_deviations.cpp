// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat get_normalized_counts(arma::sp_mat counts, const arma::vec expectation,
                                const arma::rowvec fragments_per_sample){
  arma::sp_mat::iterator start = counts.begin();
  arma::sp_mat::iterator end = counts.end();

  for(arma::sp_mat::iterator it = start; it != end; ++it)
  {
    (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
  }
  return counts;
}

