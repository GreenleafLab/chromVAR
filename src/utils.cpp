// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector col_sds( arma::mat& X, bool na_rm = false) {
  arma::vec res(X.n_cols);
  if (na_rm){
    for( int j=0; j < X.n_cols; j++ ) {
      arma::vec tmp = X.col(j);
      arma::uvec keep = arma::find_finite(tmp);
      if (keep.size() > 1){
        res(j) = arma::stddev(tmp(keep));        
      } else{
        res(j) = arma::datum::nan;
      }
    }
  } else{
    res = arma::stddev(X,0,0).t();
  }
  return NumericVector(res.begin(),res.end()); 
 }

// [[Rcpp::export]]
NumericVector row_sds( arma::mat& X, bool na_rm = false) {
  arma::vec res(X.n_rows);
  if (na_rm){
    for( int j=0; j < X.n_rows; j++ ) {
      arma::rowvec tmp = X.row(j);
      arma::uvec keep = arma::find_finite(tmp);
      if (keep.size() > 1){
        res(j) = arma::stddev(tmp(keep));
      } else{
        res(j) = arma::datum::nan;
      }
    }
  } else{
    res = arma::stddev(X,0,1);
  }
  return NumericVector(res.begin(),res.end()); 
 }

// [[Rcpp::export]]
NumericVector row_sds_perm( arma::mat& X, bool na_rm = false) {
  arma::uvec tmp = arma::linspace<arma::uvec>(0, X.n_cols - 1, X.n_cols);
  arma::uvec ix = RcppArmadillo::sample(tmp, X.n_cols, true);
  arma::mat shuffled = X.cols(ix);
  return row_sds(shuffled, na_rm); 
 }
 
 
 
 
 