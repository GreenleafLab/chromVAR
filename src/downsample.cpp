#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix downsample_dense(NumericMatrix X, double p){
  NumericMatrix out(X.nrow(),X.ncol());
  for (size_t i = 0; i < X.nrow(); i++){
    for (size_t j = 0; j < X.ncol(); j++){
      out(i,j) = sum(runif(X(i,j)) < p);
    }
  }
  return out;
}
