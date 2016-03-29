// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double pwm_euclidean(arma::mat mat1, arma::mat mat2){
  arma::mat diff_mat = arma::square(mat1 - mat2);
  int width = mat1.n_cols;
  double raw_dist = arma::sum(arma::sqrt((arma::rowvec)arma::sum(diff_mat,0)));
  double norm_dist = raw_dist / sqrt(2) / width;
  return norm_dist;
}


// [[Rcpp::export]]
double pwm_dist_single(arma::mat mat1, arma::mat mat2){
  
  double out;
  if (mat1.n_cols > mat2.n_cols){
    for (arma::uword i = 0; i <= mat1.n_cols - mat2.n_cols; i++){
      double tmp = pwm_euclidean((arma::mat)mat1.cols(i, i + mat2.n_cols - 1), mat2);
      if ( i == 0  or tmp < out){
        out = tmp;
      }
    }    
  } else if (mat1.n_cols < mat2.n_cols){
    for (arma::uword i = 0; i <= mat2.n_cols - mat1.n_cols; i++){
      double tmp = pwm_euclidean(mat1, (arma::mat)mat2.cols(i, i + mat1.n_cols - 1));
      if ( i == 0  or tmp < out){
        out = tmp;
      }
    }  
  } else{
    out = pwm_euclidean(mat1, mat2);
  }
  return(out);
}


// [[Rcpp::export]]
arma::mat pwm_dist(List pwms){
  
  arma::uword n = pwms.size();
  arma::mat out(n, n);
  for (arma::uword i = 1; i < n; i++){
    arma::mat mat1 = as<arma::mat>(pwms[i]);
    for (arma::uword j = 0; j < i; j++){
      arma::mat mat2 = as<arma::mat>(pwms[j]);
      out(i,j) = pwm_dist_single(mat1, mat2);
      out(j,i) = out(i,j);
    }
  } 
  return(out);
}
