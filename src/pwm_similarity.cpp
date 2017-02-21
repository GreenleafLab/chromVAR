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
arma::vec pwm_dist_single(arma::mat mat1, arma::mat mat2, 
                          arma::uword min_overlap){

  arma::uword n1 = mat1.n_cols;
  arma::uword n2 = mat2.n_cols;

  arma::uword m = n1 + n2 - 1;

  arma::vec out(2);
  double tmp;
  //arma::uword min_overlap = 4;
  for (arma::uword i = min_overlap - 1; i < m-min_overlap + 1; i++){
    arma::mat mat1s,mat2s;
    if (i < n2){
      if (i < n1){
        mat1s = mat1.cols(0, i);
        mat2s = mat2.cols(n2 - 1 - i, n2 - 1);
      } else{
        mat1s = mat1.cols(0,  n1 - 1);
        mat2s = mat2.cols(n2 - 1 - i, m - i - 1);
      }

    } else{
      if (i < n1){
        mat1s = mat1.cols(i - n2 + 1, i);
        mat2s = mat2.cols(0, n2 - 1);
      } else{
        mat1s = mat1.cols(i - n2 + 1, n1 -1);
        mat2s = mat2.cols(0, m - i - 1);
      }
    }

    tmp = pwm_euclidean(mat1s,mat2s);

    if (i == min_overlap - 1){
      out(1) = (double)i - (double)n2 + 1;
      out(0) = tmp;
    } else if (tmp < out(0)){
      out(1) = (double)i - (double)n2 + 1;
      out(0) = tmp;
    }
  }
  return out;
}



// [[Rcpp::export]]
List compute_pwm_dist(List pwms, arma::uword min_overlap){
  arma::uword n = pwms.size();
  arma::mat out1(n, n, arma::fill::zeros);
  arma::imat out2(n, n, arma::fill::zeros);
  CharacterMatrix out3(n,n);
  arma::vec res, res_rc;
  for (arma::uword i = 0; i < n; i++){
    arma::mat mat1 = as<arma::mat>(pwms[i]);
    for (arma::uword j = 0; j <= i; j++){
      arma::mat mat2 = as<arma::mat>(pwms[j]);
      res = pwm_dist_single(mat1, mat2, min_overlap);
      //rc
      res_rc = pwm_dist_single(mat1, arma::fliplr(arma::flipud(mat2)), 
                               min_overlap);
      if (res(0) <= res_rc(0)){
        out1(i,j) = res(0);
        out1(j,i) = res(0);
        out2(i,j) = (int)res(1);
        out2(j,i) = -(int)res(1);
        out3(i,j) = "+";
        out3(j,i) = "+";
      } else{
        out1(i,j) = res_rc(0);
        out1(j,i) = res_rc(0);
        out2(i,j) = (int)res_rc(1);
        out2(j,i) = -(int)res_rc(1);
        out3(i,j) = "-";
        out3(j,i) = "-";
      }
    }
  }
  return List::create(Rcpp::Named("dist") = out1, Rcpp::Named("offset") = out2, 
                      Rcpp::Named("strand") = out3);
}


// [[Rcpp::export]]
List compute_pwm_dist2(List pwms, List pwms2, arma::uword min_overlap){
  arma::uword n1 = pwms.size();
  arma::uword n2 = pwms2.size();
  arma::mat out1(n1, n2, arma::fill::zeros);
  arma::imat out2(n1, n2, arma::fill::zeros);
  CharacterMatrix out3(n1,n2);
  arma::vec res, res_rc;
  for (arma::uword i = 0; i < n1; i++){
    arma::mat mat1 = as<arma::mat>(pwms[i]);
    for (arma::uword j = 0; j < n2; j++){
      arma::mat mat2 = as<arma::mat>(pwms2[j]);
      res = pwm_dist_single(mat1, mat2, min_overlap);
      //rc
      res_rc = pwm_dist_single(mat1, arma::fliplr(arma::flipud(mat2)), 
                               min_overlap);
      if (res(0) <= res_rc(0)){
        out1(i,j) = res(0);
        out2(i,j) = (int)res(1);
        out3(i,j) = "+";
      } else{
        out1(i,j) = res_rc(0);
        out2(i,j) = (int)res_rc(1);
        out3(i,j) = "-";
      }
    }
  }
  return List::create(Rcpp::Named("dist") = out1, Rcpp::Named("offset") = out2,
                      Rcpp::Named("strand") = out3);
}





