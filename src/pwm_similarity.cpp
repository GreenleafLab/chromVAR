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
arma::vec pwm_dist_single(arma::mat mat1, arma::mat mat2){
  
  arma::uword n1 = mat1.n_cols;
  arma::uword n2 = mat2.n_cols;
  
  arma::mat pad1(4,n2-1, arma::fill::zeros);
  arma::mat pad2(4,n1-1, arma::fill::zeros);
  
  arma::mat mat1p = arma::join_horiz(arma::join_horiz(pad1, mat1),pad1);
  arma::mat mat2p = arma::join_horiz(arma::join_horiz(pad2, mat2),pad2);
  
  arma::uword m = n1 + n2- 1;
  
  arma::vec out(2);
  double tmp;
  for (arma::uword i = 0; i < m; i++){
    arma::mat mat1s,mat2s;
    if (i < n2){
      if (i < n1){
        mat1s = mat1p.cols(i, n2 - 1 + n1 - 1);
        mat2s = mat2p.cols(n1 - 1, n1 * 2 + n2 - 3 - i);
      } else{
        mat1s = mat1p.cols(i,  i + n2 - 1);
        mat2s = mat2p.cols(n1-1, n1 * 2 + n2 - 3 - n1 + 1);
      }

    } else{
      if (i < n1){
        mat1s = mat1p.cols(n2 - 1, n2 - 1 + n1 - 1);
        mat2s = mat2p.cols(n1 - 2 - (i - n2), n1 - 2 - (i - n2) + n1 -1);       
      } else{
        mat1s = mat1p.cols(n2 - 1, n2 - 1 + i);
        mat2s = mat2p.cols(n1 - 2 - (i - n2), n1 - 1 + n2 - 1);        
      }  
    }
    //Rcout << "i is " << i << std::endl;    
    //Rcout << "mat1s " << mat1s << std::endl;
    //Rcout << "mat2s " << mat2s << std::endl;   
    tmp = pwm_euclidean(mat1s,mat2s);     
    //Rcout << "score is " << tmp << std::endl;

    if (i == 0){
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
List compute_pwm_dist(List pwms){
  
  arma::uword n = pwms.size();
  arma::mat out1(n, n, arma::fill::zeros);
  arma::imat out2(n, n, arma::fill::zeros);
  CharacterMatrix out3(n,n);
  arma::vec res, res_rc;
  for (arma::uword i = 0; i < n; i++){
    arma::mat mat1 = as<arma::mat>(pwms[i]);
    for (arma::uword j = 0; j <= i; j++){
      arma::mat mat2 = as<arma::mat>(pwms[j]);
      res = pwm_dist_single(mat1, mat2);
      //rc
      res_rc = pwm_dist_single(mat1, arma::fliplr(arma::flipud(mat2)));
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
  return List::create(Rcpp::Named("dist") = out1, Rcpp::Named("offset") = out2, Rcpp::Named("strand") = out3);
}







