#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//adapted from RcppArmadillo sample function
// [[Rcpp::export]]
arma::urowvec ProbSampleReplace(int size, arma::vec prob){
  int nOrig = prob.size();
  arma::urowvec index(size);
  double rU;
  int ii, jj, kk; // indices, ii for loops
  // index tables, fill with zeros
  arma::vec HL_dat(nOrig);
  arma::vec alias_tab(nOrig); 
  arma::vec::iterator H, L, H0, L0;
  //HL0 = HL_dat.begin();
  H0 = H = HL_dat.begin();
  L0 = L = HL_dat.end();
  //prob *= nOrig; // scale probability table
  // fill HL_dat from beginning (small prob) and end (large prob) with indices
  for (ii = 0; ii < nOrig; ii++) {
    prob[ii] *= nOrig;
    if( prob[ii] < 1.0) {
      *(H++) = ii;
    } else {
      *(--L) = ii;
    }
  }
  
  // some of both large and small
  if ( (H > H0) && (L < L0) ) {
    for (kk = 0; kk < nOrig; kk++) {
      ii = HL_dat[kk];
      jj = *L;
      alias_tab[ii] = jj;
      prob[jj] += (prob[ii] - 1);
      if (prob[jj] < 1.) L++;
      if(L == L0) break; // now all prob >= 1
    }
  }
  for (ii = 0; ii < nOrig; ii++)  prob[ii] += ii;
  /* generate sample */
  for (ii = 0; ii < size; ii++) {
    rU = unif_rand() * nOrig;
    kk = (int) rU;
    index[ii] = (rU < prob[kk]) ? kk : alias_tab[kk];
  }
  return index + 1;
}
  


// [[Rcpp::export]]
double maha_density(size_t ix, arma::vec X1, arma::vec X2, arma::mat s, double w = 0.1) {
   arma::vec maha_dist = arma::square(X1 - X1[ix-1]) * s(0,0) + (X1 - X1[ix-1]) % (X2 - X2[ix-1]) * (s(1,0) + s(0,1))  + arma::square(X2 - X2[ix-1]) * s(1,1);
  return arma::mean(arma::exp(-1 * arma::square(maha_dist) / (2 * pow(w,2))) / (sqrt(2 * arma::datum::pi) * w));
}

//// [[Rcpp::export]]
//double maha_density(size_t ix, NumericVector X1, NumericVector X2, NumericMatrix s, double w = 0.1) {
//  NumericVector maha_dist = (pow(X1 - X1[ix-1],2) * s(0,0) + (X1 - X1[ix-1])*(X2 - X2[ix-1])*(s(1,0) + s(0,1))  + pow(X2 - X2[ix-1],2)*s(1,1));
//  return mean(dnorm(maha_dist, 0, w));
//}

// [[Rcpp::export]]
arma::urowvec bg_sample(size_t ix, arma::vec X1, arma::vec X2, arma::mat s, arma::vec dens, int n, double w = 0.1) {
  arma::vec maha_dist = arma::square(X1 - X1[ix-1]) * s(0,0) + (X1 - X1[ix-1]) % (X2 - X2[ix-1]) * (s(1,0) + s(0,1))  + arma::square(X2 - X2[ix-1]) * s(1,1);
  arma::vec p = (arma::exp(-1 * arma::square(maha_dist) / (2 * pow(w,2))) / (sqrt(2 * arma::datum::pi) * w)) / dens;
  p /= arma::sum(p);
  arma::urowvec res = ProbSampleReplace( n, p);
  return res;
}

//
//// [[Rcpp::export]]
//NumericVector bg_sample(size_t ix, NumericVector X1, NumericVector X2, NumericMatrix s, NumericVector dens, int n, double w = 0.1) {
//  NumericVector maha_dist = (pow(X1 - X1[ix-1],2) * s(0,0) + (X1 - X1[ix-1])*(X2 - X2[ix-1])*(s(1,0) + s(0,1))  + pow(X2 - X2[ix-1],2)*s(1,1));
//  NumericVector p = dnorm(maha_dist, 0, w) / dens;
//  IntegerVector indices = seq_along(p);
//  NumericVector res = RcppArmadillo::sample(as<NumericVector>(indices), n, TRUE, p);
//  return res;
//}
//

// [[Rcpp::export]]
arma::vec test_sub(arma::vec x, arma::uvec y){
  arma::vec out = x(y);
  return out;
}

// [[Rcpp::export]]
arma::umat bg_sample_helper(arma::uvec bin_membership, arma::mat bin_p, arma::vec bin_density, arma::uword niterations){
  arma::uword n = bin_membership.n_elem;
  arma::umat out = arma::umat(n, niterations);
  for (arma::uword i = 0; i < bin_density.n_elem; i++){
    arma::uvec ix = find(bin_membership == i);
    arma::vec p_tmp = bin_p.col(i);
    arma::vec p = p_tmp(bin_membership ) / bin_density(bin_membership );
    p /= arma::sum(p);
    arma::urowvec sampled = ProbSampleReplace(niterations * ix.n_elem, p);
    for (arma::uword j = 0; j < ix.n_elem; j++){
      out.row(ix(j)) = sampled(arma::span(j*niterations,((j+1)*niterations) - 1));
    }                                                             
  }
  return out;
}

