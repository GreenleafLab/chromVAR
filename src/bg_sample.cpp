#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//adapted from RcppArmadillo sample function
//Method for sampling WITH REPLACEMENT with large vector of probabilities
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
  H0 = H = HL_dat.begin();
  L0 = L = HL_dat.end();
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
  



//Helper Function for background sampler
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





//computes euclidean distance between rows of matrix x
// [[Rcpp::export]]
arma::mat euc_dist(arma::mat x) {
  arma::uword n = x.n_rows;
  arma::mat out(n, n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; i++){
    for (arma::uword j = i + 1; j < n; j ++){
      out(i,j) = sqrt(arma::sum(arma::square(x.row(i) - x.row(j))));
      out(j,i) = out(i,j);
    }
  }
  return out;
}


//Functions for slow calculation of background
// [[Rcpp::export]]
double maha_density(size_t ix, arma::vec X1, arma::vec X2, arma::mat s, double w = 0.1) {
   arma::vec maha_dist = arma::square(X1 - X1[ix-1]) * s(0,0) + (X1 - X1[ix-1]) % (X2 - X2[ix-1]) * (s(1,0) + s(0,1))  + arma::square(X2 - X2[ix-1]) * s(1,1);
  return arma::mean(arma::exp(-1 * arma::square(maha_dist) / (2 * pow(w,2))) / (sqrt(2 * arma::datum::pi) * w));
}


// [[Rcpp::export]]
arma::urowvec bg_sample(size_t ix, arma::vec X1, arma::vec X2, arma::mat s, arma::vec dens, int n, double w = 0.1) {
  arma::vec maha_dist = arma::square(X1 - X1[ix-1]) * s(0,0) + (X1 - X1[ix-1]) % (X2 - X2[ix-1]) * (s(1,0) + s(0,1))  + arma::square(X2 - X2[ix-1]) * s(1,1);
  arma::vec p = (arma::exp(-1 * arma::square(maha_dist) / (2 * pow(w,2))) / (sqrt(2 * arma::datum::pi) * w)) / dens;
  p /= arma::sum(p);
  arma::urowvec res = ProbSampleReplace( n, p);
  return res;
}