// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec col_means( arma::mat& X ) {
  
  arma::uword nCols = X.n_cols;
  arma::vec out(nCols);
  
  for( int j=0; j < nCols; j++ ) {
    out(j) = arma::mean(X.col(j));
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::vec col_sds( arma::mat& X ) {
  
  arma::uword nCols = X.n_cols;
  arma::vec out(nCols);
  
  for( int j=0; j < nCols; j++ ) {
    out(j) = arma::stddev(X.col(j));
  }
  
  return out;
  
}

// [[Rcpp::export]]
NumericVector compute_deviations_single_dense(const arma::urowvec peak_set, 
                                const arma::mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  arma::mat expected = expectation * fragments_per_sample;
  arma::mat exp_mat = (counts - expected) / sqrt(expected);
  
  arma::vec observed = (tf_mat * exp_mat).t();

  arma::uword n = tf_count * niterations;
  arma::urowvec i(n);
  
  for (arma::uword ii = 0; ii < niterations; ii++ ){
    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
  }
  
  arma::urowvec j(n);
  
  arma::uword z;
  for (arma::uword k = 0; k < niterations; k++){
    for (arma::uword l = 0; l < tf_count; l++){    
      z = (k * tf_count) + l;
      j(z) = background_peaks(peak_set[l],k) - 1;
    }
  }  
 
  arma::umat locs2(2, i.n_cols);
  locs2.row(0) = i;
  locs2.row(1) = j;
  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);

  arma::mat sampled = sample_mat * exp_mat;

  arma::vec mean_sampled = col_means(sampled);
  arma::vec sd_sampled = col_sds(sampled);
  
  arma::vec res = (observed - mean_sampled) / sd_sampled;
  
  return NumericVector(res.begin(), res.end());
}


// [[Rcpp::export]]
List compute_deviations_single_dense_with_intermediates(const arma::urowvec peak_set, 
                                const arma::mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  arma::mat expected = expectation * fragments_per_sample;
  arma::mat exp_mat = (counts - expected) / sqrt(expected);
  
  arma::vec observed = (tf_mat * exp_mat).t();

  arma::uword n = tf_count * niterations;
  arma::urowvec i(n);
  
  for (arma::uword ii = 0; ii < niterations; ii++ ){
    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
  }
  
  arma::urowvec j(n);
  
  arma::uword z;
  for (arma::uword k = 0; k < niterations; k++){
    for (arma::uword l = 0; l < tf_count; l++){    
      z = (k * tf_count) + l;
      j(z) = background_peaks(peak_set[l],k) - 1;
    }
  }  
 
  arma::umat locs2(2, i.n_cols);
  locs2.row(0) = i;
  locs2.row(1) = j;
  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);

  arma::mat sampled = sample_mat * exp_mat;

  arma::vec mean_sampled = col_means(sampled);
  arma::vec sd_sampled = col_sds(sampled);
  
  arma::vec res = (observed - mean_sampled) / sd_sampled;
  
  return List::create(Named("deviations") =  NumericVector(res.begin(), res.end()),
                      Named("observed") =  NumericVector(observed.begin(), observed.end()),
                      Named("sampled") = sampled);
}


// [[Rcpp::export]]
NumericVector compute_deviations_single_sparse(const arma::urowvec peak_set, 
                                arma::sp_mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  //normalize counts by sqrt of expectation
  arma::sp_mat::iterator start = counts.begin();
  arma::sp_mat::iterator end = counts.end();

  for(arma::sp_mat::iterator it = start; it != end; ++it)
  {
    (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
  }

  //make peak set mat
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  //get observed
  arma::vec observed = ((tf_mat * counts) - (tf_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)))).t();
 
  //get sampled peak set mat
  arma::uword n = tf_count * niterations; 
  arma::urowvec i(n);
  
  for (arma::uword ii = 0; ii < niterations; ii++ ){
    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
  }

  arma::urowvec j(n);
  
  arma::uword z;
  for (arma::uword k = 0; k < niterations; k++){
    for (arma::uword l = 0; l < tf_count; l++){    
      z = (k * tf_count) + l;
      j(z) = background_peaks(peak_set[l],k) - 1;
    }
  }  
 
  arma::umat locs2(2, i.n_cols);
  locs2.row(0) = i;
  locs2.row(1) = j;
  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);

  //get sampled deviations
  arma::mat sampled = (sample_mat * counts) - (sample_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)));

  //get z-score
  arma::vec mean_sampled = col_means(sampled);
  arma::vec sd_sampled = col_sds(sampled);
  
  arma::vec res = (observed - mean_sampled) / sd_sampled;
  
  return  NumericVector(res.begin(), res.end());
}


// [[Rcpp::export]]
List compute_deviations_single_sparse_with_intermediate(const arma::urowvec peak_set, 
                                arma::sp_mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info){

  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  //normalize counts by sqrt of expectation
  arma::sp_mat::iterator start = counts.begin();
  arma::sp_mat::iterator end = counts.end();

  for(arma::sp_mat::iterator it = start; it != end; ++it)
  {
    (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
  }

  //make peak set mat
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  //get observed
  arma::vec observed = ((tf_mat * counts) - (tf_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)))).t();
 
  //get sampled peak set mat
  arma::uword n = tf_count * niterations; 
  arma::urowvec i(n);
  
  for (arma::uword ii = 0; ii < niterations; ii++ ){
    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
  }

  arma::urowvec j(n);
  
  arma::uword z;
  for (arma::uword k = 0; k < niterations; k++){
    for (arma::uword l = 0; l < tf_count; l++){    
      z = (k * tf_count) + l;
      j(z) = background_peaks(peak_set[l],k) - 1;
    }
  }  
 
  arma::umat locs2(2, i.n_cols);
  locs2.row(0) = i;
  locs2.row(1) = j;
  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);

  //get sampled deviations
  arma::mat sampled = (sample_mat * counts) - (sample_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)));

  //get z-score
  arma::vec mean_sampled = col_means(sampled);
  arma::vec sd_sampled = col_sds(sampled);
  
  arma::vec res = (observed - mean_sampled) / sd_sampled;

  return List::create(Named("deviations") =  NumericVector(res.begin(), res.end()),
                      Named("observed") =  NumericVector(observed.begin(), observed.end()),
                      Named("sampled") = sampled);
                                  
}



