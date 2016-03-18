//// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
//using namespace Rcpp ;
//
//
//// [[Rcpp::export]]
//arma::sp_mat make_mat(arma::umat locations, arma::uword n_rows, arma::uword n_cols){
//  arma::colvec values = arma::ones<arma::vec>(arma::size(locations)(1));
//  arma::sp_mat out = arma::sp_mat(true, locations, values, n_rows, n_cols);
//  return out;
//}
//
//
//// [[Rcpp::export]]
//arma::mat compute_expected(arma::sp_mat counts, arma::colvec expectation, 
//                            arma::rowvec fragments_per_sample) {
//  arma::mat expected = expectation * fragments_per_sample;
//  arma::mat out = (counts - expected) / sqrt(expected);
//  return out;
//}
//
//
//// [[Rcpp::export]]
//NumericVector rep_example( NumericVector x, int n){
//  NumericVector res = rep_each(x, n) ;
//  return res ;
//}
//
//
//// [[Rcpp::export]]
//arma::urowvec flatten_rows(arma::uvec ix, arma::umat mat){
//  
//  arma::uword nset = ix.n_rows;
//  arma::uword ncol = mat.n_cols;
//  arma::uword n = nset * ncol;
//  arma::urowvec out(n);
//  
//  arma::uword k;
//  for (arma::uword j = 0; j < ncol; j++){
//    for (arma::uword i = 0; i < nset; i++){    
//      k = (j * nset) + i;
//      out(k) = mat(ix[i],j) - 1;
//    }
//  }  
//  return out;
//}
//
//// [[Rcpp::export]]
//arma::urowvec rep_seq(arma::uword n, arma::uword each){
//  IntegerVector tmp = seq_len(n);
//  IntegerVector tmp_i = rep_each(tmp, each) - 1;
//  arma::urowvec out = as<arma::urowvec>(tmp_i);
//  return out;
//}
//
//
//
//// [[Rcpp::export]]
//arma::mat sample_bgp(arma::sp_mat counts, arma::uvec peak_set, arma::umat background_peaks, arma::colvec expectation, 
//                            arma::rowvec fragments_per_sample){
//  
//  arma::uword niterations = background_peaks.n_cols;
//  arma::uword npeaks = background_peaks.n_rows;
//  arma::uword setlen = peak_set.n_rows;
//  
//  IntegerVector tmp = seq_len(niterations);
//  IntegerVector tmp_i = rep_each(tmp, setlen) - 1;
//  arma::urowvec i = as<arma::urowvec>(tmp_i);
//  
//  arma::uword n = setlen * niterations;
//  arma::urowvec j(n);
//  
//  arma::uword z;
//  for (arma::uword k = 0; k < niterations; k++){
//    for (arma::uword l = 0; l < setlen; l++){    
//      z = (k * setlen) + l;
//      j(z) = background_peaks(peak_set[l],k) - 1;
//    }
//  }  
// 
//  arma::umat locs(2, i.n_cols);
//  locs.row(0) = i;
//  locs.row(1) = j;
//  arma::sp_mat sample_mat = make_mat(locs, niterations, npeaks);
//  
//  arma::mat expected = expectation * fragments_per_sample;
//  arma::mat exp_mat = (counts - expected) / sqrt(expected);
//  
//  arma::mat out = sample_mat * exp_mat;
//  
//  return out.t();
//} 
//
//// [[Rcpp::export]]
//arma::mat get_observed(arma::sp_mat counts, arma::urowvec peak_set, arma::colvec expectation, 
//                            arma::rowvec fragments_per_sample){
//        
//  arma::uword npeak = peak_set.n_cols;
//  arma::umat locs(2, npeak);
//  locs.row(0) = arma::zeros<arma::urowvec>(npeak);
//  locs.row(1) = peak_set;
//  arma::sp_mat sample_mat = make_mat(locs, 1, counts.n_rows);     
//        
//  arma::mat expected = expectation * fragments_per_sample;
//  arma::mat exp_mat = (counts - expected) / sqrt(expected);
//  
//  arma::mat out = sample_mat * exp_mat;
//  return out;
//}
//  
//// [[Rcpp::export]]
//arma::urowvec test_i(arma::uword niterations, arma::uword tf_count){
//  IntegerVector tmp = seq_len(niterations);
//  IntegerVector tmp_i = rep_each(tmp, tf_count) - 1;
//  arma::urowvec i = as<arma::urowvec>(tmp_i);  
//  return i;
//}
//
//// [[Rcpp::export]]
//arma::urowvec test_i2(arma::uword niterations, arma::uword tf_count){
//  
//  arma::urowvec i(tf_count*niterations);
//  
//  for (arma::uword ii = 0; ii < niterations; ii++ ){
//    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
//  }
//  return i;
//}  
//  
//  
//// [[Rcpp::export]]
//arma::vec col_means2( arma::mat& X ) {
//  
//  int nCols = X.n_cols;
//  arma::vec out(nCols);
//  
//  for( int j=0; j < nCols; j++ ) {
//    out(j) = mean(X.col(j));
//  }
//  
//  return out;
//  
//}
//
//// [[Rcpp::export]]
//arma::vec col_means3( arma::mat& X ) {
//  
//  int nCols = X.n_cols;
//  arma::vec out(nCols);
//  
//  for( int j=0; j < nCols; j++ ) {
//    out(j) = arma::mean(X.col(j));
//  }
//  
//  return out;
//  
//}
//
//// [[Rcpp::export]]
//arma::vec col_sds2( arma::mat& X ) {
//  
//  int nCols = X.n_cols;
//  arma::vec out(nCols);
//  
//  for( int j=0; j < nCols; j++ ) {
//    out(j) = arma::stddev(X.col(j));
//  }
//  
//  return out;
//  
//}
//
//
//
//
//// [[Rcpp::export]]
//NumericVector col_means( arma::mat& X ) {
//  
//  int nCols = X.n_cols;
//  NumericVector out = no_init(nCols);
//  
//  for( int j=0; j < nCols; j++ ) {
//    NumericVector tmp = NumericVector(X.begin_col(j), X.end_col(j));
//    out[j] = mean(tmp);
//  }
//  
//  return out;
//  
//}
//
//
//// [[Rcpp::export]]
//NumericVector col_sds( arma::mat& X ) {
//  
//  int nCols = X.n_cols;
//  NumericVector out = no_init(nCols);
//  
//  for( int j=0; j < nCols; j++ ) {
//    NumericVector tmp = NumericVector(X.begin_col(j), X.end_col(j));
//    out[j] = sd(tmp);
//  }
//  
//  return out;
//  
//}
//
//
//// [[Rcpp::export]]
//NumericVector compute_deviations_single(arma::urowvec peak_set, 
//                                arma::sp_mat counts, 
//                                arma::umat background_peaks,
//                                arma::vec expectation,
//                                List counts_info){
//  
//  arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
//  arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
//  arma::uword tf_count =  peak_set.n_cols;
//  arma::uword niterations = background_peaks.n_cols;  
//  
//  arma::umat locs1(2, tf_count);
//  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
//  locs1.row(1) = peak_set;
//  
//  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
//  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
//  
//  arma::mat expected = expectation * fragments_per_sample;
//  arma::mat exp_mat = (counts - expected) / sqrt(expected);
//  
//  arma::mat observed_tmp = tf_mat * exp_mat;
//  NumericVector observed = NumericVector(observed_tmp.begin(), observed_tmp.end());
//
//  arma::uword n = tf_count * niterations;
//  arma::urowvec i(n);
//  
//  for (arma::uword ii = 0; ii < niterations; ii++ ){
//    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
//  }
//  
//  arma::urowvec j(n);
//  
//  arma::uword z;
//  for (arma::uword k = 0; k < niterations; k++){
//    for (arma::uword l = 0; l < tf_count; l++){    
//      z = (k * tf_count) + l;
//      j(z) = background_peaks(peak_set[l],k) - 1;
//    }
//  }  
// 
//  arma::umat locs2(2, i.n_cols);
//  locs2.row(0) = i;
//  locs2.row(1) = j;
//  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
//  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);
//
//  arma::mat sampled = sample_mat * exp_mat;
//
//  NumericVector mean_sampled = col_means(sampled);
//  NumericVector sd_sampled = col_sds(sampled);
//  
//  NumericVector res = (observed - mean_sampled) / sd_sampled;
//  
//  return res;
//}
//
//// [[Rcpp::export]]
//arma::mat test_repmat(arma::vec expectation, arma::uword n){
//  arma::mat tmp_mat1 = arma::repmat(sqrt(expectation), 1, n);
//  return tmp_mat1;
//}
//
//// [[Rcpp::export]]
//arma::mat test_repmat2(arma::rowvec fragments_per_sample, arma::uword n){
//  arma::mat tmp_mat2 = arma::repmat(1 / sqrt(fragments_per_sample), n, 1);
//  return tmp_mat2;
//}
//
//// [[Rcpp::export]]
//arma::sp_mat test_norm2(arma::sp_mat counts, arma::vec expectation, arma::rowvec fragments_per_sample){
//  
//  arma::sp_mat::iterator start = counts.begin();
//  arma::sp_mat::iterator end = counts.end();
//
//  arma::uword i;
//  arma::uword j;
//
//  for(arma::sp_mat::iterator it = start; it != end; ++it)
//  {
//    i = it.row();
//    j = it.col();
//    (*it) /= sqrt(expectation(i)*fragments_per_sample(j));
//  }
//  
//  return counts;
//}
//
//
//
//// [[Rcpp::export]]
//arma::sp_mat test_norm(arma::sp_mat counts, arma::rowvec fragments_per_sample){
//  arma::mat tmp_mat2 = arma::repmat(1 / sqrt(fragments_per_sample), counts.n_rows, 1);
//  arma::sp_mat out = counts % tmp_mat2; 
//  return out;
//}
//
//
//// [[Rcpp::export]]
//arma::sp_mat norm_o(arma::sp_mat counts, arma::vec expectation, arma::rowvec fragments_per_sample){
//  arma::mat tmp_mat1 = arma::repmat(sqrt(expectation), 1, counts.n_cols);
//  arma::mat tmp_mat2 = arma::repmat(sqrt(fragments_per_sample), counts.n_rows, 1);
//  arma::sp_mat out =  (counts / tmp_mat1) / tmp_mat2; 
//  return out;
//}
//
//
//
//// [[Rcpp::export]]
//NumericVector compute_deviations_single2(arma::urowvec peak_set, 
//                                arma::sp_mat counts, 
//                                arma::umat background_peaks,
//                                arma::vec expectation,
//                                List counts_info){
//  
//  arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
//  arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
//  arma::uword tf_count =  peak_set.n_cols;
//  arma::uword niterations = background_peaks.n_cols;  
//  
//  //normalize counts by sqrt of expectation
//  arma::sp_mat::iterator start = counts.begin();
//  arma::sp_mat::iterator end = counts.end();
//
//  for(arma::sp_mat::iterator it = start; it != end; ++it)
//  {
//    (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
//  }
//
//  //make peak set mat
//  arma::umat locs1(2, tf_count);
//  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
//  locs1.row(1) = peak_set;
//  
//  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
//  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
//  
//  //get observed
//  arma::mat observed_tmp = (tf_mat * counts) - (tf_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)));
//  NumericVector observed = NumericVector(observed_tmp.begin(), observed_tmp.end());
// 
//  arma::uword n = tf_count * niterations; 
//  arma::urowvec i(n);
//  
//  for (arma::uword ii = 0; ii < niterations; ii++ ){
//    i.subvec(ii*tf_count,(ii+1)*tf_count -1).fill(ii);
//  }
//
//  arma::urowvec j(n);
//  
//  arma::uword z;
//  for (arma::uword k = 0; k < niterations; k++){
//    for (arma::uword l = 0; l < tf_count; l++){    
//      z = (k * tf_count) + l;
//      j(z) = background_peaks(peak_set[l],k) - 1;
//    }
//  }  
// 
//  arma::umat locs2(2, i.n_cols);
//  locs2.row(0) = i;
//  locs2.row(1) = j;
//  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
//  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);
//
//  arma::mat sampled = (sample_mat * counts) - (sample_mat * (expectation/sqrt(expectation)) * (fragments_per_sample/sqrt(fragments_per_sample)));
//
//  NumericVector mean_sampled = col_means(sampled);
//  NumericVector sd_sampled = col_sds(sampled);
//  
//  NumericVector res = (observed - mean_sampled) / sd_sampled;
//  
//  return res;
//}



