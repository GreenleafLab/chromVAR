// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
List compute_deviations_single_dense(const arma::urowvec peak_set, 
                                arma::mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info,
                                bool intermediates,
                                bool norm){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]);
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  arma::mat expected_counts = expectation * fragments_per_sample;
  
  arma::vec observed;
  arma::vec expected;
  if (norm){
    observed = (tf_mat * (counts/arma::sqrt(expected_counts))).t();    
    expected = (tf_mat * (expected_counts/arma::sqrt(expected_counts))).t();
  } else{
    observed = (tf_mat * counts).t();    
    expected = (tf_mat * expected_counts).t();
  }
  arma::vec observed_deviation = observed - expected;

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
  
  arma::mat sampled, sampled_expected;
  if (norm){
    sampled = sample_mat * (counts/arma::sqrt(expected_counts));
    sampled_expected = sample_mat * (expected_counts/arma::sqrt(expected_counts));
  } else{
    sampled = sample_mat * counts;
    sampled_expected = sample_mat * expected_counts;   
  }
  arma::mat sampled_deviation = sampled - sampled_expected;

  //get log2 fold change 
  arma::vec logfc = arma::log2(observed/expected) - arma::log2(arma::mean(sampled/sampled_expected,0).t());
   
  //get z-score
  arma::vec mean_sampled_deviation = arma::mean(sampled_deviation,0).t();
  arma::vec sd_sampled_deviation = arma::stddev(sampled_deviation,0).t();
    
  arma::vec zscore = (observed_deviation - mean_sampled_deviation) / sd_sampled_deviation;
  
  List out;
  if (intermediates){
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()),
                      Named("observed") =  NumericVector(observed.begin(), observed.end()),
                      Named("sampled") = sampled,
                      Named("expected") =  NumericVector(expected.begin(), expected.end()),
                      Named("sampled_expected") = sampled_expected,
                      Named("observed_deviation") =  NumericVector(observed_deviation.begin(), observed_deviation.end()),
                      Named("sampled_deviation") = sampled_deviation);
  } else{
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()));
  }  
  return out;
}



// [[Rcpp::export]]
List compute_deviations_single_sparse(const arma::urowvec peak_set, 
                                arma::sp_mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info,
                                bool intermediates,
                                bool norm){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]); 
  const arma::uword tf_count =  peak_set.n_cols;
  const arma::uword niterations = background_peaks.n_cols;  
  
  //normalize counts by sqrt of expectation if norm 
  if (norm){
    arma::sp_mat::iterator start = counts.begin();
    arma::sp_mat::iterator end = counts.end();

    for(arma::sp_mat::iterator it = start; it != end; ++it)
    {
      (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
    }    
  }

  //make peak set mat
  arma::umat locs1(2, tf_count);
  locs1.row(0) = arma::zeros<arma::urowvec>(tf_count);
  locs1.row(1) = peak_set;
  
  arma::colvec values1 = arma::ones<arma::vec>(tf_count);
  arma::sp_mat tf_mat = arma::sp_mat(true, locs1, values1, 1, npeaks);
  
  //get observed
  arma::vec observed = ((arma::mat)(tf_mat * counts)).t();
  
  arma::vec expected;
  if (norm){
    expected = (tf_mat * (expectation/arma::sqrt(expectation)) * (fragments_per_sample/arma::sqrt(fragments_per_sample))).t();
  } else{
    expected = (tf_mat * expectation * fragments_per_sample).t();
  }
  arma::vec observed_deviation = observed - expected;
 
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
  arma::mat sampled = (arma::mat)(sample_mat * counts);
  arma::mat sampled_expected;
  if (norm){
    sampled_expected = sample_mat * (expectation/arma::sqrt(expectation)) * (fragments_per_sample/arma::sqrt(fragments_per_sample));
  } else{
    sampled_expected = sample_mat * expectation * fragments_per_sample;    
  }
  arma::mat sampled_deviation = sampled - sampled_expected;

  //get log2 fold change
  arma::vec logfc = arma::log2(observed/expected) - arma::log2(arma::mean(sampled/sampled_expected,0).t());
   
  //get z-score
  arma::vec mean_sampled_deviation = arma::mean(sampled_deviation,0).t();
  arma::vec sd_sampled_deviation = arma::stddev(sampled_deviation,0).t();
  
  arma::vec zscore = (observed_deviation - mean_sampled_deviation) / sd_sampled_deviation;
  
  List out;
  if (intermediates){
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()),
                      Named("observed") =  NumericVector(observed.begin(), observed.end()),
                      Named("sampled") = sampled,
                      Named("expected") =  NumericVector(expected.begin(), expected.end()),
                      Named("sampled_expected") = sampled_expected,
                      Named("observed_deviation") =  NumericVector(observed_deviation.begin(), observed_deviation.end()),
                      Named("sampled_deviation") = sampled_deviation);
  } else{
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()));
  }  
  return out;
}


// [[Rcpp::export]]
List compute_deviations_single_peak_sparse(const arma::uword peak, 
                                arma::sp_mat counts, 
                                const arma::umat background_peaks,
                                const arma::vec expectation,
                                const List counts_info,
                                bool intermediates,
                                bool norm){
  
  const arma::rowvec fragments_per_sample = as<arma::rowvec>(counts_info["fragments_per_sample"]);
  const arma::uword npeaks = as<arma::uword>(counts_info["npeak"]); 
  const arma::uword niterations = background_peaks.n_cols;  
  
  //normalize counts by sqrt of expectation if norm 
  if (norm){
    arma::sp_mat::iterator start = counts.begin();
    arma::sp_mat::iterator end = counts.end();

    for(arma::sp_mat::iterator it = start; it != end; ++it)
    {
      (*it) /= sqrt(expectation(it.row())*fragments_per_sample(it.col()));
    }    
  }

  arma::vec observed = ((arma::mat)counts.row(peak)).t();
  arma::vec expected;
  if (norm){
    expected = (arma::sum(expectation(peak)) * fragments_per_sample).t();
  } else{
    expected = (arma::sum(expectation(peak)/expectation(peak)) * (fragments_per_sample/arma::sqrt(fragments_per_sample))).t();    
  } 
  arma::vec observed_deviation = observed - expected;  
  
  //get sampled peak set mat
  arma::uword n = niterations; 
  arma::urowvec i = arma::linspace<arma::urowvec>(0,n-1,n);
  arma::urowvec j = background_peaks.row(peak) - 1;

  arma::umat locs2(2, i.n_cols);
  locs2.row(0) = i;
  locs2.row(1) = j;
  arma::colvec values2 = arma::ones<arma::vec>(i.n_cols);
  arma::sp_mat sample_mat = arma::sp_mat(true, locs2, values2, niterations, npeaks);

  //get sampled deviations
  arma::mat sampled = (arma::mat)(sample_mat * counts);
  arma::mat sampled_expected;
  if (norm){
    sampled_expected = sample_mat * (expectation/arma::sqrt(expectation)) * (fragments_per_sample/arma::sqrt(fragments_per_sample));
  } else{
    sampled_expected = sample_mat * expectation * fragments_per_sample;    
  }
  arma::mat sampled_deviation = sampled - sampled_expected;

  //get log2 fold change
  arma::vec logfc = arma::log2(observed/expected) - arma::log2(arma::mean(sampled/sampled_expected,0).t());
   
  //get z-score
  arma::vec mean_sampled_deviation = arma::mean(sampled_deviation,0).t();
  arma::vec sd_sampled_deviation = arma::stddev(sampled_deviation,0).t();
  
  arma::vec zscore = (observed_deviation - mean_sampled_deviation) / sd_sampled_deviation;
  
  List out;
  if (intermediates){
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()),
                      Named("observed") =  NumericVector(observed.begin(), observed.end()),
                      Named("sampled") = sampled,
                      Named("expected") =  NumericVector(expected.begin(), expected.end()),
                      Named("sampled_expected") = sampled_expected,
                      Named("observed_deviation") =  NumericVector(observed_deviation.begin(), observed_deviation.end()),
                      Named("sampled_deviation") = sampled_deviation);
  } else{
    out = List::create(Named("z") =  NumericVector(zscore.begin(), zscore.end()),
                      Named("fc") =  NumericVector(logfc.begin(), logfc.end()));
  }  
  return out;
}
