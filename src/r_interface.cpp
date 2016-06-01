// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "moods.h"
#include "moods_misc.h"
#include "moods_tools.h"
#include "match_types.h"
#include "motif.h"
#include "scanner.h"

using namespace Rcpp;

std::vector<std::vector<double> > mat_conversion(List& mats, size_t i){
  arma::mat tmp = as<arma::mat>(mats[i]);
  score_matrix out(tmp.n_rows);
    for (size_t i = 0; i < tmp.n_rows; ++i) {
        out[i] = arma::conv_to< std::vector<double> >::from(tmp.row(i));
    };
  return out;
}

// [[Rcpp::export]]
std::vector<double>  get_thresholds(List mats,
          const std::vector<double> nuc_freqs,
          const double p){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[n + i] = thresholds[i];
  }
  return thresholds;
}

// [[Rcpp::export]]
arma::sp_mat get_motif_ix(List mats, const std::vector<std::string> x,
          const std::vector<double> nuc_freqs,
          const double p,
          const size_t w){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[n + i] = thresholds[i];
  }
  MOODS::scan::Scanner scanner = MOODS::scan::Scanner(w);
  scanner.set_motifs(matrices, nuc_freqs, thresholds);
  size_t nstrings = x.size();
  std::vector<unsigned int> iloc;
  std::vector<unsigned int> jloc;
  //iloc.reserve(nstrings * n * 2);
  //jloc.reserve(nstrings * n * 2);
  for (size_t i = 0; i < nstrings; i++){
    auto results = scanner.scan(x[i]);
    for (size_t j = 0; j < n; j++){
      if (results[j].size() > 0){
        iloc.push_back(i);
        jloc.push_back(j);
      } else if (results[j + n].size() > 0){
        iloc.push_back(i);
        jloc.push_back(j);
      }
    }
  }
  arma::umat locs(2, iloc.size());
  locs.row(0) = arma::conv_to< arma::urowvec >::from(iloc);
  locs.row(1) = arma::conv_to< arma::urowvec >::from(jloc);
  arma::vec values(iloc.size(), arma::fill::ones);
  arma::sp_mat out = arma::sp_mat(locs, values, nstrings, n, true, true);
  return out;
}

// [[Rcpp::export]]
arma::sp_mat get_max_motif_score(List& mats, const std::vector<std::string> x,
          const std::vector<double> nuc_freqs,
          const double p,
          const size_t w){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[n + i] = thresholds[i];
  }
  MOODS::scan::Scanner scanner = MOODS::scan::Scanner(w);
  scanner.set_motifs(matrices, nuc_freqs, thresholds);
  size_t nstrings = x.size();
  std::vector<unsigned int> iloc;
  std::vector<unsigned int> jloc;
  std::vector<double> values;
  for (size_t i = 0; i < nstrings; i++){
    auto results = scanner.scan(x[i]);
    for (size_t j = 0; j < n; j++){
      if (results[j].size() > 0){
        double max_result = thresholds[j];
        for (size_t k = 0; k < results[j].size(); k++){
          if (results[j][k].score > max_result){
            max_result = results[j][k].score;
          }
        }
        if (results[n+j].size() > 0){
          for (size_t k = 0; k < results[n+j].size(); k++){
            if (results[n+j][k].score > max_result){
              max_result = results[n+j][k].score;
            }
          }
        }
        iloc.push_back(i);
        jloc.push_back(j);
        values.push_back(max_result);
      } else if (results[n+j].size() > 0){
        double max_result = thresholds[j];
        for (size_t k = 0; k < results[n+j].size(); k++){
          if (results[n+j][k].score > max_result){
            max_result = results[n+j][k].score;
          }
        }
        iloc.push_back(i);
        jloc.push_back(j);
        values.push_back(max_result);
      }
    }
  }
  arma::umat locs(2, iloc.size());
  locs.row(0) = arma::conv_to< arma::urowvec >::from(iloc);
  locs.row(1) = arma::conv_to< arma::urowvec >::from(jloc);
  arma::vec values2 = arma::conv_to< arma::vec>::from(values);
  arma::sp_mat out = arma::sp_mat(locs, values2, nstrings, n, true, true);
  return out;
}


// [[Rcpp::export]]
List get_motif_positions(List& mats, const std::vector<std::string> x,
          const std::vector<double> nuc_freqs,
          const double p,
          const size_t w){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[n + i] = thresholds[i];
  }
  MOODS::scan::Scanner scanner = MOODS::scan::Scanner(w);
  scanner.set_motifs(matrices, nuc_freqs, thresholds);
  size_t nstrings = x.size();
  std::vector<unsigned int> iloc;
  std::vector<unsigned int> jloc;
  std::vector<double> values;
  std::vector<char> strand;
  std::vector<size_t> pos;
  for (size_t i = 0; i < nstrings; i++){
    auto results = scanner.scan(x[i]);
    for (size_t j = 0; j < n; j++){
      if (results[j].size() > 0){
        for (size_t k = 0; k < results[j].size(); k++){
          iloc.push_back(i);
          jloc.push_back(j);
          values.push_back(results[j][k].score);
          strand.push_back('+');
          pos.push_back(results[j][k].pos);
        }
      }
      if (results[n+j].size() > 0){
        for (size_t k = 0; k < results[n+j].size(); k++){
          iloc.push_back(i);
          jloc.push_back(j);
          values.push_back(results[n+j][k].score);
          strand.push_back('-');
          pos.push_back(results[n+j][k].pos);
        }
      }
    }
  }
  return List::create(Named("motif_ix") = jloc,
                      Named("seq_ix") = iloc,
                      Named("strand") = strand,
                      Named("pos") = pos,
                      Named("score") = values);
}

