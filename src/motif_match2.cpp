// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double p_to_score2(const arma::mat pwm, const arma::vec bg, const double p)
{
    const double PVAL_DP_MULTIPLIER = 1000.0;
    
    arma::uword a = pwm.n_rows;
    arma::uword n = pwm.n_cols;

    arma::imat mat(a,n);

    int maxT = 0;
    int minV = INT_MAX;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            if (pwm(j,i) > 0.0){
                mat(j,i) = (int) ( PVAL_DP_MULTIPLIER * pwm(j,i) + 0.5 );
            }
            else {
                mat(j,i) = (int) ( PVAL_DP_MULTIPLIER * pwm(j,i) - 0.5 );
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        int max = mat(0,i);
        int min = max;
        for (size_t j = 1; j < a; ++j)
        {
            int v = mat(j,i);
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    int R = maxT - n * minV;

    arma::vec table0(R + 1, arma::fill::zeros);
    arma::vec table1(R + 1, arma::fill::zeros);

    for (arma::uword j = 0; j < a; ++j)
        table0[mat(j,0) - minV] += bg[j];

    for (arma::uword i = 1; i < n; ++i)
    {
        for (arma::uword j = 0; j < a; ++j)
        {
            int s = mat(j,i) - minV;
            for (int r = s; r <= R; ++r)
                table1[r] += bg[j] * table0[r - s];
        }
        for (int r = 0; r <= R; ++r)
        {
            table0[r] = table1[r];
            table1[r] = 0.0;
        }
    }

    double sum = 0.0;

    for (int r = R; r >= 0; --r)
    {
        sum += table0[r];
        if (sum > p)
        {
            return (double) ((r + n * minV + 1) / PVAL_DP_MULTIPLIER);
        }
    }

    return (double) ((n * minV) / PVAL_DP_MULTIPLIER);
}


// [[Rcpp::export]]
arma::field<arma::mat> list_to_field(List mats, bool rc = false){
  arma::uword n = mats.size(); 
  arma::field<arma::mat> tmp(n);

  for(arma::uword i = 0; i < n; i++) {
    if (rc){
      tmp(i) = arma::flipud(arma::fliplr(as<arma::mat>(mats(i))));
    } else{
      tmp(i) = as<arma::mat>(mats(i));      
    }
  }
  return tmp;
}

arma::uvec get_matlen(arma::field<arma::mat> tmp){
  arma::uword n = tmp.n_elem;
  arma::uvec matlen(n);
  
  for(arma::uword i = 0; i < n; i++) {
    matlen(i) = tmp(i).n_cols;
  }
  return matlen;
}

//// [[Rcpp::export]]
//arma::cube list_to_pwm_cube(arma::field<arma::mat> tmp, arma::uvec matlen){
//
//  arma::uword maxlen = arma::max(matlen);
//  arma::uword n = tmp
//  
//  arma::cube out(n,maxlen,5, arma::fill::zeros);
//  
//  for (arma::uword i = 0; i < mats.size(); i++){
//    for (arma::uword j = 0; j < matlen(i); j++){
//      out(i,j,0) = tmp(i)(0,j);
//      out(i,j,1) = tmp(i)(1,j);
//      out(i,j,2) = tmp(i)(2,j);
//      out(i,j,3) = tmp(i)(3,j);
//    }
//  }
//  return out;
//}

//arma::vec seq_to_ix(std::string s){
//  std::vector<char> v(s.begin(), s.end());
//  arma::vec out(v.size());
//  for (size_t j; j < v.size(); j++){
//    char letter = v[j];
//    if (letter == 'A'){
//      out(j) = 0;          
//    } else if (letter =='C'){
//      out(j) = 1;
//    } else if (letter =='G'){
//      out(j) = 2;
//    } else if (letter =='T'){
//      out(j) = 3;      
//    } else{
//      out(j) = 4;
//    }
//  }
//  return(out);
//}

// [[Rcpp::export]]
arma::field<arma::mat> seq_to_mats(std::vector< std::string > x){
  
  arma::field<arma::mat> out(x.size());
  for (size_t i = 0; i < x.size(); i++){
    std::vector<char> v(x[i].begin(), x[i].end());
    arma::mat tmp(4, v.size(), arma::fill::zeros);
    for (size_t j = 0; j < v.size(); j++){
      char letter = v[j];
      if (letter == 'A'){
        tmp(0,j) = 1;          
      } else if (letter =='C'){
        tmp(1,j) = 1;
      } else if (letter =='G'){
        tmp(2,j) = 1;
      } else if (letter =='T'){
        tmp(3,j) = 1;      
      }
    }
    out(i) = tmp;
  }
  return(out);
}




// [[Rcpp::export]]
arma::umat motif_match3(const List mats, 
                            const std::vector< std::string > x,
                            const arma::vec nuc_freqs,
                            const double p){
  
  arma::field<arma::mat> pwm_field = list_to_field(mats, false);
  arma::field<arma::mat> pwm_rc_field = list_to_field(mats, true);
  
  arma::uvec matlen = get_matlen(pwm_field);
  
  arma::uword nstrings = x.size();
  arma::uword nmotifs = mats.size();
  
  arma::vec minScores(nmotifs);
  
  for(arma::uword i = 0; i < nmotifs; i++) {
    minScores(i) = p_to_score2(as<arma::mat>(mats(i)), nuc_freqs, p);
  }
  
  arma::field<arma::mat> seq_mats = seq_to_mats(x);

  arma::uword max_motif_len = arma::max(matlen);
  arma::uword min_motif_len = arma::min(matlen);
  
  std::vector<unsigned int> iloc;
  std::vector<unsigned int> jloc;
  iloc.reserve(nstrings * nmotifs);
  jloc.reserve(nstrings * nmotifs);

  for (arma::uword y=0; y < nstrings; y++){
    for (arma::uword j = 0; j < nmotifs; j++){
      int num_substr = x[y].length() - matlen(j) + 1;
      for (arma::uword i = 0; i < num_substr; i++){
        double tmp_score = arma::dot(pwm_field(j), seq_mats(y).cols(arma::span(i,i + matlen(j) - 1))); 
        double tmp_score_rc = arma::dot(pwm_rc_field(j), seq_mats(y).cols(arma::span(i,i + matlen(j) - 1))); 
        if (tmp_score > minScores(j) || tmp_score_rc > minScores(j)){
          iloc.push_back(y);
          jloc.push_back(j);
          break;
        }
      }
    }
  }
  arma::umat locs(2, iloc.size());
  locs.row(0) = arma::conv_to< arma::urowvec >::from(iloc);
  locs.row(1) = arma::conv_to< arma::urowvec >::from(jloc);
  //arma::vec values(iloc.size(), arma::fill::ones);
  //arma::sp_mat out = arma::sp_mat(locs, values, nstrings, nmotifs, true, true);
  
  return locs;
}
//
//
//
//// [[Rcpp::export]]
//NumericVector motif_match_score(const NumericMatrix mat, const std::vector< std::string > x){
//
//  size_t motif_len = mat.ncol();
//  size_t nstrings = x.size();
//
//  NumericVector out(nstrings);
//
//  for (size_t y=0; y < nstrings; y++){
//    int num_substr = x[y].length() - motif_len + 1;
//    std::vector<char> v(x[y].begin(), x[y].end());
//    double max_score = -1000;
//    for (size_t i = 0; i < num_substr; ++i){
//      double score = 0;
//      double rc_score = 0;
//      for (size_t j = 0; j < motif_len; ++j){
//        char letter = v[i+j];
//        int j_rc = motif_len - 1 - j;
//        if (letter == 'A'){
//          score += mat(0,j);  
//          rc_score += mat(3,j_rc);
//        } else if (letter =='C'){
//          score += mat(1,j);
//          rc_score += mat(2,j_rc);
//        } else if (letter =='G'){
//          score += mat(2,j);
//          rc_score += mat(1,j_rc);
//        } else if (letter =='T'){
//          score += mat(3,j);
//          rc_score += mat(0,j_rc);
//        }
//      }
//      if (score > max_score){
//        max_score = score;
//      }
//      if (rc_score > max_score){
//        max_score = score;
//      }
//    }
//    out[y] = max_score;
//  }
//  return out;
//}