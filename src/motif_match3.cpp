// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double p_to_score(const arma::mat pwm, const arma::vec bg, const double p)
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
arma::mat dynamic_thresholds(const arma::mat mat, double threshold){
  arma::uword motif_len = mat.n_cols;
  arma::mat out(4, motif_len);
  out(0,motif_len - 1) = threshold;
  out(1,motif_len - 1) = threshold;
  out(2,motif_len - 1) = threshold;
  out(3,motif_len - 1) = threshold;
  for (arma::sword i = motif_len - 2; i >= 0; i--){
    out(0, i ) = out(0, i + 1) - arma::max(mat.col(i + 1));
    out(1, i ) = out(1, i + 1) - arma::min(mat.col(i + 1));
    out(2, i ) = out(2, i + 1) - arma::max(mat.col(motif_len - i - 2));
    out(3, i ) = out(3, i + 1) - arma::min(mat.col(motif_len - i - 2));    
  }
  return out;
}

// [[Rcpp::export]]
arma::vec seq_to_ix(std::string s){
  std::vector<char> v(s.begin(), s.end());
  arma::vec out(v.size());
  for (size_t j = 0; j < v.size(); j++){
    char letter = v[j];
    if (letter == 'A'){
      out(j) = 0;          
    } else if (letter =='C'){
      out(j) = 1;
    } else if (letter =='G'){
      out(j) = 2;
    } else if (letter =='T'){
      out(j) = 3;      
    } else{
      out(j) = 4;
    }
  }
  return(out);
}


// [[Rcpp::export]]
std::vector<int> motif_match2(const arma::mat mat, 
                            const std::vector< std::string > x,
                            const arma::vec nuc_freqs,
                            const double p){
  
  double minScore = p_to_score2(mat, nuc_freqs, p);
  arma::mat control_mat = dynamic_thresholds(mat, minScore);

  size_t motif_len = mat.n_cols;
  size_t nstrings = x.size();
  
  std::vector<int> out;
  out.reserve(nstrings);
  
  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].length() - motif_len + 1;
    std::vector<char> v(x[y].begin(), x[y].end());
    bool matched = false;
    for (size_t i = 0; i < num_substr; ++i){
      double score = 0;
      double rc_score = 0;
      for (size_t j = 0; j < motif_len; ++j){
        char letter = v[i+j];
        int j_rc = motif_len - 1 - j;
        if (letter == 'A'){
          score += mat(0,j);  
          rc_score += mat(3,j_rc);
        } else if (letter =='C'){
          score += mat(1,j);
          rc_score += mat(2,j_rc);
        } else if (letter =='G'){
          score += mat(2,j);
          rc_score += mat(1,j_rc);
        } else if (letter =='T'){
          score += mat(3,j);
          rc_score += mat(0,j_rc);
        }
        if (score >= control_mat(1,j) || rc_score >= control_mat(3,j)){
          out.push_back(y+1);
          matched = true;
          break;
        } else if (score < control_mat(0,j) && rc_score < control_mat(2,j)){
          break;
        }         
      }
      if (matched) break;
    }
  }  
  return out;
}

// [[Rcpp::export]]
std::vector<int> motif_match4(const arma::mat mat, 
                            const std::vector< std::string > x,
                            const arma::vec nuc_freqs,
                            const double p){
  
  double minScore = p_to_score2(mat, nuc_freqs, p);
  arma::mat control_mat = dynamic_thresholds(mat, minScore);
  
  size_t motif_len = mat.n_cols;
  size_t nstrings = x.size();  
  
  arma::mat rc_mat = arma::fliplr(arma::flipud(mat));
  arma::mat pad(1, motif_len, arma::fill::zeros);
  arma::mat pad_mat = arma::join_vert(mat, pad);
  arma::mat pad_rc_mat = arma::join_vert(rc_mat, pad);

  std::vector<int> out;
  out.reserve(nstrings);
  
  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].length() - motif_len + 1;
    arma::vec seq_ix = seq_to_ix(x[y]);
    bool matched = false;
    for (size_t i = 0; i < num_substr; ++i){
      double score = 0;
      double rc_score = 0;
      for (size_t j = 0; j < motif_len; ++j){
        score += pad_mat(seq_ix(i+j),j);
        rc_score += pad_rc_mat(seq_ix(i+j),j);
        if (score >= control_mat(1,j) || rc_score >= control_mat(3,j)){
          out.push_back(y+1);
          matched = true;
          break;
        } else if (score < control_mat(0,j) && rc_score < control_mat(2,j)){
          break;
        }         
      }
      if (matched) break;
    }
  }  
  return out;
}

// [[Rcpp::export]]
std::vector< std::vector<int> > seq_to_ix_list(std::vector<std::string> x){
  std::vector< std::vector<int> > out(x.size());
  for (size_t i = 0; i < x.size(); i++){
    std::vector<char> v(x[i].begin(), x[i].end());
    std::vector<int> tmp(v.size());
    for (size_t j = 0; j < v.size(); j++){
      char letter = v[j];
      if (letter == 'A'){
        tmp[j] = 0;          
      } else if (letter =='C'){
        tmp[j] = 1;
      } else if (letter =='G'){
        tmp[j] = 2;
      } else if (letter =='T'){
        tmp[j] = 3;      
      } else{
        tmp[j] = 4;
      }
    }
    out[i] = tmp;
  }
  return(out);
}

// [[Rcpp::export]]
std::vector<int> motif_match5(const arma::mat mat, 
                            const std::vector< std::vector<int> > x,
                            const arma::vec nuc_freqs,
                            const double p){
  
  double minScore = p_to_score2(mat, nuc_freqs, p);
  arma::mat control_mat = dynamic_thresholds(mat, minScore);
  
  size_t motif_len = mat.n_cols;
  size_t nstrings = x.size();  
  
  arma::mat rc_mat = arma::fliplr(arma::flipud(mat));
  arma::mat pad(1, motif_len, arma::fill::zeros);
  arma::mat pad_mat = arma::join_vert(mat, pad);
  arma::mat pad_rc_mat = arma::join_vert(rc_mat, pad);

  std::vector<int> out;
  out.reserve(nstrings);
  
  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].size() - motif_len + 1;
    bool matched = false;
    for (size_t i = 0; i < num_substr; ++i){
      double score = 0;
      double rc_score = 0;
      for (size_t j = 0; j < motif_len; ++j){
        score += pad_mat(x[y][i+j],j);
        rc_score += pad_rc_mat(x[y][i+j],j);
        if (score >= control_mat(1,j) || rc_score >= control_mat(3,j)){
          out.push_back(y+1);
          matched = true;
          break;
        } else if (score < control_mat(0,j) && rc_score < control_mat(2,j)){
          break;
        }         
      }
      if (matched) break;
    }
  }  
  return out;
}

// [[Rcpp::export]]
std::vector< std::vector<double> > motif_match_score2(const arma::mat mat, 
                            const std::vector< std::vector<int> > x,
                            const arma::vec nuc_freqs,
                            const double p){

  double minScore = p_to_score2(mat, nuc_freqs, p);
  arma::mat control_mat = dynamic_thresholds(mat, minScore);
  
  size_t motif_len = mat.n_cols;
  size_t nstrings = x.size();  
  
  arma::mat rc_mat = arma::fliplr(arma::flipud(mat));
  arma::mat pad(1, motif_len, arma::fill::zeros);
  arma::mat pad_mat = arma::join_vert(mat, pad);
  arma::mat pad_rc_mat = arma::join_vert(rc_mat, pad);

  std::vector<double> out1;
  out1.reserve(nstrings);
  
  std::vector<double> out2;
  out2.reserve(nstrings);
  
  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].size() - motif_len + 1;
    bool matched = false;
    double max_score = minScore - 1;
    for (size_t i = 0; i < num_substr; ++i){
      double score = 0;
      double rc_score = 0;
      for (size_t j = 0; j < motif_len; ++j){
        score += pad_mat(x[y][i+j],j);
        rc_score += pad_rc_mat(x[y][i+j],j);
        if (score < control_mat(0,j) && rc_score < control_mat(2,j)){
          break;
        }         
      }
      if (score > max_score){
        max_score = score;
      }
      if (rc_score > max_score){
        max_score = rc_score;
      }
    }
    if (max_score >= minScore){
      out1.push_back(y);
      out2.push_back(max_score);
    }
  }
  std::vector< std::vector<double> > out(2);
  out[0] = out1;
  out[1] = out2;
  return out;
}  




 