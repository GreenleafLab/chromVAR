#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double p_to_score(const NumericMatrix pwm, const NumericVector bg, const double p)
{
    const double PVAL_DP_MULTIPLIER = 1000.0;
    
    size_t a = pwm.nrow();
    size_t n = pwm.ncol();

    NumericMatrix mat(a,n);

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

    std::vector<double> table0(R + 1, 0.0);
    std::vector<double> table1(R + 1, 0.0);

    for (size_t j = 0; j < a; ++j)
        table0[mat(j,0) - minV] += bg[j];

    for (size_t i = 1; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
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
std::vector<int> motif_match(const NumericMatrix mat, 
                            const std::vector< std::string > x,
                            const NumericVector nuc_freqs,
                            const double p){
  
  double minScore = p_to_score(mat, nuc_freqs, p);

  size_t motif_len = mat.ncol();
  size_t nstrings = x.size();
  
  std::vector<int> out;
  out.reserve(nstrings);
  
  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].length() - motif_len + 1;
    std::vector<char> v(x[y].begin(), x[y].end());
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
      }
      if (score > minScore || rc_score > minScore){
        out.push_back(y+1);
        break;
      }
    }
  }
  
  return out;
}




// [[Rcpp::export]]
List motif_multi_match(List pwms, std::vector< std::string > x, NumericVector nuc_freqs, const double p){
  int npwm = pwms.size();
  CharacterVector pwm_names = pwms.names();
  List out(npwm);
  for (int i = 0; i < npwm; ++i){
    NumericMatrix pwm = as<NumericMatrix>(pwms[i]);
    out[i] = motif_match(pwm, x, nuc_freqs, p);
  }
  out.names() = pwm_names;
  return out;}




// [[Rcpp::export]]
NumericVector motif_match_score(NumericMatrix mat, std::vector< std::string > x){

  size_t motif_len = mat.ncol();
  size_t nstrings = x.size();

  NumericVector out(nstrings);

  for (size_t y=0; y < nstrings; y++){
    int num_substr = x[y].length() - motif_len + 1;
    std::vector<char> v(x[y].begin(), x[y].end());
    double max_score = -1000;
    for (size_t i = 0; i < num_substr; ++i){
      double score = 0;
      double rc_score = 0;
      for (size_t j = 0; j < motif_len; ++j){
        char letter = v[i+j];
        if (letter == 'A'){
          score += mat(0,j);  
          rc_score += mat(3,motif_len-1-j);
        } else if (letter =='C'){
          score += mat(1,j);
          rc_score += mat(1,motif_len-1-j);
        } else if (letter =='G'){
          score += mat(2,j);
          rc_score += mat(1,motif_len-1-j);
        } else if (letter =='T'){
          score += mat(3,j);
          rc_score += mat(0,motif_len-1-j);
        }
      }
      if (score > max_score){
        max_score = score;
      }
      if (rc_score > max_score){
        max_score = score;
      }
    }
    out[y] = max_score;
  }
  return out;
}