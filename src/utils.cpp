// [[Rcpp::depends(RcppArmadillo)]]
#include <type_traits>
#include <RcppArmadillo.h>
#include <cstdint>
#include <cmath>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector row_sds( arma::mat& X, bool na_rm = false) {
  arma::vec res(X.n_rows);
  if (na_rm){
    for( int j=0; j < X.n_rows; j++ ) {
      arma::rowvec tmp = X.row(j);
      arma::uvec keep = arma::find_finite(tmp);
      if (keep.size() > 1){
        res(j) = arma::stddev(tmp(keep));
      } else{
        res(j) = arma::datum::nan;
      }
    }
  } else{
    res = arma::stddev(X,0,1);
  }
  return NumericVector(res.begin(),res.end()); 
 }

// [[Rcpp::export]]
NumericVector row_sds_perm( arma::mat& X, bool na_rm = false) {
  arma::uvec tmp = arma::linspace<arma::uvec>(0, X.n_cols - 1, X.n_cols);
  arma::uvec ix = RcppArmadillo::sample(tmp, X.n_cols, true);
  arma::mat shuffled = X.cols(ix);
  return row_sds(shuffled, na_rm); 
 }
 
// Same bootstrap as row_sds_perm, from the transpose of X (cells x rows): each row is one contiguous column.
// [[Rcpp::export]]
NumericVector row_sds_perm_t( arma::mat& Xt, bool na_rm = false) {
  const arma::uword ncells = Xt.n_rows;
  arma::uvec tmp = arma::linspace<arma::uvec>(0, ncells - 1, ncells);
  arma::uvec ix = RcppArmadillo::sample(tmp, ncells, true);
  arma::vec res(Xt.n_cols);
  arma::vec buf(ncells);
  for (arma::uword j = 0; j < Xt.n_cols; j++) {
    const double* col = Xt.colptr(j);
    arma::uword m = 0;
    if (na_rm) {
      for (arma::uword c = 0; c < ncells; c++) {
        const double v = col[ix[c]];
        if (arma::arma_isfinite(v)) buf[m++] = v;
      }
      res(j) = (m > 1) ? arma::stddev(arma::vec(buf.memptr(), m, false, true)) : arma::datum::nan;
    } else {
      for (arma::uword c = 0; c < ncells; c++) buf[c] = col[ix[c]];
      res(j) = arma::stddev(buf);
    }
  }
  return NumericVector(res.begin(), res.end());
}

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
arma::umat bg_sample_helper(arma::uvec bin_membership, arma::mat bin_p, 
                            arma::vec bin_density, arma::uword niterations){
  arma::uword n = bin_membership.n_elem;
  arma::umat out = arma::umat(n, niterations);
  for (arma::uword i = 0; i < bin_density.n_elem; i++){
    arma::uvec ix = find(bin_membership == i);
	if (ix.n_elem == 0) continue;
    arma::vec p_tmp = bin_p.col(i);
    arma::vec p = p_tmp(bin_membership ) / bin_density(bin_membership );
    p /= arma::sum(p);
    arma::urowvec sampled = ProbSampleReplace(niterations * ix.n_elem, p);
    for (arma::uword j = 0; j < ix.n_elem; j++){
      out.row(ix(j)) = sampled(arma::span(j*niterations,
                 ((j+1)*niterations) - 1));
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


 
 
 
// Sums of the counts rows drawn by one peak set and by each of its background draws.
// Rows are visited in ascending peak index, adding multiplicity * value, the order of a CSC sparse product.
template <typename J, typename X, typename A = double>
static List deviations_kernel_impl(const IntegerVector& peak_set,
                       const IntegerMatrix& background_peaks,
                       const IntegerVector& csr_p,
                       const J* csr_j,
                       const X* csr_x,
                       const NumericVector& expectation,
                       int ncells) {
  const int n = expectation.size();
  const int niter = background_peaks.ncol();
  const int tf_count = peak_set.size();

  std::vector<int> member(n, 0);
  std::vector<int> sorted(peak_set.begin(), peak_set.end());
  for (int t = 0; t < tf_count; t++) member[peak_set[t] - 1]++;
  std::sort(sorted.begin(), sorted.end());

  NumericVector observed(ncells);
  double expected = 0;
  for (int t = 0; t < tf_count; ) {
    int k = sorted[t] - 1;
    int m = 1;
    while (t + m < tf_count && sorted[t + m] - 1 == k) m++;
    const double mult = m;
    for (int p = csr_p[k]; p < csr_p[k + 1]; p++) observed[csr_j[p]] += mult * (double) csr_x[p];
    expected += mult * expectation[k];
    t += m;
  }

  NumericMatrix sampled(niter, ncells);
  NumericVector sampled_expected(niter);
  NumericVector overlap(niter);

  // Bucket the (iteration, draw) pairs by drawn row with a counting sort; within a bucket iterations ascend.
  std::vector<int> start(n + 1, 0);
  for (int i = 0; i < niter; i++)
    for (int t = 0; t < tf_count; t++) start[background_peaks(peak_set[t] - 1, i)]++;  // 1-based draw d counts into slot d = (k = d - 1) + 1
  for (int k = 0; k < n; k++) start[k + 1] += start[k];
  std::vector<int> fill(start.begin(), start.end());
  std::vector<int> entry((size_t) niter * tf_count);
  for (int i = 0; i < niter; i++)
    for (int t = 0; t < tf_count; t++) entry[fill[background_peaks(peak_set[t] - 1, i) - 1]++] = i;

  // A = double accumulates into the output; an integer A accumulates in a scratch buffer, converted at the end.
  std::vector<A> acc_store;
  A* acc;
  if (std::is_same<A, double>::value) acc = reinterpret_cast<A*>(sampled.begin());
  else { acc_store.assign((size_t) niter * ncells, 0); acc = acc_store.data(); }
  std::vector<int> its;
  std::vector<A> mults;
  for (int k = 0; k < n; k++) {
    const int lo = start[k], hi = start[k + 1];
    if (lo == hi) continue;
    its.clear();
    mults.clear();
    for (int q = lo; q < hi; ) {
      int m = 1;
      while (q + m < hi && entry[q + m] == entry[q]) m++;
      its.push_back(entry[q]);
      mults.push_back((A) m);
      q += m;
    }
    const int nd = (int) its.size();
    for (int d = 0; d < nd; d++) {
      sampled_expected[its[d]] += (double) mults[d] * expectation[k];
      overlap[its[d]] += (double) mults[d] * member[k];
    }
    for (int p = csr_p[k]; p < csr_p[k + 1]; p++) {
      A* col = acc + (size_t) csr_j[p] * niter;
      const A x = (A) csr_x[p];
      for (int d = 0; d < nd; d++) col[its[d]] += mults[d] * x;
    }
  }
  if (!std::is_same<A, double>::value) {
    double* out = sampled.begin();
    const size_t tot = (size_t) niter * ncells;
    for (size_t q = 0; q < tot; q++) out[q] = (double) acc[q];
  }
  return List::create(_["observed"] = observed, _["expected"] = expected,
                      _["sampled"] = sampled,
                      _["sampled_expected"] = sampled_expected,
                      _["overlap"] = overlap);
}


// Rows of the CSR copy as 2-byte column indices and values held in raw vectors: 1-byte (uint8) when every value is an
// integer in [0, 255], else 4-byte float; NULL when the values are not all integers below 2^24 in magnitude or
// ncells > 65535, where float or uint16 would change a sum. The value vector carries the attribute "maxx".
// [[Rcpp::export]]
SEXP compact_csr(const IntegerVector& csr_j, const NumericVector& csr_x, int ncells) {
  const R_xlen_t nnz = csr_x.size();
  if (ncells > 65535 || csr_j.size() != nnz) return R_NilValue;
  double lo = 0, hi = 0;
  for (R_xlen_t p = 0; p < nnz; p++) {
    const double v = csr_x[p];
    if (!(v == std::floor(v)) || std::fabs(v) >= 16777216.0) return R_NilValue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  const bool bytes = lo >= 0 && hi <= 255.0;
  RawVector j((R_xlen_t) (nnz * sizeof(uint16_t))), x((R_xlen_t) (nnz * (bytes ? 1 : sizeof(float))));
  uint16_t* jp = reinterpret_cast<uint16_t*>(RAW(j));
  for (R_xlen_t p = 0; p < nnz; p++) jp[p] = (uint16_t) csr_j[p];
  if (bytes) {
    uint8_t* xp = RAW(x);
    for (R_xlen_t p = 0; p < nnz; p++) xp[p] = (uint8_t) csr_x[p];
  } else {
    float* xp = reinterpret_cast<float*>(RAW(x));
    for (R_xlen_t p = 0; p < nnz; p++) xp[p] = (float) csr_x[p];
  }
  x.attr("maxx") = std::max(hi, -lo);
  return List::create(_["j"] = j, _["x"] = x);
}

// deviations_kernel on the compact copy from compact_csr; same sums, bit for bit. The value vector's length says
// whether it is uint8 (nnz bytes) or float; uint8 values add into a uint32 accumulator when the largest possible
// sum, tf_count * maxx, stays below 2^32 (integer sums below 2^53 are exact in double too).
// [[Rcpp::export]]
List deviations_kernel_compact(const IntegerVector& peak_set,
                               const IntegerMatrix& background_peaks,
                               const IntegerVector& csr_p,
                               const RawVector& csr_j,
                               const RawVector& csr_x,
                               const NumericVector& expectation,
                               int ncells) {
  const uint16_t* jp = reinterpret_cast<const uint16_t*>(RAW(csr_j));
  const R_xlen_t nnz = csr_p[csr_p.size() - 1];
  if (csr_x.size() == nnz) {
    const double maxx = as<double>(csr_x.attr("maxx"));
    if ((double) peak_set.size() * maxx < 4294967296.0)
      return deviations_kernel_impl<uint16_t, uint8_t, uint32_t>(peak_set, background_peaks, csr_p, jp,
          RAW(csr_x), expectation, ncells);
    return deviations_kernel_impl<uint16_t, uint8_t>(peak_set, background_peaks, csr_p, jp,
        RAW(csr_x), expectation, ncells);
  }
  return deviations_kernel_impl<uint16_t, float>(peak_set, background_peaks, csr_p, jp,
      reinterpret_cast<const float*>(RAW(csr_x)), expectation, ncells);
}

// [[Rcpp::export]]
List deviations_kernel(const IntegerVector& peak_set,
                       const IntegerMatrix& background_peaks,
                       const IntegerVector& csr_p,
                       const IntegerVector& csr_j,
                       const NumericVector& csr_x,
                       const NumericVector& expectation,
                       int ncells) {
  return deviations_kernel_impl<int, double>(peak_set, background_peaks, csr_p,
      csr_j.begin(), csr_x.begin(), expectation, ncells);
}


// Column standard deviations of a numeric matrix, following the two-pass long double arithmetic of R's var (cov.c).
// A column holding NA or NaN gives NA_real_, as sd(x) does.
// [[Rcpp::export]]
NumericVector col_sds_r(const NumericMatrix& X) {
	const int n = X.nrow(), nc = X.ncol();
	NumericVector out(nc);
	for (int c = 0; c < nc; c++) {
		const double* xx = &X[0] + (R_xlen_t) c * n;
		bool has_na = false;
		for (int k = 0; k < n; k++) if (ISNAN(xx[k])) { has_na = true; break; }
		if (has_na || n < 2) { out[c] = NA_REAL; continue; }
		long double sum = 0.0;
		for (int k = 0; k < n; k++) sum += xx[k];
		long double tmp = sum / n;
		if (R_FINITE((double) tmp)) {
			sum = 0.0;
			for (int k = 0; k < n; k++) sum += (xx[k] - tmp);
			tmp = tmp + sum / n;
		}
		const long double xm = (long double) (double) tmp;
		sum = 0.0;
		for (int k = 0; k < n; k++) sum += (xx[k] - xm) * (xx[k] - xm);
		out[c] = std::sqrt((double) (sum / (n - 1)));
	}
	return out;
}
