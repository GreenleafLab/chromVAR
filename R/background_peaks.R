#' add_gc_bias
#'
#' Computes GC content for peaks
#' @param counts_mat SummarizedExperiment
#' @param peaks GenomicRanges
#' @param genome BSgenome object
#' @export
add_gc_bias <- function(counts_mat,
                   genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  peaks <- rowRanges(counts_mat)
  seqs = Biostrings::getSeq(genome, peaks)
  nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
  gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
  rowData(counts_mat)$bias <- gc
  return(counts_mat)
}


#' get_background_peaks
#'
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat SummarizedExperiment
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in
#' GC content and # of fragments across samples using the Mahalanobis distance.
#' The w paramter controls how similar background peaks should be.  The bs parameter
#' controls the precision with which the similarity is computed; increasing bs will
#' make the function run slower. Sensible default parameters are chosen for both.
#' @import SummarizedExperiment
#' @export
get_background_peaks <- function(object, niterations = 50, w = 0.1, bs = 50){

  bias = rowData(object)$bias

  if (is.null(bias)){
    stop("bias column in rowData must be set!")
  }

  fragments_per_peak <- get_fragments_per_peak(object)
  if (min(fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")

  intensity = log10(fragments_per_peak)
  norm_mat = matrix(c(intensity, bias), ncol = 2, byrow = FALSE)

  chol_cov_mat <- chol(cov(norm_mat))
  trans_norm_mat <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))

  #make bins
  bins1 = seq(min(trans_norm_mat[,1]), max(trans_norm_mat[,1]), length.out = bs)
  bins2 = seq(min(trans_norm_mat[,2]), max(trans_norm_mat[,2]), length.out = bs)

  bin_data = do.call(rbind,lapply(1:50, function(x) matrix(c(rep(bins1[x], bs),bins2), ncol = 2, byrow = FALSE)))

  bin_dist = euc_dist(bin_data)
  bin_p = dnorm(bin_dist, 0, w)

  bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx

  bin_density <- tabulate2(bin_membership, min_val = 1, max_val = bs**2)

  background_peaks <- bg_sample_helper(bin_membership-1, bin_p, bin_density, niterations)

  return(background_peaks)
}

#' get_permuted_data
#'
#' Function to get permuted data while maintaining biases
#' @param counts_mat SummarizedExperiment
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in
#' GC content and # of fragments across samples using the Mahalanobis distance.
#' The w paramter controls how similar background peaks should be.  The bs parameter
#' controls the precision with which the similarity is computed; increasing bs will
#' make the function run slower. Sensible default parameters are chosen for both.
#' @import SummarizedExperiment
#' @export
get_permuted_data <- function(object, niterations = 10, w = 0.1, bs = 50){

  out <- BiocParallel::bplapply(1:niterations, function(x) get_permuted_data_helper(object, w, bs))
  names(out) <- paste("perm",1:niterations,sep="_") 
  
  assays(object) <- c(assays(object), out)
  
  return(object)
}


get_permuted_data_helper <- function(object, w, bs){
  bgpeaks <- get_background_peaks(object, ncol(object), w, bs)
  reorder_columns(assays(object)$counts, bgpeaks)
}

reorder_columns <- function(mat, colixmat){
  reordered <- lapply(1:ncol(mat), function(x) mat[colixmat[,x],x,drop=FALSE])
  return(do.call(cBind,reordered))
}




