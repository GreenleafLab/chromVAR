

# Filter samples based on number of reads in peaks -----------------------------

#' filter_samples
#' 
#' function to get indices of samples that pass filtters
#' @param counts_mat matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param depths vector of sequencing depth per samples, as computed by
#' \code{\link{getSampleDepths}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param plot show plot?
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated based on data.
#' min_in_peaks is set to 0.5 times the median proportion of fragments in peaks.  min_depth is 
#' set to the maximum of 500 or 10% of the median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_peaks}}
#' @export          
filter_samples <- function(counts_mat, depths, min_in_peaks = NULL, min_depth = NULL, plot = TRUE){
  stopifnot(length(depths) == ncol(counts_mat))
  stopifnot(all_true(names(depths) %in% colnames(counts_mat)))
  fragments_per_sample = colSums(counts_mat)
  if (is.null(min_in_peaks)){
    min_in_peaks = round(median(fragments_per_sample/depths)*0.5, digits = 3)
    message(paste("min_in_peaks set to ",min_in_peaks,sep="",collapse=""))
  } 
  if (is.null(min_depth)){
    min_depth = max(500,median(depths)*0.1)
    message(paste("min_depth set to ",min_depth,sep="",collapse=""))
  }
  keep_samples <- intersect(which(depths >= min_depth), 
                            which(fragments_per_sample/depths >= min_in_peaks))   
  tmp_df = data.frame(x = depths, y= fragments_per_sample/depths, z = ((1:length(fragments_per_sample)) %in% keep_samples))
  p = ggplot(tmp_df, aes_string(x="x", y="y",col="z")) + geom_point() +
    xlab("Number of fragments") + ylab("Proportion of fragments in peaks") + scale_x_log10() + annotation_logticks(sides="b") + 
    scale_y_continuous(expand =c(0,0), limits =c(0, min(1,max(tmp_df$y)*1.2)))+
    scale_color_manual(name = "Pass filters?",values = c("gray","black"),breaks = c(TRUE,FALSE), labels = c("Yes","No")) +
    chromVAR_theme()
  p = p + geom_hline(yintercept = min_in_peaks, col="red", lty = 2) + geom_vline(xintercept = min_depth, col="red", lty=2)
  if (plot) print(p)
  return(keep_samples)    
}


# Filter samples based on biases  ----------------------------------------------


bias_skew <- function(counts_mat,
                      bias,
                      nbins = 10,
                      expectation = NULL,
                      norm = TRUE){
  
  if (inherits(counts_mat,"matrix")){
    counts_mat <- Matrix(counts_mat)
  }
  stopifnot(inherits(counts_mat,"Matrix"))
  
  if (inherits(bias, "GenomicRanges")){
    bias <- mcols(bias)$gc
  }
  
  keep  <- which(rowSums(counts) > 0)
  counts_mat <- counts_mat[keep,]
  bias <- bias[keep]
  
  counts_info <- counts_summary(counts_mat)
  
  bias_quantiles = quantile(bias, seq(0,1,1/nbins))
  bias_cut = cut(bias, breaks = bias_quantiles)
  bias_bins = split(1:counts_info$npeak, bias_cut)
  
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) == nrow(counts_mat))
  }
  
  sample_names <- colnames(counts_mat)
  
  if (norm){
    if (inherits(counts_mat,"dgCMatrix")){
      counts_mat <- get_normalized_counts(counts_mat,expectation, counts_info$fragments_per_sample)
    } else{
      counts_mat <- counts_mat / outer(sqrt(expectation),sqrt(counts_info$fragments_per_sample))
    }
  }
  
  bias_mat <- sparseMatrix(j = unlist(bias_bins), 
                           i = unlist(lapply(1:nbins,function(x) rep(x, length(bias_bins[[x]])))),
                           x = 1,
                           dims = c(nbins, counts_info$npeak))
  
  observed <- as.matrix(bias_mat %*% counts_mat)
  
  if (norm){
    expected <- as.matrix(bias_mat %*% (expectation/sqrt(expectation)) %*% (counts_info$fragments_per_sample/sqrt(counts_info$fragments_per_sample)))
  } else{
    expected <- as.vector(bias_mat %*% expectation %*% counts_info$fragments_per_sample)
  }
  
  out <- observed - expected
  out <- out / matrix(sapply(bias_bins, length),nrow = nrow(out), ncol = ncol(out), byrow = FALSE)
  return(out)
}

upper_bias_limit_helper <- function(x, k){
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  return(q[2] + k * (q[2] - q[1]))
}

#' bias_filtering
#' 
#' function to identify samples that show strong biases
#' @param counts counts
#' @param bias either vector with bias values for peaks, or GenomicRanges object with column named "gc"
#' @param norm count normalization
#' @param what filter based on gc bias, enrichment bias, or both? default is gc
#' @param k parameter controlling stringency of filtering. see details
#' @details Samples are scored based on the accessibility deviations for sets of peaks with a characteristic 
#' average count or gc content.  The sum of the absolute values of those deviations
#' is computed for each sample and considered the bias score.  Samples with a bias score
#' greater than Q3 + k * (Q3 - Q1) are rejected.   
#' @return vector of indices of peaks that pass bias filter
#' @export
bias_filtering <- function(counts, bias, norm = TRUE, what = c("bias","count","both"), k = 1.5){
  
  what = match.arg(what)
  
  if (inherits(bias, "GenomicRanges")){
    bias = mcols(bias)$gc
  }
  
  if (what %in% c("both","bias")){
    gc_skew <- colSums(abs(bias_skew(counts, bias, norm = norm)),na.rm = TRUE)
    gc_pass <- which(gc_skew < upper_bias_limit_helper(gc_skew, k))              
  }
  if (what %in% c("both","count")){
    count_skew <- colSums(abs(bias_skew(counts, rowSums(counts), norm = norm)),na.rm = TRUE)
    count_pass <- which(count_skew < upper_bias_limit_helper(count_skew, k))
  }
  
  if (what == "bias"){
    return(gc_pass)
  } else if (what == "count"){
    return(count_pass)
  } else{
    return(intersect(gc_pass, count_pass))    
  }
}

# Filter peaks based on counts -------------------------------------------------

#' filter_peaks
#' 
#' function to get indices of peaks that pass filters
#' @param counts_mat matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param peaks GenomicRanges object
#' \code{\link{read_peaks}}
#' @param min_fragments_per_peak minimum number of fragmints in peaks across all samples 
#' @param non_overlapping reduce peak set to non-overlapping peaks, see details
#' @details if non_overlapping is set to true, when peaks overlap the overlapping peak with lower counts is removed
#' @return vector of indices, representing peaks that should be kept
#' @seealso \code{\link{get_peaks}},  \code{\link{get_inputs}}, \code{\link{filter_samples}},
#' \code{\link{get_counts}}
#' @export          
filter_peaks <- function(counts_mat, peaks, min_fragments_per_peak = 1, non_overlapping = TRUE){
  fragments_per_peak = rowSums(counts_mat)
  keep_peaks <- which(fragments_per_peak >= min_fragments_per_peak)
  if (non_overlapping){
    strand(peaks) <- "*"    
    if (!isTRUE(all.equal(peaks, sort(peaks)))){
      stop("peaks must be sorted to be able to filter non-overlapping peaks!")
    }
    while (!(isDisjoint(peaks[keep_peaks]))){
      chr_names = seqnames(peaks[keep_peaks])
      starts = start(peaks[keep_peaks])
      ends = end(peaks[keep_peaks])      
      overlap_next = intersect(which(chr_names[1:(length(keep_peaks) -1)] == chr_names[2:(length(keep_peaks))]),
                               which(ends[1:(length(keep_peaks) -1)] >= starts[2:(length(keep_peaks))]))
      overlap_previous = overlap_next + 1
      overlap_comparison = fragments_per_peak[keep_peaks[overlap_previous]] > fragments_per_peak[keep_peaks[overlap_next]]
      discard = keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
      keep_peaks = keep_peaks[keep_peaks %ni% discard]
    }   
  }
  return(keep_peaks)
}
