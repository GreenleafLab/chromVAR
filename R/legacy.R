
compute_deviations_legacy <- function(counts_mat, 
                                       background_peaks,
                                       peak_indices = NULL, 
                                       expectation = NULL){
  
  
  stopifnot(inherits(counts_mat,"Matrix") || inherits(counts_mat,"matrix"))
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  
  if (is.null(peak_indices)){
    peak_indices <- lapply(1:counts_info$npeak, function(x) x)
  } else if (!is.list(peak_indices) && is.vector(peak_indices)){
    peak_indices = list(peak_indices)
  }
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) != nrow(counts_mat))
  }
  
  stopifnot(inherits(peak_indices,"list"))
  
  counts_info <- counts_summary(counts_mat)
  
  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names =F)
  if (is.null(tmp) ||
        !(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_info$npeak ||
        min(tmp) < 1){
    stop("peak_indices are not valid")
  }
  
  if(is.null(names(peak_indices))){
    names(peak_indices) = 1:length(peak_indices)
  }
  # remove sets of length 0
  peak_indices <- peak_indices[which(sapply(peak_indices,length)>0)]
  
  results <- BiocParallel::bplapply(peak_indices,
                                      compute_deviations_single_legacy,
                                      counts_mat, 
                                      background_peaks,
                                      expectation,
                                      counts_info)
    

  deviations <- t(vapply(results, function(x) x[["deviations"]], rep(0,counts_info$nsample)))
  variability <- sapply(results, function(x) x[["variability"]])
  results <- list(deviations = deviations, variability = variability)    
  
  return(results)
}



compute_deviations_single_legacy <- function(peak_set, 
                                      counts_mat, 
                                      background_peaks,
                                      expectation,
                                      counts_info = NULL,
                                      intermediate_results = FALSE){
    
  if (is.null(counts_info)){
    counts_info = counts_summary(counts_mat)
  }
  
  tf_count <- length(peak_set)
  
  if (tf_count == 1){
    observed <- as.vector(counts_mat[peak_set,])
    names(observed) <- colnames(counts_mat)
    sampled_counts <-  t(as.matrix(counts_mat[background_peaks[peak_set,],]))
  }
  else {
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                           dims = c(1, counts_info$npeak))
    
    observed <- as.vector(tf_vec %*% counts_mat)
    names(observed) <- colnames(counts_mat)
    
    sampled_counts <- sample_background_peaks_legacy(counts_mat = counts_mat,
                                              background_peaks = background_peaks,
                                              peak_set = peak_set,
                                              counts_info = counts_info)
  } 
  

  expected <-  counts_info$fragments_per_sample * sum(expectation[peak_set])
  expected_sampled_counts <- outer(apply(background_peaks, 2, function(x) sum(expectation[x[peak_set]])),
                                   counts_info$fragments_per_sample)
  

  observed_deviation <- observed - expected
  sampled_deviation <- sampled_counts - expected_sampled_counts
  rms_sampled_deviation <- apply(sampled_deviation, 2, rms)
  normdev <- observed_deviation / rms_sampled_deviation
  normvar <- sqrt(sum(observed_deviation**2)/mean(rowSums(sampled_deviation**2)))
  if (intermediate_results){
      out = list(deviations = normdev, variability = normvar, observed = observed,
                 sampled = sampled, expected = expected, expected_sampled = expected_sampled_counts)
  } else{
    out = list(deviations = normdev, variability = normvar)
  }
  return(out)  
}


sample_background_peaks_legacy <- function(counts_mat, peak_set, background_peaks, counts_info){
  
  niterations = ncol(background_peaks)
  sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]), 
                            i = rep(1:niterations, each = length(peak_set)), 
                            x=1, 
                            dims = c(ncol(background_peaks), counts_info$npeak))
  
  sampled_counts =  as.matrix(sample_mat %*% counts_mat)
  
  return(sampled_counts)
}



