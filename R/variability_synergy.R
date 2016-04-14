# Variability synergy ----------------------------------------------------------


get_variability_boost_helper <- function(index, 
                                         counts_mat, 
                                         background_peaks, 
                                         peak_indices, 
                                         expectation,
                                         norm,
                                         counts_info,
                                         nbg){
  
  tmpsets <- remove_nonoverlap(peak_indices, peak_indices[[index]])
  tmpvar <- do.call(c,BiocParallel::bplapply(tmpsets,
                                                 compute_variability_single,
                                                 counts_mat  = counts_mat,
                                                 background_peaks = background_peaks,
                                                 expectation = expectation,
                                                 counts_info = counts_info,
                                                 norm = norm))
  
  setlen = sapply(tmpsets,length)
  
  bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(peak_indices[[index]], x, replace=FALSE))), recursive = F)
  
  bgvar <- do.call(c,BiocParallel::bplapply(bgsets,
                                                compute_variability_single,
                                                counts_mat  = counts_mat,
                                                background_peaks = background_peaks,
                                                expectation = expectation,
                                                counts_info = counts_info,
                                                norm = norm))
    
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/ 
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
  return(var_boost)
}

get_variability_boost <- function(index, 
                                  counts_mat, 
                                  background_peaks, 
                                  peak_indices, 
                                  expectation = NULL,
                                  norm = TRUE,
                                  nbg = 50){
  
  if (inherits(counts_mat,"matrix")){
    counts_mat = Matrix(counts_mat)
  }
  stopifnot(inherits(counts_mat,"Matrix"))
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  
  counts_info <- counts_summary(counts_mat)
  if (min(counts_info$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  
  if (is.null(peak_indices)){
    peak_indices <- lapply(1:counts_info$npeak, function(x) x)
  } else if (inherits(peak_indices, "Matrix") || inherits(peak_indices, "matrix")){
    peak_indices <- convert_to_ix_list(peak_indices)
  } else if (!is.list(peak_indices) && is.vector(peak_indices)){
    peak_indices = list(peak_indices)
  }
  stopifnot(inherits(peak_indices,"list"))  
  
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) == nrow(counts_mat))
  }
  
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
  
  if (norm){
    if (inherits(counts_mat,"dgCMatrix")){
      counts_mat <- get_normalized_counts(counts_mat,expectation, counts_info$fragments_per_sample)
    } else{
      counts_mat <- counts_mat / outer(sqrt(expectation),sqrt(counts_info$fragments_per_sample))
    }
  }
  
  var_boost <- get_variability_boost_helper(index,
                                            counts_mat,
                                            background_peaks, 
                                            peak_indices,
                                            expectation,
                                            norm,
                                            counts_info,
                                            nbg)
  
  return(var_boost)
}


get_variability_synergy <- function(counts_mat, 
                                    background_peaks, 
                                    peak_indices, 
                                    variabilities = NULL,
                                    expectation = NULL,
                                    norm = TRUE,
                                    nbg = 50){
  
  if (inherits(counts_mat,"matrix")){
    counts_mat = Matrix(counts_mat)
  }
  stopifnot(inherits(counts_mat,"Matrix"))
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  
  counts_info <- counts_summary(counts_mat)
  if (min(counts_info$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  
  if (inherits(peak_indices, "Matrix") || inherits(peak_indices, "matrix")){
    peak_indices <- convert_to_ix_list(peak_indices)
  } else if (!is.list(peak_indices) && is.vector(peak_indices)){
    peak_indices = list(peak_indices)
  }
  stopifnot(inherits(peak_indices,"list"))  
  
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) == nrow(counts_mat))
  }
  
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
  
  if (norm){
    if (inherits(counts_mat,"dgCMatrix")){
      counts_mat <- get_normalized_counts(counts_mat,expectation, counts_info$fragments_per_sample)
    } else{
      counts_mat <- counts_mat / outer(sqrt(expectation),sqrt(counts_info$fragments_per_sample))
    }
  }
  
  if (is.null(variabilities)){
    variabilities = do.call(c,BiocParallel::bplapply(peak_indices,
                                                     compute_variability_single,
                                                     counts_mat  = counts_mat,
                                                     background_peaks = background_peaks,
                                                     expectation = expectation,
                                                     counts_info = counts_info,
                                                     norm = norm))
  } else if (is.data.frame(variabilities)){
    variabilities = variabilities$variability
  }
  
  stopifnot(length(variabilities) == length(peak_indices))
  
  
  l = length(peak_indices)
  outmat = matrix(nrow=l,ncol=l)
  var_order = order(variabilities, decreasing = TRUE)
  for (i in 1:(l-1)){
    outmat[i,i:l] <- get_variability_boost_helper(1,
                                                      counts_mat,
                                                      background_peaks, 
                                                      peak_indices[var_order[i:l]],
                                                      expectation,
                                                      norm,
                                                      counts_info,
                                                      nbg)
  }
  colnames(outmat) = rownames(outmat) = names(peak_indices)[var_order]
  outmat
}

