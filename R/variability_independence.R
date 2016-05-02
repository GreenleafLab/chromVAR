



get_independent_variability_helper <- function(index, 
                                               counts_mat, 
                                               background_peaks, 
                                               peak_indices, 
                                               expectation,
                                               norm,
                                               counts_info){
    
  tmpsets <- remove_overlap(peak_indices, peak_indices[[index]])
  tmpvar <- do.call(c,BiocParallel::bplapply(tmpsets,
                                             compute_variability_single,
                                             counts_mat  = counts_mat,
                                             background_peaks = background_peaks,
                                             expectation = expectation,
                                             counts_info = counts_info,
                                             norm = norm))

  return(tmpvar)
}


get_independent_variability <-  function(index, 
                                         counts_mat, 
                                         background_peaks, 
                                         peak_indices, 
                                         expectation = NULL,
                                         norm = TRUE){
  
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
  
  indep_var <- get_independent_variability_helper(index,
                                            counts_mat,
                                            background_peaks, 
                                            peak_indices,
                                            expectation,
                                            norm,
                                            counts_info)
  
  return(indep_var)
}


#' @export
get_top_sets <- function(counts_mat, 
                         background_peaks, 
                         peak_indices, 
                         variabilities = NULL,
                         expectation = NULL,
                         norm = TRUE,
                         max_iter = NULL,
                         var_threshold = NULL){
  
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
  
  
  if (is.null(var_threshold)){
    var_threshold = max(sqrt(qchisq(0.01/length(peak_indices), 
                               df = (ncol(counts)-1), 
                               lower.tail = FALSE) / (ncol(counts) - 1)),min(variabilities, na.rm=TRUE))
    message(paste("var_threshold not provided, set to ",var_threshold,sep="",collapse=""))
  }
  
  if (is.null(max_iter)){
    max_iter = length(variabilities)
  }
  
  #get max variable
  max_var <- which(variabilities == max(variabilities, na.rm=TRUE))

  #initialize output
  out <- c(max_var)
  
  #Iterate...
  iter <- 2
  candidates <- which(variabilities >= var_threshold)
  while( length(candidates) > 1 && iter <= max_iter){
    message(paste("Iteration: ",iter,", remaining candidates: ",length(candidates) - 1,sep="",collapse=""))
    tmpvar <- get_independent_variability_helper(which(candidates == max_var),
                                                 counts_mat,
                                                 background_peaks, 
                                                 peak_indices[candidates],
                                                 expectation,
                                                 norm,
                                                 counts_info)
    if (max(tmpvar, na.rm=TRUE) < var_threshold){
      break
    } else{
      candidates <- candidates[which(tmpvar >= var_threshold)]      
      max_var <- candidates[which.max(variabilities[candidates])]
      out <- c(out, max_var)
      iter <- iter + 1
    }
  }
  return(out)
}


get_unique_pwms <- function(pwms, variability, similarity_cutoff = 1){
  ds = chromVAR:::compute_pwm_dist(pwms = lapply(pwms, as.matrix))
  out = c()
  candidates = order(variability$variability, decreasing = TRUE)
  while (length(candidates) > 0){
    out = c(out, candidates[1])
    candidates = candidates[which(ds$dist[candidates[1],candidates] > 1)]
  }
  return(out)
}



    