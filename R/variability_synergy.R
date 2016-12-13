
#' variability_synergy
#' 
#' @param object 
#' @param annotations 
#' @param background_peaks 
#' @param variabilities 
#' @param expectation 
#' @param nbg 
#' @details should only be run on small number of motifs/kmers/peaksets (very slow!)
#' 
#' @return synergy matrix
#'
#'@export
variability_synergy <- function(object, 
                                    annotations, 
                                    background_peaks = NULL, 
                                    variabilities = NULL,
                                    expectation = NULL,
                                    nbg = 25){
  
  
  stopifnot(inherits(object, "SummarizedExperiment"))
  
  if (is.null(assays(object)$counts)) stop("No counts slot")
  
  if (inherits(assays(object)$counts,"matrix")){
    assays(object)$counts = Matrix(counts_mat)
  }
  stopifnot(inherits(assays(object)$counts,"Matrix"))
  
  if (min(get_fragments_per_peak(object)) <= 0) stop("All peaks must have at least one fragment in one sample")
  
  if (is.null(background_peaks)){
    background_peaks <- get_background_peaks(object)
  }
  
  stopifnot(nrow(object) == nrow(background_peaks))
  
  if (inherits(annotations, "SummarizedExperiment")){
    peak_indices <- convert_to_ix_list(assays(annotations)$match)
  } else {
    stop("peak_indices must be given as SummarizedExperiment with a 'match' slot")
  }
  
  if (is.null(expectation)){
    expectation <- compute_expectations(object)
  } else{
    stopifnot(length(expectation) == nrow(object))
  }
  
  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names =F)
  if (is.null(tmp) ||
      !(all.equal(tmp, as.integer(tmp))) ||
      max(tmp) > nrow(object) ||
      min(tmp) < 1){
    stop("peak_indices are not valid")
  }
  
  if(is.null(names(peak_indices))){
    names(peak_indices) = 1:length(peak_indices)
  }
  
  if (is.null(variabilities)){
    variabilities = do.call(c,BiocParallel::bplapply(peak_indices,
                                                     compute_variability_single,
                                                     counts_mat  = assays(object)$counts,
                                                     background_peaks = background_peaks,
                                                     expectation = expectation))
  } else if (is.data.frame(variabilities)){
    variabilities = variabilities$variability
  } else{
    stop("Incorrect variabilities input")
  }
  
  stopifnot(length(variabilities) == length(peak_indices))
  
  l = length(peak_indices)
  outmat = matrix(nrow=l,ncol=l)
  var_order = order(variabilities, decreasing = TRUE)
  for (i in 1:(l-1)){
    outmat[i,i:l] <- get_variability_boost_helper(1,
                                                  assays(object)$counts,
                                                  background_peaks, 
                                                  peak_indices[var_order[i:l]],
                                                  expectation,
                                                  nbg)
    outmat[i:l,i] <- outmat[i,i:l]
  }
  colnames(outmat) = rownames(outmat) = names(peak_indices)[var_order]
  outmat
}


# Variability synergy ----------------------------------------------------------

get_variability_boost_helper <- function(index, 
                                         counts_mat, 
                                         background_peaks, 
                                         peak_indices, 
                                         expectation,
                                         nbg){
  
  tmpsets <- remove_nonoverlap(peak_indices, peak_indices[[index]])
  tmpvar <- do.call(c,BiocParallel::bplapply(tmpsets,
                                                 compute_variability_single,
                                                 counts_mat  = counts_mat,
                                                 background_peaks = background_peaks,
                                                 expectation = expectation))
  
  setlen = sapply(tmpsets,length)
  
  bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(peak_indices[[index]], x, replace=FALSE))), recursive = F)
  
  bgvar <- do.call(c,BiocParallel::bplapply(bgsets,
                                            compute_variability_single,
                                            counts_mat  = counts_mat,
                                            background_peaks = background_peaks,
                                            expectation = expectation))
    
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/ 
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
  return(var_boost)
}



