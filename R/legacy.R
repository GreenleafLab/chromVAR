### Old method for computing background, and associated functions...


# deviationBackgroundParameters class and associated methods--------------------

#' deviationBackgroundParameters
#'
#' A legacy class designed to facilitate computation of sampled background counts for a
#' set of annotations as in Buenrostro et al. (2015).
#' @slot counts_sort
#' @slot counts_sort_rev
#' @slot bias_sort
#' @slot bias_sort_rev
#' @slot window
#' @slot npeaks
#'
deviationBackgroundParameters <- setClass("deviationBackgroundParameters",
                                          slots = c(counts_sort = 'vector',
                                                    counts_sort_rev = 'vector',
                                                    bias_sort = 'vector',
                                                    bias_sort_rev = 'vector',
                                                    window = 'numeric',
                                                    npeaks = 'numeric'))

#' @export
setBackgroundParameters <- function(counts_mat, window = 2500){
  
  #sort counts
  counts = counts_mat@fragments_per_peak
  counts_sort = sort(counts, index.return=T)$ix 
  counts_sort_rev = sort(counts_sort, index.return=T)$ix
  
  #Get sorted bias
  bias_sort = sort(counts_mat@peaks$bias, index.return=T)$ix
  bias_sort_rev = sort(bias_sort, index.return=T)$ix
  
  #make object
  out = deviationBackgroundParameters(counts_sort = counts_sort,
                                      counts_sort_rev = counts_sort_rev,
                                      bias_sort = bias_sort,
                                      bias_sort_rev = bias_sort_rev,
                                      window = window,
                                      npeaks = length(counts_mat@peaks))
  return(out)
}

setMethod("show",
          signature="deviationBackgroundParameters",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("window = ", object@window, "\n", sep="")
            cat("number of peaks = ", object@npeaks, "\n", sep ="")
            invisible(NULL)
          })

#' @export
sampleBackgroundPeaks0 <- function(object, annotation_vector, counts_mat, reads, niterations = 50){
            
  stopifnot(inherits(object, "deviationBackgroundParameters"))
            
  annotation_count = sum(annotation_vector)
            
  sampling_vec = mean_smooth(annotation_vector[object@counts_sort], object@window)
  sampling_vec = sampling_vec[object@counts_sort_rev] / sum(sampling_vec)
            
  sampled_peaks = Matrix::sparseMatrix(i= rep(1:niterations, each=annotation_count),
                               j = sample(1:object@npeaks, size = annotation_count * niterations,
                                          prob = sampling_vec, replace = TRUE),
                               dims = c(niterations, object@npeaks), 
                               x = 1)
            
  bias_vec = mean_smooth(annotation_vector[object@bias_sort], object@window)
  bias_vec = bias_vec[object@bias_sort_rev] * annotation_count / sum(bias_vec)
            
  bias_peak_vec = mean_smooth(bias_vec[object@counts_sort], object@window)
  bias_peak_vec = bias_peak_vec[object@counts_sort_rev] * annotation_count / sum(bias_peak_vec)
            
  bias_term = bias_vec %*% counts_mat@counts
            
  corr_term = bias_peak_vec %*% counts_mat@counts
            
  correction = (bias_term - corr_term) * reads/sum(bias_term)
            
  sampled_counts = t(sampled_peaks %*% counts_mat@counts) + 
    matrix(correction, byrow=F, nrow = ncol(counts_mat@counts), ncol=niterations)
            
  return(as.matrix(sampled_counts))
}


# variability calculation using deviationBackgroundParameters class  -----------
#' @export
compute_variability0 <- function(motif_indices, 
                                 counts_mat, 
                                 window = 2500,
                                 niterations = 50,
                                 metric = c("z-score","old"),
                                 BPPARAM = BiocParallel::bpparam()){
  
  metric = match.arg(metric)
  
  bg_param = setBackgroundParameters(counts_mat = counts_mat, window = window)
  
  results = deviationResultSet(BiocParallel::bplapply(motif_indices, 
                                                                compute_deviations0, 
                                                                counts_mat, 
                                                                bg_param, 
                                                                niterations, 
                                                                metric, 
                                                                BPPARAM = BPPARAM))
      
  return(results)
}


#' @export
compute_deviations0 <- function(peak_set, counts_mat, bg_param, niterations = 50,
                                metric = c("z-score","old")){
  
  metric = match.arg(metric)
  
  tf_count = length(peak_set)
  tf_vec = Matrix::sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                        dims = c(1, length(counts_mat@peaks)))
  
  observed = as.matrix(tf_vec %*% counts_mat@counts)
  
  sampled_counts = sampleBackgroundPeaks0(object = bg_param, 
                                         annotation_vector = tf_vec, 
                                         counts_mat = counts_mat,
                                         reads = sum(observed), 
                                         niterations = niterations)
  
  res <- compute_var_metrics(observed, sampled_counts, counts_mat, metric = metric)
  
  return(res)
}

