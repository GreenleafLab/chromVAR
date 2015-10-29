

getBackgroundPeakSets <- function(peak_counts, bias, niterations = 51, window = 2500, BPPARAM = BiocParallel::bpparam()){
  #Standardise peak counts and bias
  peak_counts_norm = (log(peak_counts) - mean(log(peak_counts)))/sd(log(peak_counts))
  bias_norm = (bias - mean(bias)) / sd(bias)
  norm_mat = cbind(bias_norm, peak_counts_norm)
  sample_peaks <- function(indices){
    nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
    s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
    return(s)
  }
  n = length(peak_counts)
  chunks = split(1:n, cut(1:n, breaks = seq(0, n+1, 1000) ))
  sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
  return(sampled_peaks)
}


getBackgroundPeakSets2 <- function(peak_counts, bias, niterations = 51, window = 2500, BPPARAM = BiocParallel::bpparam()){
  #Standardise peak counts and bias
  peak_counts_norm = rank(peak_counts, ties.method="average")
  bias_norm = rank(bias, ties.method="average")
  norm_mat = cbind(bias_norm, peak_counts_norm)
  sample_peaks <- function(indices){
    nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
    s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
    return(s)
  }
  n = length(peak_counts)
  chunks = split(1:n, cut(1:n, breaks = seq(0, n+1, 1000) ))
  sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
  return(sampled_peaks)
}





setBackgroundParameters <- function(counts_mat, bias, window = 2500){
  
  #sort counts
  counts = rowSums(counts_mat)
  counts_sort = sort(counts, index.return=T)$ix 
  counts_sort_rev = sort(counts_sort, index.return=T)$ix
  
  #Get sorted bias
  bias_sort = sort(bias, index.return=T)$ix
  bias_sort_rev = sort(bias_sort, index.return=T)$ix
  
  #make object
  out = deviationBackgroundParameters(counts_sort = counts_sort,
                                      counts_sort_rev = counts_sort_rev,
                                      bias_sort = bias_sort,
                                      bias_sort_rev = bias_sort_rev,
                                      window = window,
                                      npeaks = length(bias))
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



setGeneric("sampleBackgroundPeaks", function(object, annotation_vector, counts_mat, niterations, ...) standardGeneric("sampleBackgroundPeaks"))

setMethod("sampleBackgroundPeaks", "deviationBackgroundParameters",
          function(object, annotation_vector, counts_mat, reads, niterations = 50){
            
            annotation_count = sum(annotation_vector)
            
            sampling_vec = mean_smooth(annotation_vector[object@counts_sort], object@window)
            sampling_vec = sampling_vec[object@counts_sort_rev] / sum(sampling_vec)
            
            sampled_peaks = sparseMatrix(i= rep(1:niterations, each=annotation_count),
                                         j = sample(1:object@npeaks, size = annotation_count * niterations, prob = sampling_vec, replace = TRUE),
                                         dims = c(niterations, object@npeaks), 
                                         x = 1)
            
            bias_vec = mean_smooth(annotation_vector[object@bias_sort], object@window)
            bias_vec = bias_vec[object@bias_sort_rev] * annotation_count / sum(bias_vec)
            
            bias_peak_vec = mean_smooth(bias_vec[object@counts_sort], object@window)
            bias_peak_vec = bias_peak_vec[object@counts_sort_rev] * annotation_count / sum(bias_peak_vec)
            
            bias_term = bias_vec %*% counts_mat
            
            corr_term = bias_peak_vec %*% counts_mat
            
            correction = (bias_term - corr_term) * reads/sum(bias_term)
            
            sampled_counts = t(sampled_peaks %*% counts_mat) + 
              matrix(correction, byrow=F, nrow = ncol(counts_mat), ncol=niterations)
            
            return(sampled_counts)
          })




