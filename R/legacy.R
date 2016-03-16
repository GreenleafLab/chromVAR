# Old method for computing background, and associated functions...

# deviationBackgroundParameters class and associated methods--------------------

# A legacy class designed to facilitate computation of sampled background counts for a
# set of annotations as in Buenrostro et al. (2015).
deviationBackgroundParameters <- setClass("deviationBackgroundParameters",
                                          slots = c(counts_sort = 'vector',
                                                    counts_sort_rev = 'vector',
                                                    bias_sort = 'vector',
                                                    bias_sort_rev = 'vector',
                                                    window = 'numeric',
                                                    npeaks = 'numeric'))


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


# Variability calculation with legacy options ----------------------------------


compute_variability0 <- function(motif_indices, 
                                counts_mat, 
                                bg_method = c("new",old),
                                niterations = 50,
                                metric = c("z-score","old"),
                                count = NULL){
  
  bg_method <- match.arg(bg_method)
  metric <- match.arg(metric)
  
  #check class of arguments to make sure they are correct
  stopifnot(inherits(motif_indices,"list"))
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  validObject(counts_mat)
  
  if (bg_method == "new"){
    stopifnot(ncol(counts_mat@background_peaks) >= niterations)    
  } else if (bg_method == "old"){
    bg_param = setBackgroundParameters(counts_mat = counts_mat, window = window)
  }

  # check that indices fall within appropriate bounds
  tmp <- unlist(motif_indices, use.names =F)
  if (is.null(tmp) ||
        !(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_mat@npeak ||
        min(tmp) < 1){
    stop("motif_indices are not valid")
  }
  
  # check whether count data
  if (is.null(count)){
    count = all_whole(counts_mat@counts@x)
  }
  stopifnot(is.logical(count))
  
  # remove sets of length 0
  motif_indices <- motif_indices[which(sapply(motif_indices,length)>0)]
  
  if (is.installed("BiocParallel")){  
    if (bg_method == "new"){
        results <- deviationResultSet(BiocParallel::bplapply(motif_indices,
                                                         compute_deviations0,
                                                         counts_mat, 
                                                         niterations = niterations, 
                                                         metric = metric,
                                                         count = count))
    } else if (bg_method == "old"){
        results <- deviationResultSet(BiocParallel::bplapply(motif_indices, 
                                                          compute_deviations0, 
                                                          counts_mat, 
                                                          bg_param = bg_param, 
                                                          niterations = niterations, 
                                                          metric = metric,
                                                          count = count))
    }
  
    
  } else{
    if (bg_method == "new"){
      results <- deviationResultSet(lapply(motif_indices,
                                         compute_deviations0,
                                         counts_mat, 
                                         niterations = niterations, 
                                         metric = metric,
                                         count = count))
    } else if (bg_method =="old"){
      results <- deviationResultSet(lapply(motif_indices,
                                           compute_deviations0, 
                                           counts_mat, 
                                           bg_param = bg_param, 
                                           niterations = niterations, 
                                           metric = metric,
                                           count = count))
    }
  }
  
  
  return(results)
}


compute_deviations0 <- function(peak_set, 
                               counts_mat, 
                               bg_param = NULL,
                               niterations = 50,
                               metric = c("z-score","old"),
                               intermediate_results = FALSE,
                               count = NULL){
  
  metric <- match.arg(metric)
  
  tf_count <- length(peak_set)
  tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                         dims = c(1, counts_mat@npeak))
  
  observed <- as.vector(tf_vec %*% counts_mat@counts)
  names(observed) = colnames(counts_mat@counts)
  
  if (is.null(bg_param)){
     sampled_counts <- sampleBackgroundPeaks(object = counts_mat,
                                          peak_set = peak_set,
                                          niterations = niterations) 
  } else{
    sampled_counts = sampleBackgroundPeaks0(object = bg_param, 
                                          annotation_vector = tf_vec, 
                                          counts_mat = counts_mat,
                                          reads = sum(observed), 
                                          niterations = niterations)
  }
 
  res <- compute_var_metrics0(observed, sampled_counts, counts_mat,
                             metric = metric, count = count)
  if (intermediate_results){
    res@intermediate_results = list(observed = observed,
                                  sampled_counts = sampled_counts)    
  }
  
  return(res)
  
}


# helper function, not exported
compute_var_metrics0 <- function(observed, sampled_counts, counts_mat,
                                metric = c("z-score","old"), count = NULL){
  
  metric <- match.arg(metric)
  if (is.null(count)){
    count = all_whole(observed)
  }
  
  if (metric == 'z-score'){
    
    expected_prob <- sum(observed)/counts_mat@total_fragments
    expected_sampled_prob <- colSums(sampled_counts)/counts_mat@total_fragments
    
    expected <-  counts_mat@fragments_per_sample * expected_prob
    expected_sampled_counts <-  outer(counts_mat@fragments_per_sample,
                                      expected_sampled_prob)
    
    if (count){
      raw_deviation <- (observed - expected) / sqrt(expected)
      sampled_deviation <- (sampled_counts - expected_sampled_counts) / 
        sqrt(expected_sampled_counts)
    } else{
      raw_deviation <- observed - expected
      sampled_deviation <- sampled_counts - expected_sampled_counts
    }
    
    mean_sampled_deviation <- rowMeans(sampled_deviation)
    sd_sampled_deviation <- apply(sampled_deviation, 1, sd)
    
    normdev <- (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation
    
    pvals <- pnorm(normdev)
    pvals <- ifelse(pvals > 0.5, (1-pvals)*2, pvals*2)
    
    sd_normdev <- sd(normdev)
    
    bootstrap_indexes <- sample(seq_along(normdev), 
                                length(normdev)*1000,
                                replace=TRUE)
    bootstrap_sds <- sapply(1:1000, function(x) 
      sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                     (x*length(normdev))]]))
    sd_error <- quantile(bootstrap_sds, c(0.025, 0.975))
    
    p_sd <- pchisq((length(normdev)-1) * (sd_normdev**2), 
                   df = (length(normdev)-1), 
                   lower.tail = FALSE)
    
    res <- deviationResult(deviations = normdev, 
                           foldchange = log2(observed / rowMeans(sampled_counts) * rowMeans(expected_sampled_counts)/expected),
                           variability = sd_normdev, 
                           variability_bounds = sd_error,
                           p_deviations = pvals,
                           p_variability = p_sd,
                           metric = "z-score")
    
  } else if (metric == 'old'){
    
    expected <-  sum(observed)*counts_mat@fragments_per_sample / 
      counts_mat@total_fragments
    expected_sampled_counts <-  outer(counts_mat@fragments_per_sample / 
                                        counts_mat@total_fragments, 
                                      colSums(sampled_counts))
    
    raw_deviation <- observed - expected
    sampled_deviation <- sampled_counts-expected_sampled_counts
    rms_sampled_deviation <- apply(sampled_counts - expected_sampled_counts, 
                                   1, 
                                   rms)
    
    normdev <- raw_deviation / rms_sampled_deviation
    
    normvar_func <- function(raw_deviation, sampled_deviation){
      sqrt(sum(raw_deviation**2)/sum(rowMeans(sampled_deviation**2)))
    }
    
    normvar <- normvar_func(raw_deviation, sampled_deviation)
    
    bootstrap_indexes <- sample(seq_along(normdev),length(normdev)*1000,replace=T)
    bootstrap_normvars <- sapply(1:1000, function(x)
      normvar_func(raw_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                                     (x*length(normdev))]],
                   sampled_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                                         (x*length(normdev))],]))
    normvar_error <- quantile(bootstrap_normvars, c(0.025, 0.975))
    
    res <- deviationResult(deviations = normdev, 
                           foldchange = log2(observed / rowMeans(sampled_counts) * 
                                               rowMeans(expected_sampled_counts)/expected),
                           variability = normvar, 
                           variability_bounds = normvar_error,
                           metric = "old")
    
  } 
  return(res)
}








