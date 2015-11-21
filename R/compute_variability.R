
remove_overlap <- function(indices_list, index){
  toremove = indices_list[[index]]
  indices_list = indices_list[-index]
  indices_list = lapply(indices_list, function(x) x[which(x %ni% toremove)])
  return(indices_list)
} 

remove_nonoverlap <- function(indices_list, index){
  if (is.character(index)){
    index = which(names(indices_list) == index)
  }
  tokeep = indices_list[[index]]
  indices_list = indices_list[-index]
  indices_list = lapply(indices_list, function(x) x[which(x %in% tokeep)])
  return(indices_list)
} 


get_top_sets <- function(results, sets, counts_mat, bg_peaks, 
                         niterations = 50,
                         p_cutoff = 0.01,
                         max_iter = 25,
                         BPPARAM = BiocParallel::bpparam()){
  
  #Get only significant sets
  results <- subset_by_variability(results, cutoff = p_cutoff, adjusted = TRUE)
  sets <- sets[names(results)]
  p <- get_pvalues(results, adjust = TRUE)
  
  #get max variable
  max_var <- which(variability(results) == max(variability(results)))
  min_p <- p[max_var]
  
  #initialize output
  out <- results[max_var]

  #Iterate...
  iter = 2
  while( min_p < p_cutoff && iter < max_iter){
    sets <- remove_overlap(sets, max_var)
    tmpresults <- compute_variability(sets, counts_mat, bg_peaks, niterations = niterations,
                                      BPPARAM = BPPARAM)
    tmpresults <- set_nresult(tmpresults, results@nresult)
    # get most variable...
    max_var <- which(variability(tmpresults) == max(variability(tmpresults)))
    max_var_name <-  names(tmpresults)[max_var]
    min_p <- min(get_pvalues(tmpresults, adjust = TRUE))
    
    out <- c(out, tmpresults[max_var])
    
    iter <- iter + 1
  }

  out <- set_nresult(out, results@nresult)
  
  return(out)
}



peak_deviations <- function(counts_mat, window = 500, BPPARAM = BiocParallel::bpparam()){
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  #if bias not available...
  if ("bias" %ni% colnames(S4Vectors::mcols(counts_mat@peaks))){
    stop("Peaks must have metadata column named 'bias'. Compute using compute_bias.")
  }
  
  expected_probs = counts_mat@fragments_per_peak / counts_mat@total_fragments
  
  expected = outer(expected_probs,
                   counts_mat@fragments_per_sample)
  
  dev = as.matrix(counts_mat@counts - expected) / sqrt(matrix( 1- expected_probs,
                                                         nrow = counts_mat@npeak,
                                                         ncol = counts_mat@nsample,
                                                         byrow=F) * 
                                                     expected)
  

  norm_mat <- cbind(log10(counts_mat@fragments_per_peak + 1 ), counts_mat@peaks$bias)
  reflected <- add_reflections(norm_mat, window = window)
  cov_mat = cov(norm_mat)
  
  norm_dev <- function(index){
    if (!all_true(!is.na(dev[index,]))){ return(rep(NA,ncol(dev)))}
    
    mdist <- mahalanobis(reflected$data, norm_mat[index,], cov = cov_mat)
    mdist[reflected$identify(index)] = Inf
    closest <- which(mdist <= quantile(mdist, window/nrow(norm_mat), type=1))
    closest <- reflected$replace(closest)
    out <- (dev[index,] - apply(dev[closest,],2,mean)) / apply(dev[closest,],2,sd)
    
    return(out)
  }
  
  normed = t(simplify2array(BiocParallel::bplapply(1:nrow(dev), norm_dev, BPPARAM = BPPARAM)))
  return(normed)
  
}




#' compute_variability
#' 
#' Computes variability across sets of annotations
#' @param motif_indices list of indices representing different sets of peaks
#' @param counts_mat fragmentCounts object
#' @param bg_peaks backgroundPeaks object
#' @param niterations number of background sets to sample
#' @param metric which metric to use? default is z-score
#' @param BPPARAM multi-processing argument to pass to \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{deviationResultSet}}
#' @seealso \code{\link{deviationResultSet}}, \code{\link{variability}}, \code{\link{plot_variability}}
#' @export
compute_variability <- function(motif_indices, counts_mat, bg_peaks, 
                                niterations = 50,
                                metric = c("z-score","old"),
                                BPPARAM = BiocParallel::bpparam()){

  metric = match.arg(metric)
  
  #check class of arguments to make sure they are correct
  stopifnot(inherits(motif_indices,"list"))
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  stopifnot(inherits(bg_peaks,"backgroundPeaks"))
  
  #check that backgroundPeaks based on same peaks as fragmentCounts
  stopifnot(all.equal(bg_peaks@peaks, fc@peaks))
  
  #check that indices fall within appropriate bounds
  tmp = unlist(motif_indices, use.names =F)
  if (!(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_mat@npeak ||
        min(tmp) < 1){
    stop("motif_indices are not valid")
  }
  
  results = deviationResultSet(BiocParallel::bplapply(motif_indices,
                                                      compute_deviations,
                                                      counts_mat, bg_peaks,
                                                      niterations, metric,
                                                      BPPARAM = BPPARAM))

  return(results)
}


#' compute_deviations
#' 
#' Computes deviations across samples for a set of annotations
#' @param peak_set vector of indices for peaks in set
#' @param counts_mat fragmentCounts object
#' @param bg_peaks backgroundPeaks object
#' @param niterations number of background sets to sample
#' @param metric which metric to use?  default is z-score
#' @return  \code{\link{deviationResult}}
#' @seealso \code{\link{deviationResult}}, \code{\link{compute_variability}}
#' @export
compute_deviations <- function(peak_set, counts_mat, bg_peaks, niterations = 50,
                               metric = c("z-score","old")){

  metric = match.arg(metric)

  tf_count = length(peak_set)
  tf_vec = sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                        dims = c(1, length(counts_mat@peaks)))

  observed = as.matrix(tf_vec %*% counts_mat@counts)

  sampled_counts = sampleBackgroundPeaks(object = bg_peaks, peak_set = peak_set, 
                                         counts_mat = counts_mat,
                                                  niterations = niterations)

  res <- compute_var_metrics(observed, sampled_counts, counts_mat,
                             metric = metric)

  return(res)

}


# not exported
compute_var_metrics <- function(observed, sampled_counts, counts_mat,
                                metric = c("z-score","old")){

  metric = match.arg(metric)

  if (metric == 'z-score'){

    expected_prob = sum(observed)/counts_mat@total_fragments
    expected_sampled_prob = colSums(sampled_counts)/counts_mat@total_fragments

    expected =  counts_mat@fragments_per_sample * expected_prob
    expected_sampled_counts =  outer(counts_mat@fragments_per_sample,
                                     expected_sampled_prob)

    raw_deviation = (observed - expected) / sqrt(expected * (1 - expected_prob))
    sampled_deviation = (sampled_counts - expected_sampled_counts) / 
      sqrt(expected_sampled_counts * 
             (1 - matrix(expected_sampled_prob, nrow = nrow(sampled_counts),
                         ncol = ncol(sampled_counts),byrow = T)))

    mean_sampled_deviation = rowMeans(sampled_deviation)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation
    
    pvals = pnorm(normdev)
    pvals = ifelse(pvals > 0.5, (1-pvals)*2, pvals*2)

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(seq_along(normdev), 
                               length(normdev)*1000,
                               replace=TRUE)
    bootstrap_sds = sapply(1:1000, function(x) 
      sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                     (x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    p_sd = pchisq((length(normdev)-1) * (sd_normdev**2), 
                  df = (length(normdev)-1), 
                  lower.tail = FALSE)
    
    res = deviationResult(deviations = as.numeric(normdev), 
                          variability = sd_normdev, 
                          variability_bounds = sd_error,
                          p_deviations = as.numeric(pvals),
                          p_variability = p_sd,
                          metric = "z-score")

  } else if (metric == 'old'){

    expected =  sum(observed)*counts_mat@fragments_per_sample / 
      counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_sample / 
                                       counts_mat@total_fragments, 
                                     colSums(sampled_counts))

    raw_deviation = observed - expected
    sampled_deviation = sampled_counts-expected_sampled_counts
    rms_sampled_deviation = apply(sampled_counts - expected_sampled_counts, 
                                  1, 
                                  rms)

    normdev = raw_deviation / rms_sampled_deviation

    normvar_func <- function(raw_deviation, sampled_deviation){
      sqrt(sum(raw_deviation**2)/sum(rowMeans(sampled_deviation**2)))
    }

    normvar = normvar_func(raw_deviation, sampled_deviation)

    bootstrap_indexes = sample(seq_along(normdev),length(normdev)*1000,replace=T)
    bootstrap_normvars = sapply(1:1000, function(x)
      normvar_func(raw_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                                     (x*length(normdev))]],
                   sampled_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):
                                                         (x*length(normdev))],]))
    normvar_error = quantile(bootstrap_normvars, c(0.025, 0.975))

    res = deviationResult(deviations = as.numeric(normdev), 
                          variability = normvar, 
                          variability_bounds = normvar_error,
                          metric = "old")

  } 
  return(res)
}
