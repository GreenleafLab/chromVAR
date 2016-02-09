

#' compute_variability
#' 
#' Computes variability across sets of annotations
#' @param motif_indices list of indices representing different sets of peaks
#' @param counts_mat fragmentCounts object
#' @param niterations number of background sets to sample
#' @param metric which metric to use? default is z-score
#' @param BPPARAM multi-processing argument to pass to \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{deviationResultSet}}
#' @seealso \code{\link{deviationResultSet}}, \code{\link{variability}}, \code{\link{plot_variability}}
#' @export
compute_variability <- function(motif_indices, 
                                counts_mat, 
                                niterations = 50,
                                metric = c("z-score","old"),
                                BPPARAM = BiocParallel::bpparam()){

  metric <- match.arg(metric)
  
  #check class of arguments to make sure they are correct
  stopifnot(inherits(motif_indices,"list"))
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  validObject(counts_mat)
  stopifnot(ncol(counts_mat@background_peaks)>= niterations)
    
  # check that indices fall within appropriate bounds
  tmp <- unlist(motif_indices, use.names =F)
  if (!(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_mat@npeak ||
        min(tmp) < 1){
    stop("motif_indices are not valid")
  }
  
  # remove sets of length 0
  motif_indices <- motif_indices[which(sapply(motif_indices,length)>0)]
  
  results <- deviationResultSet(BiocParallel::bplapply(motif_indices,
                                                      compute_deviations,
                                                      counts_mat, 
                                                      niterations, 
                                                      metric,
                                                      BPPARAM = BPPARAM))

  return(results)
}


#' compute_deviations
#' 
#' Computes deviations across samples for a set of annotations
#' @param peak_set vector of indices for peaks in set
#' @param counts_mat fragmentCounts object
#' @param niterations number of background sets to sample
#' @param metric which metric to use?  default is z-score
#' @return  \code{\link{deviationResult}}
#' @seealso \code{\link{deviationResult}}, \code{\link{compute_variability}}
#' @export
compute_deviations <- function(peak_set, 
                               counts_mat, 
                               niterations = 50,
                               metric = c("z-score","old"),
                               intermediate_results = FALSE){

  metric <- match.arg(metric)

  tf_count <- length(peak_set)
  tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                        dims = c(1, counts_mat@npeak))

  observed <- as.matrix(tf_vec %*% counts_mat@counts)

  sampled_counts <- sampleBackgroundPeaks(object = counts_mat,
                                          peak_set = peak_set,
                                          niterations = niterations)

  res <- compute_var_metrics(observed, sampled_counts, counts_mat,
                             metric = metric)
  
  res@intermediate_results = list(observed = observed,
                                  sampled_counts = sampled_counts)

  return(res)

}


# helper function, not exported
compute_var_metrics <- function(observed, sampled_counts, counts_mat,
                                metric = c("z-score","old")){

  metric <- match.arg(metric)

  if (metric == 'z-score'){

    expected_prob <- sum(observed)/counts_mat@total_fragments
    expected_sampled_prob <- colSums(sampled_counts)/counts_mat@total_fragments

    expected <-  counts_mat@fragments_per_sample * expected_prob
    expected_sampled_counts <-  outer(counts_mat@fragments_per_sample,
                                     expected_sampled_prob)

    raw_deviation <- (observed - expected) / sqrt(expected * (1 - expected_prob))
    sampled_deviation <- (sampled_counts - expected_sampled_counts) / 
      sqrt(expected_sampled_counts * 
             (1 - matrix(expected_sampled_prob, nrow = nrow(sampled_counts),
                         ncol = ncol(sampled_counts),byrow = T)))

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
    
    res <- deviationResult(deviations = as.numeric(normdev), 
                          variability = sd_normdev, 
                          variability_bounds = sd_error,
                          p_deviations = as.numeric(pvals),
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

    res <- deviationResult(deviations = as.numeric(normdev), 
                          variability = normvar, 
                          variability_bounds = normvar_error,
                          metric = "old")

  } 
  return(res)
}
