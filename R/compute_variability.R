# 
# peak_deviation_helper2 <- function(normdev){
#   pvals <- pnorm(normdev)
#   pvals <- ifelse(pvals > 0.5, (1-pvals)*2, pvals*2)
#   
#   sd_normdev <- sd(normdev)
#   
#   bootstrap_indexes <- sample(seq_along(normdev), 
#                               length(normdev)*1000,
#                               replace=TRUE)
#   bootstrap_sds <- sapply(1:1000, function(x) 
#     sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):
#                                    (x*length(normdev))]]))
#   sd_error <- quantile(bootstrap_sds, c(0.025, 0.975))
#   
#   p_sd <- pchisq((length(normdev)-1) * (sd_normdev**2), 
#                  df = (length(normdev)-1), 
#                  lower.tail = FALSE)
#   
#   res <- deviationResult(deviations = normdev, 
#                          variability = sd_normdev, 
#                          variability_bounds = sd_error,
#                          p_deviations = pvals,
#                          p_variability = p_sd,
#                          metric = "z-score")
# }
# 
# peak_deviations_helper1 <- function(grp, counts_mat, niterations){
#   observed <- (counts_mat@counts[grp,] - outer(counts_mat@fragments_per_peak[grp], counts_mat@fragments_per_sample/counts_mat@total_fragments)) / 
#     outer(sqrt(counts_mat@fragments_per_peak[grp]), sqrt(counts_mat@fragments_per_sample/counts_mat@total_fragments))
#   sampled <- sapply(1:niterations, 
#                     function(x) as.matrix((counts_mat@counts[counts_mat@background_peaks[grp,x],] - 
#                                              outer(counts_mat@fragments_per_peak[counts_mat@background_peaks[grp,x]], 
#                                                    counts_mat@fragments_per_sample/counts_mat@total_fragments)) / 
#                                             outer(sqrt(counts_mat@fragments_per_peak[counts_mat@background_peaks[grp,x]]), 
#                                                   sqrt(counts_mat@fragments_per_sample/counts_mat@total_fragments))),
#                     simplify = "array")
#   z <- (observed - apply(sampled,1:2, mean)) / apply(sampled,1:2,sd)
#   
#   
# }
# 
# peak_deviations <- function(counts_mat, niterations = 50, count = NULL){
#   grpsize <- counts_mat@npeak / niterations
#   grps <- lapply(1:(counts_mat@npeak %/% grpsize + ((counts_mat@npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,counts_mat@npeak)))
#   
#   
#   observed = (counts_mat@counts - outer(counts_mat@fragments_per_peak, counts_mat@fragments_per_sample/counts_mat@total_fragments)) / 
#      outer(sqrt(counts_mat@fragments_per_peak), sqrt(counts_mat@fragments_per_sample/counts_mat@total_fragments))
#   sampled = sapply(1:niterations, function(x) as.matrix((counts_mat@counts[counts_mat@background_peaks[,x],] - 
#                                                  outer(counts_mat@fragments_per_peak[counts_mat@background_peaks[,x]], 
#                                                        counts_mat@fragments_per_sample/counts_mat@total_fragments)) / 
#                      outer(sqrt(counts_mat@fragments_per_peak[counts_mat@background_peaks[,x]]), sqrt(counts_mat@fragments_per_sample/counts_mat@total_fragments))), simplify = "array")
#   z = 
#   
#   
# }

peak_deviations <- function(counts_mat, niterations = ncol(fc@background_peaks)){
  peak_ix <- lapply(1:counts_mat@npeak, function(x) x)
  compute_variability2(peak_ix, counts_mat, niterations)
}

compute_expectations <- function(counts_mat, 
                                 by = c("all","annotation","reference"), 
                                 norm = FALSE,
                                 annotation = NULL,
                                 reference = NULL){
  
  stopifnot(inherits(counts_mat, "fragmentCounts"))
  by = match.arg(by)
  
  if (by == "all"){
    if (norm){
      expectation = rowSums(counts_mat@counts / 
                                         matrix(counts_mat@fragments_per_sample, 
                                                nrow = counts_mat@npeak, 
                                                ncol = counts_mat@nsample,
                                                byrow = TRUE))
    } else{
      expectation = counts_mat@fragments_per_peak / counts_mat@total_fragments
    } 
  } else if (by == "annotation"){
    anno = as.factor(counts_mat@sample_meta[,annotation])
    n_anno = length(levels(anno))
    mat = matrix(nrow = counts_mat@npeak, ncol = n_anno)
    if (norm){
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat@counts[,ix] / 
                            matrix(counts_mat@fragments_per_sample[ix], 
                                   nrow = counts_mat@npeak, 
                                   ncol = length(ix),
                                   byrow = TRUE))
      }
    } else{
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat@counts[,ix]) / sum(counts_mat@counts[,ix])
      }
    }
    expectation = rowMeans(mat)
  } else if (by =="reference"){
    stopifnot(!is.null(reference))
      if (!is.null(annotation)){
        ix = which(counts_mat@sample_meta[,annotation] == reference)
      } else{
        ix = reference
      }    
    if (norm){
      expectation = rowSums(counts_mat@counts[,ix] / 
                              matrix(counts_mat@fragments_per_sample[ix], 
                                     nrow = counts_mat@npeak, 
                                     ncol = length(ix),
                                     byrow = TRUE))
    } else{
      expectation = rowSums(counts_mat@counts[,ix]) / sum(counts_mat@counts[,ix])
    }
  }
  return(expectation)
}


#' compute_variability
#' 
#' Computes variability across sets of annotations
#' @param motif_indices list of indices representing different sets of peaks
#' @param counts_mat fragmentCounts object
#' @param niterations number of background sets to sample
#' @param count is input data count data?  will determine based on data if not given
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{deviationResultSet}}
#' @seealso \code{\link{deviationResultSet}}, \code{\link{variability}}, \code{\link{plot_variability}}
#' @export
compute_variability <- function(motif_indices, 
                                counts_mat, 
                                niterations = 50,
                                count = NULL){

  
  #check class of arguments to make sure they are correct
  stopifnot(inherits(motif_indices,"list"))
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  validObject(counts_mat)
  stopifnot(ncol(counts_mat@background_peaks)>= niterations)
    
  # check that indices fall within appropriate bounds
  tmp <- unlist(motif_indices, use.names =F)
  if (is.null(tmp) ||
        !(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_mat@npeak ||
        min(tmp) < 1){
    stop("motif_indices are not valid")
  }
  if (is.null(count)){
    count = all_whole(counts_mat@counts@x)
  }
  stopifnot(is.logical(count))
  
  if(is.null(names(motif_indices))){
    names(motif_indices) = 1:length(motif_indices)
  }
  # remove sets of length 0
  motif_indices <- motif_indices[which(sapply(motif_indices,length)>0)]
  
  
  if (is.installed("BiocParallel")){
      results <- deviationResultSet(BiocParallel::bplapply(motif_indices,
                                                      compute_deviations,
                                                      counts_mat, 
                                                      niterations = niterations, 
                                                      count = count))
  } else{
    results <- deviationResultSet(lapply(motif_indices,
                                         compute_deviations,
                                         counts_mat, 
                                         niterations = niterations, 
                                         count = count))
  }


  return(results)
}


#' compute_deviations
#' 
#' Computes deviations across samples for a set of annotations
#' @param peak_set vector of indices for peaks in set
#' @param counts_mat fragmentCounts object
#' @param niterations number of background sets to sample
#' @param intermediate_results should intermediate results be comptuted?
#' @param count is data count data? 
#' @return  \code{\link{deviationResult}}
#' @seealso \code{\link{deviationResult}}, \code{\link{compute_variability}}
#' @export
compute_deviations <- function(peak_set, 
                               counts_mat, 
                               niterations = 50,
                               intermediate_results = FALSE,
                               count = NULL){
  
  tf_count <- length(peak_set)
  tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                        dims = c(1, counts_mat@npeak))

  observed <- as.vector(tf_vec %*% counts_mat@counts)
  names(observed) = colnames(counts_mat@counts)

  sampled_counts <- sampleBackgroundPeaks(object = counts_mat,
                                          peak_set = peak_set,
                                          niterations = niterations)

  res <- compute_var_metrics(observed, sampled_counts, counts_mat,
                             count = count)
  
  if (intermediate_results){
    res@intermediate_results = list(observed = observed,
                                  sampled_counts = sampled_counts)
    
  }

  return(res)
}


# helper function, not exported
compute_var_metrics <- function(observed, sampled_counts, counts_mat, count = NULL){

  if (is.null(count)){
    count = all_whole(observed)
  }
  
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
                         foldchange = log2(observed / rowMeans(sampled_counts) * 
                                             rowMeans(expected_sampled_counts)/expected),
                         variability = sd_normdev, 
                         variability_bounds = sd_error,
                         p_deviations = pvals,
                         p_variability = p_sd,
                         metric = "z-score")
  return(res)
}

#' @export
compute_variability2 <- function(counts_mat, 
                                 motif_indices = NULL, 
                                 expectation = NULL,
                                 niterations = 50,
                                 count = NULL){
  
  if (is.null(motif_indices)){
    motif_indices <- lapply(1:counts_mat@npeak, function(x) x)
  }
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  }
  
  #check class of arguments to make sure they are correct
  stopifnot(inherits(motif_indices,"list"))
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  validObject(counts_mat)
  stopifnot(ncol(counts_mat@background_peaks)>= niterations)
  
  # check that indices fall within appropriate bounds
  tmp <- unlist(motif_indices, use.names =F)
  if (is.null(tmp) ||
        !(all.equal(tmp, as.integer(tmp))) ||
        max(tmp) > counts_mat@npeak ||
        min(tmp) < 1){
    stop("motif_indices are not valid")
  }
  if (is.null(count)){
    count = all_whole(counts_mat@counts@x)
  }
  stopifnot(is.logical(count))
  
  if(is.null(names(motif_indices))){
    names(motif_indices) = 1:length(motif_indices)
  }
  # remove sets of length 0
  motif_indices <- motif_indices[which(sapply(motif_indices,length)>0)]
  
  if (is.installed("BiocParallel")){
    results <- deviationResultSet(BiocParallel::bplapply(motif_indices,
                                                         compute_deviations2,
                                                         counts_mat, 
                                                         expectation,
                                                         niterations = niterations, 
                                                         count = count))
  } else{
    results <- deviationResultSet(lapply(motif_indices,
                                         compute_deviations2,
                                         counts_mat, 
                                         expectation,
                                         niterations = niterations, 
                                         count = count))
  }
  
  
  return(results)
}


compute_deviations2 <- function(peak_set, 
                               counts_mat, 
                               expectation,
                               niterations = 50,
                               intermediate_results = FALSE,
                               count = NULL){
  
  tf_count <- length(peak_set)
  
  if (tf_count == 1){
    observed <- (as.vector(counts_mat@counts[peak_set,]) - (expectation[peak_set] * counts_mat@fragments_per_sample)) /
      sqrt(expectation[peak_set] * counts_mat@fragments_per_sample)
    names(observed) <- colnames(counts_mat@counts)
    sampled_counts <-  t((as.matrix(counts_mat@counts[counts_mat@background_peaks[peak_set,],]) -
                          outer(expectation[counts_mat@background_peaks[peak_set,]],counts_mat@fragments_per_sample)) /    
      outer(sqrt(expectation[counts_mat@background_peaks[peak_set,]]),
            sqrt(counts_mat@fragments_per_sample)))
  }
  else {
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, 
                         dims = c(1, counts_mat@npeak))
  
    observed <- as.vector(tf_vec %*% ((counts_mat@counts - outer(expectation, counts_mat@fragments_per_sample)) / 
                                  outer(sqrt(expectation), sqrt(counts_mat@fragments_per_sample))))
    names(observed) <- colnames(counts_mat@counts)
    
    sampled_counts <- sampleBackgroundPeaks2(object = counts_mat,
                                          peak_set = peak_set,
                                          niterations = niterations)
  } 
  
  res <- compute_var_metrics2(observed, sampled_counts, counts_mat,
                             count = count)
  
  if (intermediate_results){
    res@intermediate_results = list(observed = observed,
                                    sampled_counts = sampled_counts)
    
  }
  
  return(res)
}

compute_var_metrics2 <- function(observed, sampled_counts, counts_mat, count = NULL){
  
  raw_deviation = observed
  sampled_deviation = sampled_counts
    
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
                         variability = sd_normdev, 
                         variability_bounds = sd_error,
                         p_deviations = pvals,
                         p_variability = p_sd,
                         metric = "z-score")
  return(res)
}


















