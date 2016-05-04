
#' compute_expectations
#' 
#' @param counts_mat matrix of counts
#' @param norm weight all samples equally?
#' @param annotation an annotation vector, optional
#' @details By default, this function will compute the expected fraction of reads
#' per peak as the the total fragments per peak across all samples divided by total
#' reads in peaks in all samples. Optionally, norm can be set to TRUE and then the 
#' expectation will be the average fraction of reads in a peak across the cells.  
#' This is not recommended for single cell applications as cells with very few reads
#' will have a large impact.  Another option is to give a vector of annotations, in which 
#' case the expectation will be the average fraction of reads per peak within each annotation.
#' If annotation vector is provided and norm is set to TRUE then within each annotation
#' the fraction of reads per peak is the average fraction of reads per peak in each 
#' sample.  Otherwise, the within annotation fraction of reads per peak is based on the
#'reads per peak within the sample divided by the total reads within each sample.
#' @return vector with expected fraction of reads per peak
#' @export
compute_expectations <- function(counts_mat,
                                 norm = FALSE,
                                 annotation = NULL){

  counts_info = counts_summary(counts_mat)

  if (is.null(annotation)){
    if (norm){
      expectation = rowSums(counts_mat /
                              matrix(counts_info$fragments_per_sample,
                                     nrow = counts_info$npeak,
                                     ncol = counts_info$nsample,
                                     byrow = TRUE))
    } else{
      expectation = counts_info$fragments_per_peak / counts_info$total_fragments
    }
  } else {
    anno = as.factor(annotation)
    n_anno = length(levels(anno))
    mat = matrix(nrow = counts_info$npeak, ncol = n_anno)
    if (norm){
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat[,ix] /
                            matrix(counts_info$fragments_per_sample[ix],
                                   nrow = counts_info$npeak,
                                   ncol = length(ix),
                                   byrow = TRUE))
      }
    } else{
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(counts_mat[,ix]) / sum(counts_mat[,ix])
      }
    }
    expectation = rowMeans(mat)
  }
  return(expectation)
}



#' compute_deviations
#'
#' Computes deviations across sets of annotations
#' @param peak_indices list of indices representing different sets of peaks
#' @param counts_mat matrix with counts
#' @param background_peaks background peaks matrix
#' @param expectation expectations computed using \code{\link{compute_expectations}}
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  list with two elements: 1) z: matrix with deviation z-scores,
#' 2) dev: matrix with normalized deviations
#' @seealso  \code{\link{compute_variability}}, \code{\link{plot_variability}}
#' @export
compute_deviations <- function(counts_mat,
                               background_peaks,
                               peak_indices = NULL,
                               expectation = NULL){

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

  sample_names <- colnames(counts_mat)


    if (inherits(counts_mat,"dgCMatrix")){
      counts_mat_norm <- get_normalized_counts(counts_mat,expectation, counts_info$fragments_per_sample)
    } else{
      counts_mat_norm <- counts_mat / outer(sqrt(expectation),sqrt(counts_info$fragments_per_sample))
    }

  results <- BiocParallel::bplapply(peak_indices,
                                      compute_deviations_single,
                                    counts_mat_norm = counts_mat_norm,
                                      background_peaks = background_peaks)

  out <- list()
  out$z <- t(vapply(results, function(x) x[["z"]], rep(0,counts_info$nsample)))
  out$dev <- t(vapply(results, function(x) x[["dev"]], rep(0,counts_info$nsample)))

  colnames(out$z) = sample_names
  colnames(out$dev) = sample_names

  return(out)
}

# Helper Functions -------------------------------------------------------------


compute_deviations_single <- function(peak_set,
                                      counts_mat_norm,
                                      background_peaks,
                                      intermediate_results = FALSE){
  
  #require Matrix (for some multiprocessing options)
  suppressPackageStartupMessages(library(Matrix, quietly = TRUE, warn.conflicts = FALSE))
  
  if (length(peak_set) == 0){
    return(list(z = rep(NA, ncol(counts_mat_norm)), dev = rep(NA,ncol(counts_mat_norm))))
  }

  ### counts_mat should already be normed!
  tf_count <- length(peak_set)

  if (tf_count == 1){
    observed <- as.vector(counts_mat_norm[peak_set,])
    expected <- rep(1, length(observed))
    observed_deviation <- observed - expected
    sampled <- counts_mat_norm[background_peaks[peak_set,],]
    sampled_expected <- matrix(1, nrow = nrow(sampled), ncol = ncol(sampled))
    sampled_deviation <- sampled - sampled_expected
  } else{
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1,
                           dims = c(1, nrow(counts_mat_norm)))
    
    observed <- as.vector(tf_vec %*% counts_mat_norm)
    expected <- rep(tf_count, length(observed))
    observed_deviation = (observed - expected) / tf_count
    
    niterations = ncol(background_peaks)
    sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]),
                              i = rep(1:niterations, each = tf_count),
                              x=1,
                              dims = c(niterations, nrow(counts_mat_norm)))
    
    sampled = as.matrix(sample_mat %*% counts_mat_norm)
    sampled_expected = matrix(tf_count, nrow = nrow(sampled), ncol = ncol(sampled)) 
    sampled_deviation = (sampled - sampled_expected) / tf_count    
  }
  

    
  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)

  normdev <- (observed_deviation - mean_sampled_deviation) 
  z <-  normdev / sd_sampled_deviation
  

  if (intermediate_results){
    out = list(z = z,
               dev = normdev,
               observed = observed,
               sampled = sampled,
               expected = expected,
               sampled_expected = sampled_expected,
               observed_deviation = observed_deviation ,
               sampled_deviation = sampled_deviation,
               )
  } else{
    out = list(z = z, 
               dev = normdev)
  }
  return(out)
}







