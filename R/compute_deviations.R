
#' compute_expectations
#'
#' @param object SummarizedExperiment
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
compute_expectations <- function(object,
                                 norm = FALSE,
                                 annotation = NULL){

  if (is.null(annotation)){
    if (norm){
      expectation = rowSums(assays(object)$counts /
                              matrix(get_fragments_per_sample(object),
                                     nrow = nrow(object),
                                     ncol = ncol(object),
                                     byrow = TRUE))
    } else{
      expectation = get_fragments_per_peak(object) / get_total_fragments(object)
    }
  } else {
    if (length(annotation) == 1 && annotation %in% colnames(colData(object))){
      anno = as.factor(colData(object)[[annotation]])
    } else if (length(annotation) == ncol(object)){
      anno = as.factor(annotation)
    } else{
      stop("annotation must be vector of length of columns of object, or a character vector referring to name of column in colData of object")
    }
    n_anno = length(levels(anno))
    mat = matrix(nrow = counts_info$npeak, ncol = n_anno)
    if (norm){
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(assays(object)$counts[,ix] /
                            matrix(get_fragments_per_sample(object)[ix],
                                   nrow = nrow(object),
                                   ncol = length(ix),
                                   byrow = TRUE))
      }
    } else{
      for (i in 1:n_anno){
        ix = which(anno = levels(anno)[i])
        mat[,i] = rowSums(assays(object)$counts[,ix]) / sum(assays(object)$counts[,ix])
      }
    }
    expectation = rowMeans(mat)
  }
  return(expectation)
}



#' compute_deviations
#'
#' Computes deviations across sets of annotations
#' @param object SummarizedExperiment
#' @param background_peaks background peaks matrix
#' @param annotations SummarizedExperiment
#' @param expectation expectations computed using \code{\link{compute_expectations}}
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  list with two elements: 1) z: matrix with deviation z-scores,
#' 2) dev: matrix with normalized deviations
#' @seealso  \code{\link{compute_variability}}, \code{\link{plot_variability}}
#' @export
compute_deviations <- function(object,
                               annotations = NULL,
                               background_peaks = NULL,
                               expectation = NULL){

  stopifnot(inherits(object, "SummarizedExperiment"))

  if (is.null(assays(object)$counts)) stop("No counts slot")

  if (inherits(assays(object)$counts,"matrix")){
    assays(object)$counts = Matrix(assays(object)$counts)
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
    stop("annotations must be given as SummarizedExperiment with a 'match' slot")
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

  sample_names <- colnames(object)

  results <- BiocParallel::bplapply(peak_indices,
                                      compute_deviations_single,
                                    counts_mat = assays(object)$counts,
                                      background_peaks = background_peaks,
                                    expectation = expectation)

  z <- t(vapply(results, function(x) x[["z"]], rep(0,ncol(object))))
  dev <- t(vapply(results, function(x) x[["dev"]], rep(0,ncol(object))))

  colnames(z) = colnames(dev) = sample_names

  out <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                              colData = colData(object),
                              rowData = colData(annotations))

  return(out)
}

# Helper Functions -------------------------------------------------------------


compute_deviations_single <- function(peak_set,
                                      counts_mat,
                                      background_peaks,
                                      expectation = NULL,
                                      intermediate_results = FALSE,
                                      threshold = 1){

  #require Matrix (for some multiprocessing options)
  suppressPackageStartupMessages(library(Matrix, quietly = TRUE, warn.conflicts = FALSE))

  if (length(peak_set) == 0){
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA,ncol(counts_mat))))
  }

  fragments_per_sample = colSums(counts_mat)

  ### counts_mat should already be normed!
  tf_count <- length(peak_set)

  if (tf_count == 1){
    observed <- as.vector(counts_mat[peak_set,])
    expected <- expectation[peak_set] * fragments_per_sample
    observed_deviation <- (observed - expected) / expected

    sampled <- counts_mat[background_peaks[peak_set,],]
    sampled_expected <- outer(expectation[background_peaks[peak_set,]], fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected) / sampled_expected
  } else{
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1,
                           dims = c(1, nrow(counts_mat)))

    observed <- as.vector(tf_vec %*% counts_mat)

    expected <- as.vector(tf_vec %*% expectation %*% fragments_per_sample)
    observed_deviation = (observed - expected) / expected

    niterations = ncol(background_peaks)
    sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]),
                              i = rep(1:niterations, each = tf_count),
                              x=1,
                              dims = c(niterations, nrow(counts_mat)))

    sampled = as.matrix(sample_mat %*% counts_mat)
    sampled_expected = as.matrix(sample_mat %*% expectation %*% fragments_per_sample)
    sampled_deviation = (sampled - sampled_expected) / sampled_expected

  }

  fail_filter <- which(expected < threshold)

  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)

  normdev <- (observed_deviation - mean_sampled_deviation)
  z <-  normdev / sd_sampled_deviation

  if (length(fail_filter) > 0){
    z[fail_filter] = NA
    normdev[fail_filter] = NA
  }

  if (intermediate_results){
    out = list(z = z,
               dev = normdev,
               observed = observed,
               sampled = sampled,
               expected = expected,
               sampled_expected = sampled_expected,
               observed_deviation = observed_deviation ,
               sampled_deviation = sampled_deviation
               )
  } else{
    out = list(z = z,
               dev = normdev)
  }
  return(out)
}





compute_deviations_single_no_bg <- function(peak_set,
                                      counts_mat,
                                      expectation = NULL,
                                      threshold = 1){

  #require Matrix (for some multiprocessing options)
  suppressPackageStartupMessages(library(Matrix, quietly = TRUE, warn.conflicts = FALSE))

  if (length(peak_set) == 0){
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA,ncol(counts_mat))))
  }

  fragments_per_sample = colSums(counts_mat)

  ### counts_mat should already be normed!
  tf_count <- length(peak_set)

  if (tf_count == 1){
    observed <- as.vector(counts_mat[peak_set,])
    expected <- expectation[peak_set] * fragments_per_sample
    observed_deviation <- (observed - expected) / expected

  } else{
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1,
                           dims = c(1, nrow(counts_mat)))

    observed <- as.vector(tf_vec %*% counts_mat)

    expected <- as.vector(tf_vec %*% expectation %*% fragments_per_sample)
    observed_deviation = (observed - expected) / expected

  }

  return(observed_deviation)
}

compute_deviations_no_bg <- function(object,
                                     annotations = NULL,
                                     background_peaks = NULL,
                                     expectation = NULL){

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

  sample_names <- colnames(object)

  results <- t(simplify2array(BiocParallel::bplapply(peak_indices,
                                    compute_deviations_single_no_bg,
                                    counts_mat = assays(object)$counts,
                                    expectation = expectation)))

  colnames(results) = sample_names

  out <- SummarizedExperiment(assays = list(deviations = results),
                              colData = colData(object),
                              rowData = colData(annotations))

  return(out)
}


