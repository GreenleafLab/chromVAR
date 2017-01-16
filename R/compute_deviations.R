# Generic Functions ------------------------------------------------------------

#' compute_expectations
#'
#' @param object SummarizedExperiment
#' @param norm weight all samples equally?
#' @param group an group vector, optional
#' @param ... additional arguments
#' @details By default, this function will compute the expected fraction of reads
#' per peak as the the total fragments per peak across all samples divided by total
#' reads in peaks in all samples. Optionally, norm can be set to TRUE and then the
#' expectation will be the average fraction of reads in a peak across the cells.
#' This is not recommended for single cell applications as cells with very few reads
#' will have a large impact.  Another option is to give a vector of groups, in which
#' case the expectation will be the average fraction of reads per peak within each group.
#' If group vector is provided and norm is set to TRUE then within each group
#' the fraction of reads per peak is the average fraction of reads per peak in each
#' sample.  Otherwise, the within group fraction of reads per peak is based on the
#'reads per peak within the sample divided by the total reads within each sample.
#' @return vector with expected fraction of reads per peak
#' @export
#' @export
setGeneric("compute_expectations", 
           function(object, ...) standardGeneric("compute_expectations"))

#' compute_deviations
#'
#' Computes deviations in chromatin accessibility across sets of annotations
#' @param object chromVARCounts object 
#' @param annotations chromVARAnnotations object
#' @param background_peaks (optional) background peaks matrix computed using 
#' \code{\link{get_background_peaks}}, computed internally with default
#' paramaters if not provided
#' @param expectation (optional) expectations computed using
#' \code{\link{compute_expectations}}, computed automatically if not provided
#' @param ... additional arguments
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{chromVARDeviations-class}}, which inherits from 
#' SummarizedExperiment, and has two assays: deviations and deviation scores.
#' @seealso  \code{\link{compute_variability}}, \code{\link{plot_variability}}
#' @export
setGeneric("compute_deviations", 
           function(object, annotations, ...) standardGeneric("compute_deviations"))


#' deviations
#' 
#' Accessor for bias corrected deviations from \code{\link{chromVARDeviations-class}} object
#' @rdname deviations
#' @param object chromVARDeviations object
setGeneric("deviations", function(object) standardGeneric("deviations"))

#' deviation_scores
#' 
#' Accessor for deviation Z-scores from \code{\link{chromVARDeviations-class}} object
#' @rdname deviation_scores
#' @param object chromVARDeviations object
setGeneric("deviation_scores", 
           function(object) standardGeneric("deviation_scores"))

# Accessors --------------------------------------------------------------------

#' @describeIn chromVARDeviations access deviation Z-scores matrix
#' @param object chromVARDeviations object
#' @export
setMethod("deviations", c(object = "chromVARDeviations"),
          function(object) {
            assays(object)$deviations
          })

#' @describeIn chromVARDeviations access deviation Z-scores matrix
#' @export
setMethod("deviation_scores", c(object = "chromVARDeviations"), 
          function(object) {
            assays(object)$z
          })

#' rbind method chromVARDeviations 
#'
#' Concatenates chromVARDeviations results for different sets of annotations
#' @param ... chromVARDeviations object to be combined
#' @param deparse.level See ?base::rbind for a description of this argument.
#' @export
setMethod("rbind", "chromVARDeviations",
          function(..., deparse.level=1){
            inputs = list(...)
            
            all_rowdata_colnames <- unique(do.call(c,lapply(inputs, 
                                                 function(x) 
                                                   colnames(rowData(x)))))
            
            common_colnames <- all_rowdata_colnames[sapply(all_rowdata_colnames,
                                      function(x){
                                        all(sapply(inputs, 
                                                   function(y) {
                                                     x %in% colnames(rowData(y))
                                                   }))})]
            inputs = lapply(inputs, function(x){
              rowData(x) <- rowData(x)[,common_colnames, drop = FALSE]
              x
            })
            out <- SummarizedExperiment(assays = 
                                          list(z = do.call(rbind, 
                                                           lapply(inputs,
                                                                  function(x) 
                                                                    assays(x)$z)),
                                               deviations = do.call(rbind, 
                                                                    lapply(inputs,
                                                                           function(x) 
                                                                             assays(x)$deviations))), 
                                        colData = colData(inputs[[1]]), 
                                        rowData = do.call(rbind, lapply(inputs,
                                                                        function(x) 
                                                                          rowData(x))) )
            
            return(new("chromVARDeviations", out))
          })

#' cbind method for chromVARDeviations 
#'
#' cbind returns an error when applied to chromVARDeviations because results for
#' all cells or samples should originate from same compute_deviations computation
#' @param ... chromVARDeviations object to be combined
#' @param deparse.level See ?base::rbind for a description of this argument.
#' @export
#' @export
setMethod("cbind", "chromVARDeviations",  
          function(..., deparse.level=1){
            stop("Can't concatenate chromVARDeviations horizontally")
          })

#' @describeIn compute_expectations method for Matrix or matrix
#' @export
setMethod("compute_expectations", c(object = "MatrixOrmatrix"), 
          function(object, 
                   norm = FALSE, group = NULL) {
            compute_expectations_core(object, norm = norm, group = group)
          })


#' @describeIn compute_expectations method for SummarizedExperiment with counts 
#' slot
#' @export
setMethod("compute_expectations", c(object = "SummarizedExperiment"), 
          function(object, 
                   norm = FALSE, 
                   group = NULL) {
            object <- counts_check(object)
            compute_expectations_core(counts(object), norm = norm, group = group)
          })


#' @describeIn compute_deviations object and annotations are SummarizedExperiment 
#' @export
setMethod("compute_deviations", c(object = "SummarizedExperiment", 
                                  annotations = "SummarizedExperiment"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object)) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(assays(annotations)$matches)
            compute_deviations_core(counts(object), 
                                    peak_indices,
                                    background_peaks, 
                                    expectation, 
                                    colData = colData(object),
                                    rowData = colData(annotations))
          })

#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("compute_deviations", c(object = "SummarizedExperiment", 
                                  annotations = "MatrixOrmatrix"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object)) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            compute_deviations_core(counts(object), 
                                    peak_indices, 
                                    background_peaks, 
                                    expectation, 
                                    colData = colData(object))
          })


#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("compute_deviations", c(object = "SummarizedExperiment",
                                  annotations = "list"), 
          function(object, 
                   annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object)) {
            object <- counts_check(object)
            compute_deviations_core(counts(object),
                                    annotations,
                                    background_peaks, 
                                    expectation, 
                                    colData = colData(object))
          })

#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are missing
#' @export
setMethod("compute_deviations", c(object = "SummarizedExperiment", 
                                  annotations = "missingOrNULL"), 
          function(object, 
                   annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object)) {
            message(paste0("Annotations not provided, ",
                           "so chromVAR being run on individual peaks..."))
            object <- counts_check(object)
            peak_indices <- split(seq_len(nrow(object)), seq_len(nrow(object)))
            compute_deviations_core(counts(object), 
                                    peak_indices, 
                                    background_peaks, 
                                    expectation, 
                                    colData = colData(object))
          })

#' @describeIn compute_deviations object and annotations are SummarizedExperiment 
#' @export
setMethod("compute_deviations", c(object = "MatrixOrmatrix", 
                                  annotations = "SummarizedExperiment"), 
          function(object, 
                   annotations, 
                   background_peaks, 
                   expectation = compute_expectations(object)) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(assays(annotations)$matches)
            compute_deviations_core(object, peak_indices, background_peaks, 
                                    expectation, 
                                    rowData = colData(annotations))
          })

#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("compute_deviations", c(object = "MatrixOrmatrix", 
                                  annotations = "MatrixOrmatrix"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object)) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            compute_deviations_core(object, 
                                    peak_indices, 
                                    background_peaks,
                                    expectation)
          })


#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("compute_deviations", c(object = "MatrixOrmatrix", 
                                  annotations = "list"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object)) {
            compute_deviations_core(object, annotations, 
                                    background_peaks, expectation)
          })

#' @describeIn compute_deviations object is SummarizedExperiment, 
#' annotations are missing
#' @export
setMethod("compute_deviations", c(object = "MatrixOrmatrix", 
                                  annotations = "missingOrNULL"), 
          function(object, annotations,
                   background_peaks, 
                   expectation = compute_expectations(object)) {
            message(paste0("Annotations not provided, ", 
                           "so chromVAR being run on individual peaks..."))
            peak_indices <- split(seq_len(nrow(object)), seq_len(nrow(object)))
            compute_deviations_core(object, peak_indices, background_peaks, 
                                    expectation)
          })

compute_deviations_core <- function(counts_mat,
                                    peak_indices, 
                                    background_peaks,
                                    expectation, 
                                    rowData = NULL, 
                                    colData = NULL) {
  
  if (min(get_fragments_per_peak(counts_mat)) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))
  
  if (is.null(colData)) 
    colData <- DataFrame(seq_len(ncol(counts_mat)), 
                         row.names = colnames(counts_mat))[,FALSE]
  
  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names = FALSE)
  if (is.null(tmp) || 
      !(all.equal(tmp, as.integer(tmp))) || 
      max(tmp) > nrow(counts_mat) || 
      min(tmp) < 1) {
    stop("Annotations are not valid")
  }
  
  if (is.null(names(peak_indices))) {
    names(peak_indices) <- seq_along(peak_indices)
  }
  
  sample_names <- colnames(counts_mat)
  
  results <- bplapply(peak_indices,
                      compute_deviations_single, 
                      counts_mat = counts_mat, 
                      background_peaks = background_peaks, 
                      expectation = expectation)
  
  z <- t(vapply(results, function(x) x[["z"]], rep(0, ncol(counts_mat))))
  dev <- t(vapply(results, function(x) x[["dev"]], rep(0, ncol(counts_mat))))
  
  colnames(z) <- colnames(dev) <- sample_names
  
  out <- SummarizedExperiment(assays = list(deviations = dev, z = z), 
                              colData = colData, 
                              rowData = rowData)
  
  return(new("chromVARDeviations", out))
}




# Helper Functions -------------------------------------------------------------

compute_expectations_core <- function(object, norm = FALSE, group = NULL) {
  
  if (is.null(group)) {
    if (norm) {
      expectation <- rowSums(object/matrix(get_fragments_per_sample(object), 
                                           nrow = nrow(object), 
                                           ncol = ncol(object), 
                                           byrow = TRUE))
    } else {
      expectation <- get_fragments_per_peak(object) / 
        get_total_fragments(object)
    }
  } else {
    if (length(group) == 1 && group %in% colnames(colData(object))) {
      anno <- as.factor(colData(object)[[group]])
    } else if (length(group) == ncol(object)) {
      anno <- as.factor(group)
    } else {
      stop(paste0("Group must be vector of length of columns of object, ",
                  "or a character vector referring to name of column in",
                  "colData of object"))
    }
    n_anno <- length(levels(anno))
    mat <- matrix(nrow = nrow(object), ncol = n_anno)
    if (norm) {
      for (i in seq_along(n_anno)) {
        ix <- which(anno == levels(anno)[i])
        mat[, i] <- rowSums(object[, ix] / 
                              matrix(get_fragments_per_sample(object)[ix], 
                                     nrow = nrow(object), 
                                     ncol = length(ix), 
                                     byrow = TRUE))
      }
    } else {
      for (i in seq_along(n_anno)) {
        ix <- which(anno == levels(anno)[i])
        mat[, i] <- rowSums(object[, ix])/sum(object[, ix])
      }
    }
    expectation <- rowMeans(mat)
  }
  return(expectation)
}


compute_deviations_single <- function(peak_set, 
                                      counts_mat, 
                                      background_peaks, 
                                      expectation = NULL, 
                                      intermediate_results = FALSE, 
                                      threshold = 1) {
  
  if (length(peak_set) == 0) {
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA, ncol(counts_mat))))
  }
  
  fragments_per_sample <- colSums(counts_mat)
  
  ### counts_mat should already be normed!
  tf_count <- length(peak_set)
  
  if (tf_count == 1) {
    observed <- as.vector(counts_mat[peak_set, ])
    expected <- expectation[peak_set] * fragments_per_sample
    observed_deviation <- (observed - expected)/expected
    
    sampled <- counts_mat[background_peaks[peak_set, ], ]
    sampled_expected <- outer(expectation[background_peaks[peak_set, ]], 
                              fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected)/sampled_expected
  } else {
    tf_vec <- sparseMatrix(j = peak_set, 
                           i = rep(1, tf_count),
                           x = 1, 
                           dims = c(1, 
                                    nrow(counts_mat)))
    
    observed <- as.vector(tf_vec %*% counts_mat)
    
    expected <- as.vector(tf_vec %*% expectation %*% fragments_per_sample)
    observed_deviation <- (observed - expected)/expected
    
    niterations <- ncol(background_peaks)
    sample_mat <- sparseMatrix(j = as.vector(background_peaks[peak_set, seq_len(niterations)]), 
                               i = rep(seq_len(niterations), each = tf_count), 
                               x = 1,
                               dims = c(niterations, 
                                        nrow(counts_mat)))
    
    sampled <- as.matrix(sample_mat %*% counts_mat)
    sampled_expected <- as.matrix(sample_mat %*% expectation %*% fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected)/sampled_expected
    
  }
  
  fail_filter <- which(expected < threshold)
  
  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)
  
  normdev <- (observed_deviation - mean_sampled_deviation)
  z <- normdev/sd_sampled_deviation
  
  if (length(fail_filter) > 0) {
    z[fail_filter] <- NA
    normdev[fail_filter] <- NA
  }
  
  if (intermediate_results) {
    out <- list(z = z, dev = normdev, 
                observed = observed, 
                sampled = sampled, 
                expected = expected, 
                sampled_expected = sampled_expected, 
                observed_deviation = observed_deviation, 
                sampled_deviation = sampled_deviation)
  } else {
    out <- list(z = z, dev = normdev)
  }
  return(out)
}




setAs("SummarizedExperiment", "chromVARDeviations", function(from) {
  out <- new("chromVARDeviations", from)
  validObject(out)
  out
})
