# Generic Functions ------------------------------------------------------------

#' computeExpectations
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
#' reads per peak within the sample divided by the total reads within each sample.
#' The group can also be given by a length 1 character vector representing the
#' name of a column in the colData of the input object if the input is a
#' SummarizedExperiment
#' @return vector with expected fraction of reads per peak.
#' @export
#' @author Alicia Schep
#' @examples
#'
#' # First get some data
#' data(mini_counts, package = "chromVAR")
#'
#' # Compute expectations
#' expectations <- computeExpectations(mini_counts)
#'
setGeneric("computeExpectations",
           function(object, ...) standardGeneric("computeExpectations"))

#' computeDeviations
#'
#' Computes deviations in chromatin accessibility across sets of annotations
#' @param object chromVARCounts object
#' @param annotations chromVARAnnotations object
#' @param background_peaks (optional) background peaks matrix computed using
#' \code{\link{getBackgroundPeaks}}, computed internally with default
#' paramaters if not provided
#' @param expectation (optional) expectations computed using
#' \code{\link{computeExpectations}}, computed automatically if not provided
#' @param ... additional arguments
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{chromVARDeviations-class}}, which inherits from
#' SummarizedExperiment, and has two assays: deviations and deviation scores.
#' @seealso  \code{\link{computeVariability}}, \code{\link{plotVariability}}
#' @export
#' @author Alicia Schep
#' @examples
#' # Register BiocParallel
#' BiocParallel::register(BiocParallel::SerialParam())
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' # load annotation matrix; result from matchMotifs
#' data(mini_ix, package = "chromVAR")
#'
#' # computing deviations
#' dev <- computeDeviations(object = mini_counts,
#'                          annotations = mini_ix)
setGeneric("computeDeviations",
           function(object, annotations, ...) standardGeneric("computeDeviations"))


#' deviations
#'
#' Accessor for bias corrected deviations from \code{\link{chromVARDeviations-class}} object
#' @rdname deviations
#' @name deviations
#' @aliases deviations,chromVARDeviations-method
#' @param object chromVARDeviations object
#' @return matrix of bias corrected deviations
#' @author Alicia Schep
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' bias_corrected_deviations <- deviations(mini_dev)
setGeneric("deviations", function(object) standardGeneric("deviations"))

#' deviationScores
#'
#' Accessor for deviation Z-scores from \code{\link{chromVARDeviations-class}} object
#' @rdname deviationScores
#' @name deviationScores
#' @aliases deviationScores,chromVARDeviations-method
#' @param object chromVARDeviations object
#' @return The deviationScores and deviations accessors both return matrices.
#' @return matrix of deviation Z-scores
#' @author Alicia Schep
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' scores <- deviationScores(mini_dev)
setGeneric("deviationScores",
           function(object) standardGeneric("deviationScores"))

# Accessors --------------------------------------------------------------------

#' @rdname deviations
#' @export
setMethod("deviations", c(object = "chromVARDeviations"),
          function(object) {
            assays(object)$deviations
          })

#' @rdname deviationScores
#' @export
setMethod("deviationScores", c(object = "chromVARDeviations"),
          function(object) {
            assays(object)$z
          })

#' rbind method chromVARDeviations
#'
#' Concatenates chromVARDeviations results for different sets of annotations
#' @param ... chromVARDeviations object to be combined
#' @param deparse.level See ?base::rbind for a description of this argument.
#' @export
#' @return chromVARDeviations object
#' @author Alicia Schep
#' @seealso \code{\link{chromVARDeviations-class}}
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' doubledev <- rbind(mini_dev, mini_dev) #concatenate two of the same tother
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
#' all cells or samples should originate from same computeDeviations computation
#' @param ... chromVARDeviations object to be combined
#' @param deparse.level See ?base::rbind for a description of this argument.
#' @return chromVARDeviations object
#' @author Alicia Schep
#' @seealso \code{\link{chromVARDeviations-class}}
#' @export
setMethod("cbind", "chromVARDeviations",
          function(..., deparse.level=1){
            stop("Can't concatenate chromVARDeviations horizontally")
          })

#' @describeIn computeExpectations method for Matrix or matrix
#' @export
setMethod("computeExpectations", c(object = "MatrixOrmatrix"),
          function(object,
                   norm = FALSE, group = NULL) {
            if (!is.null(group)){
              if (length(group) == ncol(object)) {
                group <- as.factor(group)
              } else {
                stop("Group must be vector of length of columns of object")
              }
            }
            compute_expectations_core(object, norm = norm, group = group)
          })


#' @describeIn computeExpectations method for SummarizedExperiment with counts
#' slot
#' @export
setMethod("computeExpectations", c(object = "SummarizedExperiment"),
          function(object,
                   norm = FALSE,
                   group = NULL) {
            object <- counts_check(object)
            if (!is.null(group)){
              if (length(group) == 1 && group %in% colnames(colData(object))) {
                group <- as.factor(colData(object)[[group]])
              } else if (length(group) == ncol(object)) {
                group <- as.factor(group)
              } else {
                stop("Group must be vector of length of columns of object, ",
                     "or a character vector referring to name of column in",
                     "colData of object")
              }
            }
            compute_expectations_core(counts(object), norm = norm, group = group)
          })


#' @describeIn computeDeviations object and annotations are SummarizedExperiment
#' @export
setMethod("computeDeviations", c(object = "SummarizedExperiment",
                                  annotations = "SummarizedExperiment"),
          function(object, annotations,
                   background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            compute_deviations_core(counts(object),
                                    peak_indices,
                                    background_peaks,
                                    expectation,
                                    colData = colData(object),
                                    rowData = colData(annotations))
          })

#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are Matrix
#' @export
setMethod("computeDeviations", c(object = "SummarizedExperiment",
                                  annotations = "MatrixOrmatrix"),
          function(object, annotations,
                   background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
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


#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are list
#' @export
setMethod("computeDeviations", c(object = "SummarizedExperiment",
                                  annotations = "list"),
          function(object,
                   annotations,
                   background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
            object <- counts_check(object)
            compute_deviations_core(counts(object),
                                    annotations,
                                    background_peaks,
                                    expectation,
                                    colData = colData(object))
          })

#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are missing
#' @export
setMethod("computeDeviations", c(object = "SummarizedExperiment",
                                  annotations = "missingOrNULL"),
          function(object,
                   annotations,
                   background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
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

#' @describeIn computeDeviations object and annotations are SummarizedExperiment
#' @export
setMethod("computeDeviations", c(object = "MatrixOrmatrix",
                                  annotations = "SummarizedExperiment"),
          function(object,
                   annotations,
                   background_peaks,
                   expectation = computeExpectations(object)) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            compute_deviations_core(object, peak_indices, background_peaks,
                                    expectation,
                                    rowData = colData(annotations))
          })

#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are Matrix
#' @export
setMethod("computeDeviations", c(object = "MatrixOrmatrix",
                                  annotations = "MatrixOrmatrix"),
          function(object, annotations, background_peaks,
                   expectation = computeExpectations(object)) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            compute_deviations_core(object,
                                    peak_indices,
                                    background_peaks,
                                    expectation)
          })


#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are list
#' @export
setMethod("computeDeviations", c(object = "MatrixOrmatrix",
                                  annotations = "list"),
          function(object, annotations, background_peaks,
                   expectation = computeExpectations(object)) {
            compute_deviations_core(object, annotations,
                                    background_peaks, expectation)
          })

#' @describeIn computeDeviations object is SummarizedExperiment,
#' annotations are missing
#' @export
setMethod("computeDeviations", c(object = "MatrixOrmatrix",
                                  annotations = "missingOrNULL"),
          function(object, annotations,
                   background_peaks,
                   expectation = computeExpectations(object)) {
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

  if (min(getFragmentsPerPeak(counts_mat)) <= 0)
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

  rowData$fractionMatches <- vapply(results, function(x) x[["matches"]], 0)
  rowData$fractionBackgroundOverlap <- vapply(results, function(x) x[["overlap"]], 
                                              0)
  
  out <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                              colData = colData,
                              rowData = rowData)

  return(new("chromVARDeviations", out))
}




# Helper Functions -------------------------------------------------------------

compute_expectations_core <- function(object, norm = FALSE, group = NULL) {

  if (is.null(group)) {
    if (norm) {
      expectation <- rowSums(object/matrix(getFragmentsPerSample(object),
                                           nrow = nrow(object),
                                           ncol = ncol(object),
                                           byrow = TRUE))
    } else {
      expectation <- getFragmentsPerPeak(object) /
        getTotalFragments(object)
    }
  } else {
    n_anno <- length(levels(group))
    mat <- matrix(nrow = nrow(object), ncol = n_anno)
    if (norm) {
      for (i in seq_len(n_anno)) {
        ix <- which(group == levels(group)[i])
        mat[, i] <- rowSums(object[, ix] /
                              matrix(getFragmentsPerSample(object)[ix],
                                     nrow = nrow(object),
                                     ncol = length(ix),
                                     byrow = TRUE))
      }
    } else {
      for (i in seq_len(n_anno)) {
        ix <- which(group == levels(group)[i])
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
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA, ncol(counts_mat)),
                matches = 0, overlap = NA))
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
    bg_overlap <- sum(background_peaks[peak_set,] == peak_set[1]) / 
      ncol(background_peaks)
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
    sample_mat <-
      sparseMatrix(j = 
                     as.vector(background_peaks[peak_set, 
                                                seq_len(niterations)]),
                   i = rep(seq_len(niterations), each = tf_count),
                   x = 1,
                   dims = c(niterations,
                            nrow(counts_mat)))
    
    sampled <- as.matrix(sample_mat %*% counts_mat)
    sampled_expected <- as.matrix(sample_mat %*% expectation %*% 
                                    fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected)/sampled_expected

    bg_overlap <- mean(sample_mat %*% t(tf_vec)) / tf_count
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
                sampled_deviation = sampled_deviation,
                matches = tf_count / nrow(counts_mat),
                overlap = bg_overlap)
  } else {
    out <- list(z = z, dev = normdev, matches = tf_count / nrow(counts_mat),
                overlap = bg_overlap)
  }
  return(out)
}




setAs("SummarizedExperiment", "chromVARDeviations", function(from) {
  out <- new("chromVARDeviations", from)
  validObject(out)
  out
})
