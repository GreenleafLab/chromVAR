
#' getAnnotationSynergy
#'
#' @param object result from computeDeviations
#' @param annotations SummarizedExperiment of annotation matches
#' @param background_peaks optional, matrix of background peaks
#' @param variabilities optional, variabilities computed from computeVariability
#' @param expectation optional, expected fraction of reads per peak, as computed
#' by computeExpectations
#' @param nbg number of background iterations
#' @param ... additional arguments
#' @details should only be run on small number of motifs/kmers/peaksets 
#' (very slow!)
#'
#' @return synergy matrix
#'
#' @export
setGeneric("getAnnotationSynergy", 
           function(object, annotations, ...) 
             standardGeneric("getAnnotationSynergy"))



#' @describeIn getAnnotationSynergy object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("getAnnotationSynergy", c(object = "SummarizedExperiment", 
                                      annotations = "SummarizedExperiment"), 
          function(object, annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object), 
                   variabilities = NULL, nbg = 25) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            getAnnotationSynergy_core(counts(object), 
                                        peak_indices,
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })

#' @describeIn getAnnotationSynergy object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("getAnnotationSynergy", c(object = "SummarizedExperiment", 
                                      annotations = "MatrixOrmatrix"), 
          function(object, annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object), 
                   variabilities = NULL, nbg = 25) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            getAnnotationSynergy_core(counts(object), 
                                        peak_indices, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn getAnnotationSynergy object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("getAnnotationSynergy", c(object = "SummarizedExperiment",
                                      annotations = "list"), 
          function(object, 
                   annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object),
                   variabilities = NULL, nbg = 25) {
            object <- counts_check(object)
            getAnnotationSynergy_core(counts(object),
                                        annotations,
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn getAnnotationSynergy object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("getAnnotationSynergy", c(object = "MatrixOrmatrix", 
                                      annotations = "SummarizedExperiment"), 
          function(object, 
                   annotations, 
                   background_peaks, 
                   expectation = computeExpectations(object),
                   variabilities = NULL, 
                   nbg = 25) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            getAnnotationSynergy_core(object, 
                                        peak_indices, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })

#' @describeIn getAnnotationSynergy object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("getAnnotationSynergy", c(object = "MatrixOrmatrix", 
                                      annotations = "MatrixOrmatrix"), 
          function(object, annotations, background_peaks, 
                   expectation = computeExpectations(object), 
                   variabilities = NULL, nbg = 25) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            getAnnotationSynergy_core(object, 
                                        peak_indices, 
                                        background_peaks,
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn getAnnotationSynergy object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("getAnnotationSynergy", c(object = "MatrixOrmatrix", 
                                      annotations = "list"), 
          function(object, annotations, background_peaks, 
                   expectation = computeExpectations(object),
                   variabilities = NULL, nbg = 25) {
            getAnnotationSynergy_core(object, annotations, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })



getAnnotationSynergy_core <- function(counts_mat,
                                        peak_indices, 
                                        background_peaks, 
                                        expectation , 
                                        variabilities = NULL,
                                        nbg = 25) {
  
  
  if (min(getFragmentsPerPeak(counts_mat)) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))
  
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
  
  
  if (is.null(variabilities)) {
    variabilities <- do.call(c, bplapply(peak_indices, 
                                         compute_variability_single, 
                                         counts_mat = counts_mat, 
                                         background_peaks = background_peaks, 
                                         expectation = expectation))
  } else if (is.data.frame(variabilities)) {
    variabilities <- variabilities$variability
  } else {
    stop("Incorrect variabilities input")
  }
  
  
  stopifnot(length(variabilities) == length(peak_indices))
  
  l <- length(peak_indices)
  stopifnot(l >= 2)
  outmat <- matrix(nrow = l, ncol = l)
  var_order <- order(variabilities, decreasing = TRUE)
  for (i in seq_len(l - 1)) {
    outmat[i, (i + 1):l] <- 
      get_variability_boost_helper(peak_indices[[var_order[i]]], 
                                   counts_mat, 
                                   background_peaks, 
                                   peak_indices[var_order[(i + 1):l]], 
                                   expectation, 
                                   nbg)
    outmat[(i + 1):l, i] <- outmat[i, (i + 1):l]
  }
  diag(outmat) <- 0
  colnames(outmat) <- rownames(outmat) <- names(peak_indices)[var_order]
  outmat
}

#' getAnnotationCorrelation
#'
#' @param object result from computeDeviations
#' @param annotations SummarizedExperiment of annotation matches
#' @param background_peaks optional, matrix of background peaks
#' @param variabilities optional, variabilities computed from 
#' computeVariability
#' @param expectation optional, expected fraction of reads per peak, as computed
#' by computeExpectations
#' @param ... additional arguments
#' @details should only be run on small number of motifs/kmers/peaksets 
#' (very slow!)
#'
#' @return correlation matrix
#'
#'@export
setGeneric("getAnnotationCorrelation", 
           function(object, annotations, ...) 
             standardGeneric("getAnnotationCorrelation"))



#' @describeIn getAnnotationCorrelation object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("getAnnotationCorrelation", 
          c(object = "SummarizedExperiment", 
            annotations = "SummarizedExperiment"), 
          function(object, annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            getAnnotationCorrelation_core(counts(object), 
                                            peak_indices,
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })

#' @describeIn getAnnotationCorrelation object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("getAnnotationCorrelation", c(object = "SummarizedExperiment", 
                                          annotations = "MatrixOrmatrix"), 
          function(object, annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            getAnnotationCorrelation_core(counts(object), 
                                            peak_indices, 
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn getAnnotationCorrelation object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("getAnnotationCorrelation", c(object = "SummarizedExperiment",
                                          annotations = "list"), 
          function(object, 
                   annotations, 
                   background_peaks = getBackgroundPeaks(object), 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            object <- counts_check(object)
            getAnnotationCorrelation_core(counts(object),
                                            annotations,
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn getAnnotationCorrelation object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("getAnnotationCorrelation", 
          c(object = "MatrixOrmatrix", 
            annotations = "SummarizedExperiment"), 
          function(object, 
                   annotations, 
                   background_peaks, 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            getAnnotationCorrelation_core(object, 
                                            peak_indices, 
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })

#' @describeIn getAnnotationCorrelation object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("getAnnotationCorrelation", c(object = "MatrixOrmatrix", 
                                          annotations = "MatrixOrmatrix"), 
          function(object, annotations, background_peaks, 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            getAnnotationCorrelation_core(object, 
                                            peak_indices, 
                                            background_peaks,
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn getAnnotationCorrelation object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("getAnnotationCorrelation", c(object = "MatrixOrmatrix", 
                                          annotations = "list"), 
          function(object, annotations, background_peaks, 
                   expectation = computeExpectations(object), 
                   variabilities = NULL) {
            getAnnotationCorrelation_core(object, annotations, 
                                            background_peaks, expectation, 
                                            variabilities = variabilities)
          })


getAnnotationCorrelation_core <- function(counts_mat, 
                                            peak_indices, 
                                            background_peaks,
                                            expectation, 
                                            variabilities) {
  
  
  if (min(getFragmentsPerPeak(counts_mat)) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))
  
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
  
  
  if (is.null(variabilities)) {
    variabilities <- do.call(c, 
                             bplapply(peak_indices, 
                                      compute_variability_single, 
                                      counts_mat = counts_mat, 
                                      background_peaks = background_peaks, 
                                      expectation = expectation))
  } else if (is.data.frame(variabilities)) {
    variabilities <- variabilities$variability
  } else {
    stop("Incorrect variabilities input")
  }
  
  stopifnot(length(variabilities) == length(peak_indices))
  
  l <- length(peak_indices)
  stopifnot(l >= 2)
  outmat <- matrix(1, nrow = l, ncol = l)
  var_order <- order(variabilities, decreasing = TRUE)
  ixs <- which(upper.tri(outmat), arr.ind = TRUE)
  tmp <- lapply(seq_len(nrow(ixs)), function(x) unname(ixs[x, ]))
  outmat[ixs] <- outmat[ixs[, c(2:1)]] <- 
    do.call(c, bplapply(tmp, 
                        get_cor_helper, 
                        peak_indices[var_order], 
                        counts_mat,
                        background_peaks, 
                        expectation))
  
  colnames(outmat) <- rownames(outmat) <- names(peak_indices)[var_order]
  outmat
}



# Variability synergy ----------------------------------------------------------

get_variability_boost_helper <- function(peak_set, 
                                         counts_mat,
                                         background_peaks, 
                                         peak_indices, 
                                         expectation, 
                                         nbg) {
  
  tmpsets <- remove_nonoverlap(peak_indices, peak_set)
  tmpvar <- do.call(c, bplapply(tmpsets, 
                                compute_variability_single, 
                                counts_mat = counts_mat, 
                                background_peaks = background_peaks, 
                                expectation = expectation))
  
  setlen <- vapply(tmpsets, length, 0)
  
  bgsets <- unlist(lapply(setlen, 
                          function(x) 
                            lapply(seq_len(nbg), 
                                   function(y) sample(peak_set, 
                                                      x, 
                                                      replace = FALSE))), 
                   recursive = FALSE)
  
  bgvar <- do.call(c, bplapply(bgsets, 
                               compute_variability_single, 
                               counts_mat = counts_mat, 
                               background_peaks = background_peaks, 
                               expectation = expectation))
  
  var_boost <- 
    vapply(seq_along(tmpvar), 
           function(x) {
             (tmpvar[x] - mean(bgvar[((x - 1) * nbg + 1):(x * nbg)])) /
               sd(bgvar[((x - 1) * nbg + 1):(x * nbg)])
             },
           0)
  return(var_boost)
}


get_cor_helper <- function(ixs, peak_indices, counts_mat, background_peaks, 
                           expectation) {
  
  peak_set1 <- peak_indices[[ixs[1]]][which(peak_indices[[ixs[1]]] %ni% 
                                              peak_indices[[ixs[2]]])]
  dev1 <- compute_deviations_single(peak_set1, counts_mat, 
                                    background_peaks, 
                                    expectation)
  peak_set2 <- peak_indices[[ixs[2]]][which(peak_indices[[ixs[2]]] %ni% 
                                              peak_indices[[ixs[1]]])]
  dev2 <- compute_deviations_single(peak_set2, counts_mat, 
                                    background_peaks, 
                                    expectation)
  cor(dev1$dev, dev2$dev)
}


