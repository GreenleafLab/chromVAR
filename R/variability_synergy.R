
#' get_annotation_synergy
#'
#' @param object result from compute_deviations
#' @param annotations SummarizedExperiment of annotation matches
#' @param background_peaks optional, matrix of background peaks
#' @param variabilities optional, variabilities computed from compute_variability
#' @param expectation optional, expected fraction of reads per peak, as computed
#' by compute_expectations
#' @param nbg number of background iterations
#' @param ... additional arguments
#' @details should only be run on small number of motifs/kmers/peaksets 
#' (very slow!)
#'
#' @return synergy matrix
#'
#' @export
setGeneric("get_annotation_synergy", 
           function(object, annotations, ...) standardGeneric("get_annotation_synergy"))



#' @describeIn get_annotation_synergy object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("get_annotation_synergy", c(object = "SummarizedExperiment", 
                                      annotations = "SummarizedExperiment"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object), 
                   variabilities = NULL, nbg = 25) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotation_matches(annotations))
            get_annotation_synergy_core(counts(object), 
                                        peak_indices,
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })

#' @describeIn get_annotation_synergy object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("get_annotation_synergy", c(object = "SummarizedExperiment", 
                                      annotations = "MatrixOrmatrix"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object), 
                   variabilities = NULL, nbg = 25) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            get_annotation_synergy_core(counts(object), 
                                        peak_indices, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn get_annotation_synergy object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("get_annotation_synergy", c(object = "SummarizedExperiment",
                                      annotations = "list"), 
          function(object, 
                   annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object),
                   variabilities = NULL, nbg = 25) {
            object <- counts_check(object)
            get_annotation_synergy_core(counts(object),
                                        annotations,
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn get_annotation_synergy object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("get_annotation_synergy", c(object = "MatrixOrmatrix", 
                                      annotations = "SummarizedExperiment"), 
          function(object, 
                   annotations, 
                   background_peaks, 
                   expectation = compute_expectations(object),
                   variabilities = NULL, 
                   nbg = 25) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotation_matches(annotations))
            get_annotation_synergy_core(object, 
                                        peak_indices, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })

#' @describeIn get_annotation_synergy object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("get_annotation_synergy", c(object = "MatrixOrmatrix", 
                                      annotations = "MatrixOrmatrix"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object), 
                   variabilities = NULL, nbg = 25) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            get_annotation_synergy_core(object, 
                                        peak_indices, 
                                        background_peaks,
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })


#' @describeIn get_annotation_synergy object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("get_annotation_synergy", c(object = "MatrixOrmatrix", 
                                      annotations = "list"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object),
                   variabilities = NULL, nbg = 25) {
            get_annotation_synergy_core(object, annotations, 
                                        background_peaks, 
                                        expectation, 
                                        variabilities = variabilities,
                                        nbg = nbg)
          })



get_annotation_synergy_core <- function(counts_mat,
                                        peak_indices, 
                                        background_peaks, 
                                        expectation , 
                                        variabilities = NULL,
                                        nbg = 25) {
  
  
  if (min(get_fragments_per_peak(counts_mat)) <= 0) 
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
    names(peak_indices) <- 1:length(peak_indices)
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
  for (i in 1:(l - 1)) {
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

#' get_annotation_correlation
#'
#' @param object result from compute_deviations
#' @param annotations SummarizedExperiment of annotation matches
#' @param background_peaks optional, matrix of background peaks
#' @param variabilities optional, variabilities computed from 
#' compute_variability
#' @param expectation optional, expected fraction of reads per peak, as computed
#' by compute_expectations
#' @param ... additional arguments
#' @details should only be run on small number of motifs/kmers/peaksets 
#' (very slow!)
#'
#' @return correlation matrix
#'
#'@export
setGeneric("get_annotation_correlation", 
           function(object, annotations, ...) 
             standardGeneric("get_annotation_correlation"))



#' @describeIn get_annotation_correlation object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("get_annotation_correlation", c(object = "SummarizedExperiment", 
                                          annotations = "SummarizedExperiment"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotation_matches(annotations))
            get_annotation_correlation_core(counts(object), 
                                            peak_indices,
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })

#' @describeIn get_annotation_correlation object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("get_annotation_correlation", c(object = "SummarizedExperiment", 
                                          annotations = "MatrixOrmatrix"), 
          function(object, annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            get_annotation_correlation_core(counts(object), 
                                            peak_indices, 
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn get_annotation_correlation object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("get_annotation_correlation", c(object = "SummarizedExperiment",
                                          annotations = "list"), 
          function(object, 
                   annotations, 
                   background_peaks = get_background_peaks(object), 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            object <- counts_check(object)
            get_annotation_correlation_core(counts(object),
                                            annotations,
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn get_annotation_correlation object and annotations are 
#' SummarizedExperiment 
#' @export
setMethod("get_annotation_correlation", c(object = "MatrixOrmatrix", 
                                          annotations = "SummarizedExperiment"), 
          function(object, 
                   annotations, 
                   background_peaks, 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotation_matches(annotations))
            get_annotation_correlation_core(object, 
                                            peak_indices, 
                                            background_peaks, 
                                            expectation, 
                                            variabilities = variabilities)
          })

#' @describeIn get_annotation_correlation object is SummarizedExperiment, 
#' annotations are Matrix
#' @export
setMethod("get_annotation_correlation", c(object = "MatrixOrmatrix", 
                                          annotations = "MatrixOrmatrix"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            stopifnot(canCoerce(annotations, "lMatrix"))
            annotations <- as(annotations, "lMatrix")
            peak_indices <- convert_to_ix_list(annotations)
            get_annotation_correlation_core(object, 
                                            peak_indices, 
                                            background_peaks,
                                            expectation, 
                                            variabilities = variabilities)
          })


#' @describeIn get_annotation_correlation object is SummarizedExperiment, 
#' annotations are list
#' @export
setMethod("get_annotation_correlation", c(object = "MatrixOrmatrix", 
                                          annotations = "list"), 
          function(object, annotations, background_peaks, 
                   expectation = compute_expectations(object), 
                   variabilities = NULL) {
            get_annotation_correlation_core(object, annotations, 
                                            background_peaks, expectation, 
                                            variabilities = variabilities)
          })


get_annotation_correlation_core <- function(counts_mat, 
                                            peak_indices, 
                                            background_peaks,
                                            expectation, 
                                            variabilities) {
  
  
  if (min(get_fragments_per_peak(counts_mat)) <= 0) 
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
    names(peak_indices) <- 1:length(peak_indices)
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
  
  setlen <- sapply(tmpsets, length)
  
  bgsets <- unlist(lapply(setlen, 
                          function(x) lapply(1:nbg, 
                                             function(y) sample(peak_set, 
                                                                x, 
                                                                replace = FALSE))), 
                   recursive = FALSE)
  
  bgvar <- do.call(c, bplapply(bgsets, 
                                             compute_variability_single, 
                                             counts_mat = counts_mat, 
                                             background_peaks = background_peaks, 
                                             expectation = expectation))
  
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) 
                        (tmpvar[x] - mean(bgvar[((x -1) * nbg + 1):(x * nbg)])) /
                        sd(bgvar[((x - 1) * nbg + 1):(x * nbg)]))
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


