#' makeBiasBins
#'
#' Makes bins based on fragment counts
#' @param object fragment counts stored as RangedSummarizedExperiment, 
#' SummarizedExperiment, matrix, or Matrix
#' @param bias vector of some bias signal (usually gc content) for each row of 
#' object
#' @param nbins number of bins for each category, see Details
#' @param frac fraction of peaks within given bin to select randomly
#' @param ... addtional arguments
#' @return SummarizedExperiment storing bias bins annotation
#' @details Will create nbins * 3 annotations based on sampling from peaks with 
#' a certain fragment count, fragment count, or fragment count & bias.
#' @export
#' @author Alicia Schep
#' @examples
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' bb <- makeBiasBins(mini_counts)
setGeneric("makeBiasBins", 
           function(object, ...) standardGeneric("makeBiasBins"))

#' @describeIn makeBiasBins method for SummarizedExperiment
#' @export
setMethod(makeBiasBins, c(object = "SummarizedExperiment"), 
          function(object, 
                   bias = rowData(object)$bias, nbins = 25, frac = 0.3) {
            object <- counts_check(object)
            makeBiasBins_core(object, bias, nbins, frac)
          })

#' @describeIn makeBiasBins method for RangedSummarizedExperiment
#' @export
setMethod(makeBiasBins, c(object = "RangedSummarizedExperiment"), 
          function(object, 
                   bias = rowRanges(object)$bias, nbins = 25, frac = 0.3) {
            object <- counts_check(object)
            makeBiasBins_core(object, bias, nbins, frac)
          })

#' @describeIn makeBiasBins method for Matrix or matrix
#' @export
setMethod(makeBiasBins, c(object = "MatrixOrmatrix"), 
          function(object, 
                   bias, nbins = 25, frac = 0.3) {
            makeBiasBins_core(object, bias, nbins, frac)
          })


makeBiasBins_core <- function(object, bias, nbins = 25, frac = 0.3) {
  npeaks <- nrow(object)
  fragments_per_peak <- getFragmentsPerPeak(object)
  # make bias bins
  bias_quantiles <- quantile(bias, seq(0, 1, 1/nbins))
  bias_cut <- cut(bias, breaks = bias_quantiles)
  bias_bins <- split(seq_len(npeaks), bias_cut)
  names(bias_bins) <- vapply(seq_len(nbins), 
                             function(x) paste("bias_bin_", x, sep = "", 
                                               collapse = ""),
                             "")
  # make count bins
  pseudo_counts <- fragments_per_peak + runif(npeaks, min = 0, max = 0.1)
  count_quantiles <- quantile(pseudo_counts, seq(0, 1, 1/nbins))
  count_cut <- cut(pseudo_counts, breaks = count_quantiles)
  count_bins <- split(seq_len(npeaks), count_cut)
  names(count_bins) <- vapply(seq_len(nbins), 
                              function(x) paste("count_bin_", x, sep = "", 
                                                collapse = ""),
                              "")
  # make bias / count bins
  nbins <- round(sqrt(nbins))
  bias_quantiles <- quantile(bias, seq(0, 1, 1/nbins))
  bias_cut <- cut(bias, breaks = bias_quantiles)
  tmp_bias_bins <- split(seq_len(npeaks), bias_cut)
  count_quantiles <- quantile(pseudo_counts, seq(0, 1, 1/nbins))
  count_cut <- cut(pseudo_counts, breaks = count_quantiles)
  tmp_count_bins <- split(seq_len(npeaks), count_cut)
  bias_count_bins <- unlist(lapply(seq_len(nbins), 
                            function(x) 
                              lapply(seq_len(nbins), 
                                     function(y)
                                       intersect(tmp_bias_bins[[y]], 
                                                 tmp_count_bins[[x]]))),
                            recursive = FALSE)
  names(bias_count_bins) <- vapply(seq_len(nbins), 
                                   function(x) 
                                     vapply(seq_len(nbins), 
                                            function(y) 
                                              paste("bias_count_bin_", 
                                                    x, "_", y, sep = "", 
                                                    collapse = ""),
                                            ""),
                                   rep("",nbins))
  tmp <- c(bias_bins, count_bins, bias_count_bins)
  sets <- lapply(tmp, function(x) sample(x, size = length(x) * frac))
  mean_bias <- vapply(sets, function(x) mean(bias[x]),0)
  mean_counts <- vapply(sets, function(x) mean(fragments_per_peak[x]),0)
  rd <- DataFrame(bias = mean_bias, counts = mean_counts)
  out <- SummarizedExperiment(assays = 
                                list(annotationMatches = 
                                       convert_from_ix_list(sets, 
                                                            nrow(object))),
                              colData = rd)
  return(out)
}


#' makePermutedSets
#'
#' Makes annotations sets with similar bias to input sets
#' @param object fragment counts stored as RangedSummarizedExperiment, 
#' SummarizedExperiment, matrix, or Matrix
#' @param annotations annotations as SummarizedExperiment, matrix, or list
#' @param bias vector of some bias signal (usually gc content) for each row of 
#' object
#' @param window number of nearest neighbors to consider
#' @param ... additional arguments
#' @return SummarizedExperiment storing bias bins annotation
#' @details Will create nbins * 3 annotations based on sampling from peaks with 
#' a certain fragment count, fragment count, or fragment count & bias.
#' @author Alicia Schep
#' @export
#' @examples
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' data(example_motifs, package = "motifmatchr")
#' library(motifmatchr)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' motif_ix <- matchMotifs(example_motifs, mini_counts, 
#'                          genome = BSgenome.Hsapiens.UCSC.hg19)
#'
#' perm_sets <- makePermutedSets(mini_counts, motif_ix)
setGeneric("makePermutedSets", 
           function(object, annotations, ...) 
             standardGeneric("makePermutedSets"))

#' @describeIn makePermutedSets method for SummarizedExperiment and 
#' SummarizedExperiment
#' @export
setMethod(makePermutedSets, c(object = "SummarizedExperiment", 
                              annotations = "SummarizedExperiment"), 
          function(object,  annotations,
                   bias = rowData(object)$bias, window = 10) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            makePermutedSets_core(object, peak_indices, bias, window = window, 
                                    colData = colData(annotations))
          })

#' @describeIn makePermutedSets method for RangedSummarizedExperiment and 
#' SummarizedExperiment
#' @export
setMethod(makePermutedSets, c(object = "RangedSummarizedExperiment", 
                              annotations = "SummarizedExperiment"), 
          function(object,  annotations,
                   bias = rowRanges(object)$bias, window = 10) {
            object <- counts_check(object)
            annotations <- matches_check(annotations)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            makePermutedSets_core(object, peak_indices, bias, window = window, 
                                    colData = colData(annotations))
          })

#' @describeIn makePermutedSets method for Matrix or matrix and 
#' SummarizedExperiment
#' @export
setMethod(makePermutedSets, c(object = "MatrixOrmatrix", 
                              annotation = "SummarizedExperiment"), 
          function(object, annotations,
                   bias, window = 10) {
            annotations <- matches_check(annotations)
            peak_indices <- 
              convert_to_ix_list(annotationMatches(assays(annotations)))
            makePermutedSets_core(object, peak_indices, bias, window = window, 
                                    colData = colData(annotations))
          })


#' @describeIn makePermutedSets method for SummarizedExperiment and 
#' MatrixOrmatrix
#' @export
setMethod(makePermutedSets, c(object = "SummarizedExperiment", 
                              annotations = "MatrixOrmatrix"), 
          function(object,  annotations,
                   bias = rowData(object)$bias, window = 10) {
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            makePermutedSets_core(object, peak_indices, bias, window = window)
          })

#' @describeIn makePermutedSets method for RangedSummarizedExperiment and 
#' MatrixOrmatrix
#' @export
setMethod(makePermutedSets, c(object = "RangedSummarizedExperiment", 
                              annotations = "MatrixOrmatrix"), 
          function(object,  annotations,
                   bias = rowRanges(object)$bias, window = 10) {
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotations)
            makePermutedSets_core(object, peak_indices, bias, window = window)
          })

#' @describeIn makePermutedSets method for Matrix/matrix and Matrix/matrix
#' @export
setMethod(makePermutedSets, c(object = "MatrixOrmatrix", 
                              annotation = "MatrixOrmatrix"), 
          function(object, annotations,
                   bias, window = 10) {
            peak_indices <- convert_to_ix_list(annotations)
            makePermutedSets_core(object, peak_indices, bias, window = window)
          })


#' @describeIn makePermutedSets method for SummarizedExperiment and list
#' @export
setMethod(makePermutedSets, c(object = "SummarizedExperiment", 
                              annotations = "list"), 
          function(object,  annotations,
                   bias = rowData(object)$bias, window = 10) {
            object <- counts_check(object)
            makePermutedSets_core(object, annotations, bias, window = window)
          })

#' @describeIn makePermutedSets method for RangedSummarizedExperiment and list
#' @export
setMethod(makePermutedSets, c(object = "RangedSummarizedExperiment",
                              annotations = "list"), 
          function(object,  annotations,
                   bias = rowRanges(object)$bias, window = 10) {
            object <- counts_check(object)
            makePermutedSets_core(object, annotations, bias, window = window)
          })

#' @describeIn makePermutedSets method for Matrix or matrix and list
#' @export
setMethod(makePermutedSets, c(object = "MatrixOrmatrix",
                              annotation = "list"), 
          function(object, annotations,
                   bias, window = 10) {
            makePermutedSets_core(object, annotations, bias, window = window)
          })



makePermutedSets_core <- function(object, peak_indices, bias, window = 10, 
                                  colData = NULL) {
  
  
  
  bg <- get_background_peaks_alternative(object, bias, niterations = 1, 
                                         window = window)
  sets <- lapply(seq_along(peak_indices), function(x) bg[peak_indices[[x]], 1])
  names(sets) <- vapply(names(peak_indices), 
                        function(x) 
                          paste("permuted.", x, collapse = "", sep = ""),
                        "")
  
  if (is.null(colData)) {
    colData <- DataFrame(name = names(sets))
  } else{
    if ("name" %in% colnames(colData)) {
      colData$name <- vapply(colData$name, 
                             function(x) 
                               paste("permuted.", x, collapse = "", sep = ""),
                             "")
    } else {
      colData$name <- names(sets)
    }}

  rownames(colData) <- names(sets)

  out <- SummarizedExperiment(assays = 
                                list(matches = 
                                       convert_from_ix_list(sets, 
                                                            nrow(object))), 
                              colData = colData)
  
  return(out)
}

get_background_peaks_alternative <- function(object, bias, niterations = 50,
                                             window = 500) {
  
  fragments_per_peak <- getFragmentsPerPeak(object)
  if (min(fragments_per_peak) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  npeak <- nrow(object)
  
  norm_mat <- cbind(fragments_per_peak, bias)
  
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat), t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(seq_len(npeak%/%grpsize + ((npeak%%grpsize) != 0)), 
                 function(x) ((x - 
                                 1) * grpsize + 1):(min(x * grpsize, npeak)))
  
  bghelper <- function(grp, tmp_vals, niterations) {
    tmp_nns <- nabor::knn(tmp_vals, tmp_vals[grp, ], window + 1, eps = 0)$nn.idx
    if (niterations == 1) {
      return(matrix(vapply(seq_len(nrow(tmp_nns)), 
                           function(x){
                             sample(tmp_nns[x, ][tmp_nns[x,] != grp[x]],
                                    niterations, replace = TRUE)}, 
                           0),
                    ncol = 1))
    } else {
      return(t(vapply(seq_len(nrow(tmp_nns)), 
                      function(x) {
                        sample(tmp_nns[x, ][tmp_nns[x, ] != grp[x]], 
                               niterations, replace = TRUE)},
                      rep(0, niterations))))
    }
  }
  
  background_peaks <- do.call(rbind, bplapply(grps, bghelper, tmp_vals, 
                                              niterations))
  
  return(background_peaks)
}
