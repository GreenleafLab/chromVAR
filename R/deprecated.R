#' Deprecated functions in chromVAR
#' 
#' chromVAR has moved functions and methods to camelCase from snake_case.
#' The following functions have been deprecated and replaced with a different
#' name:
#' \itemize{
#'   \item get_annotation_synergy is now \code{\link{getAnnotationSynergy}}
#'   \item make_bias_bins is now \code{\link{makeBiasBins}}
#'   \item make_permuted_sets is now \code{\link{makePermutedSets}}
#'   \item pwm_distance is now \code{\link{pwmDistance}}
#'   \item plot_variability is now \code{\link{plotVariability}} 
#'   \item plot_deviations_tsne is now \code{\link{plotDeviationsTsne}}
#'   \item get_jaspar_motifs is now \code{\link{getJasparMotifs}}
#'   \item match_kmers is now \code{\link{matchKmers}}
#'   \item deviations_covariability is now \code{\link{deviationsCovariability}}
#'   \item plot_kmer_mismatch is now \code{\link{plotKmerMismatch}}
#'   \item get_total_fragments is now \code{\link{getTotalFragments}}
#'   \item get_fragments_per_peak is now \code{\link{getFragmentsPerPeak}}
#'   \item get_fragments_per_sample is now \code{\link{getFragmentsPerSample}}
#'   \item get_peaks is now \code{\link{getPeaks}}
#'   \item get_counts is now \code{\link{getCounts}}
#'   \item get_sample_depths is now \code{\link{getSampleDepths}}
#'   \item read_macs2_narrowpeaks is now \code{\link{readNarrowpeaks}}
#'   \item get_annotations is now \code{\link{getAnnotations}}
#'   \item filter_samples is now \code{\link{filterSamples}}
#'   \item filter_samples_plot is now \code{\link{filterSamplesPlot}}
#'   \item filter_peaks is now \code{\link{filterPeaks}}
#'   \item annotation_matches is now \code{\link{annotationMatches}}
#'   \item deviations_tsne is now \code{\link{deviationsTsne}}
#'   \item differential_deviations is now \code{\link{differentialDeviations}}
#'   \item differential_variability is now \code{\link{differentialVariability}}
#'   \item compute_variability is now \code{\link{computeVariability}}
#'   \item compute_expectations is now \code{\link{computeExpectations}}
#'   \item computeDeviations is now \code{\link{computeDeviations}}
#'   \item deviation_scores is now \code{\link{deviationScores}}
#'   \item get_sample_distance is now \code{\link{getSampleDistance}}
#'   \item get_sample_correlation is now \code{\link{getSampleCorrelation}}
#'   \item get_cis_groups is now \code{\link{getCisGroups}}
#'   \item add_gc_bias is now \code{\link{addGCBias}}
#'   \item get_background_peaks is now \code{\link{getBackgroundPeaks}}
#'   \item get_permuted_data is now \code{\link{getPermutedData}}
#'   \item assemble_kmers is now \code{\link{assembleKmers}}
#' }
#' @author Alicia Schep
#' @rdname chromVAR_deprecated
#' @param ... arguments passed to new function
#' @return calls the replacement functions
#' @name chromVAR_deprecated
NULL

#' @rdname chromVAR_deprecated
#' @export
assemble_kmers <- function(...)
{
  .Deprecated("assembleKmers")
  assembleKmers(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_annotation_synergy <- function(...)
{
  .Deprecated("getAnnotationSynergy")
  getAnnotationSynergy(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_annotation_correlation <- function(...)
{
  .Deprecated("getAnnotationCorrelation")
  getAnnotationCorrelation(...)
}

#' @rdname chromVAR_deprecated
#' @export
make_bias_bins <- function(...)
{
  .Deprecated("makeBiasBins")
  makeBiasBins(...)
}

#' @rdname chromVAR_deprecated
#' @export
make_permuted_sets <- function(...)
{
  .Deprecated("makePermutedSets")
  makeBiasBins(...)
}

#' @rdname chromVAR_deprecated
#' @export
pwm_distance <- function(...)
{
  .Deprecated("pwmDistance")
  makeBiasBins(...)
}

#' @rdname chromVAR_deprecated
#' @export
plot_variability <- function(...)
{
  .Deprecated("plotVariability")
  plotVariability(...)
}

#' @rdname chromVAR_deprecated
#' @export
plot_deviations_tsne <- function(...)
{
  .Deprecated("plotDeviationsTsne")
  plotDeviationsTsne(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_jaspar_motifs <- function(...)
{
  .Deprecated("getJasparMotifs")
  getJasparMotifs(...)
}

#' @rdname chromVAR_deprecated
#' @export
match_kmers <- function(...)
{
  .Deprecated("matchKmers")
  matchKmers(...)
}

#' @rdname chromVAR_deprecated
#' @export
deviations_covariability <- function(...)
{
  .Deprecated("deviationsCovariability")
  deviationsCovariability(...)
}

#' @rdname chromVAR_deprecated
#' @export
plot_kmer_mismatch <- function(...)
{
  .Deprecated("plotKmerMismatch")
  plotKmerMismatch(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_total_fragments <- function(...)
{
  .Deprecated("getTotalFragments")
  getTotalFragments(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_fragments_per_peak <- function(...)
{
  .Deprecated("getFragmentsPerPeak")
  getFragmentsPerPeak(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_fragments_per_sample <- function(...)
{
  .Deprecated("getFragmentsPerSample")
  getFragmentsPerSample(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_peaks <- function(...)
{
  .Deprecated("getPeaks")
  getPeaks(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_counts <- function(...)
{
  .Deprecated("getCounts")
  getCounts(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_sample_depths <- function(...)
{
  .Deprecated("getSampleDepths")
  getSampleDepths(...)
}

#' @rdname chromVAR_deprecated
#' @export
read_macs2_narrowpeaks <- function(...)
{
  .Deprecated("readNarrowpeaks")
  readNarrowpeaks(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_annotations <- function(...)
{
  .Deprecated("getAnnotations")
  getAnnotations(...)
}

#' @rdname chromVAR_deprecated
#' @export
filter_samples <- function(...)
{
  .Deprecated("filterSamples")
  filterSamples(...)
}

#' @rdname chromVAR_deprecated
#' @export
filter_samples_plot <- function(...)
{
  .Deprecated("filterSamplesPlot")
  filterSamplesPlot(...)
}

#' @rdname chromVAR_deprecated
#' @export
filter_peaks <- function(...)
{
  .Deprecated("filterPeaks")
  filterPeaks(...)
}

#' @rdname chromVAR_deprecated
#' @export
annotation_matches <- function(...)
{
  .Deprecated("annotationMatches")
  annotationMatches(...)
}

#' @rdname chromVAR_deprecated
#' @export
deviations_tsne <- function(...)
{
  .Deprecated("deviationsTsne")
  deviationsTsne(...)
}

#' @rdname chromVAR_deprecated
#' @export
differential_deviations <- function(...)
{
  .Deprecated("differentialDeviations")
  differentialDeviations(...)
}

#' @rdname chromVAR_deprecated
#' @export
differential_variability <- function(...)
{
  .Deprecated("differentialVariability")
  differentialVariability(...)
}

#' @rdname chromVAR_deprecated
#' @export
compute_variability <- function(...)
{
  .Deprecated("computeVariability")
  computeVariability(...)
}

#' @rdname chromVAR_deprecated
#' @export
compute_expectations <- function(...)
{
  .Deprecated("computeExpectations")
  computeExpectations(...)
}

#' @rdname chromVAR_deprecated
#' @export
compute_deviations <- function(...)
{
  .Deprecated("computeDeviations")
  computeDeviations(...)
}

#' @rdname chromVAR_deprecated
#' @export
deviation_scores <- function(...)
{
  .Deprecated("deviationScores")
  deviationScores(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_sample_distance <- function(...)
{
  .Deprecated("getSampleDistance")
  getSampleDistance(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_sample_correlation <- function(...)
{
  .Deprecated("getSampleCorrelation")
  getSampleCorrelation(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_cis_groups <- function(...)
{
  .Deprecated("getCisGroups")
  getCisGroups(...)
}

#' @rdname chromVAR_deprecated
#' @export
add_gc_bias <- function(...)
{
  .Deprecated("addGCBias")
  addGCBias(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_background_peaks <- function(...)
{
  .Deprecated("getBackgroundPeaks")
  getBackgroundPeaks(...)
}

#' @rdname chromVAR_deprecated
#' @export
get_permuted_data <- function(...)
{
  .Deprecated("getPermutedData")
  getPermutedData(...)
}




















