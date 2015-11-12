#' An S4 class the chromatin accessibility deviation for a set of motifs
#'
#' @slot counts_sort
#' @slot counts_sort_rev
#' @slot bias_sort
#' @slot bias_sort_rev
#' @slot window
#' @slot npeaks
#'
deviationBackgroundParameters <- setClass("deviationBackgroundParameters",
                               slots = c(counts_sort = 'vector',
                                         counts_sort_rev = 'vector',
                                         bias_sort = 'vector',
                                         bias_sort_rev = 'vector',
                                         window = 'numeric',
                                         npeaks = 'numeric'))


#' An S4 class that holds fragment counts across peak regions
#'
#' @slot counts
#' @slot peaks
#' @slot total_fragments
#' @slot fragments_per_cell
#' @slot fragments_per_peak
#'
fragmentCounts <- setClass("fragmentCounts",
                                          slots = c(counts = 'Matrix',
                                                    peaks = 'GenomicRanges',
                                                    total_fragments = 'numeric',
                                                    fragments_per_cell = 'numeric',
                                                    fragments_per_peak = 'numeric'
                                                    ))


#' An S4 class that holds background peak set
#'
#' @slot peaks
#' @slot background_peaks
#'
backgroundPeaks <- setClass("backgroundPeaks",
                           slots = c(peaks = 'GenomicRanges',
                                     background_peaks = 'matrix'
                           ))




#' An S4 class the chromatin accessibility deviation for a set of motifs
#'
#' @slot results a list of objects (of type deviationResult), one for each motif
#' @slot adjusted_p_values a vector of p values that have been corrected for multiple testing
#' @slot sim_adjusted_p_values a vector of p values for a permuted set of matched peaks that have been corrected for multiple testing
#' @slot tfs a vector with the names for the motifs
deviationResultSet <- setClass("deviationResultSet",
                            slots = c(results = 'list',
                                      adjusted_p_values = 'vector',
                                      tfs = 'vector'))




#' An S4 class to store stuff
#'
#'
deviationResult <- setClass("deviationResult",
                            slots = c(normalized_deviations = 'numeric',
                                      variability = 'numeric',
                                      variability_error = 'vector',
                                      #normalized_deviations_error = 'matrix',
                                      name = 'character'))




#' An S4 class for mapping kmers to wildcard kmers
#'
#' @slot k length of kmers
#' @slot m number of consecutive mismatches allowed
#' @slot mapping dgCMatrix with mapping
#
# kmerMapping <- setClass("kmerMapping",
#                             slots = c(mapping = 'Matrix::dgCMatrix',
#                                       colMapping = 'vector',
#                                       k = 'numeric',
#                                       m = 'numeric',
#                                       l = 'numeric'))
#
# setMethod("show",
#           signature="kmerMapping",
#           definition = function(object){
#             cat("An object of class ", class(object), "\n", sep = "")
#             cat("Mapping of kmers with k = ", object@k, " to kmers of the same length with up to m = ",
#                 object@m, " possible consecutive mismatches. \n", sep = "")
#             invisible(NULL)
#           })
#
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "missing"),
#           definition = function(x, i, j) {
#             x@mapping[i, colnames(x@mapping)]
#             })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "character"),
#           definition = function(x, i ,j ) {
#             stopifnot(is.character(i), is.character(j))
#             x@mapping[i, x@colMapping[j]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "character"),
#           definition = function(x, i, j ) {
#             x@mapping[as.character(i), x@colMapping[j]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "missing"),
#           definition = function(x, i ,j ) {
#             x@mapping[as.character(i), colnames(x@mapping)]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "DNAStringSet"),
#           definition = function(x, i ,j) {
#             x@mapping[as.character(i), x@colMapping[as.character(j)]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "DNAStringSet"),
#           definition = function(x, i, j) {
#             x@mapping[i, x@colMapping[as.character(j)]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "ANY"),
#           definition = function(x, i) {
#            stop(paste("Must supply list of kmers as first argument.",
#                       " kmers list must be either character vector or DNAStringSet", sep="\n", collapse=""), call. = F)
#           })
#
#












