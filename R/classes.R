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
                                                    peaks = 'GRanges',
                                                    total_fragments = 'numeric',
                                                    fragments_per_cell = 'numeric',
                                                    fragments_per_peak = 'numeric'))



#' An S4 class the chromatin accessibility deviation for a set of motifs
#'
#' @slot results a list of objects (of type deviationResult), one for each motif
#' @slot adjusted_p_values a vector of p values that have been corrected for multiple testing
#' @slot sim_adjusted_p_values a vector of p values for a permuted set of matched peaks that have been corrected for multiple testing
#' @slot tfs a vector with the names for the motifs
deviationResultSet <- setClass("deviationResultSet",
                            slots = c(results = 'list',
                                      adjusted_p_values = 'vector',
                                      sim_adjusted_p_values = 'vector',
                                      tfs = 'vector'))




#' An S4 class the chromatin accessibility deviation for a motif
#'
#' @slot sampled_deviations a matrix with the deviation (observed-expected) for background sets of peaks
#' @slot observed_deviations a vector with the deviation (oberved-expected) for # of fragments in peaks with motif for each cell
#' @slot extra_deviations  a vector with the deviation (observed-expected) for # of fragments for an extra set of background peaks
#' @slot normalized_deviations a vector with the observed_deviations normalized by the RMS of the sampled_deviations
#' @slot normalized_variability a numeric metric for the square root of the sum of the observed squared deviations divided by the sum of the mean background squared deviations
#' @slot error_variability a numeric metric for the standard deviation of the variability metric calculated using different background sets
#' @slot z_scores Z scores for the observed_deviations given the sampled_deviations
#' @slot z_sd the standard deviation of the Z scores
#' @slot error_z_sd bootstrapped 95% confidence interval for the z_sd
#' @slot p_value p value for the standard deviation being greater than 1, based on chi-square test for variance
#' @slot sim_normalized_deviations a vector with the extra_deviations normalized by the RMS of the sampled_deviations
#' @slot sim_normalized_variability a numeric metric for the square root of the sum of the squared deviations of the extra permutation divided by the sum of the mean background squared deviations
#' @slot sim_error_variability a numeric metric for the standard deviation of the sim_normalized_variability calculated using different background sets
#' @slot sim_z_scores Z scores for the extra_deviations given the sampled_deviations
#' @slot sim_z_sd the standard deviation of the Z scores of the extra permutation
#' @slot sim_error_z_sd bootstrapped 95% confidence interval for the z_sd
#' @slot sim_p_value p value for the standard deviation being greater than 1 for the sim_z_sd, based on chi-square test for variance
#' 
deviationResult <- setClass("deviationResult",
                            slots = c(sampled_deviations = 'matrix',
                                      observed_deviations = 'vector',
                                      extra_deviations = 'vector',
                                      normalized_deviations = 'vector',
                                      normalized_variability = 'numeric',
                                      error_variability = 'numeric',
                                      z_scores = 'vector',
                                      z_sd = 'numeric',
                                      error_z_sd = 'numeric',
                                      p_value = 'numeric',
                                      sim_normalized_deviations = 'vector',
                                      sim_normalized_variability = 'numeric',
                                      sim_error_variability = 'numeric',
                                      sim_z_scores = 'vector',
                                      sim_z_sd = 'numeric',
                                      sim_error_z_sd = 'numeric',
                                      sim_p_value = 'numeric',                                      
                                      tf = 'character'))

setMethod("show",
          signature="deviationResult",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Deviations for ", object@tf, " for ",
                  length(object@observed_deviations), " cells/samples. \n", sep = "")
            invisible(NULL)
          })

setValidity("deviationResult", function(object) {
  msg <- NULL
  valid <- TRUE
  #Check dimensions of sampled_deviations matched observed_deviations
  if (nrow(object@sampled_deviations) != length(object@observed_deviations)){
    valid <- FALSE
    msg <- c(msg, "Number of rows of sampled_deviations must be same as length of observed_deviations")
  }
  #Check for minimum iterations
  if (ncol(object@sampled_deviations) < 10){
    valid <- FALSE
    msg <- c(msg, "Must have at least 10 columns in sampled_deviations (iterations)")
  }
  if (valid) TRUE else msg})



#' @param object An object
#' @param ... Extra named arguments passed to FUN
#' @rdname sampled_deviations
#' @export 
setGeneric("sampled_deviations", function(object, ...) standardGeneric("sampled_deviations"))

#' @describeIn deviationResult return deviationResult
setMethod("sampled_deviations", "deviationResult",
          function(object) object@sampled_deviations)

setGeneric("observed_deviations", function(object, ...) standardGeneric("observed_deviations"))

setMethod("observed_deviations", "deviationResult",
          function(object) object@observed_deviations)

setGeneric("extra_deviations", function(object, ...) standardGeneric("extra_deviations"))

setMethod("extra_deviations", "deviationResult",
          function(object) object@extra_deviations)


setGeneric("compute_z_score", function(object) standardGeneric("compute_z_score"))

setMethod("compute_z_score","deviationResult",
          function(object){            
            #Compute for data
            object@z_scores = (object@observed_deviations - rowMeans(object@sampled_deviations)) / apply(object@sampled_deviations,1,sd)
            object@z_sd = sd(object@z_scores)
            object@p_value = pchisq((length(object@z_scores)-1)*(object@z_sd**2), length(object@z_scores)-1, lower.tail = F)
            object@error_z_sd = quantile(sapply(1:100, function(x){
              s = sample(1:length(object@z_scores), size = length(object@z_scores), replace=T)
              return(sd(object@z_scores[s]))
            }),c(0.025,0.975))
            #Compute for extra iteration
            object@sim_z_scores = (object@extra_deviations - rowMeans(object@sampled_deviations)) / apply(object@sampled_deviations,1,sd)
            object@sim_z_sd = sd(object@sim_z_scores)
            object@sim_p_value = pchisq((length(object@sim_z_scores)-1)*(object@sim_z_sd**2), length(object@sim_z_scores)-1, lower.tail = F)
            object@sim_error_z_sd = quantile(sapply(1:100, function(s){
              s = sample(1:length(object@sim_z_scores), size = length(object@sim_z_scores), replace=T)
              return(sd(object@sim_z_scores[s]))
            }),c(0.025,0.975))
            #return object
            return(object)
          })


setGeneric("compute_variability", function(object) standardGeneric("compute_variability"))

setMethod("compute_variability","deviationResult",
          function(object){            
            #Compute for data
            object@normalized_deviations = object@observed_deviations / sqrt(rowMeans(object@sampled_deviations**2))
            object@normalized_variability = sqrt(sum(object@observed_deviations**2) / sum(rowMeans(object@sampled_deviations**2)))
            object@error_variability = sd( sqrt(sum(object@observed_deviations**2) / colSums(object@sampled_deviations**2)))
            #Compute for extra iteration
            object@sim_normalized_deviations = object@extra_deviations / sqrt(rowMeans(object@sampled_deviations**2))
            object@sim_normalized_variability = sqrt(sum(object@extra_deviations**2) / sum(rowMeans(object@sampled_deviations**2)))
            object@sim_error_variability = sd( sqrt(sum(object@extra_deviations**2) / colSums(object@sampled_deviations**2)))
            #return object
            return(object)
          })



#' An S4 class for mapping kmers to wildcard kmers
#'
#' @slot k length of kmers
#' @slot m number of consecutive mismatches allowed
#' @slot mapping dgCMatrix with mapping

kmerMapping <- setClass("kmerMapping",
                            slots = c(mapping = 'dgCMatrix',
                                      colMapping = 'vector',
                                      k = 'numeric',
                                      m = 'numeric',
                                      l = 'numeric'))

setMethod("show",
          signature="kmerMapping",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Mapping of kmers with k = ", object@k, " to kmers of the same length with up to m = ",
                object@m, " possible consecutive mismatches. \n", sep = "")
            invisible(NULL)
          })


setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "missing"),
          definition = function(x, i, j) {
            x@mapping[i, colnames(x@mapping)]
            })

setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "character"),
          definition = function(x, i ,j ) {
            stopifnot(is.character(i), is.character(j))
            x@mapping[i, x@colMapping[j]]
          })

setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "character"),
          definition = function(x, i, j ) {
            x@mapping[as.character(i), x@colMapping[j]]
          })

setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "missing"),
          definition = function(x, i ,j ) {
            x@mapping[as.character(i), colnames(x@mapping)]
          })

setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "DNAStringSet"),
          definition = function(x, i ,j) {
            x@mapping[as.character(i), x@colMapping[as.character(j)]]
          })

setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "DNAStringSet"),
          definition = function(x, i, j) {
            x@mapping[i, x@colMapping[as.character(j)]]
          })

setMethod("[", signature = signature(x = "kmerMapping", i = "ANY"),
          definition = function(x, i) {
           stop(paste("Must supply list of kmers as first argument.",
                      " kmers list must be either character vector or DNAStringSet", sep="\n", collapse=""), call. = F)
          })






            
            






                                      