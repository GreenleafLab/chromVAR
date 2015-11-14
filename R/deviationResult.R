
#' An S4 class to store chromVAAR result
#' 
#' @slot deviations 
#' @slot variability
#' @slot variability_bounds
#' @slot name
#' @slot metric
#' @slot pvalues
#'
#' @aliases deviations, variability, variability_bounds, name, metric, pvalues
deviationResult <- setClass("deviationResult",
                            slots = c(deviations = 'numeric',
                                      variability = 'numeric',
                                      variability_bounds = 'numeric',
                                      name = 'character',
                                      metric = 'character',
                                      pvalues = 'numeric'))



setMethod("show",
          signature="deviationResult",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Deviations for ", object@tf, " for ",
                length(object@deviations), " cells/samples. \n", sep = "")
            invisible(NULL)
          })

#' metric
#' 
#' Accessor for which variability metric was used.
#' @param object either deviationResult or deviationResultSet object
#' @seealso \code{\link{variability}} 
setGeneric("metric", function(object, ...) standardGeneric("metric"))

#' @describeIn metric  returns name of metric used
#' @export
setMethod("metric", "deviationResult",
          function(object){object@metric})


#' variability
#' 
#' Accessor for variability metric in deviationResult or deviationResultSet.
#' Use \code{\link{metric}} to find metric
#' @param object either deviationResult or deviationResultSet object
#' @seealso \code{\link{metric}}, \code{\link{variability_bounds}} 
setGeneric("variability", function(object, ...) standardGeneric("variability"))

#' @describeIn variability  returns single numeric value representing the variability metric for that result
#' @export
setMethod("variability", "deviationResult",
          function(object){object@variability})

#' variability_bounds
#' 
#' Accessor for variability error in deviationResult or deviationResultSet.
#' variability error is computed by bootstrapping the samples and determining a 95% confidence interval
#' @param object either deviationResult or deviationResultSet object
#' @param lower return lowerbounds?  only applicable if object is deviationResultSet
#' @seealso \code{\link{metric}}, \code{\link{variability}} 
setGeneric("variability_bounds", function(object, ...) standardGeneric("variability_bounds"))

#' @describeIn variability_bounds returns numeric vector of length two, with first  
#' value the lower bound and second value the upper bound
#' @export
setMethod("variability_bounds", "deviationResult",
          function(object){object@variability_bounds})



setGeneric("differentialSamples", function(object, fdr = 0.1) standardGeneric("differentialSamples"))

setMethod("differentialSamples", "deviationResult",
          function(object, fdr = 0.1){
            fdrs = deviationFDR(object)
            sig = which(fdrs < fdr)
            up = intersect(sig, which(object@deviations > 0))
            down = intersect(sig, which(object@deviations < 0))
            return(list(up = up, down = down))
          })








