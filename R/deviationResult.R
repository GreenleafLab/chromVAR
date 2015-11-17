
#' An S4 class to store chromVAAR result
#' 
#' @slot deviations 
#' @slot variability
#' @slot variability_bounds
#' @slot name
#' @slot metric
#' @slot p_deviations
#' @slot p_variability
#'
#' @aliases deviations, variability, variability_bounds, name, metric, pvalues
deviationResult <- setClass("deviationResult",
                            slots = c(deviations = 'numeric',
                                      variability = 'numeric',
                                      variability_bounds = 'numeric',
                                      name = 'character',
                                      metric = 'character',
                                      p_deviations = 'numeric',
                                      p_variability = 'numeric'))



setMethod("show",
          signature="deviationResult",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            if (nchar(object) >= 1){
              cat("Deviations for ", object@name, " for ",
                length(object@deviations), " cells/samples. \n", sep = "")
            } else{
              cat("Deviations for ",
                  length(object@deviations), " cells/samples. \n", sep = "")
            }
            invisible(NULL)
          })

#' metric
#' 
#' Accessor for which variability metric was used.
#' @param object either \code{\link{deviationResult}} or \code{\link{deviationResultSet}} object
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
#' @param object either \code{\link{deviationResult}} or \code{\link{deviationResultSet}} object
#' @seealso \code{\link{metric}}, \code{\link{variability_bounds}} 
setGeneric("variability", function(object, ...) standardGeneric("variability"))

#' @describeIn variability  returns single numeric value representing the variability metric for that result
#' @export
setMethod("variability", "deviationResult",
          function(object){object@variability})

#' variability_bounds
#' 
#' Accessor for variability bounds in deviationResult or deviationResultSet.
#' The bounds on the variability metric are computed by bootstrapping the 
#' samples and determining a 95 percent confidence interval
#' @param object either \code{\link{deviationResult}} or \code{\link{deviationResultSet}} object
#' @seealso \code{\link{metric}}, \code{\link{variability}} 
setGeneric("variability_bounds", function(object, ...) standardGeneric("variability_bounds"))

#' @describeIn variability_bounds returns numeric vector of length two, with first  
#' value the lower bound and second value the upper bound
#' @export
setMethod("variability_bounds", "deviationResult",
          function(object){object@variability_bounds})

#' get_pvalues
#' 
#' Accessor for p values
#' @param object either \code{\link{deviationResult}} or \code{\link{deviationResultSet}} object
#' @seealso \code{\link{metric}}, \code{\link{variability}}, \code{\link{deviations}}
setGeneric("get_pvalues", function(object, ...) standardGeneric("get_pvalues"))

#' @param what return p values for deviations or variability?
#' @describeIn get_pvalues returns p values for either deviations or variability
#' @export
setMethod("get_pvalues", "deviationResult",
          function(object, what = c("deviations","variability")){
            stopifnot(metric(object)!="old")
            what = match.arg(what)
            if (what == "deviations"){
              return(object@p_deviations)
            } else if (what == "variability"){
              return(object@p_variability)
            }
          })

#' differential_samples
#' 
#' get differential_samples
#' @param object \code{\link{deviationResult}} object
#' @seealso \code{\link{get_pvalues}}, \code{\link{metric}},
#' \code{\link{variability}}, \code{\link{deviations}}
setGeneric("differential_samples", function(object, ...) standardGeneric("differential_samples"))

#' @param cutoff p-value cutoff for differential samples
#' @rdname differential_samples
#' @export
setMethod("differential_samples", "deviationResult",
          function(object, cutoff = 0.05){
            fdrs = get_pvalues(object, what = "deviations")
            sig = which(fdrs <= cutoff)
            up = intersect(sig, which(object@deviations > 0))
            down = intersect(sig, which(object@deviations < 0))
            return(list(up = up, down = down))
          })








