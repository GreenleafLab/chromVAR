#' @include deviationResult.R
NULL

#' A class to store deviationResult objects
#' @seealso \code{\link{metric}}, \code{\link{variability}}, \code{\link{variability_bounds}} 
deviationResultSet <- setClass("deviationResultSet", contains = "list")

setMethod("show",
          signature="deviationResultSet",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("deviationResults for ", length(object), " annotation sets. \n ", sep = "")
            invisible(NULL)
          })

#' @describeIn metric  returns name of metric used
#' @export
setMethod("metric", "deviationResultSet",
          function(object){object[[1]]@metric})


#' @describeIn variability  returns numeric vector of variability metric for all results in set
#' @export
setMethod("variability", "deviationResultSet",
          function(object){sapply(object, function(x) x@variability)})

#' @param lower return lower bound?  if FALSE, returns upper bound
#' @describeIn variability_bounds returns only lower or upper bound (based on lower argument) for all results in set
#' @export
setMethod("variability_bounds", "deviationResultSet",
          function(object, lower = T){
            if (lower){
              out = sapply(object, function(x) x@variability_bounds[1])
            } else{
              out = sapply(object, function(x) x@variability_bounds[2])
            }
            return(out)
          })

#' plot_variability
#' 
#' Plots variability for deviationResultSet
#' @param object deviationResultSet object
#' @param top number of top results to label
#' @param xlab label for x-axis (default is "sorted bins")
#' @return ggplot 
#' @seealso \code{\link{metric}}, \code{\link{variability}} 
setGeneric("plot_variability", function(object, ...) standardGeneric("plot_variability"))

#' @rdname plot_variability
#' @import ggplot2
#' @export
setMethod("plot_variability", "deviationResultSet",
          function(object, top = 3, xlab = "Sorted TFs"){

            res_df = data.frame(var = variability(object),
                                min = variability_bounds(object, lower = T),
                                max = variability_bounds(object, lower = F),
                                tf = names(object),
                                ranks = rank(-1 * variability(object)))

            ylab = ifelse(metric(object) == "z-score","SD of Z-scores", "Normalized Variability")
            
            out = ggplot2::ggplot(res_df, ggplot2::aes_string(x = "ranks", y = "var", min = "min", max = "max")) + ggplot2::geom_point()+ ggplot2::geom_errorbar() +
              ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))

            if (!is.na(top) && top > 0){
              top_df = res_df[res_df$ranks <= top,]
              out = out + ggplot2::geom_text(data = top_df, ggplot2::aes_string(x = "ranks", y = "var", label = "tf"), size = 2, hjust=-0.15,col = "Black")
            }

            return(out)

          })

