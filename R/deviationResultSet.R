#' @include deviationResult.R
NULL

#' deviationResultSet
#' 
#' A class to store deviationResult objects
#' @seealso \code{\link{metric}}, \code{\link{variability}}, 
#' \code{\link{variability_bounds}} 
deviationResultSet <- setClass("deviationResultSet", contains = "list")

setMethod("show",
          signature="deviationResultSet",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("deviationResults for ", length(object), " annotation sets. \n ", 
                sep = "")
            invisible(NULL)
          })

#' @rdname deviationResultSet-class
setMethod("[", signature = signature(x = "deviationResultSet"),
          definition = function(x, ...) {
            tmp <- callNextMethod()
            out <- deviationResultSet(tmp)
            return(out)
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


#' @describeIn get_pvalues returns p values for variability for each set
#' @param adjust adjust p-values for multiple testing
#' @export
setMethod("get_pvalues", "deviationResultSet",
          function(object, adjust = TRUE){
            stopifnot(metric(object)!="old")
            out <- sapply(object, get_pvalues, what = "variability")
            if (adjust){
              out <- p.adjust(out, method = "BH")
            }
            return(out)
          })

            

#' plot_variability
#' 
#' Plots variability for deviationResultSet
#' @return ggplot 
#' @seealso \code{\link{metric}}, \code{\link{variability}} 
setGeneric("plot_variability", function(object, ...) standardGeneric("plot_variability"))

#' @param object \code{\link{deviationResultSet}} object
#' @param xlab label for x-axis (default is "sorted bins")
#' @param label_type which points to label? See Details.
#' @param label_option See Details.
#' @param only FALSE, p-value cutoff, or number of points.  See Details.
#' @details Labelling of points is controlled by label_type and label_options.  
#' label_type can be either "top", "names", or "none".  If "top", the top k most 
#' variable points are labelled with their names. k can be set by passing an integer
#' value to label_option.  By default, k will be 3.  Specific points chosen by 
#' the user can be labelled by setting label_type to "name".  The names of the points 
#' to label should then be passed to label_option.  Either character vectors of the 
#' names can be given, or a numeric vector of the indices. If no input is given 
#' to label_option when label_type is "name", all points are labelled. If 
#' label_type is "none", no labelling is performed.\cr \cr
#' Rather than plotting all results, you can also plot
#' only significant results, a certain number of top results, or results by index
#' or name.  For plotting only significant results, pass a numeric greater than
#' 0 and less than 1 to the argument only.  For plotting only the top points, 
#' pass a numeric value greater than 1 to the argument only.  For plotting results
#' by index or name, give a vector of indices or names.  If only is FALSE 
#' (the default), then all results are plotted.  
#' @rdname plot_variability
#' @import ggplot2
#' @export
setMethod("plot_variability", "deviationResultSet",
          function(object, xlab = "Sorted TFs", 
                   label_type = c("top","names","none"), 
                   label_option = NULL, 
                   only = FALSE){
            
            label_type = match.arg(label_type)

            if (only != FALSE){
              if (length(only) > 1){
                object <- object[only]
              } else if (is.numeric(only) && only > 0 && only < 1){
                object <- object[which(get_pvalues(object) <= only)]
              } else if (is.numeric(only) && only > 1){
                object <- object[which(rank(-1 * variability(object),
                                            ties.method="random") <= only)]
              } else{
                message("Invalid option to only, plotting all results. See 
                        documentation for details on usage of only")
              }
            }
            
            res_df = data.frame(var = variability(object),
                                min = variability_bounds(object, lower = TRUE),
                                max = variability_bounds(object, lower = FALSE),
                                tf = names(object),
                                ranks = rank(-1 * variability(object),
                                             ties.method="random"))

            ylab = ifelse(metric(object) == "z-score",
                          "SD of Z-scores", 
                          "Normalized Variability")
            
            out = ggplot2::ggplot(res_df, ggplot2::aes_string(x = "ranks", 
                                                              y = "var",
                                                              min = "min", 
                                                              max = "max")) + 
              ggplot2::geom_point()+ ggplot2::geom_errorbar() +
              ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + 
              ggplot2::scale_y_continuous(expand=c(0,0),
                                          limits=c(0,max(res_df$max)*1.05))

            if (label_type =="top"){
              if (!is.numeric(label_option)){
                label_option = 3
              }
              top_df = res_df[res_df$ranks <= label_option,]
              out = out + ggplot2::geom_text(data = top_df, 
                                             ggplot2::aes_string(x = "ranks", 
                                                                 y = "var", 
                                                                 label = "tf"),
                                             size = 2, hjust=-0.15,col = "Black")
            } else if (label_type =="name"){
              if (is.null(label_option)){
                label_option = 1:length(object)
              } else if (is.character(label_option)){
                label_option = which(names(object) %in% label_option)
              } else if (!is.numeric(label_option)){
                stop("Labels must be numeric or character vector")
              }
              label_df = res_df[label_option,]
              out = out + ggplot2::geom_text(data = label_df, 
                                             ggplot2::aes_string(x = "ranks", 
                                                                 y = "var",
                                                                 label = "tf"), 
                                             size = 2, hjust=-0.15,col = "Black")
              
            }

            return(out)

          })

