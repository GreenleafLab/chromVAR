#' @include deviationResult.R
NULL

# Class Definition and basic methods -------------------------------------------

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

#' @describeIn deviations returns matrix of deviations, with columns representing 
#' samples and rows representing annotation sets
#' @export
setMethod("deviations", "deviationResultSet",
          function(object){
            out <- t(sapply(object, deviations))
            colnames(out) = names(object[[1]]@deviations)
            out
            })

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

# Plotting ---------------------------------------------------------------------            

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
                label_option = seq_along(object)
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



#' subset_by_variability
#' 
#' function to isolate results with high variability
#' @param object \code{\link{deviationResultSet}} object
#' @param by choice of "pvalue", "top", "metric", or "bound".  See Details.
#' @param cutoff -- cutoff for top sites. See Details.
#' @param adjusted if by is set to "pvalue", adjust for multiple comparisons?
#' @return \code{\link{deviationResultSet}} object
#' @details Subsetting by variability can be done in sevarl ways, determined by 
#' the options "by" and "cutoff".  If "by" is set to "pvalue" (the default),
#' then the "cutoff" is an upper bound cutoff for the pvalue.  If "by" is set to 
#' "top", then the top K most variable are returned, where K is set by "cutoff".  
#' If "by" is "metric" then subsetting is based on the value of the variability 
#' metric, and "cutoff" is a lower bound on that metric.  If "by" is "bound" then
#' subsetting is based on the value of the lower bound for the variability metric,
#' and "cutoff" is the lower bound for the lower bounds of the variability metric.
#' @export
subset_by_variability <- function(object, by = c("pvalue","top","metric","bound"), cutoff = 0.05, adjusted = TRUE){
  
  #Perform checks on arguments
  stopifnot(inherits(object,"deviationResultSet"))
  stopifnot(is.numeric(cutoff))
  by = match.arg(by)

  if (by =="pvalue"){
    stopifnot(cutoff > 0 && cutoff < 1)
    pvals = get_pvalues(object, adjust = adjusted)
    object <- object[which(pvals <= cutoff)]
  } else if (by == "top"){
    stopifnot(cutoff > 2)
    object <- object[which(rank(-1 * variability(object),
                                ties.method="random") <= cutoff)]
  } else if (by == "metric"){
    object <- object[which(variability(object) >= cutoff)]
  } else if (by == "bound"){
    object <- object[which(variability_bound(object, lower = TRUE) >= cutoff)]
  }
  
  return(object)}
  


# Clustering Sets --------------------------------------------------------------

get_cluster_number <-function(object, 
                              max_k = 15, 
                              out = c("plot","k"),
                              method = "globalSEmax"){
  
  out = match.arg(out)
  
  #Perform checks on arguments
  stopifnot(inherits(object,"deviationResultSet"))
  stopifnot(as.integer(max_k) == max_k && max_k > 1)

  mat = deviations(object)
  
  tmpfun <- function(x, k){
    d = as.dist(1-cor(t(x)))
    h = hclust(d)
    cl = cutree(h, k)
    return(list("cluster" = cl))
  }
  
  g = cluster::clusGap(mat, tmpfun, K.max = max_k)
  rec = cluster::maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"], method=method) 
  
  if (out == "plot"){
    
    tmp = data.frame(k = 1:max_k,
                     Gap = g$Tab[,"gap"],
                     ymax =  g$Tab[,"gap"] +  g$Tab[,"SE.sim"],
                     ymin = g$Tab[,"gap"] -  g$Tab[,"SE.sim"],
                     rec = (1:max_k == rec))
    
    tmp2 = data.frame(label = paste0("k = ",rec), k = tmp[tmp$rec,"k"],
                      Gap = tmp[tmp$rec,"ymax"] + 0.1 *(max(tmp$ymax) - min(tmp$ymin)),
                      ymax = tmp[tmp$rec,"ymax"], ymin = tmp[tmp$rec,"ymin"], rec = TRUE,
                      stringsAsFactors = FALSE)
    print(tmp2)
    p = ggplot2::ggplot(tmp,
                    ggplot2::aes_string(x = "k", y = "Gap", ymax = "ymax", 
                                        ymin = "ymin", col = "rec")) +
      ggplot2::geom_point() + ggplot2::geom_errorbar() + 
      geom_text(data = tmp2, 
                ggplot2::aes_string(x = "k", y = "Gap", label = "label"),col="red") +
      scale_color_manual(values = c("black","red")) + 
      theme(legend.position="none")
    
    print(p)
    invisible(list("plot" = p, "value" = rec))
  } else{
    return(rec)
  }  
}

#' cluster_sets
#' 
#' function to custer the top annotation sets
#' @param object \code{\link{deviationResultSet}} object
#' @param k number of clusters -- if not provided, will be estimated. See Details
#' @param max_k if k is not provided, max_k allowed
#' @return list with two elements.  "hclust" is \code{\link[stats]{hclust}} object
#' "cluster" is cluster membership vector.
#' @export
cluster_sets <- function(object, k = NULL, max_k = 15, plot = TRUE){
  
  #Perform checks on arguments
  stopifnot(inherits(object,"deviationResultSet"))
  stopifnot(as.integer(max_k) == max_k && max_k > 1)
  
  if (is.null(k)){
    cn = get_cluster_number(object, max_k = max_k,  out = ifelse(plot, "plot","value"))
    if (plot){
      k = cn$value
    } else{
      k = cn 
    }
  }
  
  mat = deviations(object)

  d = as.dist(1 - cor(t(mat)))
  h = hclust(d)
  cl = cutree(h, k)
  
  if (plot){
    out = list("hclust" = h, "cluster" = cl, "plot" = cn$plot)
  } else{
    out = list("hclust" = h, "cluster" = cl)
  }
  
  return(out)
  
}
  









