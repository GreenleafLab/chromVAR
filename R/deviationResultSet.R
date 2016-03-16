#' @include deviationResult.R
NULL

# Class Definition and basic methods -------------------------------------------

#' deviationResultSet
#' 
#' A class to store deviationResult objects
#' @seealso \code{\link{metric}}, \code{\link{variability}}, 
#' \code{\link{variability_bounds}} 
deviationResultSet <- setClass("deviationResultSet", 
                               slots = c(names = 'character',
                                         nresult = 'numeric') ,
                               contains = "list")

setValidity("deviationResultSet", function(object) {
  msg <- NULL
  valid <- TRUE
  #Check that all elements are deviationResult
  if (!all_true(sapply(object, function(x) inherits(x,"deviationResult")))){
    valid <- FALSE
    msg <- c(msg, "All elements of deviationResultSet must be of class deviationResult\n")
  } 
  if (!all_false(duplicated(names(object)))){
    valid <- FALSE
    msg <- c(msg, "Names must be unique")
  }
  if (valid) TRUE else msg})



#' @rdname deviationResultSet-class
setMethod("initialize",
          "deviationResultSet",
          function(.Object, ...){
            .Object <- callNextMethod()
            if (length(.Object@nresult) == 0){
              .Object@nresult = length(.Object)
            }
            .Object@.Data <- lapply(.Object@.Data, function(x) set_nresult(x, .Object@nresult))
            .Object
          })

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
            nresult <- x@nresult
            tmp <- callNextMethod()
            if (length(tmp) ==0){
              return(NULL)
            } else{
              out <- deviationResultSet(tmp, nresult = nresult)
              return(out)              
            }
          })

#' @section c:
#' Using 'c' to combine deviationResultSet will reset the nresult slot,
#' affecting multiple hypothesis testing...
#' @rdname deviationResultSet-class
#' @export
setMethod("c", signature = signature(x = "deviationResultSet"),
          definition = function(x, ...) {
            if (!identical(recursive, FALSE)) 
              stop("\"c\" method for deviationResultSet objects ", "does not support the 'recursive' argument")
            if (missing(x)) {
              args <- list(...)
              x <- args[[1L]]
            }
            else {
              args <- list(x, ...)
            }
            names <- unlist(lapply(args,function(obj) obj@names), recursive=FALSE)
            results <- unlist(lapply(args,function(obj) obj@.Data), recursive = FALSE)
            out <- deviationResultSet(results, names = names)
            out})   
            

setMethod("set_nresult", "deviationResultSet",
          function(object, nresult){
            object@nresult <- nresult
            object@.Data <- lapply(object@.Data, function(x) set_nresult(x, nresult))
            return(object)
          })



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

#' @describeIn foldchange returns matrix of fold changes, with columns representing 
#' samples and rows representing annotation sets
#' @export
setMethod("foldchange", "deviationResultSet",
          function(object){
            out <- t(sapply(object, foldchange))
            colnames(out) = names(object[[1]]@foldchange)
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
#' @export
setMethod("get_pvalues", "deviationResultSet",
          function(object, adjust = TRUE){
            out <- sapply(object, get_pvalues, what = "variability")
            if (adjust){
              out <- p.adjust(out, method = "BH", n = object@nresult)
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
    stopifnot(cutoff >= 1)
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
  
  g = cluster::clusGap(t(mat), tmpfun, K.max = max_k)
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

#' cluster_samples
#' 
#' function to custer samples
#' @param object \code{\link{deviationResultSet}} object
#' @param k number of clusters 
#' @return list with two elements.  "hclust" is \code{\link[stats]{hclust}} object
#' "cluster" is cluster membership vector.
#' @export
cluster_samples <- function(object, k = 5, plot = TRUE){
  
  #Perform checks on arguments
  stopifnot(inherits(object,"deviationResultSet"))
  stopifnot(k >= 1)
  
  mat = deviations(object)

  d = as.dist(1 - cor(mat))
  h = hclust(d)
  cl = cutree(h, k)
  
  out = list("hclust" = h, "cluster" = cl)
  
  return(out)
  
}
  








