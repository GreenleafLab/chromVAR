# utility classes

setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
setClassUnion("missingOrNULL", c("missing", "NULL"))


#' chromVARDeviations
#' 
#' Class for storing results from \code{\link{computeDeviations}} function.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @details This class inherits from 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and most methods 
#' for that class should work for objects of this class as well. Additionally, 
#' two accessor functions are defined for extracting bias corrected deviations 
#' (\code{\link{deviations}}) and deviation Z-scores 
#' (\code{\link{deviationScores}})
setClass("chromVARDeviations", contains = "SummarizedExperiment")


setValidity("chromVARDeviations",
            function(object) {
              if ("deviations" %ni% assayNames(object) || "z" %ni% 
                  assayNames(object)) 
                return("The assays slot must contain 'deviations' and 'z'")
              return(TRUE)
            })