#' getFragmentsPerPeak
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return vector with sum across rows of counts assay within chromVARCounts
#' @export
#' @seealso \code{\link{getFragmentsPerSample}},
#' \code{\link{getTotalFragments}}
#' @examples 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' 
#' frags_per_peak <- getFragmentsPerPeak(mini_counts)
setGeneric("getFragmentsPerPeak", 
           function(object) standardGeneric("getFragmentsPerPeak"))

#' getFragmentsPerSample
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return vector with sum across columns of counts assay within chromVARCounts
#' @export
#' @seealso \code{\link{getFragmentsPerPeak}},
#' \code{\link{getTotalFragments}}
#' @examples 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' frags_per_sample <- getFragmentsPerSample(mini_counts)
setGeneric("getFragmentsPerSample", 
           function(object) standardGeneric("getFragmentsPerSample"))

#' getTotalFragments
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return sum of all counts within object
#' @export
#' @seealso \code{\link{getFragmentsPerSample}}, 
#' \code{\link{getFragmentsPerPeak}}
#' @examples 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' total_frags <- getTotalFragments(mini_counts)
setGeneric("getTotalFragments", 
           function(object) standardGeneric("getTotalFragments"))


#' @describeIn getFragmentsPerPeak method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("getFragmentsPerPeak", c(object = "SummarizedExperiment"), 
          function(object) {
            Matrix::rowSums(counts(object))
          })

#' @describeIn getFragmentsPerPeak method for Matrix or matrix object
#' @export
setMethod("getFragmentsPerPeak", c(object = "MatrixOrmatrix"), 
          function(object) {
            Matrix::rowSums(object)
          })


#' @describeIn getFragmentsPerSample method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("getFragmentsPerSample", c(object = "SummarizedExperiment"), 
          function(object) {
            Matrix::colSums(counts(object))
          })

#' @describeIn getFragmentsPerSample method for Matrix or matrix object
#' @export
setMethod("getFragmentsPerSample", c(object = "MatrixOrmatrix"), 
          function(object) {
            Matrix::colSums(object)
          })

#' @describeIn getTotalFragments method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("getTotalFragments", c(object = "SummarizedExperiment"), 
          function(object) {
            sum(counts(object))
          })

#' @describeIn getTotalFragments method for Matrix or matrix object
#' @export
setMethod("getTotalFragments", c(object = "MatrixOrmatrix"), 
          function(object) {
            sum(object)
          })


