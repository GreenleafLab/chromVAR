#' get_fragments_per_peak
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return vector with sum across rows of counts assay within chromVARCounts
#' @export
setGeneric("get_fragments_per_peak", 
           function(object) standardGeneric("get_fragments_per_peak"))

#' get_fragments_per_sample
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return vector with sum across columns of counts assay within chromVARCounts
#' @export
setGeneric("get_fragments_per_sample", 
           function(object) standardGeneric("get_fragments_per_sample"))

#' get_total_fragments
#' 
#' @param object SummarizedExperiment, matrix, or Matrix object
#' @return sum of all counts within object
#' @export
setGeneric("get_total_fragments", 
           function(object) standardGeneric("get_total_fragments"))


#' @describeIn get_fragments_per_peak method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("get_fragments_per_peak", c(object = "SummarizedExperiment"), 
          function(object) {
            Matrix::rowSums(counts(object))
          })

#' @describeIn get_fragments_per_peak method for Matrix or matrix object
#' @export
setMethod("get_fragments_per_peak", c(object = "MatrixOrmatrix"), 
          function(object) {
            Matrix::rowSums(object)
          })


#' @describeIn get_fragments_per_sample method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("get_fragments_per_sample", c(object = "SummarizedExperiment"), 
          function(object) {
            Matrix::colSums(counts(object))
          })

#' @describeIn get_fragments_per_sample method for Matrix or matrix object
#' @export
setMethod("get_fragments_per_sample", c(object = "MatrixOrmatrix"), 
          function(object) {
            Matrix::colSums(object)
          })

#' @describeIn get_total_fragments method for SummarizedExperiment object 
#' with counts assay
#' @export
setMethod("get_total_fragments", c(object = "SummarizedExperiment"), 
          function(object) {
            sum(counts(object))
          })

#' @describeIn get_total_fragments method for Matrix or matrix object
#' @export
setMethod("get_total_fragments", c(object = "MatrixOrmatrix"), 
          function(object) {
            sum(object)
          })


