
#' get_fragments_per_peak
#' 
#' @param input SummarizedExperiment object
#' @return vector with sum across rows of counts assay within SummarizedExperiment
#' @export
get_fragments_per_peak <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  rowSums(assays(input)$counts)
}

#' get_fragments_per_peak
#' 
#' @param input SummarizedExperiment object
#' @return vector with sum across columns of counts assay within SummarizedExperiment
#' @export
get_fragments_per_sample <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  colSums(assays(input)$counts)
}

#' get_total_fragments
#' 
#' @param input SummarizedExperiment object
#' @return sum of all counts within input
#' @export
get_total_fragments <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  sum(assays(input)$counts)
}