

get_fragments_per_peak <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  rowSums(assays(input)$counts)
}

get_fragments_per_sample <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  colSums(assays(input)$counts)
}

get_total_fragments <- function(input){
  stopifnot(inherits(input,"SummarizedExperiment"))
  sum(assays(input)$counts)
}