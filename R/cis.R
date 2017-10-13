#' getCisGroups
#' 
#' Function for grouping peaks based on proximity along chromosomes
#' @param object GenomicRanges or RangedSummarizedExperiment
#' @param grpsize number of peaks to include in each grouop
#' @param stepsize number of peaks between each new set of groups
#' @param ... additional arguments
#' @export
#' @author Alicia Schep
#' @return SummarizedExperiment with annotationMatches assay storing which peaks
#' belong to which groups
#' @examples 
#' 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' mini_counts <- sort(mini_counts)
#' cisg <- getCisGroups(mini_counts)
setGeneric("getCisGroups", 
           function(object, ...) standardGeneric("getCisGroups"))


#' @describeIn getCisGroups method for RangedSummarizedExperiment
#' @export
setMethod("getCisGroups", c(object = "RangedSummarizedExperiment"), 
          function(object, grpsize = 25, stepsize = 10) {
            get_cis_groups_core(rowRanges(object), grpsize, stepsize)
          })

#' @describeIn getCisGroups method for GenomicRanges
#' @export
setMethod("getCisGroups", c(object = "GenomicRanges"), 
          function(object, grpsize = 25, stepsize = 10) {
            get_cis_groups_core(object, grpsize, stepsize)
          })

get_cis_groups_core <- function(peaks, grpsize = 25, stepsize = 10) {
  
  if (!isSorted(peaks)) {
    stop("peaks must be sorted to be able to get cis groups")
  }
  
  chrs <- seqlevels(peaks)
  out <- list()
  out <- do.call(c, bplapply(seq_along(chrs), function(i) {
    chr_ix <- which(as.character(seqnames(peaks)) == chrs[i])
    if (length(chr_ix) > stepsize) {
      tmp <- lapply(seq_len(length(chr_ix)%/%stepsize),  function(x) { 
          start_ix <- ((x - 1) * stepsize + 1)
          end_ix <- min((x - 1) * stepsize + grpsize, length(chr_ix))
          chr_ix[start_ix:end_ix]
          }
        )
      names(tmp) <- vapply(seq_along(tmp), 
                           function(x) 
                             paste(chrs[i], x, sep = "_", collapse = ""), 
                           "")
      return(tmp)
    } else {
      return(list())
    }
  }))
  out <- getAnnotations(out, rowRanges = peaks)
  return(out)
}

