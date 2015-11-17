#' @include fragmentCounts.R
NULL


#' compute_bias
#' 
#' Computes GC content for peaks
#' @param object GenomicRanges or fragmentCounts
setGeneric("compute_bias", function(object, ...) standardGeneric("compute_bias"))

#' @param genome BSgenome object, e.g. BSgenome.Hsapiens.UCSC.hg19
#' @return If object is GenomicRanges, returns GenomicRanges with additional meta-data
#'  column named "gc" with gc content
#' @rdname compute_bias
setMethod("compute_bias","GenomicRanges",
          function(object, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
            seqs = Biostrings::getSeq(genome, object)
            nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
            gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
            S4Vectors::mcols(object)$bias <- gc
            return(object)
          })

#' @inheritParams compute_bias
#' @return If object is fragmentCounts, returns fragmentCounts with additional meta-data
#'  column for peaks slot named "gc" with gc content
#' @rdname compute_bias
setMethod("compute_bias","fragmentCounts",
          function(object, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
            seqs = Biostrings::getSeq(genome, object@peaks)
            nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
            gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
            S4Vectors::mcols(object@peaks)$bias <- gc
            return(object)
          })

#' @keywords internal
#' @export
make_bias_bins <- function(counts_mat, nbins = 25){
  npeaks = length(counts_mat@peaks)
  #make bias bins
  bias_quantiles = quantile(counts_mat@peaks$bias, seq(0,1,1/nbins))
  bias_cut = cut(counts_mat@peaks$bias, breaks = bias_quantiles)
  bias_bins = split(1:npeaks, bias_cut)
  names(bias_bins) = sapply(1:nbins, function(x) paste("bias_bin_",x,sep="",collapse=""))
  #make count bins
  pseudo_counts = counts_mat@fragments_per_peak + runif(npeaks,min = 0, max = 0.1)
  count_quantiles = quantile(pseudo_counts, seq(0,1,1/nbins))
  count_cut = cut(pseudo_counts, breaks = count_quantiles)
  count_bins = split(1:npeaks, count_cut)
  names(count_bins) = sapply(1:nbins, function(x) paste("count_bin_",x,sep="",collapse=""))
  #make bias / count bins
  nbins = round(sqrt(nbins))
  bias_quantiles = quantile(counts_mat@peaks$bias, seq(0,1,1/nbins))
  bias_cut = cut(counts_mat@peaks$bias, breaks = bias_quantiles)
  tmp_bias_bins = split(1:npeaks, bias_cut)
  count_quantiles = quantile(pseudo_counts, seq(0,1,1/nbins))
  count_cut = cut(pseudo_counts, breaks = count_quantiles)
  tmp_count_bins = split(1:npeaks, count_cut)
  bias_count_bins = sapply(1:nbins, function(x) sapply(1:nbins, function(y) intersect(tmp_bias_bins[[y]], tmp_count_bins[[x]])))
  names(bias_count_bins) = sapply(1:nbins, function(x) sapply(1:nbins, function(y) paste("bias_count_bin_",x,"_",y,sep="",collapse="")))
  return(c(bias_bins, count_bins, bias_count_bins))
}

#' @keywords internal
#' @export
make_permuted_sets <- function(counts_mat, motif_indices, window = 10, BPPARAM = BPPARAM){
  bg <- getBackgroundPeakSets(counts_mat, niterations = 1, window = window, BPPARAM = BPPARAM)
  sets <- lapply(1:length(motif_indices), function(x) bg@background_peaks[motif_indices[[x]],1])
  names(sets) <- sapply(names(motif_indices), function(x) paste("permuted_",x,collapse=""))
  return(sets)
}

