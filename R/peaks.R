
#' read_peaks
#' 
#' Read in peaks from a bed file.  
#' @param filename filename of bed file
#' @param extra_cols extra columns to read in beyond first three
#' @return GenomicRanges containing peaks in bed file
#' @details As in standard definition of bed file, first column is assumed to be chromosome, 
#' second is assumed to be start of peak (0-based), and third is assumed to be end of peak (1-based).
#' Extra columns can be added as metadata or strand information if provided, but the user must 
#' indicate column index and name using named vector for extra_cols
#'  @export
read_peaks <- function(filename, extra_cols = c()){
  if (is.installed('readr')){
    bed <- readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]
  } else{
    bed <- read.delim(file = filename, col.names = FALSE, sep = "\t", stringsAsFactors = FALSE)[, c(1:3, extra_cols)]
  }
  colnames(bed) <- c("chr", "start", "end", names(extra_cols))
  bed[,"end"] <- bed[,"end"] - 1
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
#   if (!isDisjoint(bed)){
#     stop("Peaks should be non-overlapping!
#          Use read_top_peaks to read peaks & select ")
#   }
  if (sum(width(bed) == width(bed[1])) != length(bed)){
    warning('Peaks are not equal width! 
            Use resize(peaks, width = x, fix = "center") to make peaks equal in size, 
            where x is the desired size of the peaks)')
  }
  return(bed)
}


#' reduce_peaks
#' 
#' Function to eliminate overlaping peaks from counts_mat
#' @param counts_mat \code{\link{fragmentCounts}}
#' @details This function will look for overlapping peaks and remove overlapping peaks
#' that have fewer fragments than the peaks they overlap
#' @return \code{\link{fragmentCounts}} with only non-overlapping peaks
reduce_peaks <- function(counts_mat){
  strand(counts_mat@peaks) <- "*"
  if (!isTRUE(all.equal(counts_mat@peaks, sort(counts_mat@peaks)))){
    counts_mat <- counts_mat[GenomicRanges::order(counts_mat@peaks),]
  }
  while (!(isDisjoint(counts_mat@peaks))){
    chr_names = seqnames(counts_mat@peaks)
    starts = start(counts_mat@peaks)
    ends = end(counts_mat@peaks)
  
    overlap_next = intersect(which(chr_names[1:(counts_mat@npeak -1)] == chr_names[2:(counts_mat@npeak)]),
                             which(ends[1:(counts_mat@npeak -1)] > starts[2:(counts_mat@npeak)]))
    overlap_previous = overlap_next + 1
    overlap_comparison = counts_mat@fragments_per_peak[overlap_previous] > counts_mat@fragments_per_peak[overlap_next]
    discard = c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])
    keep = which(1:counts_mat@npeak %ni% discard)
    counts_mat <- counts_mat[keep,]
  }
  return(counts_mat) 
}



