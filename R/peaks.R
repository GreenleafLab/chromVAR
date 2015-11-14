
#' read_peaks
#' 
#' Read in peaks from a bed file.  
#' @param filename filename of bed file
#' @param extra_cols extra columns to read in beyond first three
#' @param extra_names extra names corresponding to extra columns
#' @return GenomicRanges containing peaks in bed file
#' @details As in standard definition of bed file, first column is assumed to be chromosome, 
#' second is assumed to be start of peak (0-based), and third is assumed to be end of peak (1-based).
#' Extra columns can be added as metadata or strand information if provided, but the user must 
#' indicate column index and name
#'  @export
read_peaks <- function(filename, extra_cols = c(), extra_names = c()){
  stopifnot(length(extra_cols) == length(extra_names))
  if (is.element('readr',installed.packages()[,1])){
    bed <- readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]
  } else{
    bed <- read.delim(file = filename, col.names = FALSE, sep = "\t", stringsAsFactors = FALSE)[, c(1:3, extra_cols)]
  }
  colnames(bed) <- c("chr", "start", "end", extra_names)
  bed[,"end"] <- bed[,"end"] - 1
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  if (sum(width(bed) == width(bed[1])) != length(bed)){
    stop("All peaks in bed file must be of equal width!")
  }
  return(bed)
}
