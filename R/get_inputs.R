################# Functions for reading in inputs ##############################

# Read in peaks ----------------------------------------------------------------

#' getPeaks
#'
#' Read in peaks from a bed file.
#' @param filename filename of bed file
#' @param extra_cols extra columns to read in beyond first three
#' @param sort_peaks sort the peaks?
#' @return \code{\link[GenomicRanges]{GenomicRanges}} containing peaks from file
#' @details As in standard definition of bed file, first column is assumed to be
#'  chromosome, second is assumed to be start of peak (0-based), and third is
#'  assumed to be end of peak (1-based). Note that in output GenomicRanges
#'  output, start and end indices are both 1-based. Extra columns can be added
#'  as metadata or strand information if provided, but the user must indicate 
#'  column index and name using named vector for extra_cols.
#' @seealso \code{\link{getCounts}}, \code{\link{filterPeaks}}, 
#' \code{\link{readNarrowpeaks}}
#' @export
#' @examples 
#' peaks_file <- system.file("extdata", "test_bed.txt", package = "chromVAR")
#' peaks <- getPeaks(peaks_file, sort = TRUE)
getPeaks <- function(filename, extra_cols = c(), sort_peaks = FALSE) {
  if (is.installed("readr")) {
    bed <- as.data.frame(
      suppressMessages(readr::read_tsv(file = filename, 
                                       col_names = FALSE)[,c(1:3, extra_cols)]))
  } else {
    bed <- read.delim(file = filename, header = FALSE, sep = "\t", 
                      stringsAsFactors = FALSE)[, c(1:3, extra_cols)]
  }
  colnames(bed) <- c("chr", "start", "end", names(extra_cols))
  bed[, "start"] <- bed[, "start"] + 1
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  if (!isDisjoint(bed)) {
    warning("Peaks are overlapping!",
      "After getting peak counts, peaks can be reduced to non-overlapping set",
      " using filterPeaks function")
  }
  if (sum(width(bed) == width(bed[1])) != length(bed)) {
    warning("Peaks are not equal width!",
      "Use resize(peaks, width = x, fix = \"center\") to make peaks equal in ",
      "size, where x is the desired size of the peaks)")
  }
  #sorted_bed <- sortSeqlevels(bed)
  #sorted_bed <- sort(sorted_bed, ignore.strand = TRUE)
  if (sort_peaks) {
    if (!isSorted(bed)) {
      message("Peaks sorted")
    }
    sorted_bed <- sortSeqlevels(bed)
    sorted_bed <- sort(sorted_bed, ignore.strand = TRUE)
    return(sorted_bed)
  } else {
    if (!isSorted(bed)) {
      warning("Peaks not sorted")
    }
    return(bed)
  }
}

#' readNarrowpeaks
#' 
#' Reads in peaks in narrowpeaks format, as output by macs2. Uses summit as 
#' center of peak, and makes peak the given 'width'. By default removes 
#' overlapping peaks to get set of peaks with no overlaps
#' @param filename filename
#' @param width desired width of peaks
#' @param non_overlapping remove overlapping peaks
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @export
readNarrowpeaks <- function(filename, 
                                   width = 500, 
                                   non_overlapping = TRUE) {
  cn <- c("chr", "start", "end", "name", "score", "strand", "fc", 
          "pval", "qval", "summit")
  if (is.installed("readr")) {
    bed <- as.data.frame(readr::read_tsv(file = filename, col_names = cn))
  } else {
    bed <- read.delim(file = filename, 
                      header = FALSE, 
                      sep = "\t", 
                      stringsAsFactors = FALSE, 
                      col.names = cn)
  }
  bed[, "summit"] <- bed[, "start"] + bed[, "summit"]
  bed <- as(bed, "DataFrame")
  bed <- makeGRangesFromDataFrame(bed[, c("chr", "summit", "score", 
                                          "qval", "name")], 
                                  start.field = "summit", end.field = "summit", 
                                  keep.extra.columns = TRUE)
  bed <- resize(bed, width = width, fix = "center")
  if (non_overlapping) {
    bed <- sortSeqlevels(bed)
    bed <- sort(bed)
    keep_peaks <- seq_along(bed)
    while (!(isDisjoint(bed[keep_peaks]))) {
      chr_names <- as.character(seqnames(bed[keep_peaks]))
      starts <- start(bed[keep_peaks])
      ends <- end(bed[keep_peaks])
      overlap_next <- 
        intersect(which(chr_names[seq_len(length(keep_peaks) - 1)] == 
                          chr_names[seq_len(length(keep_peaks) - 1) + 1]),
                  which(ends[seq_len(length(keep_peaks) - 1)] >= 
                          starts[seq_len(length(keep_peaks) - 1) + 1]))
      overlap_previous <- overlap_next + 1
      overlap_comparison <- bed[keep_peaks[overlap_previous]]$qval > 
        bed[keep_peaks[overlap_next]]$qval
      discard <- keep_peaks[c(overlap_previous[!overlap_comparison], 
                              overlap_next[overlap_comparison])]
      keep_peaks <- keep_peaks[keep_peaks %ni% discard]
    }
    bed <- bed[keep_peaks]
  }
  return(bed)
}

# Main function for reading in counts from bam ---------------------------------

#' Accessors for the 'counts' slot of a SummarizedExperiment
#'
#' @docType methods
#' @rdname counts
#' @param object SummarizedExperiment object
#' @aliases counts,SummarizedExperiment-method
#' @importMethodsFrom BiocGenerics counts
#' @export
#' @return Matrix of counts
#' @examples 
#' data(mini_counts, package = "chromVAR")
#' fragment_counts <- counts(mini_counts)
setMethod("counts", signature(object = "SummarizedExperiment"), 
          function(object) {
  stopifnot("counts" %in% assayNames(object))
  assays(object)$counts
})


#' @rdname counts
#' @param value matrix of counts
#' @aliases counts<-,SummarizedExperiment-method
#' @importMethodsFrom BiocGenerics counts<-
#' @exportMethod 'counts<-'
setReplaceMethod("counts", signature(object = "SummarizedExperiment", 
                                     value = "MatrixOrmatrix"), 
                 function(object, value) {
                   assays(object)[["counts"]] <- value
                   validObject(object)
                   object
                 })



#' getCounts
#'
#' makes matrix of fragment counts in peaks using one or multiple bam or bed 
#' files
#' @param alignment_files filenames for bam or bed files with aligned reads
#' @param peaks GRanges object with peaks
#' @param by_rg use RG tags in bam to separate groups?
#' @param paired paired end data?
#' @param format bam or bed?  default is bam
#' @param colData sample annotation DataFrame
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#'  object
#' @seealso \code{\link{getSampleDepths}},  \code{\link{getPeaks}}, 
#' \code{\link{filterSamples}}
#' @export
#' @examples 
#' 
#' # First we'll read in some peaks
#' peaks_file <- system.file("extdata", "test_bed.txt", package = "chromVAR")
#' test_peaks <- getPeaks(peaks_file, sort = TRUE)
#' 
#' # With single bam with RG tags (can also give multiple bams with RG)
#' test_rg <- system.file("extdata", "test_RG.bam", package = "chromVAR")
#' test_counts <- getCounts(test_rg, peaks = test_peaks, by_rg = TRUE, 
#'                           paired = TRUE, 
#'                           colData = S4Vectors::DataFrame(condition ="A"))
#'                           
#' 
#' # Multiple bams without RG tags
#' test_bam1 <- system.file("extdata", "test_single1.bam", package = "chromVAR")
#' test_bam2 <- system.file("extdata", "test_single2.bam", package = "chromVAR")
#' test_bam3 <- system.file("extdata", "test_single3.bam", package = "chromVAR")
#' test_counts2 <- getCounts(c(test_bam1, test_bam2,test_bam3), 
#'                            peaks = test_peaks, by_rg = FALSE, 
#'                            paired = TRUE, 
#'                            colData = S4Vectors::DataFrame(celltype = 
#'                                                           c("A","B","C"))) 
#'                            
#' # Bed file with reads (can give multiple bed files, here we will just read 1)
#' test_bed <- system.file("extdata", "test_reads.bed", package = "chromVAR")
#' test_counts3 <- getCounts(test_bed, test_peaks, by_rg = FALSE, 
#'                            paired = FALSE, 
#'                            format = "bed")                           
getCounts <- function(alignment_files, peaks, paired, by_rg = FALSE, 
                       format = c("bam", "bed"), colData = NULL) {
  
  format <- match.arg(format)
  if (format == "bam") {
    return(get_counts_from_bams(alignment_files, peaks, paired, by_rg, colData))
  } else {
    return(get_counts_from_beds(alignment_files, peaks, paired, colData))
  }
}


get_counts_from_bams <- function(bams, peaks, paired, by_rg = FALSE,
                                 sample_annotation = NULL) {
  
  if (by_rg) {
    tmp <- lapply(bams, getFragmentCountsByRG, peaks = peaks, paired = paired)
    if (!is.null(sample_annotation) && nrow(sample_annotation) == length(bams)){
      sample_annotation <- as(sample_annotation, "DataFrame")
      l <- vapply(tmp, function(x) length(x$depths), 0)
      sample_annotation <- 
        do.call(rbind, 
                lapply(seq_along(l), 
                       function(x){
                         rep(sample_annotation[x, ,drop = FALSE], l[x])}))
    }
    counts_mat <- do.call(cbind, lapply(tmp, function(x) x$counts))
    depths <- do.call(c, lapply(tmp, function(x) x$depths))
  } else {
    mat <- matrix(nrow = length(peaks), ncol = length(bams))
    depths <- vector("numeric", length(bams))
    
    for (i in seq_along(bams)) {
      message("Reading in file: ", bams[i])
      fragments <- bamToFragments(bams[i], paired = paired)
      depths[i] <- length(fragments)
      mat[, i] <- countOverlaps(peaks, fragments, type = "any", 
                                ignore.strand = TRUE)
    }
    colnames(mat) <- basename(bams)
    counts_mat <- Matrix::Matrix(mat)
  }
  if (is.null(sample_annotation)) {
    sample_annotation <- DataFrame(depth = depths)
  } else {
    sample_annotation$depth <- depths
  }
  out <- SummarizedExperiment(assays = list(counts = counts_mat), 
                              rowRanges = peaks, 
                              colData = sample_annotation)
  return(out)
}


get_counts_from_beds <- function(beds, peaks, paired, colData = NULL) {
  
  
  results <- bplapply(seq_along(beds), function(i) {
    fragments <- readAlignmentFromBed(beds[i], paired = paired)
    if (!isTRUE(all.equal(sort(seqlevels(fragments)), sort(seqlevels(peaks))))){
      merged_seq <- unique(c(seqlevels(fragments), seqlevels(peaks)))
      seqlevels(fragments) <- merged_seq
      seqlevels(peaks) <- merged_seq
    }
    return(list(counts = countOverlaps(peaks, fragments, type = "any", 
                                       ignore.strand = TRUE), 
                depth = length(fragments)))
  })
  
  mat <- vapply(results, function(x) x[["counts"]], rep(0, length(peaks)))
  depths <- vapply(results, function(x) x[["depth"]], 0)
  
  colnames(mat) <- basename(beds)
  counts_mat <- Matrix::Matrix(mat)
  if (is.null(colData)) {
    colData <- DataFrame(depth = depths)
  } else {
    colData$depth <- depths
  }
  out <- SummarizedExperiment(assays = list(counts = counts_mat), 
                              rowRanges = peaks, 
                              colData = colData)
  return(out)
}



# Helper functions for reading in counts from bam ------------------------------


readAlignmentFromBed <- function(filename, paired) {
  if (is.installed("readr")) {
    tmp <- suppressMessages(readr::read_tsv(file = filename, col_names = FALSE))
  } else {
    tmp <- read.delim(file = filename, col.names = FALSE, sep = "\t", 
                      stringsAsFactors = FALSE)
  }
  strand_col <- which(apply(tmp[seq_len(min(100, nrow(tmp))), ], 2, 
                            function(x) all(x %in%  c("+", "-", "*"))))
  if (length(strand_col) == 1) {
    tmp <- tmp[, c(1:3, strand_col)]
    colnames(tmp) <- c("chr", "start", "end", "strand")
    tmp[, "start"] <- tmp[, "start"] + 1
    tmp <- GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end), 
                             strand = tmp$strand)
  } else {
    tmp <- tmp[, 1:6]
    colnames(tmp) <- c("chr", "start", "end")
    tmp[, "start"] <- tmp[, "start"] + 1
    tmp <- GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end))
  }
  if (paired) {
    left <- resize(tmp, width = 1, fix = "start", ignore.strand = TRUE)
    right <- resize(tmp, width = 1, fix = "end", ignore.strand = TRUE)
    out <- left_right_to_grglist(left, right)
  } else {
    out <- resize(tmp, width = 1, ignore.strand = FALSE)
  }
  return(out)
}

#' @importFrom IRanges PartitioningByEnd
#' @importFrom BiocGenerics relist
left_right_to_grglist <- function(left, right) {
  stopifnot(length(left) == length(right))
  if (length(left) == 0) {
    return(GenomicRangesList())
  }
  x <- c(left, right)[as.vector(matrix(seq_len(2L * length(left)), nrow = 2L, 
                                       byrow = TRUE))]
  p <- PartitioningByEnd(cumsum(rep(2, length(x)/2)))
  out <- relist(x, p)
  return(out)
}


bamToFragmentsByRG <- function(bamfile, paired) {
  
  if (paired) {
    scanned <- scanBam(bamfile, 
                       param = ScanBamParam(
                         flag = scanBamFlag(isMinusStrand = FALSE, 
                                            isProperPair = TRUE),
                         what = c("rname", "pos", "isize"),
                         tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))
    
    out <- bplapply(RG_tags, function(RG) {
      match_RG <- which(scanned$tag$RG == RG)
      scanned_left <- GRanges(seqnames = scanned$rname[match_RG], 
                              IRanges(start = scanned$pos[match_RG], 
                                      width = 1), strand = "+")
      scanned_right <- GRanges(seqnames = scanned$rname[match_RG], 
                               IRanges(start = scanned$pos[match_RG] + 
                                         abs(scanned$isize[match_RG]) - 1, 
                                       width = 1), strand = "-")
      return(left_right_to_grglist(scanned_left, scanned_right))
    })
  } else {
    scanned <- scanBam(bamfile,
                       param = ScanBamParam(what = c("rname", 
                                                     "pos", "strand", "qwidth"),
                                            tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))
    
    out <- bplapply(RG_tags, function(RG) {
      match_RG <- which(scanned$tag$RG == RG)
      return(GRanges(seqnames = scanned$rname[match_RG], 
                     IRanges(start = ifelse(scanned$strand[match_RG] ==  "-", 
                                            scanned$pos[match_RG] + 
                                              scanned$qwidth[match_RG] - 1, 
                                            scanned$pos[match_RG]), 
                             width = 1)))
    })
  }
  
  names(out) <- RG_tags
  
  return(out)
}


bamToFragments <- function(bamfile, paired) {
  if (paired) {
    scanned <- scanBam(bamfile, 
                       param = 
                         ScanBamParam(flag = 
                                        scanBamFlag(isMinusStrand = FALSE, 
                                                    isProperPair = TRUE),
                                      what = c("rname", "pos", "isize")))[[1]]
    scanned_left <- GRanges(seqnames = scanned$rname, 
                            IRanges(start = scanned$pos, 
                                    width = 1), strand = "+")
    scanned_right <- GRanges(seqnames = scanned$rname, 
                             IRanges(start = scanned$pos + 
                                       abs(scanned$isize) - 1, width = 1),
                             strand = "-")
    out <- left_right_to_grglist(scanned_left, scanned_right)
  } else {
    scanned <- scanBam(bamfile,
                       param = ScanBamParam(what = c("rname", 
                                                     "pos", 
                                                     "strand", 
                                                     "qwidth")))[[1]]
    out <- GRanges(seqnames = scanned$rname, 
                   IRanges(start = ifelse(scanned$strand == "-", 
                                          scanned$pos + scanned$qwidth - 1, 
                                          scanned$pos),
                           width = 1))
  }
  return(out)
  
}

getFragmentCountsByRG <- function(bam, peaks, paired) {
  message("Reading in file: ", bam)
  rg_fragments <- bamToFragmentsByRG(bam, paired)
  
  tmpfun <- function(frags) {
    overlaps <- as.data.frame(findOverlaps(peaks, 
                                           frags, 
                                           type = "any", 
                                           ignore.strand = TRUE))
    return(overlaps)
  }
  
  all_overlaps <- bplapply(rg_fragments, tmpfun)
  counts_mat <- sparseMatrix(
    i = do.call(rbind, all_overlaps)$queryHits,
    j = unlist(lapply(seq_along(all_overlaps), 
                      function(y) rep(y, 
                                      nrow(all_overlaps[[y]]))), 
               use.names = FALSE),
    x = 1, 
    dims = c(length(peaks), length(rg_fragments)), 
    dimnames = list(NULL, names(rg_fragments)))
  
  return(list(counts = counts_mat, 
              depths = vapply(rg_fragments, length, 0)))
}

# Read in depths from bam ------------------------------------------------------

#' getSampleDepths
#'
#' makes vector of read depths in bam files or RG groups within bam files
#' @param alignment_files filenames for bam or bed file(s) with aligned reads
#' @param by_rg use RG tags to separate groups?
#' @param paired paired end data?
#' @param format bam or bed format? default is bam
#' @return numeric vector
#' @seealso \code{\link{getCounts}},  \code{\link{filterSamples}}
#' @export
#' @examples 
#' 
#' # With single bam with RG tags (can also give multiple bams with RG)
#' test_rg <- system.file("extdata", "test_RG.bam", package = "chromVAR")
#' test_counts <- getSampleDepths(test_rg,  by_rg = TRUE, 
#'                           paired = TRUE)
#'                           
#' 
#' # Multiple bams without RG tags
#' test_bam1 <- system.file("extdata", "test_single1.bam", package = "chromVAR")
#' test_bam2 <- system.file("extdata", "test_single2.bam", package = "chromVAR")
#' test_bam3 <- system.file("extdata", "test_single3.bam", package = "chromVAR")
#' test_counts2 <- getSampleDepths(c(test_bam1, test_bam2,test_bam3), 
#'                            by_rg = FALSE, 
#'                            paired = TRUE) 
#'                            
getSampleDepths <- function(alignment_files, 
                              paired = TRUE, 
                              by_rg = FALSE, 
                              format = c("bam", "bed")) {
  format <- match.arg(format)
  if (format == "bam") {
    return(get_sample_depths_from_bams(alignment_files, paired, by_rg))
  } else {
    return(get_sample_depths_from_beds(alignment_files))
  }
}


get_sample_depths_from_bams <- function(bams, paired = TRUE, by_rg = FALSE) {
  if (by_rg) {
    out <- do.call(c, lapply(bams, getSampleDepthsByRG, paired = paired))
  } else {
    if (paired) {
      out <- vapply(bams, 
                    function(x) 
                      countBam(x,
                               param = 
                                 ScanBamParam(
                                   scanBamFlag(isMinusStrand = FALSE, 
                                               isProperPair = TRUE)))$records, 
                    0)
    } else {
      out <- vapply(bams, function(x) countBam(x)$records, 0)
    }
    names(out) <- vapply(bams, basename, "")
  }
  return(out)
}

get_sample_depths_from_beds <- function(beds) {
  if (is.installed("readr")) {
    out <- 
      do.call(c, 
              bplapply(beds,function(filename){ 
                nrow(suppressMessages(readr::read_tsv(file = filename, 
                                                      col_names = FALSE)))}))
  } else {
    out <- do.call(c, 
                   bplapply(beds, 
                            function(filename) 
                              nrow(read.delim(file = filename, 
                                              header = FALSE, 
                                              sep = "\t", 
                                              stringsAsFactors = FALSE))))
  }
  names(out) <- vapply(beds, basename,"")
  return(out)
}

getSampleDepthsByRG <- function(bamfile, paired = TRUE) {
  if (paired) {
    tags <- scanBam(bamfile, 
                    param = ScanBamParam(flag = 
                                           scanBamFlag(isMinusStrand = FALSE, 
                                                       isProperPair = TRUE),
                                         tag = "RG"))[[1]]$tag$RG
  } else {
    tags <- scanBam(bamfile, param = ScanBamParam(tag = "RG"))[[1]]$tag$RG
  }
  
  RG_tags <- mxsort(unique(tags))
  out <- tabulate(factor(tags, levels = RG_tags, ordered = TRUE))
  names(out) <- RG_tags
  return(out)
}



