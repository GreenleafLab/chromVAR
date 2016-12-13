################# Functions for reading in inputs ##############################

# Read in peaks ----------------------------------------------------------------

#' get_peaks
#'
#' Read in peaks from a bed file.
#' @param filename filename of bed file
#' @param extra_cols extra columns to read in beyond first three
#' @return \code{\link[GenomicRanges]{GenomicRanges}} containing peaks in bed file
#' @details As in standard definition of bed file, first column is assumed to be chromosome,
#' second is assumed to be start of peak (0-based), and third is assumed to be end of peak (1-based).
#' Note that in output GenomicRanges output, start and end indices are both 1-based.
#' Extra columns can be added as metadata or strand information if provided, but the user must
#' indicate column index and name using named vector for extra_cols.
#' @seealso \code{\link{get_counts}}, \code{\link{get_inputs}}, \code{\link{filter_peaks}}
#' @export
get_peaks <- function(filename, extra_cols = c(), sort_peaks = TRUE){
  if (is.installed('readr')){
    bed <- as.data.frame(suppressMessages(readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]))
  } else{
    bed <- read.delim(file = filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, c(1:3, extra_cols)]
  }
  colnames(bed) <- c("chr", "start", "end", names(extra_cols))
  bed[,"start"] <- bed[,"start"] + 1
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  if (!isDisjoint(bed)){
       warning("Peaks are overlapping!
            After getting peak counts, peaks can be reduced to non-overlapping set
              using filter_peaks function")
     }
  if (sum(width(bed) == width(bed[1])) != length(bed)){
    warning('Peaks are not equal width!
            Use resize(peaks, width = x, fix = "center") to make peaks equal in size,
            where x is the desired size of the peaks)')
  }
  bed <- GenomeInfoDb::sortSeqlevels(bed)
  sorted_bed = sort(bed, ignore.strand = TRUE)
  if (sort_peaks){
    if (!isTRUE(all.equal(sorted_bed, bed))){
      message("Peaks sorted")
    }
    return(sorted_bed)
  } else{
    if (!isTRUE(all.equal(sorted_bed, bed))){
      warning("Peaks not sorted")
    }
    return(bed)
  }
}

#'@export
read_macs2_narrowpeaks <- function(filename, width = 500, non_overlapping = TRUE){
  if (is.installed('readr')){
    bed <- as.data.frame(readr::read_tsv(file = filename,
                                         col_names = c("chr","start","end","name","score","strand","fc","pval","qval","summit")))
  } else{
    bed <- read.delim(file = filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names= c("chr","start","end","name","score","strand","fc","pval","qval","summit"))
  }
  bed[,"summit"] <- bed[,"start"] + bed[,"summit"]
  bed <- as(bed, "DataFrame")
  bed <- makeGRangesFromDataFrame(bed[,c("chr","summit","score","qval","name")],
                                  start.field = "summit",
                                  end.field = "summit",
                                  keep.extra.columns = TRUE)
  bed <- resize(bed, width = width, fix = "center")
  bed <- GenomeInfoDb::sortSeqlevels(bed)
  bed <- sort(bed)
  if (non_overlapping){
    keep_peaks = 1:length(bed)
    while (!(isDisjoint(bed[keep_peaks]))){
      chr_names = as.character(seqnames(bed[keep_peaks]))
      starts = start(bed[keep_peaks])
      ends = end(bed[keep_peaks])
      overlap_next = intersect(which(chr_names[1:(length(keep_peaks) -1)] == chr_names[2:(length(keep_peaks))]),
                               which(ends[1:(length(keep_peaks) -1)] >= starts[2:(length(keep_peaks))]))
      overlap_previous = overlap_next + 1
      overlap_comparison = bed[keep_peaks[overlap_previous]]$qval > bed[keep_peaks[overlap_next]]$qval
      discard = keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
      keep_peaks = keep_peaks[keep_peaks %ni% discard]
    }
    bed = bed[keep_peaks]
  }
  return(bed)
}

# Main function for reading in counts from bam ---------------------------------

#' get_counts
#'
#' makes matrix of fragment counts in peaks using one or multiple bam or bed files
#' @param alignment_files filenames for bam or bed files with aligned reads
#' @param peaks GRanges object with peaks
#' @param by_rg use RG tags in bam to separate groups?
#' @param paired paired end data?
#' @param format bam or bed?  default is bam
#' @return \code{\link[Matrix]{Matrix}} object
#' @seealso \code{\link{get_sample_depths}},  \code{\link{get_inputs}}, \code{\link{filter_samples}}
#' @export
get_counts<- function(alignment_files, peaks, paired, by_rg = FALSE, format = c("bam","bed"), colData = NULL){

  format = match.arg(format)
  if (format == "bam"){
    return(get_counts_from_bams(alignment_files, peaks, paired, by_rg, colData))
  } else{
    return(get_counts_from_beds(alignment_files, peaks, paired, colData))
  }
}

#' @import SummarizedExperiment
get_counts_from_bams <- function(bams, peaks, paired, by_rg = FALSE, sample_annotation = NULL){

  if (by_rg){
    tmp <- lapply(bams, getFragmentCountsByRG, peaks = peaks, paired = paired)
    if (!is.null(sample_annotation) && nrow(sample_annotation) == length(bams)){
      sample_annotation <- as(sample_annotation,"DataFrame")
      l <- sapply(tmp, function(x) length(x$depths))
      sample_annotation <- do.call(rbind, lapply(seq_along(l), function(x) rep(sample_annotation[x,,drop=FALSE],l[x])))
    }
    counts_mat <- do.call(cBind, lapply(tmp, function(x) x$counts))
    depths <- do.call(c, lapply(tmp, function(x) x$depths))
  } else{
    mat = matrix(nrow = length(peaks), ncol = length(bams))
    depths = vector("numeric", length(bams))

    for (i in seq_along(bams)){
      message(paste("Reading in file: ",bams[i], sep="",collapse=""))
      fragments = bamToFragments(bams[i], paired = paired)
      depths[i] = length(fragments)
      mat[,i] = GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=T)
    }
    colnames(mat) = basename(bams)
    counts_mat = Matrix::Matrix(mat)
  }
  if (is.null(sample_annotation)){
    sample_annotation <- DataFrame(depth = depths)
  } else{
    sample_annotation$depth = depths
  }
  out <- SummarizedExperiment(assays = list(counts = counts_mat), rowRanges = peaks, colData = sample_annotation)
  return(out)
}


get_counts_from_beds <- function(beds, peaks, paired, colData = NULL){


  results <- BiocParallel::bplapply(seq_along(beds), function(i){
    fragments = readAlignmentFromBed(beds[i], paired = paired)
    if (!isTRUE(all.equal(sort(seqlevels(fragments)), sort(seqlevels(peaks))))){
      merged_seq <- unique(c(seqlevels(fragments), seqlevels(peaks)))
      seqlevels(fragments) <- merged_seq
      seqlevels(peaks) <- merged_seq
    }
    return(list(counts = GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=TRUE),
                depth = length(fragments)))
  })

  mat <-  vapply(results, function(x) x[["counts"]], rep(0, length(peaks)))
  depths <- vapply(results, function(x) x[["depth"]], 0)

  colnames(mat) = basename(beds)
  counts_mat = Matrix::Matrix(mat)
  if (is.null(colData)){
    colData <- DataFrame(depth = depths)
  } else{
    colData$depth = depths
  }
  out <- SummarizedExperiment(assays = list(counts = counts_mat), rowRanges = peaks, colData = colData)
  return(out)
}



# Helper functions for reading in counts from bam ------------------------------

#' @export
readAlignmentFromBed <- function(filename, paired){
  if (is.installed('readr')){
    tmp <- suppressMessages(readr::read_tsv(file = filename, col_names = FALSE))
  } else{
    tmp <- read.delim(file = filename, col.names = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
   strand_col <- which(apply(tmp[1:min(100, nrow(tmp)),], 2, function(x) all_true(x %in% c("+","-","*"))))
   if (length(strand_col) == 1){
     tmp <- tmp[,c(1:3,strand_col)]
     colnames(tmp) <- c("chr", "start", "end", "strand")
     tmp[,"start"] <- tmp[,"start"] + 1
     tmp <- with(tmp, GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end), strand = strand))
   } else{
      tmp <- tmp[,1:6]
      colnames(tmp) <- c("chr", "start", "end")
      tmp[,"start"] <- tmp[,"start"] + 1
      tmp <- with(tmp, GRanges(tmp$chr, ranges = IRanges(tmp$start, tmp$end)))
   }
  if (paired){
    left <- resize(tmp, width = 1, fix = "start", ignore.strand = TRUE)
    right <- resize(tmp, width = 1, fix = "end", ignore.strand = TRUE)
    out <- left_right_to_grglist(left,right)
  } else{
    out <- resize(tmp, width = 1, ignore.strand = FALSE)
  }
  return(out)
}

left_right_to_grglist <- function(left, right){
  stopifnot(length(left) == length(right))
  if (length(left) == 0){
    return(GenomicRangesList())
  }
  x = c(left,right)[as.vector(matrix(seq_len(2L * length(left)), nrow=2L, byrow=TRUE))]
  p = IRanges::PartitioningByEnd(cumsum(rep(2,length(x)/2)))
  out = BiocGenerics::relist(x,p)
  return(out)
}

#' @export
bamToFragmentsByRG <- function(bamfile, paired){

  if (paired){
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                                                isProperPair = TRUE),
                                                                  what = c("rname","pos","isize"),
                                                                  tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))

    out <- BiocParallel::bplapply(RG_tags, function(RG){
      match_RG = which(scanned$tag$RG == RG)
      scanned_left <- GRanges(seqnames = scanned$rname[match_RG],
                              IRanges::IRanges(start = scanned$pos[match_RG],
                                               width = 1),
                              strand = "+")
      scanned_right <- GRanges(seqnames = scanned$rname[match_RG],
                               IRanges::IRanges(start = scanned$pos[match_RG] +
                                                  abs(scanned$isize[match_RG]) - 1,
                                                width = 1),
                               strand = "-")
      return(left_right_to_grglist(scanned_left, scanned_right))
    })
  } else{
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(what = c("rname","pos","strand","qwidth"),
                                                                  tag = "RG"))[[1]]
    RG_tags <- mxsort(unique(scanned$tag$RG))

    out <- BiocParallel::bplapply(RG_tags, function(RG){
      match_RG = which(scanned$tag$RG == RG)
      return(GRanges(seqnames = scanned$rname[match_RG],
                     IRanges::IRanges(start = ifelse(scanned$strand[match_RG] == "-", scanned$pos[match_RG] + scanned$qwidth[match_RG]-1, scanned$pos[match_RG]),
                                      width = 1)))
    })
  }

  names(out) = RG_tags

  return(out)
}

#' @export
bamToFragments <- function(bamfile, paired){
  if (paired){
    scanned <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE), what = c("rname","pos","isize")))[[1]]
    scanned_left <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos, width = 1), strand = "+")
    scanned_right <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos + abs(scanned$isize) - 1, width = 1), strand = "-")
    out <- left_right_to_grglist(scanned_left, scanned_right)
  } else{
    scanned <- Rsamtools::scanBam(bamfile,
                                  param = Rsamtools::ScanBamParam(what = c("rname","pos","strand","qwidth")))[[1]]
    out <- GRanges(seqnames = scanned$rname,
                   IRanges::IRanges(start = ifelse(scanned$strand == "-", scanned$pos + scanned$qwidth-1, scanned$pos),
                                    width = 1))
  }
  return(out)

}

getFragmentCountsByRG <- function(bam, peaks, paired){
  message(paste("Reading in file: ",bam, sep="",collapse=""))
  rg_fragments <- bamToFragmentsByRG(bam, paired)

  tmpfun <- function(frags){
    overlaps = as.data.frame(GenomicRanges::findOverlaps(peaks, frags, type="any", ignore.strand=TRUE))
    return(overlaps)
  }

  all_overlaps <-  BiocParallel::bplapply(rg_fragments,tmpfun)
  counts_mat <- sparseMatrix(i = do.call(rbind,all_overlaps)$queryHits,
                             j = unlist(lapply(seq_along(all_overlaps),function(y) rep(y,nrow(all_overlaps[[y]]))),use.names=FALSE),
                             x = 1, dims = c(length(peaks),length(rg_fragments)), dimnames = list(NULL,names(rg_fragments)))

  return(list(counts = counts_mat, depths = sapply(rg_fragments, length)))
}

# Read in depths from bam ------------------------------------------------------

#' get_sample_depths
#'
#' makes vector of read depths in bam files or RG groups within bam files
#' @param alignment_files filenames for bam or bed file(s) with aligned reads
#' @param by_rg use RG tags to separate groups?
#' @param paired paired end data?
#' @param format bam or bed format? default is bam
#' @return numeric vector
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_samples}}
#' @export
get_sample_depths <- function(alignment_files, paired = TRUE, by_rg = FALSE, format = c("bam","bed")){
  format = match.arg(format)
  if (format == "bam"){
    return(get_sample_depths_from_bams(alignment_files, paired, by_rg))
  } else{
    return(get_sample_depths_from_beds(alignment_files))
  }
}
get_sample_depths_from_bams <- function(bams, paired = TRUE, by_rg = FALSE){
  if (by_rg){
    out <- do.call(c, lapply(bams, getSampleDepthsByRG, paired = paired))
  } else{
    if (paired){
      out <- sapply(bams,
                    Rsamtools::countBam,
                    param = Rsamtools::ScanBamParam(Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                     isProperPair = TRUE)))
    } else{
      out <- sapply(bams, Rsamtools::countBam)
    }
    names(out) <- sapply(bams, basename)
  }
  return(out)
}

get_sample_depths_from_beds <- function(beds){
  if (is.installed('readr')){
    out = do.call(c, BiocParallel::bplapply(beds, function(filename) nrow(suppressMessages(readr::read_tsv(file = filename, col_names = FALSE)))))
  } else{
    out = do.call(c, BiocParallel::bplapply(beds, function(filename) nrow(read.delim(file = filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE))))
  }
  names(out) <- sapply(beds, basename)
  return(out)
}

getSampleDepthsByRG <- function(bamfile, paired = TRUE){
  if (paired){
    tags <- Rsamtools::scanBam(bamfile,
                               param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE,
                                                                                             isProperPair = TRUE),
                                                               tag = "RG"))[[1]]$tag$RG
  } else{
    tags <- Rsamtools::scanBam(bamfile,
                               param = Rsamtools::ScanBamParam(tag = "RG"))[[1]]$tag$RG
  }

  RG_tags <- mxsort(unique(tags))
  out <- tabulate(factor(tags, levels = RG_tags, ordered = TRUE))
  names(out) <- RG_tags
  return(out)
}



