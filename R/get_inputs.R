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
#'  @export
get_peaks <- function(filename, extra_cols = c()){
  if (is.installed('readr')){
    bed <- readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]
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
get_counts<- function(alignment_files, peaks, paired, by_rg = FALSE, format = c("bam","bed")){
  
  format = match.arg(format)
  if (format == "bam"){
    return(get_counts_from_bams(alignment_files, peaks, paired, by_rg))
  } else{
    return(get_counts_from_beds(alignment_files, peaks, paired))
  }
}


get_counts_from_bams <- function(bams, peaks, paired, by_rg = FALSE){
  
  if (by_rg){
    counts_mat <- do.call(cBind, lapply(bams, getFragmentCountsByRG, peaks = peaks, paired = paired))
  } else{
    mat = matrix(nrow = length(peaks), ncol = length(bams))
    
    for (i in seq_along(bams)){
      message(paste("Reading in file: ",bams[i], sep="",collapse=""))
      fragments = bamToFragments(bams[i], paired = paired)
      mat[,i] = GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=T)
    }
    colnames(mat) = basename(bams)
    counts_mat = Matrix(mat)    
  }
  
  return(counts_mat)
}


get_counts_from_beds <- function(beds, peaks, paired){
  
  mat = simplify2array(BiocParallel::bplapply(seq_along(beds), function(i){
    fragments = readAlignmentFromBed(beds[i], paired = paired)
    return(GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=TRUE))
  }))
  colnames(mat) = basename(beds)
  counts_mat = Matrix(mat)   
  
  return(counts_mat)
}



# Helper functions for reading in counts from bam ------------------------------

readAlignmentFromBed <- function(filename, paired){
  if (is.installed('readr')){
    tmp <- readr::read_tsv(file = filename, col_names = FALSE)
  } else{
    tmp <- read.delim(file = filename, col.names = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
    tmp <- tmp[,1:3]
    colnames(tmp) <- c("chr", "start", "end")
    tmp[,"start"] <- tmp[,"start"] + 1
    tmp <- makeGRangesFromDataFrame(tmp)  
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
    overlaps = as.data.frame(GenomicRanges::findOverlaps(peaks, frags, type="any", ignore.strand=T))
    return(overlaps)
  }
  
  all_overlaps <-  BiocParallel::bplapply(rg_fragments,tmpfun)
  counts_mat <- sparseMatrix(i = do.call(rbind,all_overlaps)$queryHits, 
                             j = unlist(lapply(seq_along(all_overlaps),function(y) rep(y,nrow(all_overlaps[[y]]))),use.names=F),
                             x = 1, dims = c(length(peaks),length(rg_fragments)), dimnames = list(NULL,names(rg_fragments)))
  
  return(counts_mat)
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
    out = do.call(c, BiocParallel::bplapply(beds, function(filename) nrow(readr::read_tsv(file = filename, col_names = FALSE))))
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


# Filter samples based on number of reads in peaks -----------------------------

#' filter_samples
#' 
#' function to get indices of samples that pass filtters
#' @param counts_mat matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param depths vector of sequencing depth per samples, as computed by
#' \code{\link{getSampleDepths}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param plot show plot?
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated based on data.
#' min_in_peaks is set to 0.5 times the median proportion of fragments in peaks.  min_depth is 
#' set to the maximum of 500 or 10% of the median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_peaks}}
#' @export          
filter_samples <- function(counts_mat, depths, min_in_peaks = NULL, min_depth = NULL, plot = TRUE){
  stopifnot(length(depths) == ncol(counts_mat))
  stopifnot(all_true(names(depths) %in% colnames(counts_mat)))
  fragments_per_sample = colSums(counts_mat)
  if (is.null(min_in_peaks)){
    min_in_peaks = round(median(fragments_per_sample/depths)*0.5, digits = 3)
    message(paste("min_in_peaks set to ",min_in_peaks,sep="",collapse=""))
  } 
  if (is.null(min_depth)){
    min_depth = max(500,median(depths)*0.1)
    message(paste("min_depth set to ",min_depth,sep="",collapse=""))
  }
  keep_samples <- intersect(which(depths >= min_depth), 
                            which(fragments_per_sample/depths >= min_in_peaks))   
  tmp_df = data.frame(x = depths, y= fragments_per_sample/depths, z = ((1:length(fragments_per_sample)) %in% keep_samples))
  p = ggplot(tmp_df, aes_string(x="x", y="y",col="z")) + geom_point() +
    xlab("Library size") + ylab("Proportion of fragments in peaks") + scale_x_log10() + annotation_logticks(sides="b") + 
    scale_y_continuous(expand =c(0,0), limits =c(0, min(1,max(tmp_df$y)*1.2)))+
    scale_color_manual(name = "Pass filters?",values = c("gray","black"),breaks = c(TRUE,FALSE), labels = c("Yes","No")) +
                           chromVAR_theme()
  p = p + geom_hline(yintercept = min_in_peaks, col="red", lty = 2) + geom_vline(xintercept = min_depth, col="red", lty=2)
  if (plot) print(p)
  return(keep_samples)    
}

# Filter peaks based on counts -------------------------------------------------

#' filter_peaks
#' 
#' function to get indices of peaks that pass filters
#' @param counts_mat matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param peaks GenomicRanges object
#' \code{\link{read_peaks}}
#' @param min_fragments_per_peak minimum number of fragmints in peaks across all samples 
#' @param non_overlapping reduce peak set to non-overlapping peaks, see details
#' @details if non_overlapping is set to true, when peaks overlap the overlapping peak with lower counts is removed
#' @return vector of indices, representing peaks that should be kept
#' @seealso \code{\link{get_peaks}},  \code{\link{get_inputs}}, \code{\link{filter_samples}},
#' \code{\link{get_counts}}
#' @export          
filter_peaks <- function(counts_mat, peaks, min_fragments_per_peak = 1, non_overlapping = TRUE){
  fragments_per_peak = rowSums(counts_mat)
  keep_peaks <- which(fragments_per_peak >= min_fragments_per_peak)
  if (non_overlapping){
    strand(peaks) <- "*"
    while (!(isDisjoint(peaks[keep_peaks]))){
      chr_names = seqnames(peaks[keep_peaks])
      starts = start(peaks[keep_peaks])
      ends = end(counts_mat@peaks[keep_peaks])
      
      overlap_next = intersect(which(chr_names[1:(length(keep_peaks) -1)] == chr_names[2:(length(keep_peaks))]),
                               which(ends[1:(length(keep_peaks) -1)] > starts[2:(length(keep_peaks))]))
      overlap_previous = overlap_next + 1
      overlap_comparison = fragments_per_pea[keep_peaks[overlap_previous]] > fragments_per_peak[keep_peaks[overlap_previous]]
      discard = keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
      keep_peaks = keep_peaks[keep_peaks %ni% discard]
    }   
  }
  return(keep_peaks)
}



