#' fragmentCounts
#' 
#' fragmentCounts is a class to store fragment counts for genomic regions across samples. 
#' The class is specifically designed to be able to handle sparse data, as would 
#' arise from single-cell ATAC-seq experiment.  The counts slot holds the count data as
#' a Matrix (from package Matrix); the specific class of the Matrix is determined by
#' sparsity of the data.  The fragmentCounts class also stores the genomic regions that
#' correspond to the rows of the counts matrix. Optionally, the class can store sample meta data
#' as a data.frame with rows corresponding to samples.   
#' @section Accessors:
#' Accessors are available for counts, peaks, and meta.
#' @slot counts Matrix of fragment counts, with each row representing a peak, each column a sample
#' @slot peaks GRanges object with peaks for each row of counts
#' @slot total_fragments total fragments within peaks
#' @slot fragments_per_sample total fragments within peaks per sample
#' @slot fragments_per_peak total fragments within each peak across samples
#' @slot nsample number of samples
#' @slot npeak number of peaks
#' @slot sample_meta optional data.frame with meta data for samples, with rows corresponding to samples
#' @slot depth total read depth for each sample
#' @export
fragmentCounts <- setClass("fragmentCounts",
                           slots = c(counts = 'Matrix',
                                     peaks = 'GenomicRanges',
                                     total_fragments = 'numeric',
                                     fragments_per_sample = 'numeric',
                                     fragments_per_peak = 'numeric',
                                     nsample = 'numeric',
                                     npeak = 'numeric',
                                     sample_meta = 'data.frame',
                                     depth = 'numeric'
                           ))


#' @export 
combine_samples <- function(samples, meta_name = "sample_group"){
  stopifnot(all_true(sapply(samples, inherits, 'fragmentCounts')))
  stopifnot(all_true(sapply(samples, function(x) all.equal(x@peaks,samples[[1]]@peaks))))
  for (i in seq_along(samples)){
    sample = samples[[i]]
    if (ncol(sample@sample_meta)==0){
      sample@sample_meta = data.frame(meta_name = rep(names(samples)[i], sample@nsample))
    } else{
      sample@sample_meta[,meta_name] = rep(names(samples)[i], sample@nsample)
    }    
  } 
  out <- fragmentCounts(counts = do.call(cBind, lapply(samples, function(x) x@counts)),
                        peaks = samples[[1]]@peaks,
                        sample_meta = do.call(rbind, lapply(samples, function(x) x@sample_meta)),
                        depth = do.call(c,sapply(samples, function(x) x@depth)))
  
  return(out)
}


#' @rdname fragmentCounts-class
#' @export
setMethod("initialize",
          "fragmentCounts",
          function(.Object, ...){
            .Object <- callNextMethod()
            .Object@total_fragments <- sum(.Object@counts)
            .Object@fragments_per_sample <- colSums(.Object@counts)
            .Object@fragments_per_peak <- rowSums(.Object@counts)
            .Object@nsample <- ncol(.Object@counts)
            .Object@npeak <- nrow(.Object@counts)
            .Object
          })

#' @rdname fragmentCounts-class
#' @export
setMethod("show",
          signature="fragmentCounts",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Counts for ", ncol(object@counts), " cells/samples and ",
                nrow(object@counts), " peaks. \n", sep = "")
            invisible(NULL)
          })

setValidity("fragmentCounts", function(object) {
  msg <- NULL
  valid <- TRUE
  #Check dimensions
  if (nrow(object@counts) != length(object@peaks)){
    valid <- FALSE
    msg <- c(msg, "Number of rows of counts matrix must be same as length of peaks")
  }
  if (valid) TRUE else msg})


#' @rdname fragmentCounts-class
#' @export
setMethod("[", signature = signature(x = "fragmentCounts"),
          definition = function(x, i ,j ) {
            if (missing(i)){
              i = 1:x@npeak
            }
            if (missing(j)){
              j = 1:x@nsample
            }
            x@counts = x@counts[i,j]
            x@fragments_per_sample = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@peaks = x@peaks[i]
            x@total_fragments = sum(x@counts)
            x@npeak = length(x@peaks)
            x@nsample = ncol(x@counts)
            x@sample_meta = x@sample_meta[j,]
            x@depth = x@depth[j]
            return(x)
          })


#' getFragmentCountsByRG
#' 
#' make fragmentCounts object using bam file with RG tags 
#' @param bam filename for bam file with aligned reads
#' @param peaks GRanges object with peaks 
#' @return \code{\link{fragmentCounts}} object
#' @seealso \code{\link{getFragmentCounts}}, \code{\link{fragmentCounts}}, \code{\link{filterFragmentCounts}}
#' @export
getFragmentCountsByRG <- function(bam, peaks, BPPARAM = BiocParallel::bpparam()){

  rg_fragments <- bamToFragmentsByRG(bam, BPPARAM)
  depth <- sapply(rg_fragments, length)
  tmpfun <- function(frags){
    overlaps = as.data.frame(GenomicRanges::findOverlaps(peaks, frags, type="any", ignore.strand=T))
    return(overlaps)
  }

  all_overlaps <-  BiocParallel::bplapply(rg_fragments,tmpfun,
                                        BPPARAM = BPPARAM)
  counts_mat <- sparseMatrix(i = do.call(rbind,all_overlaps)$queryHits, 
                                   j = unlist(lapply(seq_along(all_overlaps),function(y) rep(y,nrow(all_overlaps[[y]]))),use.names=F),
                                   x = 1, dims = c(length(peaks),length(rg_fragments)), dimnames = list(NULL,names(rg_fragments)))

  counts <- fragmentCounts(counts = counts_mat,
                                    peaks = peaks,
                                    depth = depth)

  return(counts)
}


#' getFragmentCounts
#' 
#' make fragmentCounts object using multiple bam files
#' @param bams filenames for bam file with aligned reads
#' @param peaks GRanges object with peaks 
#' @return \code{\link{fragmentCounts}} object
#' @seealso \code{\link{getFragmentCountsByRG}}, \code{\link{fragmentCounts}}, \code{\link{filterFragmentCounts}}
#' @export
getFragmentCounts <- function(bams, peaks){

  depth = rep(0, length(bams))
  mat = matrix(nrow = length(peaks), ncol = length(bams))
  
  for (i in seq_along(bams)){
    fragments = bamToFragments(bams[i])
    depth[i] = length(fragments)
    mat[,i] = countOverlaps(peaks, fragments, type="any", ignore.strand=T)
  }
  
  counts_mat = Matrix(mat)
  
  counts <- fragmentCounts(counts = counts_mat,
                           peaks = peaks,
                           depth = depth)

  return(counts)
}

#not exported
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
bamToFragmentsByRG <- function(bamfile, BPPARAM = BiocParallel::bpparam()){

  scanned <- Rsamtools::scanBam(bamfile,
                                param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, 
                                                                                              isProperPair = TRUE),
                                                                         what = c("rname","pos","isize"),
                                                                         tag = "RG"))[[1]]
  RG_tags <- gtools::mixedsort(unique(scanned$tag$RG), decreasing = TRUE)

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
  }, BPPARAM = BPPARAM)

  names(out) = RG_tags

  return(out)
}

#' @export
bamToFragments <- function(bamfile){

  scanned <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE), what = c("rname","pos","isize")))[[1]]
  scanned_left <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos, width = 1), strand = "+")
  scanned_right <- GRanges(seqnames = scanned$rname, IRanges::IRanges(start = scanned$pos + abs(scanned$isize) - 1, width = 1), strand = "-")
  out <- left_right_to_grglist(scanned_left, scanned_right)

  return(out)

}


#' filterFragmentCounts
#' 
#' Filter Fragment Counts by quality
#' @param counts_mat fragmentCounts object
#' @param min_in_peaks minimum fraction of fragments within peaks for a sample (default = 0.25)
#' @param min_fragments minimum fragments in peaks for a sample(default = 2500)
#' @param min_fragments_per_peak minimum fragments across samples that fall into a peak
#' @details Filtering of samples based on min_in_peaks and min_fragments arguments is performed first.
#' Then peaks are filtered by min_fragments_per_peak argument (so only fragments in samples that survived
#' filtering of samples are counted).
#' @return \code{\link{fragmentCounts}} object
#' @seealso \code{\link{getFragmentCounts}}, \code{\link{getFragmentCountsByRG}}, \code{\link{fragmentCounts}}
#' @export          
filterFragmentCounts <- function(counts_mat, min_in_peaks = 0.25, min_fragments = 2500, min_fragments_per_peak = 3){
  stopifnot(inherits(counts_mat, "fragmentCounts"))
  stopifnot(sum(counts_mat@depth >= counts_mat@fragments_per_sample) == counts_mat@nsample)
  keep_samples <- intersect(which(counts_mat@fragments_per_sample >= min_fragments), 
                    which(counts_mat@fragments_per_sample/counts_mat@depth >= min_in_peaks))  
  if (length(keep_samples) == 0){
    stop("No samples passed filters!")
  }
  counts_mat <- counts_mat[,keep_samples]
  keep_peaks <- which(counts_mat@fragments_per_peak >= min_fragments_per_peak)
  if (length(keep_peaks) == 0){
    stop("No peaks passed filters!")
  }
  counts_mat <- counts_mat[keep_peaks,]
  return(counts_mat)
}







