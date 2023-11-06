
.importFragments <- function(fragData, chunkSize=10)
  {
    if(is.list(fragData))
    {
      if(grepl(".bam", fragData[[1]], fixed=TRUE))
      {
        print("Importing .bam files")
        param <- Rsamtools::ScanBamParam(what=c('pos', 'qwidth', 'isize'))
        fragData <- lapply(fragData, function(bamPath){
          readPairs <- GenomicAlignments::readGAlignmentPairs(bamPath, param=param)
          # get fragment coordinates from read pairs
          frags <- GRanges(seqnames(GenomicAlignments::first(readPairs)), 
                           IRanges(start=pmin(GenomicAlignments::start(GenomicAlignments::first(readPairs)), 
                                              GenomicAlignments::start(GenomicAlignments::second(readPairs))), 
                                   end=pmax(GenomicAlignments::end(GenomicAlignments::first(readPairs)), 
                                            GenomicAlignments::end(second(readPairs)))))
          frags <- granges(frags, use.mcols=TRUE)
          
          # ATAC shift
          start(frags) <- ifelse(strand(frags) == "+", start(frags) + 4, start(frags))
          end(frags) <- ifelse(strand(frags) == "-", end(frags) - 5, end(frags))
          frags <- as.data.table(frags)
          setnames(frags, c("seqnames"), c("chr"))
          frags$count <- 1
          frags
        })
      }
      else if(grepl(".bed", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, fread, select=1:3, col.names=c("chr", "start", "end"))
        fragData <- lapply(fragData, function(dt){dt$count <- 1
        dt <- dt[,c("chr", "start", "end", "count"), with=FALSE]
        dt})
      }
      else if(grepl(".rds", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, function(gr){
          gr <- readRDS(gr)
          dt <- as.data.table(gr)
          setnames(dt, c("seqnames"), c("chr"))
          dt$count <- 1
          dt <- dt[,c("chr", "start", "end", "count")]
          dt
        })
      }
      else if(is.data.table(fragData[[1]]))
      {
        if(!(grepl(paste(c("chr", "start", "end"), collapse=";"),
                   paste(colnames(fragData[[1]]),collapse=";"))))
        {stop("data.table list elements of data need to contain columns: chr, start, end")}
      }
      else stop("List elements of data need to be either data.tables, .bams-, .beds- .rds-files")
    }
    else
    {
      stop("Data needs to be a list")
    }
    
    return(fragData)
}

#' @description
#' Resize the peaks
#' 
#' @param peakRanges: a GRanges object of peak ranges
#' @param width: the re-defined size of each peak
#' @return a GRange object with resized ranges

.resizeRanges <- function(peakRanges, 
  width = 200, 
  fix = c("center", "start", "end", "summit"),
  ...) {
  
      fix <- match.arg(fix, choices = c("center", "start", "end", "summit"))
      # Sanity check
      if (!class(peakRanges) == "GRanges") {
        stop("peakRanges must be a GRanges object")
      }
      
      if (fix == "summit") {
        start(peakRanges) <- round(peakRanges$summit-width/2)
        end(peakRanges) <- start(peakRanges)+width-1
      } else {
        peakRanges <- resize(peakRanges, width = width, fix = fix)
      }
      
      return(peakRanges)
      
}

#' @param peakRanges: a GRanges object of peak ranges
#' @param genome: a BSgenome object, the corresponding genome 
#' @return a GRanges objects with an additional metadata column gc that contains
#' GC content
.getGCContent <- function(peakRanges, genome) {
    # Sanity check
    if (!class(peakRanges) == "GRanges") {
        stop("peakRanges must be a GRanges object")
    }
        
    peakSeqs <- getSeq(x = genome, peakRanges)
    mcols(peakRanges)$gc <- letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
    peakRanges
} 

#' @param peakRanges a GRange or data.table object that contains peak ranges
#' @param motifRanges a GRange or data.table object that contains motif ranges
#' @param flankSize integer, the number of nucleotides to define the buffer/flanking 
#' region near the motif instance
#' @param shiftATAC logic if shifting the ATAC-seq fragment data table
.getInsertionCounts <- function(peakRanges, 
  motifRanges,
  #rangeType = 
  flankSize = 30,
  shiftATAC = FALSE
  ) {
    # Sanity check
    if (is.data.table(motifRanges) == FALSE) {
      motifRanges <- data.table::as.data.table(motifRanges)
    }
    
    if ("seqnames" %in% colnames(motifRanges))
      colnames(motifRanges)[which(colnames(motifRanges)=="seqnames")] <- "chr"
    
    if (is.data.table(peakRanges) == FALSE) {
      peakRanges <- data.table::as.data.table(peakRanges)
    }
  
    if ("seqnames" %in% colnames(peakRanges))
      colnames(peakRanges)[which(colnames(peakRanges)=="seqnames")] <- "chr"
    
    # why?
    motifData <- data.table::copy(motifRanges)
    atacFrag <- data.table::copy(peakRanges)
    
    if (shiftATAC == TRUE) {
      atacFrag[, start := ifelse(strand == "+", start + 4, start)]
      atacFrag[, end   := ifelse(strand == "-", end   - 5, end)]
    }
    
    # set up flanking region
    motifData[, start_margin := start - flankSize]
    motifData[, end_margin   := end + flankSize]
    motifData$motifind <- 1:nrow(motifData)
    
    # returning the overlaying counts and within motif counts
    # ?fragment length is usually larger than motif
    res <- motifData[, .(motifind, 
      Start_Count = atacFrag[motifData, 
        on = .(start >= start_margin, start <= start, chr == chr), 
        .N, by = .EACHI]$N,
      End_Count   = atacFrag[motifData, 
        on = .(end >= end, end <= end_margin, chr == chr), 
        .N, by = .EACHI]$N,
      counts = atacFrag[motifData, 
        on = .(start >= start, end <= end, chr == chr), 
        .N, by = .EACHI]$N)] 

  
    res <- res[, .(motifind, 
      total_counts = Start_Count + End_Count + counts,
      flanking_counts = Start_Count + End_Count,
      within_counts = counts) ] 
}

#' @description
#' count the number of fragments in each peak region
#' 
#' @param peakRanges a GRange or data.table object that contains the peak ranges
#' @param fragRanges a list of data.table that contains the fragment ranges, 
#' each data.table represents a sample
.getCountsOverlaps <- function(peakRanges, 
  fragRanges,
  shiftATAC = FALSE,
  by = c("count", "weight")) {
  
    by <- match.arg(by, choices = c("count", "weight"))
  
    # sanity check
    if (is.data.table(peakRanges) == FALSE) 
      peakRanges <- data.table::as.data.table(peakRanges)
    
    if (is.data.table(fragRanges) == FALSE) 
      fragRanges <- data.table::as.data.table(fragRanges)
    
    if ("seqnames" %in% colnames(peakRanges))
      colnames(peakRanges)[which(colnames(peakRanges)=="seqnames")] <- "chr"
    
    if (shiftATAC == TRUE) {
      peakRanges[, start := ifelse(strand == "+", start + 4, start)]
      peakRanges[, end   := ifelse(strand == "-", end   - 5, end)]
    }  
    
    peaks <- data.table::copy(peakRanges)
    frags <- data.table::copy(fragRanges)
    peaks$peakID <- 1:nrow(peaks)
    
    if (by == "count") {
      fragCounts <- lapply(frags, \(frag) {
        res <- peaks[, .(peakID, 
          counts = frag[peaks, 
            on = .(start >= start, end <= end, chr == chr), 
            .N, by = .EACHI]$N)]$counts
      })
      counts <- do.call(cbind, fragCounts)
      rownames(counts) <- peaks$peakID
    }
    
    counts
    
    
}

#' change cuts into bin size
.getFLD <- function(dts, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(names(dts), \(sample) {
        dt <- dts[[sample]]
        dt[,width:=end-start+1]
        # check if max(cuts) covers max(width)
        if (max(dt$width) > cuts[length(cuts)]) 
            cuts[length(cuts)+1] <- max(dt$width)
        
        dt[,type:=as.numeric(cut(width, breaks = cuts))]
        if (length(unique(dt$type))==1) 
            warnings("All fragments fell into 1 type")
        
        
        dtAg <- dt[, .(sum_counts=.N), by=type]
        dtAg[,total:=sum(sum_counts)]
        dtAg[,prop:=sum_counts/total]
        dtAg[,sample_id:=sample]
    })
    names(res) <- names(dts)
    res
}

.weightFragments <- function (dts, cuts = c(0,120,300,500)) {
    fragDist <- .getFLD(dts, cuts = cuts)
    fragDT <- rbindlist(fragDist)
    #fragDT[,sample_type:=paste0(sample_id, ",", type)]
    fragDT[,mean_prop:=mean(prop), by=type]
    fragDT[,weight:=mean_prop/prop]
    res <- split(fragDT, fragDT$sample_id)
}

#' getCounts
#' 
#' @description
#' makes matrix of fragment counts in peaks using a list of bam or bed files
#' @param files a list of filenames for bam or bed files with aligned reads
#' @param paired paired end data?
#' @param resize a logical to determine if peaks should be resized
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#'  object

getCounts <- function (files,
    #motifs, 
    ranges, # motif matches or peaks, set a parameter to define
    rowType = c("peaks", "motifs"),
    paired=TRUE,
    resize=TRUE, 
    width=100,
    binsize,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    ...) {
    
    rowType <- match.arg(rowType, choices=c("peaks", "motifs"))  
    
    dts <- .importFragments(files)
    
    if(rowType=='peak' & resize){
      
      rs <- .resizeRanges(ranges, ...)
    }
    
    rs <- .getGCContent(rs)
    dts <- .weightFragments(dts,...)
    
    if(type=='peak'){
      .getOverlapCounts(fragDTs, rs, cuts, countCol)
    }
    else if(type=="motif"){
      .getInsertionCounts(fragDTs, rs, countCol)
    }
    
    
    
    # define resize
    # if (resize) {
    #   dts <- .resizeRanges(dts)
    # }
    #dts <- lapply(dts, .getGCContent, genome = genome)
    #weights <- .weightFragments(dts, cuts = c(0,10,50,80))
    #motif_pos <- lapply(dts, .getInsertionCounts, motif = motif)

}

# resize before insertionCounts
