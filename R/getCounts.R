
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

#' @param ranges: a GenomicRanges objects that contains 
#' @param size: the size of each region/peak
#' return a list of GRange with resized ranges
#' return a GR
.resizeRanges <- function(ranges, 
  size = 200, 
  fix = c("center", "start", "summit")
  #summit = FALSE,
  ...) {
      # fix <- match.args()
      
      res <- lapply(ranges, \(range) {
        # Sanity check
        if (is.data.table(range) == FALSE) {
          range <- data.table::as.data.table(range)
          if (!all(c("chr", "start", "end") %in% colnames(range))) {
            stop("The range data.table must contain chr, start and end columns")
          }
        }
        range[,width:=end-start+1]
        gr <- GRanges(seqnames = range$chr, 
          ranges = IRanges(start=range$start, 
            end=range$end,
            width=range$width))
        gr <- resize(gr, width = size, fix = fix)
        gr <- as.data.table(gr)
        colnames(gr)[which(colnames(gr) == "seqnames")] <- "chr"
        gr
       })
}

#' @param ranges: a data.table contains `chr`, `start` and `end` columns
#' @param genome: a BSgenome object, the corresponding genome 
#' @return a GRanges objects with an additional metadata column gc that contains
#' GC content
.getGCContent <- function(range, genome) {
        if (!all(c("chr", "start", "end") %in% colnames(range))) {
          stop("The range data.table must contain chr, start and end columns")
        }
        gr <- GRanges(seqnames = range$chr, 
          ranges = IRanges(start=range$start, 
            end=range$end,
            width=range$width))
        peakSeqs <- getSeq(x = genome, gr)
        mcols(gr)$gc <- letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
        gr <- as.data.table(gr)
        colnames(gr)[which(colnames(gr) == "seqnames")] <- "chr"
        gr
} 

#' add a parameter of motif or peak
#' if peak, just within counts
#' ranges is the peak ranges
#' 
.getInsertionCounts <- function(range, 
  motif,
  flankSize = 30,
  shiftATAC = FALSE
  ) {
    # Sanity check
    if (is.data.table(motif) == FALSE) {
      motif <- data.table::as.data.table(motif)
      if (!all(c("chr", "start", "end") %in% colnames(motif))) {
          stop("The motif data.table must contain chr, start and end columns")
      }
    }
    
    if (is.data.table(range) == FALSE) {
      range <- data.table::as.data.table(range)
      if (!all(c("chr", "start", "end") %in% colnames(range))) {
        stop("The range data.table must contain chr, start and end columns")
      }
    }
    
    # why?
    motifData <- data.table::copy(motif)
    atacFrag <- data.table::copy(range)
    
    if (shiftATAC == TRUE) {
      atacFrag[, start := ifelse(strand == "+", start + 4, start)]
      atacFrag[, end   := ifelse(strand == "-", end   - 5, end)]
    }
    
    # set up flanking region
    motifData[, start_margin := start - flankSize]
    motifData[, end_margin   := end + flankSize]
    motifData$motifind <- 1:nrow(motifData)
    
    # returning the overlaying counts and within motif counts
    res <- motifData[, .(motifind, 
      Start_Count = atacFrag[motifData, 
        on = .(start >= start_margin, start <= start, chr == chr), 
        .N, by = .EACHI]$N,
      End_Count   = atacFrag[motifData, 
        on = .(end >= end, end <= end_margin, chr == chr), 
        .N, by = .EACHI]$N,
      counts = atacFrag[motifData, 
        on = .(start <= start, end >= end, chr == chr), 
        .N, by = .EACHI]$N)] 

  
    res <- res[, .(motifind, 
      total_counts = Start_Count + End_Count + counts,
      flanking_counts = Start_Count + End_Count,
      within_counts = counts) ] 
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
    type = c("peaks", "motifs"),
    paired=TRUE,
    resize=TRUE, 
    width=100,
    binsize,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    ...) {
    
    type <- match.arg(type, choices=c("peaks", "motifs"))  
    
    dts <- .importFragments(files)
    
    if(type=='peak' & resiye){
      
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
