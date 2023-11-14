
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
          #setnames(frags, c("seqnames"), c("chr"))
          frags$count <- 1
          frags
        })
      }
      else if(grepl(".bed", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, fread, select=1:3, col.names=c("seqnames", "start", "end"))
        fragData <- lapply(fragData, function(dt){dt$count <- 1
        dt <- dt[,c("seqnames", "start", "end", "count"), with=FALSE]
        dt})
      }
      else if(grepl(".rds", fragData[[1]], fixed=TRUE))
      {
        fragData <- lapply(fragData, function(gr){
          gr <- readRDS(gr)
          dt <- as.data.table(gr)
          #setnames(dt, c("seqnames"), c("chr"))
          dt$count <- 1
          dt <- dt[,c("seqnames", "start", "end", "count")]
          dt
        })
      }
      else if(is.data.table(fragData[[1]]))
      {
        if(!(grepl(paste(c("seqnames", "start", "end"), collapse=";"),
                   paste(colnames(fragData[[1]]),collapse=";"))))
        {stop("data.table list elements of data need to contain columns: seqnames, start, end")}
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
.getInsertionCounts <- function(fragRanges, 
  motifRanges,
  #rangeType = 
  flankSize = 30,
  shiftATAC = FALSE
  ) {
    # Sanity check
    if (is.data.table(motifRanges) == FALSE) {
      motifRanges <- data.table::as.data.table(motifRanges)
    }
    
    if (is.data.table(fragRanges) == FALSE) {
      peakRanges <- data.table::as.data.table(fragRanges)
    }
    
    motifData <- data.table::copy(motifRanges)
    frags <- data.table::copy(fragRanges)
    
    if (shiftATAC == TRUE) {
      frags[, start := ifelse(strand == "+", start + 4, start)]
      frags[, end   := ifelse(strand == "-", end   - 5, end)]
    }
    
    # set up flanking region
    motifData[, start_margin := start - flankSize]
    motifData[, end_margin   := end + flankSize]
    motifData[, motifID := seq_len(nrow(motifData))]
    
    fragCounts <- lapply(frags, \(frag) {
      res <- motifData[, .(motifID, 
        start_count = frag[motifData, 
          on = .(start >= start_margin, start <= start, seqnames == seqnames), 
          .N, by = .EACHI]$N,
        end_count   = frag[motifData, 
          on = .(end >= end, end <= end_margin, seqnames == seqnames), 
          .N, by = .EACHI]$N,
        within_count = frag[motifData, 
          on = .(start >= start, end <= end, seqnames == seqnames), 
          .N, by = .EACHI]$N)] 
      res[,total_count:=start_count+end_count+within_count]
      res
    })
    
    cols <- names(fragCounts[[1]])[grepl("count", 
      names(fragCounts[[1]]))]
    allCounts <- lapply(cols, \(x) {
      lst <- lapply(fragCounts, \(.) data.frame(.)[,x])
      mat <- do.call(cbind, lst)
      colnames(mat) <- names(fragCounts)
      mat
    })
    names(allCounts) <- cols
    allCounts
}

#' @description
#' count the number of fragments in each peak region
#' 
#' @param peakRanges a GRange or data.table object that contains the peak ranges
#' @param fragRanges a list of data.table that contains the fragment ranges, 
#' each data.table represents a sample
#' return a list of count table; by all, nucleosome-free, mono...
.getOverlapCounts <- function(peakRanges, 
  fragRanges,
  by = c("count", "weight"),
  cuts,
  genome) {
  
    by <- match.arg(by, choices = c("count", "weight"))
  
    # sanity check
    if (is.data.table(peakRanges)) {
        peaks <- peakRanges
        peakGR <- makeGRangesListFromDataFrame(as.data.frame(peakRanges))
    } else if (class(peakRanges) == "GRanges") {
        peaks <- data.table::as.data.table(peakRanges)
        peakGR <- peakRanges
    }
    
    frags <- data.table::copy(fragRanges)
    peaks$peakID <- seq_len(nrow(peaks))
    
    if (by == "weight") {
      frags <- .weightFragments(frags, genome = genome)
    }
    
    frags <- .getType(frags, cuts = cuts)
    
    fragCounts <- lapply(frags, \(frag) {
      if (by == "weight") {
        frag[,count:=weight]
        frag[,(types) := lapply(.SD, function(x) x*weight), 
          .SDcols = types]
      }
      fragGR <- makeGRangesFromDataFrame(as.data.frame(frag))
      hits <- findOverlaps(fragGR, peakGR)
      overlaps <- cbind(frag[queryHits(hits),],
               peaks[subjectHits(hits), c("peakID")])
      types <- names(frag)[grepl("^type_", names(frag))]
      tmp <- overlaps[, c(list(counts = sum(count, na.rm = TRUE)), 
        lapply(.SD, sum, na.rm = TRUE)), 
        by = peakID, .SDcols = c(types)]
      
      res <- data.table(peakID = seq_len(nrow(peaks)))
      res <- merge(res, tmp, all =TRUE)
      res[is.na(res)] <- 0
      res
      
    })
    
    cols <- names(fragCounts[[1]])[grepl("^type_|counts", 
      names(fragCounts[[1]]))]
    allCounts <- lapply(cols, \(x) {
      lst <- lapply(fragCounts, \(.) data.frame(.)[,x])
      mat <- do.call(cbind, lst)
      colnames(mat) <- names(fragCounts)
      mat
    })
    names(allCounts) <- cols
    
    allCounts
    
    
}

#' change cuts into bin size; given the parameter of nBin;
#' first calculate intervals
.getType <- function(fragRanges, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(names(fragRanges), \(sample) {
        dt <- fragRanges[[sample]]
        dt[,width:=end-start+1]
        # check if max(cuts) covers max(width)
        if (max(dt$width) > cuts[length(cuts)]) 
            cuts[length(cuts)+1] <- max(dt$width)
        
        dt[,type:=as.numeric(cut(width, breaks = cuts))]
        if (length(unique(dt$type))==1) {
            warnings("All fragments fell into 1 type")
        } else if(length(unique(dt$type))>=6) {
            stop("Too many types!")
        }
        for (i in unique(dt$type)) {
          dt[, paste0("type_",i) := as.integer(type == i)]
        }
        dt

    })
    names(res) <- names(fragRanges)
    res
}


.getBins <- function(fragRanges, 
  nWidthBins = 10, 
  nGCBins = 10, 
  genome) {
    
    fragDts <- lapply(names(fragRanges), function(x){
      dt <- fragRanges[[x]]
      dt[,width:=end-start+1]
      gr <- makeGRangesFromDataFrame(as.data.frame(dt))
      gr <- .getGCContent(gr, genome = genome)
      dt <- as.data.table(gr)
      dt[,sample:=x]
      dt
    })
    fragDt <- rbindlist(fragDts)
    
    widthIntervals <- seq(min(fragDt$width), max(fragDt$width), 
      length.out=nWidthBins)
    GCIntervals <-  seq(0, 1, length.out=nGCBins)
    
    fragDt[,widthBin:=cut(width, 
        breaks=widthIntervals, 
        include.lowest=TRUE)]
    fragDt[,GCBin:=as.numeric(cut(gc, 
        breaks=GCIntervals, 
        include.lowest=TRUE))]
    fragDt
    
}
  
  

## bin by both fragment length and GC content


.weightFragments <- function (fragRanges, 
  #cuts = c(0,120,300,500),
  genome,
  ...) {
  
    fragDt <- .getBins(fragRanges, genome = genome)
    fragDt[, bin:=paste0(widthBin, GCBin)]
    fragDt[, bin:=as.numeric(as.factor(bin))]
    fragDt[,count_bin:=.N, by=c("sample", "bin")]
    tmp <- fragDt[,.(mean_count_bin=mean(count_bin, na.rm=TRUE)), by=c("bin")]
    fragDt <- merge(fragDt, tmp, by = "bin")
    fragDt[,weight:=mean_count_bin/count_bin]
    res <- split(fragDt, fragDt$sample)
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
    ranges, # motif matches or peaks, set a parameter to define
    genome,
    rowType = c("peaks", "motifs"),
    by = c("count", "weight"),
    paired = TRUE,
    resize = TRUE, 
    width = 200,
    nWidthBins = 20,
    nGCBins = 20,
    cuts = c(0,120,300,500),
    ...) {
    
    rowType <- match.arg(rowType, choices=c("peaks", "motifs"))  
    by <- match.arg(by, choices=c("count", "weight"))  
      
    fragRanges <- .importFragments(files)

    
    if (rowType=='peaks' & resize) {
      ranges <- .resizeRanges(peakRanges = ranges, width = width)
    }
    
    if (rowType=='peaks') {
        asy <- .getOverlapCounts(peakRanges = ranges, 
          fragRanges = fragRanges,
          by = by,
          cuts = cuts,
          genome = genome)
    }
    #else if(type=="motif"){
      #asy <- .getInsertionCounts(fragRanges, motifRanges = )
    #}
    
    SummarizedExperiment(assays = asy, rowRanges = ranges)

}

# Sys.time()
# se <- getCounts(beds, ranges = peakRanges,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   by = "weight",
#   rowType = "peaks",
#   cuts = c(40,60,80))
# Sys.time()
