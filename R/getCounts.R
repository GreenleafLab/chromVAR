
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
    
    if ("seqnames" %in% colnames(motifRanges))
      colnames(motifRanges)[which(colnames(motifRanges)=="seqnames")] <- "chr"
    
    if (is.data.table(fragRanges) == FALSE) {
      peakRanges <- data.table::as.data.table(fragRanges)
    }
  
    if ("seqnames" %in% colnames(fragRanges))
      colnames(fragRanges)[which(colnames(fragRanges)=="seqnames")] <- "chr"
    
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
          on = .(start >= start_margin, start <= start, chr == chr), 
          .N, by = .EACHI]$N,
        end_count   = frag[motifData, 
          on = .(end >= end, end <= end_margin, chr == chr), 
          .N, by = .EACHI]$N,
        within_count = frag[motifData, 
          on = .(start >= start, end <= end, chr == chr), 
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
#  shiftATAC = FALSE,
  by = c("count", "weight"),
  cuts,
  genome = genome) {
  
    by <- match.arg(by, choices = c("count", "weight"))
  
    # sanity check
    if (is.data.table(peakRanges) == FALSE) 
      peakRanges <- data.table::as.data.table(peakRanges)
    
    # if (is.data.table(fragRanges) == FALSE) 
    #   fragRanges <- data.table::as.data.table(fragRanges)
    
    if ("seqnames" %in% colnames(peakRanges))
      colnames(peakRanges)[which(colnames(peakRanges)=="seqnames")] <- "chr"
    
    # if ("seqnames" %in% colnames(fragRanges))
    #   colnames(fragRanges)[which(colnames(fragRanges)=="seqnames")] <- "chr"
    
    
    # if (shiftATAC == TRUE) {
    #   peakRanges[, start := ifelse(strand == "+", start + 4, start)]
    #   peakRanges[, end   := ifelse(strand == "-", end   - 5, end)]
    # }  
    
    peaks <- data.table::copy(peakRanges)
    frags <- data.table::copy(fragRanges)
    
    peaks$peakID <- 1:nrow(peaks)
    
    if (by == "count") {
      frags <- .getType(frags, cuts = cuts)
      fragCounts <- lapply(frags, \(frag) {
        tmp <- frag[peaks, 
               on = .(start >= start, end <= end, chr == chr)]
        
        types <- names(tmp)[grepl("^type_", names(tmp))]
        
        res <- tmp[, c(list(counts = sum(count, na.rm = TRUE)), 
          lapply(.SD, sum, na.rm = TRUE)), 
          by = peakID, .SDcols = c(types)]
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
    } else if (by == "weight") {
          frags <- .weightFragments(frags, genome = genome)
          frags <- .getType(frags, cuts = cuts)
          fragCounts <- lapply(frags, \(frag) {
            colnames(frag)[which(colnames(frag) == "seqnames")] <- "chr"
            types <- names(frag)[grepl("^type_", names(frag))]
            sel <- c("chr", "start", "end", "weight", types)
            frag <- frag[,..sel]
            frag$count <- 1
            tmp <- frag[peaks, 
              on = .(start >= start, end <= end, chr == chr)]
            tmp[,count:=count*weight]
            tmp[,(types) := lapply(.SD, function(x) x * weight), .SDcols = types]
            
            res <- tmp[, c(list(counts = sum(count, na.rm = TRUE)), 
              lapply(.SD, sum, na.rm = TRUE)), 
              by = peakID, .SDcols = c(types)]
            
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
    }
    
    # add count by type: 
    
    
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

## test
motifs <- readRDS("/Volumes/jiayiwang/chromVAR2/data/motif.rds")
mgr <- lapply(names(motifs), \(x) {
  gr <- as.data.table(motifs[[x]])
  gr$motif <- x
  gr
})

motifRanges <- rbindlist(mgr)




