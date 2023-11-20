
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
          start(frags) <- ifelse(strand(frags) == "+", 
            start(frags) + 4, start(frags))
          end(frags) <- ifelse(strand(frags) == "-", 
            end(frags) - 5, end(frags))
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
#' Sanity check to ensure the input arguments have the correct classes
#' 
#' @param atacFrag: a list of data.tables that contain the ranges of fragments
#' @param peakRanges: a GRange object that contains the ranges of peaks
#' @param motifRanges:  a GRange object that contains the ranges of motifs, 
#' check metadata columns
#' check seqnames to factor in datatable
.sanityCheck <- function(atacFrag, 
  ranges, 
  type = c("peaks", "motifs")) {
  nPeaks <- nrow(atacFrag[[1]])
  lapply(atacFrag, \(frag) {
      if (!is.data.table(frag)) {
        stop("Each element in atacFrag should be a data.table")
      }
      if(nrow(frag) != nPeaks) {
        stop("The number of peaks does not match between samples")
      }
  })
  
  if (!class(ranges)=="GRanges") {
    stop("ranges should be a GRanges object")
  }
  
  type <- match.arg(type, choices = c("peaks", "motifs"))
  if (type=="motifs") {
    if (!("motif" %in% names(mcols(ranges)))) {
      stop("There is no motif names in metadata columns")
    }
  }
    
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


.filterFrags <- function(atacFrag, min = 30, max = 2000) {
    res <- lapply(atacFrag, \(frag) {
        frag[,width:=end-start+1]
        frag[width>=min & width<=max,]
        frag
    })
    res
}

.matchSeqlevels <- function(atacFrag, ranges) {
    frags <- rbindlist(atacFrag)
    fragSeq <- unique(frags$seqnames)
    rangeSeq <- GenomicRanges::seqnames(ranges) 
    common <- intersect(fragSeq,rangeSeq)
    atacFrag <- lapply(atacFrag, \(frag) {
     frag <- frag[seqnames %in% common,]
     frag$seqnames <- factor(frag$seqnames)
     frag})
    ranges <- ranges[seqnames(ranges) %in% common,]
    # turn seqnames to factor
    list(atacFrag=atacFrag, ranges=ranges)
}

#TODO: Functionality
# - filter max and min length of fragments (30 & 2000 as defaults)
# - match seqLevels of ranges and atac frags in the beginning
# - write helper functions for sanity checks, i.e. for data.table / granges conversions etc
# - write a helper function for ATAC shifting
# - conversion seqLevels to integer

#TODO: Speed-up & memory usage
# - use dtToGr helper for data.table to GenomicRanges conversion
# - overlap with findOverlaps and use indices 
# - convert samples (!), chr & motif ids to integers / numerical factors wherever possible
# - test chunking across chromosomes: Different chunk sizes

#TODO:
# 1. unify the arguments 
# 2. filter functions
# 3. seqLevels

#' @Author: Emanuel Sonder
dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end"){
  setnames(dt, seqCol, "seqnames", skip_absent = TRUE)
  gr <- GRanges(seqnames=dt[[seqCol]], ranges=IRanges(start=dt[[startCol]], 
                                                      end=dt[[endCol]]))
  if(startCol==endCol)
  {
    gr <- GPos(seqnames=dt[[seqCol]], pos=dt[[startCol]])
  }
  
  return(gr)
}

#' @Author: Emanuel Sonder
.getInsertionProfiles <- function(atacFrag, 
                                  motifRanges,
                                  margin=100,
                                  aggFun=sum,
                                  #minWidth=30,
                                  #maxWidth=2000,
                                  chunk=TRUE){
  
  # prep motif data
  motifData <- as.data.table(motifRanges)
  #setnames(motifData, "seqnames", "chr")
  chrLevels <- unique(motifData$chr)
  motifLevels <- unique(motifData$motif_id)
  
  # convert to factors (memory usage)
  motifData[,chr:=as.integer(factor(chr, 
    levels=chrLevels, ordered=TRUE))]
  motifData[,motif_id:=as.integer(factor(motif_id, 
    levels=motifLevels, ordered=TRUE))]
  
  # determine margins
  motifData[,start_margin:=start-margin]
  motifData[,end_margin:=end+margin]
  
  # determine motif center
  motifData[,motif_center:=floor((end_margin-start_margin)/2)+start_margin]
  
  # convert to factors (memory usage)
  atacFrag <- copy(atacFrag) #TODO: take out that copy 
  atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
  nSamples <- length(unique(atacFrag$sample))
  
  if(chunk){
    setorder(motifData, chr)
    setorder(atacFrag, chr)
    motifData <- split(motifData, by="chr")
    atacFrag <- split(atacFrag, by="chr")
    
    atacInserts <- mapply(function(md,af){
      
      # convert to granges for faster overlapping
      motifMarginRanges <- dtToGr(md, startCol="start_margin", 
        endCol="end_margin")
      atacStartRanges <- dtToGr(af, startCol="start", endCol="start")
      atacEndRanges <- dtToGr(af, startCol="end", endCol="end")
      
      startHits <- findOverlaps(atacStartRanges, 
        motifMarginRanges, type="within") # check if type within faster or slower
      endHits <- findOverlaps(atacEndRanges, motifMarginRanges, type="within") 
      
      # get overlapping insertion sites
      atacStartInserts <- af[queryHits(startHits), 
        c("sample", "start"), with=FALSE]
      atacEndInserts <-af[queryHits(endHits), c("sample", "end"), with=FALSE]
      setnames(atacStartInserts, "start", "insert")
      setnames(atacEndInserts, "end", "insert")
      
      ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
                  rbindlist(list(
                    md[subjectHits(startHits), c("motif_center", "motif_id")],
                    md[subjectHits(endHits), c("motif_center", "motif_id")])))
      
      # count insertions around motif
      ai[,rel_pos:=abs(insert-motif_center)]
      ai[,.(pos_count_sample=.N), by=.(motif_id, rel_pos, sample)]}, 
      motifData, 
      atacFrag, 
      SIMPLIFY=FALSE)
    
    # combine insertion counts across chromosomes
    atacInserts <- rbindlist(atacInserts)
    atacInserts <- atacInserts[,.(pos_count_sample=sum(pos_count_sample)), 
                               by=.(motif_id, rel_pos, sample)]
  }
  else
  {
    # convert to granges for faster overlapping
    motifMarginRanges <- dtToGr(motifData, startCol="start_margin", 
      endCol="end_margin")
    atacStartRanges <- dtToGr(atacFrag, startCol="start", endCol="start")
    atacEndRanges <- dtToGr(atacFrag, startCol="end", endCol="end")
    
    startHits <- findOverlaps(atacStartRanges, motifMarginRanges, type="within") # check if type within faster or slower
    endHits <- findOverlaps(atacEndRanges, motifMarginRanges, type="within") 
    
    # get overlapping insertion sites
    atacStartInserts <- atacFrag[queryHits(startHits),
      c("sample", "start"),with=FALSE]
    atacEndInserts <-atacFrag[queryHits(endHits),
      c("sample", "end"),with=FALSE]
    setnames(atacStartInserts, "start", "insert")
    setnames(atacEndInserts, "end", "insert")
    
    atacFragInserts <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
                             rbindlist(list(
                               motifData[subjectHits(startHits),
                                 c("motif_center", "motif_id")],
                               motifData[subjectHits(endHits),
                                 c("motif_center", "motif_id")])))
    
    # count insertions around motif
    atacFragInserts[,rel_pos:=abs(insert-motif_center)]
    atacInserts <- atacFragInserts[,.(pos_count_sample=.N), 
                                   by=.(motif_id, rel_pos, sample)]
  }
  
  atacInserts[,pos_count_global:=sum(pos_count_sample), by=.(motif_id, rel_pos)]
  
  # ensure ordering of positions
  setorder(atacInserts, motif_id, rel_pos)
  
  # calculate and smooth weights
  atacInserts[,w:=smooth(pos_count_global/sum(pos_count_global),
                         twiceit=TRUE), by=motif_id]
  
  # calculate per sample motif scores
  atacInserts[,score:=w*pos_count_sample]
  
  motifAct <- atacInserts[,.(activity=aggFun(score)), by=.(sample, motif_id)]
  actMat <- data.table::dcast(motifAct, motif_id~sample, value.var="activity")
  setorder(actMat, motif_id)
  rownames(actMat) <- motifLevels
  
  se <- SummarizedExperiment(list(scores=actMat[,setdiff(colnames(actMat), 
    "motif_id"), with=FALSE]))
  S4Vectors::metadata(se)$globalProfiles <- atacInserts
  
  return(se)
}


# Comment: 
# - I removed the sanity checks as we either write a seperate function for that 
# or we do it once in the main function at the beginning
# - same for the atac shifts

.atacShift <- function(atacFrag, shifts=c(4L, 5L)){
  atacFrag[, start := ifelse(strand == "+", start + shifts[1], start)]
  atacFrag[, end   := ifelse(strand == "-", end   - shifts[2], end)]
  atacFrag
}

# something like this: might add more checks
.checkRanges <- function(ranges){
  # Sanity check
  if (!is.data.table(ranges)) {
    dt <- data.table::as.data.table(ranges)
  } else {
    dt <- data.table::copy(ranges)
  }
  
  # if ("seqnames" %in% colnames(dt))
  #   colnames(dt)[which(colnames(dt)=="seqnames")] <- "chr"
  
  return(dt)
}

#' @param atacFrag a GRange or data.table object that contains atac fragment ranges
#' @param motifRanges a GRange or data.table object that contains motif ranges
#' @param flankSize integer, the number of nucleotides to define the buffer/flanking 
#' region near the motif instance
#' @param shiftATAC logic if shifting the ATAC-seq fragment data table
#' @Author: Emanuel Sonder
.getInsertionCounts <- function(atacFrag, 
                                motifRanges,
                                mode=c("total", "weight"),
                                flankSize=30,
                                shiftATAC=FALSE,
                                weightCol="weight", 
                                #addProfile=FALSE,
                                ...) {
  
  mode <- match.arg(mode, choices=c("total", "weight"))
  
  atacFrag <- .checkRanges(atacFrag)
  motifData <- .checkRanges(motifRanges)
  
  # ATAC insertion shifts
  if (shiftATAC == TRUE) {
    atacFrag <- .atacShift(atacFrag, ...)
  }
  
  # set up flanking region
  motifData[, start_margin := start - flankSize]
  motifData[, end_margin   := end + flankSize]
  motifData$match_id <- seq_len(nrow(motifData))
  
  # Either count weights or total number of fragments
  if(mode=="total"){
    atacFrag$weight <- 1
  } else {
    setnames(atacFrag, weightCol, "weight")
  }
  
  # get number of samples and motif matches
  nSample <- length(unique(atacFrag$sample))
  nMatch <- length(unique(motifData$match_id))
  
  # maybe its faster to chunk that part across samples -------------------------
  # convert to granges for faster overlaps
  fragStarts <- dtToGr(atacFrag, startCol="start", endCol="start")
  fragEnds <- dtToGr(atacFrag, startCol="end", endCol="end")
  motifMarginRanges <- dtToGr(motifData, 
                              startCol="start_margin", endCol="end_margin")
  
  # refactor this: Try 
  # - with sparse matrix depending on zero fraction
  # - findOverlaps invert=TRUE
  # - (!) if this is run in a loop by sample also countByOverlaps could be used
  startInsertsMargin <- as.data.table(findOverlaps(motifMarginRanges, fragStarts))
  startInsertsMargin <- cbind(startInsertsMargin, 
                              atacFrag[startInsertsMargin$subjectHits, 
                                c("sample", "weight")])
  endInsertsMargin <- as.data.table(findOverlaps(motifMarginRanges, fragEnds))
  endInsertsMargin <- cbind(endInsertsMargin, 
                            atacFrag[endInsertsMargin$subjectHits, 
                              c("sample", "weight")])
  
  # get matrix with inserts within the motif and the margin
  insertsMargin <- rbind(startInsertsMargin, endInsertsMargin)
  # get zero inserts
  insertsMargin <- rbind(insertsMargin, data.table(
    queryHits=setdiff(seq_len(length(motifRanges)), 
      unique(insertsMargin$queryHits)),
    weight=0, 
    sample=insertsMargin$sample[1]), # give some dummy sample
    use.names=TRUE, fill=TRUE)
  
  insertsMargin <- dcast(insertsMargin, queryHits~sample, 
                         value.var="weight",
                         fill=0,
                         fun.aggregate=sum, drop=FALSE)
  insertsMargin$queryHits <- NULL
  
  # get the inserts within the motif
  startInsertsWithin <- as.data.table(findOverlaps(motifRanges, fragStarts))
  startInsertsWithin <- cbind(startInsertsWithin, 
                              atacFrag[startInsertsWithin$subjectHits, 
                                c("sample", "weight")])
  endInsertsWithin <- as.data.table(findOverlaps(motifRanges, fragEnds))
  endInsertsWithin <- cbind(endInsertsWithin, 
                            atacFrag[endInsertsWithin$subjectHits, 
                              c("sample", "weight")])
  
  insertsWithin <- rbind(startInsertsWithin, endInsertsWithin)
  # get missing combinations 
  insertsWithin <- rbind(insertsWithin, data.table(
    queryHits=setdiff(seq_len(length(motifRanges)), 
      unique(insertsWithin$queryHits)),
    weight=0, 
    sample=insertsWithin$sample[1]), # give some dummy sample
    use.names=TRUE, fill=TRUE)
  
  insertsWithin <- dcast(insertsWithin, queryHits~sample, 
                         value.var="weight",
                         fill=0,
                         fun.aggregate=sum, drop=FALSE)
  insertsWithin$queryHits <- NULL
  
  # ----------------------------------------------------------------------------
  
  flankingCounts <- insertsMargin - insertsWithin
  
  return(list(total_counts=insertsMargin, 
              flanking_counts=flankingCounts, 
              within_counts=insertsWithin))
}

#' @description
#' count the number of fragments in each peak region
#' 
#' @param peakRanges a GRange or data.table object that contains the peak ranges
#' @param atacFrag a list of data.table that contains the fragment ranges, 
#' each data.table represents a sample
#' return a list of count table; by all, nucleosome-free, mono...
.getOverlapCounts <- function(peakRanges, 
  atacFrag,
  mode = c("total", "weight"),
  cuts,
  genome,
  overlap = c("any", "start", "end", "within", "equal"),
  ...
  #overlaps = c("any",...)
  ) {
  
    by <- match.arg(mode, choices = c("total", "weight"))
      
    peaks <- data.table::as.data.table(peakRanges)
    peakGR <- peakRanges
    
    frags <- data.table::copy(atacFrag)
    peaks$peakID <- seq_len(nrow(peaks))
    
    if (by == "weight") {
      frags <- .weightFragments(frags, genome = genome)
    }
    
    frags <- .getType(frags, cuts = cuts)
    
    fragCounts <- lapply(frags, \(frag) {
      types <- names(frag)[grepl("^type_", names(frag))]
      if (by == "weight") {
        frag[,count:=weight]
        frag[,(types) := lapply(.SD, function(x) x*weight), 
          .SDcols = types]
      }
      #fragGR <- makeGRangesFromDataFrame(as.data.frame(frag))
      fragGR <- dtToGr(frag)
      hits <- findOverlaps(fragGR, peakGR, type = overlap)
      overlaps <- cbind(frag[queryHits(hits),],
               peaks[subjectHits(hits), c("peakID")])
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
.getType <- function(atacFrag, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(names(atacFrag), \(sample) {
        dt <- atacFrag[[sample]]
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
    names(res) <- names(atacFrag)
    res
}


.getBins <- function(atacFrag, 
  nWidthBins = 10, 
  nGCBins = 10, 
  genome) {
    
    fragDts <- lapply(names(atacFrag), function(x){
      dt <- atacFrag[[x]]
      dt[,width:=end-start+1]
      #gr <- makeGRangesFromDataFrame(as.data.frame(dt))
      gr <- dtToGr(dt)
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


.weightFragments <- function (atacFrag, 
  #cuts = c(0,120,300,500),
  genome,
  ...) {
  
    fragDt <- .getBins(atacFrag, genome = genome)
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
#' @param ranges a GRanges that contains the ranges of peaks or motifs. If 
#' rowType == "motifs", ranges must contain a mcols called `motif` that contains
#' the name of motif for each row
#' @param rowType using peaks or motifs in rows
#' @param mode if the assays return original counts or weighted counts
#' @param paired paired end data?
#' @param resize a logical to determine if peaks should be resized
#' @param width if resize is `TRUE`, the new width of each peak
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}}
#'  object

getCounts <- function (files,
    ranges, # motif matches or peaks, set a parameter to define
    genome,
    rowType = c("peaks", "motifs"),
    mode = c("total", "weight"),
    paired = TRUE,
    resize = TRUE, 
    width = 200,
    nWidthBins = 20,
    nGCBins = 20,
    cuts = c(0,120,300,500),
    minFrag = 30,
    maxFrag = 3000,
    ...) {
    
    rowType <- match.arg(rowType, choices=c("peaks", "motifs"))  
    mode <- match.arg(mode, choices=c("total", "weight"))  
    # get fragments ranges
    atacFrag <- .importFragments(files)
    
    # sanity check
    .sanityCheck(atacFrag, ranges, type = rowType)
    
    # filter too short or too long fragments
    atacFrag <- .filterFrags(atacFrag, min = minFrag, max = maxFrag)
    
    # match seqLevels
    res <- .matchSeqlevels(atacFrag, ranges)
    atacFrag <- res$atacFrag
    ranges <- res$ranges

    
    if (rowType=='peaks' & resize) {
      ranges <- .resizeRanges(peakRanges = ranges, width = width)
    }
    
    if (rowType=='peaks') {
        asy <- .getOverlapCounts(peakRanges = ranges, 
          atacFrag = atacFrag,
          mode = mode,
          cuts = cuts,
          genome = genome)
    } else if (rowType=="motifs"){
        if (mode=="weight") {
          atacFrag <- .weightFragments(atacFrag, genome)
        }
        lst <- lapply(names(atacFrag), \(x) {
            frag <- atacFrag[[x]]
            frag$sample <- x
            frag
        })
        atacFrags <- rbindlist(lst)
        asy <- .getInsertionCounts(atacFrags, 
          motifRanges = ranges, 
          mode = mode)
    }
    
    SummarizedExperiment(assays = asy, rowRanges = ranges)

}

# files <- list("ctrl1"="ctrl1.bed", "ctrl2"="ctrl2.bed",
# "treat1"="treat1.bed", "treat2"="treat2.bed")

# motifData <- lapply(names(motifRanges), \(x) {
#   gr <- motifRanges[[x]]
#   mcols(gr)$motif <- x
#   gr
# })
# 
# motifs <- bind_ranges(motifData)


# s <- Sys.time()
# se <- getCounts(files, ranges = motifs, rowType = "motifs", mode = "weight",
#     genome = BSgenome.Hsapiens.UCSC.hg38)
# e <- Sys.time()
# e-s

# compare t-values of perturbed motifs to original chromVAR
