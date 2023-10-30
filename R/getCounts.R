
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

.resizeRanges <- function(){}

.getGCCountent <- function(){} 

.getInsertionCounts <- function(){}

.getFLD <- function(dts, cuts=c(0,120,300,500)) {
    # check the min of cuts
    if (cuts[1] != 0) cuts <- c(0,cuts)
    res <- lapply(names(dts), \(sample) {
        dt <- dts[[sample]]
        dt[,width:=end-start]
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

getCounts <- function (data,
                      ranges,
                      paired=TRUE,
                      resize=TRUE, 
                      width=100,
    ...) {
    
    dt <- .importFragments(data)

}
