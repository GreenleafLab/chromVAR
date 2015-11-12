'%ni%' = Negate('%in%')

mean_smooth <- function(X, window){
  if (is.na(as.integer(window)) || length(window) != 1 || window < 2 || window >= length(X)){
    stop("window must be an integer between 2 and length(X)")
  }
  pad_left = rev(X[1:(window %/% 2)])
  pad_right = rev(X[(length(X)-((window-1) %/% 2)):length(X)])
  cx <- c(0,cumsum(c(pad_left,X,pad_right)))
  return((cx[(window+1):(length(cx)-1)] - cx[1:(length(cx)-1-window)])/window)
}

filterSeqs <- function(seqs){
  freqs <- Biostrings::alphabetFrequency(seqs)
  seqs <- seqs[which(rowSums(freqs[,-(1:4)])==0)]
  return(seqs)
}

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

rms <- function(x) sqrt(mean(x**2))

merge_lists <- function(..., by = c("order","name")){
  by <- match.arg(by)
  if (by == "order"){
    return(Map(c,...))
  } else{
    input <- list(...)
    if (length(input) ==1){
      output = input[[1]]
      if (sum(duplicated(output)) == 0){
        return(output)
      } else{
        unames = unique(names(output))
        output = lapply(lapply(unames, function(x) l[names(l) == x]), unlist, recursive= FALSE, use.names=F)
        names(output) = unames
        return(output)
      }
    } else{
      output = merge_lists(input[[1]], by = "name")
      for (j in 2:length(input)){
        tmp = merge_lists(input[[j]], by = "name")
        common_names = names(output)[which(names(output) %in% names(tmp))]
        unique_old = names(output)[which(names(output) %ni% common_names)]
        unique_new = names(tmp)[which(names(tmp) %ni% common_names)]
        output = c(Map(c,output[common_names],tmp[common_names]),output[unique_old],tmp[unique_new])
      }
      return(output)
    }
  }
}



read_peaks <- function(filename, extra_cols = c(), extra_names = c()){
  bed <- readr::read_tsv(file = filename, col_names = FALSE)[, c(1:3, extra_cols)]
  colnames(bed) <- c("chr", "start", "end", extra_names)
  bed[,"end"] <- bed[,"end"] - 1
  bed <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  if (sum(GenomicRanges::width(bed) == GenomicRanges::width(bed[1])) != length(bed)){
    stop("All peaks in bed file must be of equal width!")
  }
  return(bed)
}

compute_gc_content <- function(seqs){
  nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
  gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
  return(gc)
}
