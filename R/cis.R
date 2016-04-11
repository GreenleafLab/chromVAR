
#'@importMethodsFrom GenomeInfoDb seqlevels seqnames
#'@export
get_cis_groups <- function(peaks, grpsize = 25, stepsize = 10){
  
  if (!isTRUE(all.equal(peaks, sort(peaks)))){
    stop("peaks must be sorted to be able to get cis groups")
  }
  
  chrs = seqlevels(peaks)
  out = list()
  out = do.call(c, BiocParallel::bplapply(seq_along(chrs), function(i) {
    chr_ix = which(seqnames(peaks) == chrs[i])
    if (length(chr_ix) > stepsize){
      tmp = lapply(1:(length(chr_ix) %/% stepsize), function(x) chr_ix[((x-1)*stepsize +1):(min((x-1)*stepsize + grpsize,length(chr_ix)))])
      names(tmp) = sapply(1:length(tmp), function(x) paste(chrs[i], x,sep="_",collapse=""))
      return(tmp)
    } else{
      return(list())
    }
  }))
  return(out)
}

compute_cis_correlations <- function(counts, peaks){
  
  
  
  
}