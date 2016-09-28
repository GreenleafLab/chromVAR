




get_gapped_kmer_indices <- function(peaks, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, k = 6, m = 1){
  
  l = k + m

  lmer_indices <- get_kmer_indices(peaks, genome, k = l)
  
  all_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),
                                                               width = k))
  
  Ns = paste(rep("N",m),collapse="")
  
  gapped_kmer = Biostrings::replaceAt(all_kmer,1,value=Ns)
  
  for (i in 2:(k+1)){
    gapped_kmer = c(gapped_kmer,Biostrings::replaceAt(all_kmer,i,value=Ns))
  }
  
  gapped_kmer = remove_rc(gapped_kmer)
  
  pd = Biostrings::PDict(Biostrings::DNAStringSet(colnames(lmer_indices)))
  mapping  = Biostrings::vwhichPDict(pd,gapped_kmer,fixed = "pattern")
  mapping.rc  = Biostrings::vwhichPDict(pd,
                                        Biostrings::reverseComplement(gapped_kmer),
                                        fixed = "pattern")
  
  mapping = merge_lists(mapping, mapping.rc)
  
 map_mat <- sparseMatrix(j = unlist(lapply(seq_along(mapping),
                                        function(x) rep(x, length(mapping[[x]]))),use.names = FALSE),
                      i = unlist(mapping, use.names = FALSE),
                      x = 1,
                      dims =  c(ncol(lmer_indices), length(mapping)))
 out <- lmer_indices %*% map_mat
 colnames(out) = as.character(gapped_kmer)
 out@x = rep(1, length(out@x)) 
 
  return(out)
}

remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs),
                as.character(Biostrings::reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  Biostrings::DNAStringSet(unique(temp[,1]))
}
