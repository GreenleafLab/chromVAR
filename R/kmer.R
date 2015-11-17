#' @export
get_kmer_indices <- function(peaks, 
                             genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                             k = 6){
  
  if (k > 8){
    stop("k must be less than 8")
  }
  if (k < 5){
    stop("k must be greater than or equal to 5")
  }
  seqs <- Biostrings::getSeq(genome, peaks)
  kmers <- Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),
                                                             width = k))
  pd <- Biostrings::PDict(kmers)
  indices  <- Biostrings::vwhichPDict(pd,seqs)
  
  tmp <- data.frame(peak_ix = unlist(lapply(1:length(indices), 
                                            function(x)rep(x, length(indices[[x]])))), 
                    kmer_ix = factor(unlist(indices),
                                     levels = 1:length(kmers),
                                     ordered=T))
  out <- split(tmp$peak_ix,tmp$kmer_ix)
  names(out) <- as.character(kmers)
  
  #remove reverse complements
  tmp <- cbind(names(out), as.character(Biostrings::reverseComplement(kmers)))
  names(out)[tmp[,1] > tmp[,2]] <- tmp[tmp[,1] > tmp[,2],2]
  out <- merge_lists(out, by = "name")
  
  return(out)
}



# Gapped k-mers ----------------------------------------------------------------

#' @export
get_gapped_kmer_indices <- function(peaks, k, m){
  
  l = k + m
  lmer_indices = get_kmer_indices(peaks, k = l)

  all_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),
                                                               width = k))

  Ns = paste(rep("N",m),collapse="")

  gapped_kmer = Biostrings::replaceAt(all_kmer,1,value=Ns)

  for (i in 2:(k+1)){
    gapped_kmer = c(gapped_kmer,Biostrings::replaceAt(all_kmer,i,value=Ns))
  }

  gapped_kmer = remove_rc(gapped_kmer)

  pd = Biostrings::PDict(Biostrings::DNAStringSet(names(lmer_indices)))
  mapping  = Biostrings::vwhichPDict(pd,gapped_kmer,fixed = "pattern")
  mapping.rc  = Biostrings::vwhichPDict(pd,
                                        Biostrings::reverseComplement(gapped_kmer),
                                        fixed = "pattern")

  mapping = merge_lists(mapping, mapping.rc)
  
  out <- lapply(1:length(mapping), 
                function(x) unique(unlist(lmer_indices[mapping[[x]]],
                                          use.names=FALSE)))
  
  names(out) <- as.character(gapped_kmer)
  
  return(out)
}

#Don't export
remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs), 
                as.character(Biostrings::reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  Biostrings::DNAStringSet(unique(temp[,1]))
}

