# Functions for finding kmer indices -------------------------------------------

#' get_kmer_indices
#' 
#' For each possible kmer of size k, finds which peaks contain 1 or more copies 
#' of that kmer
#' @param peak \code{\link[GenomicRanges]{GenomicRanges}} object
#' @param genome \code{\link[BSgenome]{BSgenome-class}} object, default is hg19
#' @param k length of kmer, default is 6, must be between 5 and 8
#' @return list of vectors with indices of peaks containing kmer
#' @export
get_kmer_indices <- function(peaks,
                             genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                             k = 6){
  if (k > 8){
    stop("k must be less than or equal to 8")
  }
  if (k < 5){
    stop("k must be greater than or equal to 5")
  }
  seqs <- Biostrings::getSeq(genome, peaks)
  kmers <- Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),
                                                             width = k))
  kmers <- remove_rc(kmers)
  
  out <- match_kmers(kmers, seqs)
  
  return(out)
}

match_kmers <- function(kmers, seqs, var = FALSE){
  if (is.character(kmers)) kmers = Biostrings::DNAStringSet(kmers)
  stopifnot(inherits(kmers,"DNAStringSet"))
  if (!all_true(width(kmers) == width(kmers[1])) || var){
    indices  <- Biostrings::vwhichPDict(kmers,seqs, fixed = FALSE)
    indices_rc <- Biostrings::vwhichPDict(kmers,Biostrings::reverseComplement(seqs), 
                                          fixed = FALSE)
  } else {
    pd <- Biostrings::PDict(kmers)
    indices  <- Biostrings::vwhichPDict(pd,seqs)
    indices_rc <- Biostrings::vwhichPDict(pd,Biostrings::reverseComplement(seqs))
  } 
  indices <- merge_lists(indices, indices_rc, by = "order")
  indices <- lapply(indices, unique)
    
  out <- sparseMatrix(i = unlist(lapply(seq_along(indices),
                                        function(x) rep(x, length(indices[[x]]))),use.names = FALSE),
                      j = unlist(indices, use.names = FALSE),
                      x = 1,
                      dims = c(length(seqs), length(kmers)),
                      dimnames = list(NULL, as.character(kmers)))
  
  return(out)
}



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
