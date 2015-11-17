
get_kmer_indices <- function(peaks, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, k = 6){
  
  if (k > 8){
    stop("k must be less than 8")
  }
  if (k < 5){
    stop("k must be greater than or equal to 5")
  }
  seqs <- Biostrings::getSeq(genome, peaks)
  kmers <- Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))
  pd <- Biostrings::PDict(kmers)
  indices  <- Biostrings::vwhichPDict(pd,seqs)
  
  tmp <- data.frame(peak_ix = unlist(lapply(1:length(indices), function(x) rep(x, length(indices[[x]])))), kmer_ix = factor(unlist(indices), levels = 1:length(kmers), ordered=T))
  out <- split(tmp$peak_ix,tmp$kmer_ix)
  names(out) <- as.character(kmers)
  
  #remove reverse complements
  tmp <- cbind(names(out), as.character(Biostrings::reverseComplement(kmers)))
  names(out)[tmp[,1] > tmp[,2]] <- tmp[tmp[,1] > tmp[,2],2]
  out <- merge_lists(out, by = "name")
  
  return(out)
}

get_gapped_kmer_indices <- function(peaks, k, m){
  
  l = k + m
  lmer_indices = get_kmer_indices(peaks, k = l)

  all_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))

  Ns = paste(rep("N",m),collapse="")

  gapped_kmer = Biostrings::replaceAt(all_kmer,1,value=Ns)

  for (i in 2:(k+1)){
    gapped_kmer = c(gapped_kmer,Biostrings::replaceAt(all_kmer,i,value=Ns))
  }

  gapped_kmer = remove_rc(gapped_kmer)

  pd = Biostrings::PDict(Biostrings::DNAStringSet(names(lmer_indices)))
  mapping  = Biostrings::vwhichPDict(pd,gapped_kmer,fixed = "pattern")
  mapping.rc  = Biostrings::vwhichPDict(pd,Biostrings::reverseComplement(gapped_kmer),fixed = "pattern")

  mapping = merge_lists(mapping, mapping.rc)
  
  out <- lapply(1:length(mapping), function(x) unique(unlist(lmer_indices[mapping[[x]]],use.names=FALSE)))
  
  names(out) <- as.character(gapped_kmer)
  
  return(out)
}


remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs), as.character(Biostrings::reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  Biostrings::DNAStringSet(unique(temp[,1]))
}




# 
# 
# ###Find kmer that are over-represented
# 
# #Find expectation
# n_kmers = (width(seqs[1]) - k + 1) * length(seqs)
# nuc_freqs = colSums(Biostrings::letterFrequency(seqs,c("A","C","G","T")))
# nuc_freqs = nuc_freqs/sum(nuc_freqs)
# #nuc_freqs = c(nuc_freqs, c("N"=1)) #* (n_kmers ** (1/k))
# kmer_nuc_freqs = Biostrings::letterFrequency(all_kmer,c("A","C","G","T"))
# kmer_exp =  apply(nuc_freqs ** t(kmer_nuc_freqs),2,prod) * n_kmers
# wc_kmer_exp = as.numeric(kmer_exp %*% kmer.match.mat)
# wc_kmer_counts = colSums(seq_kmer_mat)
# wc_kmer_enrich = wc_kmer_counts/wc_kmer_exp
# 
# 
# sparse_jaccard <- function(x){
# 
#   x[which(x>1)]=1
#   i = t(x) %*% x
#   tmp = colSums(x)
#   u = matrix(tmp, nrow = nrow(i), ncol=ncol(i), byrow=T) + matrix(tmp, nrow = nrow(i), ncol=ncol(i), byrow=F) - i
#   jac = i/u
#   return(jac)
#   #return(as.dist(jac))
# 
# }
# 
# tmp = summary(seq_kmer_mat)
# seq_kmer_mat_bin = sparseMatrix(i = tmp[,1], j = tmp[,2], x = 1, dims = dim(seq_kmer_mat), dimnames = list(NULL,colnames(seq_kmer_mat)))
# 
# sparse_jaccard_help <- function(y,x){
# 
#   i = t(y) %*% x
#   u = outer(colSums(y), colSums(x), FUN = "+") - i
#   jac = i/u
#   return(jac)
#   #return(as.dist(jac))
# 
# }
# 
# 
# 
# wc_kmer_dist = vegan::vegdist(seq_kmer_mat[,which(wc_kmer_enrich>=2)], method="jaccard",binary=T)
# 
# 
# reduce_kmer <- function(){
#   max_jaccard = 0.9
# 
# 
# 
# 
# 
# }
# 
# 
# 
# 
# 
# add_wildcard_assoc <- function()
# 
# 
# 
# 
# wildcardInit <- function(seqs, k, m){
# 
# }
# 
# wildCardSub <- function(seqs, )
