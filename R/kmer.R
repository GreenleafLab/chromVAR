

compute_kmer_assoc <- function(bed, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, k = 6){
  seqs <- Biostrings::getSeq(genome,bed)
  rc_seqs <- Biostrings::reverseComplement(seqs)
  assoc <- Biostrings::oligonucleotideFrequency(seqs, k) + Biostrings::oligonucleotideFrequency(rc_seqs,k)
  kmer_indices <- lapply(1:ncol(assoc), function(x) which(assoc[,x]>0))
  names(kmer_indices) <- colnames(assoc)
  ###Reduce by removing reverse complements
  kmer = Biostrings::DNAStringSet(names(kmer_indices))
  rc_kmer = Biostrings::reverseComplement(kmer)
  dups = sapply(1:length(kmer_indices), function(i) kmer[i] %in% rc_kmer[1:(i-1)])
  kmer_indices = kmer_indices[which(!dups)]
  return(kmer_indices)
}



'%ni%' = Negate('%in%')

remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs), as.character(reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  DNAStringSet(unique(temp[,1]))
}

make_gapped_kmer <-function(k, m){
  
  all_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))
  all_lmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k+m))
  
  Ns = paste(rep("N",m),collapse="")
  
  gapped_kmer = replaceAt(all_kmer,1,value=Ns)
  
  for (i in 2:(k+1)){
    gapped_kmer = c(gapped_kmer,replaceAt(all_kmer,i,value=Ns))
  }
  
  gapped_kmer = remove_rc(gapped_kmer)
    
  pd = PDict(all_lmer)
  mapping  = vwhichPDict(pd,gapped_kmer,fixed = "pattern")
  mapping.rc  = vwhichPDict(pd,Biostrings::reverseComplement(gapped_kmer),fixed = "pattern")
  
  #rows are kmer, col are wc-kmer
  kmer.match.mat = sparseMatrix(i = c(unlist(mapping),unlist(mapping.rc)),
                                j = c(unlist(sapply(1:length(mapping), function(x) rep(x,length(mapping[[x]])))),
                                      unlist(sapply(1:length(mapping.rc), function(x) rep(x,length(mapping.rc[[x]]))))),
                                x = 1, dims = c(length(all_lmer),length(gapped_kmer)), use.last.ij = TRUE,
                                dimnames = list(as.character(all_lmer), as.character(gapped_kmer))
  )
  
  kmer.col.mat = rep(colnames(kmer.match.mat),2)
  names(kmer.col.mat) = c(colnames(kmer.match.mat), as.character(Biostrings::reverseComplement(gapped_kmer)))
  
  kmerMapping(k = k, m = m, l = k + m, mapping = kmer.match.mat, colMapping = kmer.col.mat)
  
}  



make_wildcard_kmer <-function(k, m, adjacent_n = TRUE){
  
  all_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))
  
  wc_kmer = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T","N"),width = k))
  
  #reduce the space
  wc_kmer = wc_kmer[which(letterFrequency(wc_kmer,"N")<=m)]
  wc_kmer = wc_kmer[which(letterFrequency(subseq(wc_kmer,start=1,end=1),"N")==0)]
  
  if (adjacent_n){
    for (i in 2:m){
      mult = which(letterFrequency(wc_kmer,"N") < i)
      consec = which(vcountPattern(paste(rep("N",i),collapse=""),wc_kmer)>=1)
      wc_kmer = wc_kmer[c(mult,consec)]
    } 
  }
  
  wc_kmer <- remove_rc(wc_kmer)
  
  pd = PDict(all_kmer)
  mapping  = vwhichPDict(pd,wc_kmer,fixed = "pattern")
  mapping.rc  = vwhichPDict(pd,Biostrings::reverseComplement(wc_kmer),fixed = "pattern")

  #rows are kmer, col are wc-kmer
  kmer.match.mat = sparseMatrix(i = c(unlist(mapping),unlist(mapping.rc)),
                                j = c(unlist(sapply(1:length(mapping), function(x) rep(x,length(mapping[[x]])))),
                                      unlist(sapply(1:length(mapping.rc), function(x) rep(x,length(mapping.rc[[x]]))))),
                                x = 1, dims = c(length(all_kmer),length(wc_kmer)), use.last.ij = TRUE,
                                dimnames = list(as.character(all_kmer), as.character(wc_kmer))
                                )
  
  kmer.col.mat = rep(colnames(kmer.match.mat),2)
  names(kmer.col.mat) = c(colnames(kmer.match.mat), as.character(Biostrings::reverseComplement(wc_kmer)))
  
  kmerToWildcard(k = k, m = m, mapping = kmer.match.mat, colMapping = kmer.col.mat)
  
}  
  

get_gapped_kmer_counts <- function(seqs, kmer_mapping, verbose = F){
  seqs = as.character(seqs)
  seq_kmer_mat = sparseMatrix(i={},j={}, dims = c(length(seqs), ncol(kmer_mapping@mapping)), dimnames=list(NULL,colnames(kmer_mapping@mapping)))
  add_kmer_mat = sparseMatrix(i={},j={}, dims = c(length(seqs), ncol(kmer_mapping@mapping)), dimnames=list(NULL,colnames(kmer_mapping@mapping)))
  blocks = 20
  if (verbose){print(paste("Scanning across ", nchar(seqs[1])," bases per peak for kmers",collapse=""))}
  for (i in 1:(nchar(seqs[1]) - kmer_mapping@l + 1)){
    if (i %% blocks == 0){
      if (verbose){print(paste("Scanned ",i," bases"))}
      seq_kmer_mat = seq_kmer_mat + add_kmer_mat
      add_kmer_mat = sparseMatrix(i={},j={}, dims = c(length(seqs), ncol(kmer_mapping@mapping)), dimnames=list(NULL,colnames(kmer_mapping@mapping)))
    }
    add_kmer_mat = add_kmer_mat + kmer_mapping[substr(seqs,i,i+(kmer_mapping@l-1))]
  }
  seq_kmer_mat = seq_kmer_mat + add_kmer_mat
  seq_kmer_mat
}


pd = Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))
kmer_matchs = vwhichPDict(pd,seqs)




get_enriched_kmer<- function(seqs, genome, seq_kmer_mat){
  
  
}



  wc_kmer_ac = as.character(wc_kmer) 
  #kmer.col.mat = rep(colnames(kmer.match.mat),2)
  #names(kmer.col.mat) = c(colnames(kmer.match.mat), as.character(Biostrings::reverseComplement(wc_kmer)))
  seq_kmer_mat = sparseMatrix(i={},j={}, dims = c(length(seqs), length(wc_kmer)),dimnames=list(NULL,wc_kmer_ac))
  seqs_ac = as.character(seqs)
  for (i in 1:(width(seqs[1])-k + 1)){
    seq_kmer_mat = seq_kmer_mat + kmer.match.mat[substr(seqs_ac,i,i+(k-1)),wc_kmer_ac]
  }
  seq_kmer_mat
}
  

###Find kmer that are over-represented

#Find expectation
n_kmers = (width(seqs[1]) - k + 1) * length(seqs)
nuc_freqs = colSums(Biostrings::letterFrequency(seqs,c("A","C","G","T")))
nuc_freqs = nuc_freqs/sum(nuc_freqs)
#nuc_freqs = c(nuc_freqs, c("N"=1)) #* (n_kmers ** (1/k))
kmer_nuc_freqs = Biostrings::letterFrequency(all_kmer,c("A","C","G","T"))  
kmer_exp =  apply(nuc_freqs ** t(kmer_nuc_freqs),2,prod) * n_kmers
wc_kmer_exp = as.numeric(kmer_exp %*% kmer.match.mat) 
wc_kmer_counts = colSums(seq_kmer_mat)
wc_kmer_enrich = wc_kmer_counts/wc_kmer_exp


sparse_jaccard <- function(x){
  
  x[which(x>1)]=1
  i = t(x) %*% x
  tmp = colSums(x)
  u = matrix(tmp, nrow = nrow(i), ncol=ncol(i), byrow=T) + matrix(tmp, nrow = nrow(i), ncol=ncol(i), byrow=F) - i
  jac = i/u
  return(jac)
  #return(as.dist(jac))
  
}

tmp = summary(seq_kmer_mat)
seq_kmer_mat_bin = sparseMatrix(i = tmp[,1], j = tmp[,2], x = 1, dims = dim(seq_kmer_mat), dimnames = list(NULL,colnames(seq_kmer_mat)))

sparse_jaccard_help <- function(y,x){
  
  i = t(y) %*% x
  u = outer(colSums(y), colSums(x), FUN = "+") - i
  jac = i/u
  return(jac)
  #return(as.dist(jac))
  
}



wc_kmer_dist = vegan::vegdist(seq_kmer_mat[,which(wc_kmer_enrich>=2)], method="jaccard",binary=T)


reduce_kmer <- function(){
  max_jaccard = 0.9
  
  
  
  
  
}





add_wildcard_assoc <- function()
  
  


wildcardInit <- function(seqs, k, m){
  
}

wildCardSub <- function(seqs, )