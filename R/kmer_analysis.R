
get_reverse_complement <- function(x){
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

get_single_mm_effect <- function(kmer, kmer_variability){
  
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  kmer_nucs = strsplit(kmer,"")[[1]]
  kvar <- kmer_variability[kmer,"variability"]
  
  out <- matrix(nrow = 4,ncol = k)
  rownames(out) <- nucs
  for (i in 1:k){
    for (j in nucs[nucs != kmer_nucs[i]]){
      kmod = kmer_nucs
      kmod[i] = j
      kmod = paste(kmod, sep="",collapse="")
      if (kmod %ni% row.names(kmer_variability)){
        kmod = get_reverse_complement(kmod)
      }
      out[j,i] = kmer_variability[kmod,"variability"]
    }
    #out[kmer_nucs[i],i] <- kvar
  }
  out <- (out**2 - 1)/(kvar**2 - 1)
  out[out<0]=0
  return(out)
}

plot_mm_effect <- function(mm_mat){
  
  mm_df <- data.frame(val = as.vector(mm_mat), pos = rep(1:ncol(mm_mat),each = 4),
                      Nucleotide = rep(c("A","C","G","T"),ncol(mm_mat)))
  mm_df <- mm_df[!is.na(mm_df$val),]
  
  ggplot(mm_df, aes(x = pos, y = val, col = Nucleotide)) + 
    geom_point(size = 3, position = position_jitter(height = 0, width = 0.1))  + 
    geom_hline(yintercept = 1, lty = 2, col = "gray") + ylab("% excess variance") +
    xlab("Position") + scale_x_continuous(breaks = c(1:ncol(mm_mat)))+
    chromVAR_theme()
  
}



get_kmer_group <- function(seed, counts_mat, bg, peak_indices, max_extend = 2){
  #get candidate kmers that are extension of up to max_extend bases
  
  #get variability boost for those kmers
  
  #get maximal alignment for those that boost
  
  return(out)
}








