
# kmer comparison and grouping functions ---------------------------------------

get_single_mm_var <- function(kmer, kmer_variability){
  
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  kmer_nucs = strsplit(kmer,"")[[1]]
  kvar <- kmer_variability[kmer,"variability"]
  
  out <- matrix(nrow = 4,ncol = k)
  out2 <- matrix(nrow = 4,ncol = k)
  rownames(out) <- nucs
  rownames(out2) <- nucs
  for (i in 1:k){
    for (j in nucs[nucs != kmer_nucs[i]]){
      kmod = kmer_nucs
      kmod[i] = j
      kmod = paste(kmod, sep="",collapse="")
      if (kmod %ni% row.names(kmer_variability)){
        kmod = get_reverse_complement(kmod)
      }
      out[j,i] = kmer_variability[kmod,"variability"]
      out2[j,i] = kmod
    }
  }
  return(list(var = out, kmers = out2))
}

get_single_mm_cor <- function(kmer, kmer_deviations){
  
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  kmer_nucs = strsplit(kmer,"")[[1]]
  
  out <- matrix(nrow = 4,ncol = k)
  out2 <- matrix(nrow = 4,ncol = k)
  rownames(out) <- nucs
  rownames(out2) <- nucs
  for (i in 1:k){
    for (j in nucs[nucs != kmer_nucs[i]]){
      kmod = kmer_nucs
      kmod[i] = j
      kmod = paste(kmod, sep="",collapse="")
      if (kmod %ni% row.names(kmer_deviations$z)){
        kmod = get_reverse_complement(kmod)
      }
      out[j,i] = cor(as.vector(kmer_deviations$z[kmod,]),as.vector(kmer_deviations$z[kmer,]), use = "pairwise.complete.obs")
      out2[j,i] = kmod
    }
  }
  return(list(cor = out, kmers = out2))
}


plot_mm_var_effect <- function(mm_var, k_var){
  
  mm_var <- (mm_var**2 - 1)/(k_var**2 - 1)
  mm_var[mm_var<0]=0
  
  mm_df1 <- data.frame(val = as.vector(mm_var), pos = rep(1:ncol(mm_var),each = 4),
                       Nucleotide = rep(c("A","C","G","T"),ncol(mm_var)))
  mm_df1 <- mm_df1[!is.na(mm_df1$val),]
  
  ggplot(mm_df1, aes(x = pos, y = val, col = Nucleotide)) + 
    geom_point(size = 3, position = position_jitter(height = 0, width = 0.1))  + 
    geom_hline(yintercept = 1, lty = 2, col = "gray") + ylab("% excess variance") +
    xlab("Position") + scale_x_continuous(breaks = c(1:ncol(mm_var)))+
    chromVAR_theme()
  
}


plot_mm_cor_effect <- function(mm_cor){
  
  mm_df2 <- data.frame(val = as.vector(mm_cor), pos = rep(1:ncol(mm_cor),each = 4),
                       Nucleotide = rep(c("A","C","G","T"),ncol(mm_cor)))
  mm_df2 <- mm_df2[!is.na(mm_df2$val),]
  
  ggplot(mm_df2, aes(x = pos, y = val, col = Nucleotide)) + 
    geom_point(size = 3, position = position_jitter(height = 0, width = 0.1))  + 
    geom_hline(yintercept = c(0,1,-1), lty = 2, col = "gray") + ylab("Correlation in deviations") +
    xlab("Position") + scale_x_continuous(breaks = c(1:ncol(mm_cor)))+
    chromVAR_theme() 
  
}

get_overlap_kmers <- function(kmer, max_extend, dir = "both"){
  out = c()
  k <- nchar(kmer)
  nucs <- c("A","C","G","T") 
  if (max_extend == 1){
    if (dir == "both"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)    
      return(c(out1,out2))
    } else if (dir == "left"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      return(out1)
    } else if (dir == "right"){
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)    
      return(out2)
    }   
  } else{
    if (dir == "both"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)    
      return(c(out1,out2, sapply(out1, get_overlap_kmers, max_extend = max_extend - 1, dir = "left"),
               sapply(out2, get_overlap_kmers, max_extend = max_extend - 1, dir = "right")))
    } else if (dir == "left"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      return(c(out1,sapply(out1, get_overlap_kmers, max_extend = max_extend - 1, dir = "left")))
    } else if (dir == "right"){
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)    
      return(c(out2, 
               sapply(out2, get_overlap_kmers, max_extend = max_extend - 1, dir = "right")))
    }
  }
}

get_kmer_overlap_effect <- function(kmer, kmer_variability, kmer_deviations, counts_mat, background_peaks, peak_indices, expectation = NULL, norm = TRUE, max_extend = 2){
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  kmer_nucs = strsplit(kmer,"")[[1]]
  out = list()
  ncands = sum(sapply(1:max_extend, function(x) 4^x*2))
  #get candidate kmers that are extension of up to max_extend bases
  out$kmers = get_overlap_kmers(kmer, max_extend)
  kmers2 = ifelse(out$kmers %in% colnames(peak_indices), out$kmers, get_reverse_complement((out$kmers)))
  out$shift = do.call(c, sapply(1:max_extend, function(x) c(rep(-x,4^x),rep(x,4^x))))
  #get variability boost for those kmers
  out$var_boost = chromVAR:::get_variability_boost(1, 
                                                   counts_mat, 
                                                   background_peaks, 
                                                   peak_indices[,c(kmer, kmers2)], 
                                                   expectation,
                                                   norm)[2:(length(kmers2)+1)]
  #get variability of those kmers
  out$var = kmer_variability[kmers2,"variability"]
  #get deviation correlation for those kmers
  out$cor = cor(t(kmer_deviations$z[c(kmer,kmers2),]), use ="pairwise.complete.obs")[1,2:(length(kmers2)+1)]
  return(out)
}


#'@export
make_kmer_group <- function(kmer, kmer_variability, kmer_deviations, counts_mat, background_peaks, peak_indices, expectation = NULL, norm = TRUE, max_extend = 2,
                            cor_threshold = 0.67, var_threshold = 0.67, boost_threshold = 3){
  
  mm_var = get_single_mm_var(kmer, kmer_variability)
  mm_cor = get_single_mm_cor(kmer, kmer_deviations)
  
  out = data.frame(kmer = kmer, type = "seed", var = kmer_variability[kmer,"variability"], cor = 1, boost = NA, shift = 0, stringsAsFactors = FALSE)
  mm_keep = intersect(which((mm_var$var**2 -1)/(out$var**2 -1) > var_threshold), which(mm_cor$cor > cor_threshold))
  
  if (length(mm_keep) >= 1){   
    out = rbind(out, data.frame( kmer = mm_var$kmers[mm_keep], type = "mismatch", var = mm_var$cor[mm_keep], cor = mm_cor$cor[mm_keep], boost = NA, shift =0, stringsAsFactors = FALSE))
  }
  
  overlap_effect = get_kmer_overlap_effect(kmer, kmer_variability, kmer_deviations, counts_mat, background_peaks, peak_indices, expectation, norm , max_extend = 2)
  overlap_keep = intersect(which(overlap_effect$var_boost > boost_threshold),which(overlap_effect$cor > cor_threshold))
  
  if (length(overlap_keep) >= 1){
    out = rbind(out, 
                data.frame(kmer = overlap_effect$kmers[overlap_keep], 
                           type = "overlap", 
                           var = overlap_effect$var[overlap_keep], 
                           cor = overlap_effect$cor[overlap_keep], 
                           boost = overlap_effect$var_boost[overlap_keep],
                           shift = overlap_effect$shift[overlap_keep], stringsAsFactors = FALSE))    
  }
  
  
  for (j in mm_var$kmers[mm_keep]){
    overlap_effect = get_kmer_overlap_effect(j, kmer_variability, kmer_deviations, counts_mat, background_peaks, peak_indices, expectation, norm , max_extend = 2)
    overlap_keep = intersect(which(overlap_effect$var_boost > boost_threshold),which(overlap_effect$cor > cor_threshold), which(overlap_effect %ni% out$kmer))
    if (length(overlap_keep) >= 1){
      out = rbind(out, 
                  data.frame(kmer = overlap_effect$kmers[overlap_keep], 
                             type = "mismatch_overlap", 
                             var = overlap_effect$var[overlap_keep], 
                             cor = overlap_effect$cor[overlap_keep], 
                             boost = overlap_effect$var_boost[overlap_keep],
                             shift = overlap_effect$shift[overlap_keep], stringsAsFactors = FALSE))
    }
  }
  row.names(out) = NULL
  
  return(out)
}



get_similar_motifs0 <- function(kmer_group, motif_list, top = 3){
  
  consensus = Biostrings::consensusMatrix(Biostrings::DNAStringSet(kmer_group$kmer),
                                          shift = kmer_group$shift - min(kmer_group$shift))[1:4,]
  consensus_shift = align_pwms(consensus, seq_to_pwm(kmer_group$kmer[1]), both_strands = FALSE)$pos
  helpfun <- function(x){
    x = as.matrix(x)
    #x = x + 0.1
    #x = x / matrix(colSums(x), byrow = TRUE, nrow = nrow(x), ncol= ncol(x))
    #x = log(x/0.25)
    align_pwms(x, consensus)[1,]
  }
  motif_matches <- lapply(motif_list, helpfun)
  scores <- sapply(motif_matches, function(x) x$score)
  return(motif_list[order(-scores)[1:top]])
}

#'@export
get_similar_motifs <- function(kmer_groups, motif_list, top = 3){
  if (is.data.frame(kmer_groups)){
    return(get_similar_motifs0(kmer_groups, motif_list, top))
  } else if (is.list(kmer_groups)){
    return(lapply(kmer_groups, get_similar_motifs0, motif_list, top))
  } else{
    stop("input should be list of kmer groups or single kmer group")
  }
}



