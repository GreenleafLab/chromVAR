# Utility functions ------------------------------------------------------------

get_reverse_complement <- function(x){
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

ambiguity_mapping <- lapply(Biostrings::IUPAC_CODE_MAP, function(x){
  letters = strsplit(x,"")[[1]]
  out = rep(0, length(letters))
  out[which(letters == "A")] = 1
  out[which(letters == "C")] = 2
  out[which(letters == "G")] = 3
  out[which(letters == "T")] = 4
  return(out)  
})


seq_to_pwm <- function(in_seq, mismatch = 0){
  mat <- matrix(mismatch, ncol = nchar(in_seq), nrow = 4)
  rownames(mat) <- c("A","C","G","T")
  in_seq <- strsplit(in_seq,"")[[1]]
  for (x in seq_along(in_seq)){
    mat[ambiguity_mapping[[in_seq[x]]],x] = 1
  }
  return(mat) 
}



# kmer relatedness and alignment functions -------------------------------------


align_pwms <- function(pwm1,pwm2, minimum = FALSE, both_strands = TRUE){
  w1 = ncol(pwm1)
  w2 = ncol(pwm2)
  cc = sapply(1:(w1+w2-1), function(x) sum(pwm1[,max(w1-x+1,1):min(w1,w1+w2-x)]*pwm2[,max(1,x-w1+1):min(x,w2)]))
  if (both_strands){
    pwm1_rev = pwm1[c(4,3,2,1),w1:1]
    cc_rev = sapply(1:(w1+w2-1), function(x) sum(pwm1_rev[,max(w1-x+1,1):min(w1,w1+w2-x)]*pwm2[,max(1,x-w1+1):min(x,w2)]))
  }
  out = data.frame()
  if (minimum){
    if (max(cc) >= minimum){
      wm = which(cc >= minimum)
      pos = wm - w1
      out = rbind(out, data.frame(score = cc[wm], pos = pos, strand = 1))
    }
    if (both_strands && max(cc_rev >= minimum)){
      wm = which(cc_rev >= minimum)
      pos = wm - w1
      out = rbind(out, data.frame(score = cc_rev[wm], pos = pos, strand = -1))
    }
  } else{
    if (!both_strands || (max(cc) >= max(cc_rev))){
      wm = which(cc == max(cc))
      pos = wm - w1
      out = rbind(out, data.frame(score = max(cc), pos = pos, strand = 1))
    }
    if (both_strands && (max(cc_rev) >= max(cc))){
      wm = which(cc_rev == max(cc_rev))
      pos = wm - w1
      out = rbind(out,data.frame(score = max(cc_rev), pos = pos, strand = -1))
    }
  }  
  return(out)
}


# kmer position functions ------------------------------------------------------

get_kmer_positions <- function(kmer, peaks, seqs){
  matches <- Biostrings::vmatchPattern(Biostrings::DNAString(kmer), seqs, fixed=FALSE)
  rc_matches <- Biostrings::vmatchPattern(Biostrings::reverseComplement(Biostrings::DNAString(kmer)), seqs, fixed=FALSE)
  
  tmp1 <- elementLengths(matches)
  tmp2 <-  unlist(sapply(1:length(peaks), function(x) rep(x,tmp1[x])), use.names=F)
  f_pos <- resize(shift(peaks[tmp2],shift = start(unlist(matches)))+1,width = 1)
  BiocGenerics::strand(f_pos) <- "+"
  
  tmp1 <- elementLengths(rc_matches)
  tmp2 <-  unlist(sapply(1:length(peaks), function(x) rep(x,tmp1[x])), use.names=F)
  r_pos <- resize(shift(peaks[tmp2], shift = start(unlist(rc_matches)) -1),width = 1)
  BiocGenerics::strand(r_pos) <- "-"  
  
  return(BiocGenerics::sort(c(f_pos, r_pos)))
}

get_kmer_to_kmer_dist <- function(kmer1, kmer2, peaks, seqs, all_seqs, max_dist = 25){
  
  ranges1 <- get_kmer_positions(kmer1, peaks, seqs)
  ranges2 <- get_kmer_positions(kmer2, peaks, seqs)
  
  ranges1_mod = ranges1
  BiocGenerics::strand(ranges1_mod)="*"
  close = as.data.frame(findOverlaps(ranges1_mod, ranges2, maxgap = max_dist, select="all", type="any"))
  ss = which(getstrand(ranges1[close$queryHits]) == getstrand(ranges2[close$subjectHits]) )
  ds = which(getstrand(ranges1[close$queryHits]) != getstrand(ranges2[close$subjectHits]) )
  
  dists_f = ifelse(getstrand(ranges1[close[ss,1]]) == "-", 
                   BiocGenerics::start(ranges1[close[ss,1]]) - BiocGenerics::start(ranges2[close[ss,2]]),
                   BiocGenerics::start(ranges2[close[ss,2]]) - BiocGenerics::start(ranges1[close[ss,1]]))
  
  dists_r = ifelse(getstrand(ranges1[close[ds,1]]) == "-", 
                   BiocGenerics::start(ranges1[close[ds,1]]) - BiocGenerics::start(ranges2[close[ds,2]]),
                   BiocGenerics::start(ranges2[close[ds,2]]) - BiocGenerics::start(ranges1[close[ds,1]]))
  
  
  
  expected = rep(sum(vcountPattern(kmer2, all_seqs))/ sum(width(all_seqs)-nchar(kmer2)+1) * length(ranges1), 2* max_dist +1)
  expected_rc = expected
  
  out = list(forward = tabulate2(dists_f, max_val = max_dist, min_val = -max_dist), 
             reverse_complement = tabulate2(dists_r, max_val = max_dist, min_val = -max_dist))
  
  for (i in -(nchar(kmer1)-1):(nchar(kmer1)-1)){
    if (out$forward[as.character(i)] > 0 ){
      if (i < 0){
        tmp_mer = substr(kmer2, 1, nchar(kmer2) + i)
        expected[max_dist + i + 1] = sum(vcountPattern(tmp_mer, all_seqs))/ sum(width(all_seqs)-nchar(tmp_mer)+1) * length(ranges1)
      } else if (i > 0){
        tmp_mer = substr(kmer2, nchar(kmer2) -i, nchar(kmer2))
        expected[max_dist + i + 1] = sum(vcountPattern(tmp_mer, all_seqs))/ sum(width(all_seqs)-nchar(tmp_mer)+1) * length(ranges1)       
      } else{
        expected[max_dist + 1] = length(ranges1) 
      }
    }
    if (out$reverse[as.character(i)] > 0 ){
      rev_kmer2 = get_reverse_complement(kmer2)
      if (i < 0){
        tmp_mer = substr(rev_kmer2, 1, nchar(rev_kmer2) + i)
        expected_rc[max_dist + i + 1] = sum(vcountPattern(tmp_mer, all_seqs))/ sum(width(all_seqs)-nchar(tmp_mer)+1) * length(ranges1)
      } else if (i > 0){
        tmp_mer = substr(rev_kmer2, nchar(rev_kmer2) - i, nchar(rev_kmer2))
        expected_rc[max_dist + i + 1] = sum(vcountPattern(tmp_mer, all_seqs))/ sum(width(all_seqs)-nchar(tmp_mer)+1) * length(ranges1)       
      } else{
        expected_rc[max_dist + 1] = length(ranges1) 
      }
    }  
  }
  
  out$forward = -1 * ppois(out$forward,expected, lower.tail = FALSE,log.p=TRUE)
  out$reverse_complement = -1 * ppois(out$reverse_complement,expected_rc, lower.tail = FALSE,log.p=TRUE)
  
  return(out)
}


get_sequence_flanking_kmer <- function(kmer, 
                                       peaks, 
                                       genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                                       flank = 5){
  
  seqs = Biostrings::getSeq(genome, peaks)
  
  matches <- Biostrings::vmatchPattern(Biostrings::DNAString(kmer), seqs, fixed=FALSE)
  rc_matches <- Biostrings::vmatchPattern(Biostrings::reverseComplement(Biostrings::DNAString(kmer)), seqs, fixed=FALSE)
  
  tmp1 <- elementLengths(matches)
  tmp2 <-  unlist(sapply(1:length(peaks), function(x) rep(x,tmp1[x])), use.names=F)
  f_seqs <- Biostrings::getSeq(genome,resize(shift(peaks[tmp2], 
                                                   shift = start(unlist(matches)) - flank - 1), 
                                             width = flank*2 + 1 + nchar(kmer)))
  tmp1 <- elementLengths(rc_matches)
  tmp2 <-  unlist(sapply(1:length(peaks), function(x) rep(x,tmp1[x])), use.names=F)
  r_seqs <- Biostrings::reverseComplement(Biostrings::getSeq(genome,resize(shift(peaks[tmp2], 
                                                                                 shift = start(unlist(rc_matches)) - flank - 2), 
                                                                           width = flank*2 + 1 + nchar(kmer))))
  return(c(f_seqs,r_seqs))
}


