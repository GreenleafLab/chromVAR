
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
    stop("k must be less than 8")
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
    indices_rc <- Biostrings::vwhichPDict(kmers,Biostrings::reverseComplement(seqs), fixed = FALSE)
  } else {
    pd <- Biostrings::PDict(kmers)
    indices  <- Biostrings::vwhichPDict(pd,seqs)
    indices_rc <- Biostrings::vwhichPDict(pd,Biostrings::reverseComplement(seqs))
  } 
  indices <- merge_lists(indices, indices_rc, by = "order")
  indices <- lapply(indices, unique)
  
  tmp <- data.frame(peak_ix = unlist(lapply(seq_along(indices),
                                            function(x)rep(x, length(indices[[x]])))),
                    kmer_ix = factor(unlist(indices),
                                     levels = seq_along(kmers),
                                     ordered=T))
  out <- split(tmp$peak_ix,tmp$kmer_ix)
  names(out) <- as.character(kmers)
  return(out)
}

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
  
  out <- lapply(seq_along(mapping),
                function(x) unique(unlist(lmer_indices[mapping[[x]]],
                                          use.names=FALSE)))
  
  names(out) <- as.character(gapped_kmer)
  
  return(out)
}

# kmer utilility functions -----------------------------------------------------

remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs),
                as.character(Biostrings::reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  Biostrings::DNAStringSet(unique(temp[,1]))
}

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

get_kmer_alignment <- function(kmers, reference, mismatch = -1, minimum = FALSE, 
                               both_strands = TRUE){
  out = data.frame()
  for (kmer in kmers){
    pwm1 = seq_to_pwm(kmer)
    if (is.character(reference)){
      pwm2 = seq_to_pwm(reference, mismatch = mismatch)      
    } else {
      pwm2 = reference
    }
    alignment <- align_pwms(pwm1, pwm2, minimum = minimum, 
                            both_strands = both_strands)
    if (nrow(alignment) > 0){
      out = rbind(out,
                  data.frame(kmer = ifelse(alignment$strand ==-1,
                                           as.character(Biostrings::reverseComplement(Biostrings::DNAString(kmer))),
                                           kmer),
                             shift = alignment$pos,
                             score = alignment$score,
                             input = kmer, stringsAsFactors=FALSE))
    }
  }
  return(unique(out))
}

get_similar_kmers <- function(kmer, kmers, cutoff = 4){
  a = get_kmer_alignment(kmers,kmer, mismatch = 0, minimum = cutoff)
  return(unique(a$input))
}

get_overlapping_kmers <- function(kmer, kmers, cutoff = 2){
  a = get_kmer_alignment(kmers,kmer, mismatch = -Inf, minimum = cutoff)
  return(unique(a$input))
}

# kmer grouping functions ------------------------------------------------------

make_kmer_group <- function(seed, candidates, sets, counts_mat, minScore = nchar(seed)-2){
  preliminary_alignment = get_kmer_alignment(kmers = c(seed ,candidates), reference = seed, mismatch = -100)
  candidates = unique(preliminary_alignment$input[which(preliminary_alignment$score >= minScore)])
  print(length(candidates))
  var_boost = get_variability_boost(sets[[seed]], sets[candidates], counts_mat, nbg = 25, out = "z-score")
  print(var_boost)
  boost = which(pnorm(var_boost,lower.tail = FALSE) < 0.05)
  print(length(boost))
  if (length(boost)>0){
    a = get_kmer_alignment(kmers = c(seed ,candidates[boost]), reference = seed, mismatch = -100)
  }
  else{
    a = get_kmer_alignment(kmers = seed, reference = seed, mismatch = -100)
  }
  highscore = which(a$score >= minScore)
  consensus = Biostrings::consensusString(Biostrings::DNAStringSet(a$kmer[highscore]),
                                          shift = a$shift[highscore] - min(a$shift[highscore]), 
                                          threshold = 0.1)
  rc_kmers = unique(get_reverse_complement(a$kmer))
  rc_kmers = rc_kmers[rc_kmers %ni% a$kmer]
  a2 = get_kmer_alignment(kmers = c(a$kmer, rc_kmers), reference = consensus, mismatch = -100, both_strands = FALSE)
  a2 = a2[which(a2$score >= minScore),]
  a2$shift = a2$shift - a2$shift[which(a2$kmer == seed)]
  seed_block = a2[which(a2$kmer == seed),]
  forward_block = a2[intersect(which(a2$kmer != seed),which(a2$kmer %in% a$kmer)),]
  reverse_block = a2[which(a2$kmer %in% rc_kmers),]
  out = rbind(seed_block,
              forward_block[order(forward_block$score),],
              reverse_block[order(reverse_block$score),])
  return(out)
}


group_kmers <- function(results, sets, counts_mat, max_groups = 10){
  results <- subset_by_variability(results)
  indep_results <- results
  i = 1
  out = list()
  while (i <= max_groups){
    print(i)
    top = names(subset_by_variability(indep_results, by = "top", cutoff =1))
    candidates = get_correlated_results(top, results)   
    new_grp = make_kmer_group(top, candidates, sets, counts_mat)   
    out = c(out, list(new_grp))
    candidates = candidates[which(candidates %in% names(indep_results))]
    indep_var = get_independent_variability(to.remove = sets[[top]],
                                            sets = sets[candidates],
                                            counts_mat = counts_mat)
    indep = which(get_pvalues(indep_var) < 0.01)
    keep = unique(c(names(indep_var[indep]),names(indep_results)[which(names(indep_results) %ni% c(top,candidates))]))
    if (length(keep) > 0 ){
      indep_results = indep_results[keep]
    } else{
      break
    }
    i = i + 1
  }
  return(out)
}


get_similar_motifs0 <- function(a, motif_list, top = 3){
  consensus = Biostrings::consensusMatrix(Biostrings::DNAStringSet(a$kmer),
                                          shift = a$shift - min(a$shift))[1:4,]
  consensus_shift = align_pwms(consensus, seq_to_pwm(a$kmer[1]), both_strands = FALSE)$pos
  helpfun <- function(x){
    x = as.matrix(x)
    x = x + 0.1
    x = x / matrix(colSums(x), byrow = TRUE, nrow = nrow(x), ncol= ncol(x))
    x = log(x/0.25)
    align_pwms(x, consensus)[1,]
  }
  motif_matches <- lapply(motif_list, helpfun)
  scores <- sapply(motif_matches, function(x) x$score)
  return(motif_list[order(-scores)[1:top]])
}

get_similar_motifs <- function(kmer_groups, motif_list, top = 3){
  if (is.data.frame(kmer_groups)){
    return(get_similar_motifs0(kmer_groups, motif_list, top))
  } else if (is.list(kmer_groups)){
    return(lapply(kmer_groups, get_similar_motifs0, motif_list, top))
  } else{
    stop("input should be list of kmer groups or single kmer group")
  }
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
