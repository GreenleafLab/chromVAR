
get_nuc_freqs <- function(seqs){
  #get nucleotide frequencies
  nucFreqs <- colSums(Biostrings::letterFrequency(seqs, c("A","C","G","T")))
  nucFreqs <- nucFreqs/sum(nucFreqs)
  return(nucFreqs)
}


get_motif_matches <- function(motif, seqs, nucfreqs = get_nuc_freqs(seqs), min.score = "90%"){
  
  #concatenate sequences to speed up motif matching
  seq_widths <- width(seqs)
  n_seq <- length(seqs)
  ranges <- IRanges(start = cumsum(c(1,seq_widths)[1:n_seq]), width = seq_widths)
  tmpseq <- paste(seqs, collapse="")
  
  #find motifs on forward and reverse strands
  forward <- matchPWM(log(motif/nucfreqs),tmpseq, with.score=T, min.score = min.score)
  reverse <- matchPWM(log(reverseComplement(motif)/nucfreqs),tmpseq, with.score = T, min.score = min.score)

  #map motif matches back to peaks
  forward_matches <- findOverlaps(forward,ranges,type="within",select="all")
  reverse_matches <- findOverlaps(reverse,ranges,type="within",select="all")
  
  if (length(forward_matches)>0 | length(reverse_matches)>0){
    #get max motif score per range
    score_mat <- rbind(cbind(subjectHits(forward_matches),mcols(forward[queryHits(forward_matches)])$score),
                       cbind(subjectHits(reverse_matches),mcols(reverse[queryHits(reverse_matches)])$score))
    colnames(score_mat) <- c("id","score")
    max_scores <- plyr::ddply(as.data.frame(score_mat),"id",plyr::summarize,max_score = max(score))
    return(max_scores)
  } else{
    return(data.frame("id"=NULL,"max_score"=NULL))
  }}
  #scores <- rep(0,n_seq)
  #scores[max_scores$id] <- max_scores$max_score
  
  #return(scores)}



get_motifs <- function(species = "Hsapiens", dataSource = "JASPAR"){
  
  if (length(species)!=1){
    stop("Must input species")}
  
  species_indices = grep(species, values (MotifDb::MotifDb)$organism)
  
  if (length(species_indices) < 1){
    stop("No motifs found species. Perhaps incorrect species name used?")}
  
  if (length(dataSource)==1){
    data_indices <- grep(dataSource, values (MotifDb::MotifDb)$dataSource)
  } else if (length(dataSrouce>1)){
    data_indices <- do.call(union,lapply(dataSource, function(x) grep(x, values (MotifDb::MotifDb)$dataSource)))
  } else{
    data_indices <- species_indices }
  
  if (length(data_indices) < 1){
    stop("No motifs found for given dataSource. Perhaps incorrect spelling/name?")
  }
  
  motifs <- MotifDb::MotifDb[intersect(species_indices, data_indices)]
  
  return(motifs)
}


get_motif_table <- function(bed, motifs, seqs, BPPARAM = bpparam()){

  #get motifs from MotifDb
  #motifs <- get_motifs(species, dataSource) 

  #get sequences
  #seqs <- getSeq(genome,bed)
  
  #get nucleotide frequencies
  nucFreqs <-  get_nuc_freqs(seqs)
  
  .get_motif_matches <- function(index){
    res = get_motif_matches(motifs[[index]], seqs, nucFreqs)
    if (nrow(res)==0){
      return(cbind(data.frame("motif"=NULL),res))
    } else{
      return(cbind(data.frame("motif" = names(motifs[index])),res))
    }
  }
  
  #find motif matches
  motif_df <- do.call(rbind,bplapply(1:length(motifs), .get_motif_matches, BPPARAM = BPPARAM))
  motif_mat <- sparseMatrix(i = motif_df$id, j = as.integer(motif_df$motif), x= motif_df$max_score, 
                            dims =c(length(seqs),length(levels(motif_df$motif))), dimnames = c(NULL,levels(motif_db$motif)))
  
  return(motif_mat)
}
  

filter_motif_table <- function(motif_table, min.matches = 250, max.matches = length(bed)*0.25){
  
  n_matches = apply(motif_table, 2, function(x) sum(x>0))
  
  ##for motifs with too many matches, keep only top matches
  high = which(n_matches > max.matches)
  for (i in high){
    nz = which(motif_table[,i]>0)
    cutoff = quantile( motif_table[nz,i], 1 - (max.matches/n_matches[i]))
    below_cutoff = nz[which(motif_table[nz,i]<cutoff)]
    if (length(below_cutoff) > 0){
      motif_table[below_cutoff,i] = 0
    }
  }
  #remove motifs with too few matches
  motif_table = motif_table[,which(n_matches > min.matches)]
  return(motif_table)
}


  
  
  
  