
get_nuc_freqs <- function(seqs){
  #get nucleotide frequencies
  nucFreqs <- colSums(Biostrings::letterFrequency(seqs, c("A","C","G","T")))
  nucFreqs <- nucFreqs/sum(nucFreqs)
  return(nucFreqs)
}


get_motif_matches <- function(motif, seqs, nucfreqs = get_nuc_freqs(seqs), min.score = "90%", max.matches = length(seqs)/5){

  #concatenate sequences to speed up motif matching
  seq_widths <- Biostrings::width(seqs)
  n_seq <- length(seqs)
  ranges <- IRanges::IRanges(start = cumsum(c(1,seq_widths)[1:n_seq]), width = seq_widths)
  tmpseq <- paste(seqs, collapse="")

  #find motifs on forward and reverse strands
  forward <- Biostrings::matchPWM(log(motif/nucfreqs),tmpseq, with.score=T, min.score = min.score)
  reverse <- Biostrings::matchPWM(log(Biostrings::reverseComplement(motif)/nucfreqs),tmpseq, with.score = T, min.score = min.score)

  #map motif matches back to peaks
  forward_matches <- IRanges::findOverlaps(forward,ranges,type="within",select="all")
  reverse_matches <- IRanges::findOverlaps(reverse,ranges,type="within",select="all")

  if (length(forward_matches)>0 | length(reverse_matches)>0){
    #get max motif score per range
    score_mat <- rbind(cbind(IRanges::subjectHits(forward_matches),S4Vectors::mcols(forward[IRanges::queryHits(forward_matches)])$score),
                       cbind(IRanges::subjectHits(reverse_matches),S4Vectors::mcols(reverse[IRanges::queryHits(reverse_matches)])$score))
    colnames(score_mat) <- c("id","score")
    max_scores <- plyr::ddply(as.data.frame(score_mat),"id",plyr::summarize,max_score = max(score))
    if (nrow(max_scores) > max.matches){
      max_scores <- max_scores[which(max_scores$max_score >= quantile(max_scores$max_score, 1 - (max.matches/nrow(max_scores)))),]
    }
    return(max_scores)
  } else{
    return(data.frame("id"=NULL,"max_score"=NULL))
  }}


get_motifs <- function(species = "Hsapiens", dataSource = "JASPAR"){

  if (length(species)!=1){
    stop("Must input species")}

  species_indices = grep(species, S4Vectors::values(MotifDb::MotifDb)$organism)

  if (length(species_indices) < 1){
    stop("No motifs found species. Perhaps incorrect species name used?")}

  if (length(dataSource)==1){
    data_indices <- grep(dataSource, S4Vectors::values(MotifDb::MotifDb)$dataSource)
  } else if (length(dataSource>1)){
    data_indices <- do.call(union,lapply(dataSource, function(x) grep(x, S4Vectors::values(MotifDb::MotifDb)$dataSource)))
  } else{
    data_indices <- species_indices }

  if (length(data_indices) < 1){
    stop("No motifs found for given dataSource. Perhaps incorrect spelling/name?")
  }

  motifs <- MotifDb::MotifDb[intersect(species_indices, data_indices)]

  return(motifs)
}


get_motif_indices <- function(motifs, peaks, genome =  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, BPPARAM = BiocParallel::bpparam(), min.score = "90%",
                              max.matches = length(peaks)/5, min.matches = 50){

  #get sequences
  seqs <- Biostrings::getSeq(genome,peaks)

  #get nucleotide frequencies
  nucFreqs <-  get_nuc_freqs(seqs)

  #find motif matches
  motif_dfs <- BiocParallel::bplapply(motifs, get_motif_matches, seqs, nucFreqs, min.score, max.matches, BPPARAM = BPPARAM)
  n_matches <- lapply(motif_dfs, nrow)
  keep <- which(n_matches > min.matches)

  motif_indices <- lapply(motif_dfs, function(x) x$id)[keep]
  motif_scores <- lapply(motif_dfs, function(x) x$score)[keep]

  return(list(indices = motif_indices, scores = motif_scores))
}








