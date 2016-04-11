# Not exported
get_nuc_freqs <- function(seqs){
  #get nucleotide frequencies
  nucFreqs <- colSums(Biostrings::letterFrequency(seqs, c("A","C","G","T")))
  nucFreqs <- nucFreqs/sum(nucFreqs)
  return(nucFreqs)
}

#' get_motifs
#' 
#' Function to get motifs from JASPAR database
#' @param species Which species?  use eithe jaspar code or latin name. default is "Homo sapiens"
#' @param collection Which collection to use?  default is "CORE"
#' @param ... additional arguments to opts for \code{\link[TFBSTools]{getMatrixSet}}
#' @details Simply a wrapper function for \code{\link[TFBSTools]{getMatrixSet}} that calls
#' JASPAR2014 database using \code{\link[JASPAR2014]{JASPAR2014}}
#' @return \code{\link[TFBSTools]{PFMatrixList}}
#' @export
get_motifs <- function(species = "Homo sapiens", collection = "CORE", ...){
  opts = list()
  opts['species'] = species
  opts['collection'] = collection
  opts = c(opts, list(...))
  TFBSTools::getMatrixSet(JASPAR2014::JASPAR2014, opts)
}


#' get_motif_indices
#' 
#' Function to get indices of peaks that contain motifs
#' @param motifs \code{\link[TFBSTools]{PFMatrixList}} or \code{\link[TFBSTools]{PWMatrixList}}
#' @param peaks \code{\link[GenomicRanges]{GenomicRanges}} 
#' @param genome \code{\link[BSgenome]{BSgenome}} object
#' @param p.cutoff default is 0.00005
#' @return A list with a vector of peak indices for each motif.
#'@export
get_motif_indices <- function(motifs, 
                               peaks, 
                               genome =  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                               p.cutoff = 0.00005){
  stopifnot(inherits(peaks,"GenomicRanges"))
  # If PFMatrixList or PFMatrix, convert to PWM
  if (inherits(motifs,"PFMatrixList")){
    motifs <- do.call(TFBSTools::PWMatrixList,lapply(motifs, TFBSTools::toPWM))
  } else if (inherits(motifs, "PFMatrix")){
    motifs <- TFBSTools::toPWM(motifs)
  }
  # Check class
  if (!(inherits(motifs, "PWMatrixList") | inherits(motifs, "PWMatrix"))){
    stop("motifs must be PFMatrixList, PWMatrixList, PFMatrix, or PWMatrix object")
  }
  
  seqs <- Biostrings::getSeq(genome, peaks)
  
  #get nucleotide frequencies
  nucFreqs <-  get_nuc_freqs(seqs)
  
  seqs <- as.character(seqs)
  
  motif_mats <- lapply(motifs, function(x) TFBSTools::as.matrix(x))
  
  grpsize <- BiocParallel::bpparam()
  if (grpsize > 2){
    motif_mat_grps <- lapply(1:(length(motif_mats) %/% grpsize + ((length(motif_mats) %% grpsize)!=0)), function(x) motif_mats[((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak))])
    motif_ix = do.call(c, BiocParallel::bplapply(motif_mats,
                                      motif_match_helper,
                                      seqs,
                                      nucFreqs,
                                      p.cutoff))
  } else{
    motif_ix = multi_motif_match_helper(motif_mats, seqs, nucFreqs, p.cutoff)
  }

  return(motif_ix)
}

motif_match_helper <- function(motif_mats, seqs, nuc_freqs, p.cutoff){
  out = multi_motif_match(motif_mats, seqs, nuc_freqs, p.cutoff)
  names(out) = names(motif_mats)
  return(out)
}


#' get_max_motif_scores
#' 
#' Function to get max motif score within each peak
#' @param motif \code{\link[TFBSTools]{PFMatrix}} or \code{\link[TFBSTools]{PWMatrix}}
#' @param peaks \code{\link[GenomicRanges]{GenomicRanges}} 
#' @param genome \code{\link[BSgenome]{BSgenome}} object
#' @return A numeric vector with the max motif score for each peak
#' @export
get_max_motif_scores <- function(motif, 
                                 peaks, 
                                 genome=  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                                 p = 0.00005){
  
  if (inherits(motifs, "PFMatrix")){
    motif <- TFBSTools::toPWM(motif)
  }
  # Check class
  if (!inherits(motifs, "PWMatrix")){
    stop("motif must be PFMatrix or PWMatrix object")
  }
  seqs <- Biostrings::getSeq(genome, peaks)  
  nucFreqs <-  get_nuc_freqs(seqs)  
  seqs <- as.character(seqs)
  motif_mat <- TFBSTools::as.matrix(motif)
  out <- motif_match_score(motif_mat, seq_ixs, nucFreqs, p)
  names(out) = c("ix","scores")
  return(out)
}


get_motif_positions <- function(motif, 
                                peaks,
                                matches = NULL,
                                genome=  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                                p.cutoff = 0.00005,
                                top = TRUE){
  if (inherits(motif, "PFMatrix")){
    motif <- TFBSTools::toPWM(motif)
  }
  # Check class
  if (!inherits(motif, "PWMatrix")){
    stop("motif must be PFMatrix or PWMatrix object")
  }
  seqs <- Biostrings::getSeq(genome, peaks)
  nucFreqs <-  get_nuc_freqs(seqs)
  motif_mat <- TFBSTools::as.matrix(motif)
  minScore <- p_to_score(motif_mat, nucFreqs, p.cutoff)
  if (is.null(matches)){
    matches = motif_match(motif_mat, as.character(seqs), nucFreqs, p.cutoff)
  }
  motif_pos <- do.call(c,BiocParallel::bplapply(matches, function(x){
    positions = GenomicRanges::GRanges()
    forward_matches <- Biostrings::matchPWM(motif_mat, seqs[[x]], 
                                            min.score = minScore, with.score = TRUE)
    if (length(forward_matches) > 0){
      tmp_positions = GenomicRanges::resize(GenomicRanges::shift(rep(peaks[x],length(forward_matches)), shift = start(forward_matches)-1), width(forward_matches),fix = "start")
      strand(tmp_positions) = "+"
      mcols(tmp_positions)$score = mcols(forward_matches)$score
      positions = c(positions, tmp_positions)
    }
    reverse_matches <- Biostrings::matchPWM(Biostrings::reverseComplement(motif_mat),
                                            seqs[[x]], 
                                            min.score = minScore, with.score = TRUE)
    if (length(reverse_matches) > 0){
      tmp_positions = GenomicRanges::resize(GenomicRanges::shift(rep(peaks[x],length(reverse_matches)), shift = start(reverse_matches)-1), width(reverse_matches),fix = "start")
      strand(tmp_positions) = "-"
      mcols(tmp_positions)$score = mcols(reverse_matches)$score
      positions = c(positions, tmp_positions)
    }
    if (top){
      top_match = positions[which.max(mcols(positions)$score)]
      return(top_match)
    } else{
      return(positions)      
    }
  }))
  return(motif_pos)
}



