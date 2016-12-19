
match_pwms_helper <- function(pwms, seqs, bg, p.cutoff, w, out, ranges){
  
  motif_mats <- convert_pwms(pwms, bg)
  
  if (out == "matches"){
    tmp_out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
    if (is.null(ranges)){
      out <- SummarizedExperiment(assays = list(matches = tmp_out),
                                  colData = DataFrame(pwm = pwms, name = TFBSTools::name(pwms), 
                                                      row.names = names(pwms)))
    } else{
      out <- SummarizedExperiment(assays = list(matches = tmp_out),
                                rowRanges = ranges, 
                                colData = DataFrame(pwm = pwms, name = TFBSTools::name(pwms), 
                                                    row.names = names(pwms)))
    }
  } else if (out =="scores"){
    tmp_out <- get_motif_ix_plus(motif_mats,seqs,bg,p.cutoff,w)
    if (is.null(ranges)){
      out <- SummarizedExperiment(assays = tmp_out,
                                  colData = DataFrame(pwm = pwms, name = TFBSTools::name(pwms), 
                                                      row.names = names(pwms)))
    } else{
      out <- SummarizedExperiment(assays = tmp_out,
                                  rowRanges = ranges, 
                                  colData = DataFrame(pwm = pwms, name = TFBSTools::name(pwms), 
                                                      row.names = names(pwms)))
    }
  } else{
    tmp_out <- get_motif_positions(motif_mats,seqs,bg,p.cutoff,w)
    if (is.null(ranges)){
      out <- lapply(1:length(motif_mats), function(x){
        m_ix <- which(tmp_out$motif_ix == x - 1)
        tmp <- IRanges::IRanges(start = tmp_out$pos[m_ix] + 1,
                                width = ncol(motif_mats[[x]]))
        mcols(tmp) <- DataFrame(strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])
        tmp
      })
      names(out) <- names(motif_mats)
    } else{
      out <- lapply(1:length(motif_mats), function(x){
      m_ix <- which(tmp_out$motif_ix == x - 1)
      GenomicRanges::GRanges(GenomicRanges::seqnames(ranges)[tmp_out$seq_ix[m_ix] + 1],
                             IRanges::IRanges(start = GenomicRanges::start(ranges[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                              width = ncol(motif_mats[[x]])),
                             strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])
      })
    }
    names(out) <- names(pwms)
  }
  return(out)
}  




#' match_pwms
#'
#' Find pwm matches
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}}, \code{\link[TFBSTools]{PFMatrixList}},
#' \code{\link[TFBSTools]{PWMatrix}}, \code{\link[TFBSTools]{PWMatrixList}}
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}}, \code{\link[Biostrings]{DNAStringSet}},
#' \code{\link[Biostrings]{DNAString}}, or character vector
#' @param genome BSgenome object, only used if subect is \code{\link[GenomicRanges]{GenomicRanges}}
#' @param bg background nucleotide frequencies. if not provided, computed from subject
#' @param out what to return? see details
#' @param p.cutoff p-value cutoff for returning motifs
#' @param w parameter controlling size of window for filtration; default is 7
#' @param if subject is not GenomicRanges, ranges to use when out is positions
#' @details Can either return a SummarizedExperiment with just sparse matrix with values set to 1 for a match (if return == "matches"), 
#' a SummarizedExperiment with a matches matrix as well as matrices with the maximum motif score and total motif counts (if return == "scores"), or
#' a GenomicRanges or IRanges object with all the positions of matches   
#' or \code{\link[GenomicRanges]{GenomicRanges}} if positions
#' @import Biostrings
#' @import TFBSTools
#' @import Matrix
#'@export
setGeneric("match_pwms", function(pwms, subject,...) standardGeneric("match_pwms"))

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAStringSet"),
          function(pwms, subject, genome = NULL, bg = NULL, out = c("matches","scores","positions"), p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)
            
            match_pwms_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "character"),
          function(pwms, subject, genome = NULL, bg = NULL,
                   out = c("matches","scores","positions"), p.cutoff = 0.00005, w = 7, ranges = NULL){
            
            out = match.arg(out)
            
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            match_pwms_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAString"),
          function(pwms, subject, genome = NULL,  bg = NULL,out = c("matches","scores","positions"),p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            
            seqs <- as.character(subject)
            
            match_pwms_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL, out = c("matches","scores","positions"),p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            GenomicRanges::strand(subject) <- "+"
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            
            seqs <- as.character(seqs)
            
            match_pwms_helper(pwms, seqs, bg, p.cutoff, w, out, subject)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "RangedSummarizedExperiment"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,  bg = NULL,out = c("matches","scores","positions"),p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            match_pwms(pwms, rowRanges(subject), genome, bg, out, p.cutoff, w)
          })

### PFMatrixList ---------------------------------------------------------------


#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "ANY"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL,out = c("matches","scores","positions"), p.cutoff = 0.00005, w =7, ranges = NULL){
            out <- match.arg(out)
            pwms_list <- do.call(TFBSTools::PWMatrixList, lapply(pwms, TFBSTools::toPWM))
            match_pwms(pwms_list, subject, genome = genome, bg = bg, out = out, p.cutoff = p.cutoff, w = w, ranges = ranges)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "ANY"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL,out = c("matches","scores","positions"), p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            pwms_list <- TFBSTools::PWMatrixList(pwms)
            match_pwms(pwms_list, subject, genome = genome, bg = bg, out = out, p.cutoff = p.cutoff, w = w, ranges = ranges)
          })


# Single PFM -------------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "ANY"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL,out = c("matches","scores","positions"), p.cutoff = 0.00005, w = 7, ranges = NULL){
            out = match.arg(out)
            pwms_list <- TFBSTools::PWMatrixList(TFBSTools::toPWM(pwms))
            match_pwms(pwms_list, subject, genome = genome, bg = bg, out = out, p.cutoff = p.cutoff, w = w, ranges = ranges)
          })
