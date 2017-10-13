get_reverse_complement <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

ambiguity_mapping <- lapply(Biostrings::IUPAC_CODE_MAP, function(x) {
  letters <- strsplit(x, "")[[1]]
  out <- rep(0, length(letters))
  out[which(letters == "A")] <- 1
  out[which(letters == "C")] <- 2
  out[which(letters == "G")] <- 3
  out[which(letters == "T")] <- 4
  return(out)
})


seq_to_pwm <- function(in_seq, mismatch = 0) {
  mat <- matrix(mismatch, ncol = nchar(in_seq), nrow = 4)
  rownames(mat) <- c("A", "C", "G", "T")
  in_seq <- strsplit(in_seq, "")[[1]]
  for (x in seq_along(in_seq)) {
    mat[ambiguity_mapping[[in_seq[x]]], x] <- 1
  }
  return(mat)
}



#' deviationsCovariability
#'
#' @param object deviations result
#'
#' @return 'covariability' matrix
#' @details Returns the 'covariability' between motifs/kmers/peaksets. 
#' Covariability' is defined as covariance between Z-scores divided by variance 
#' of Z-scores for one motif/kmer/peakset (the row).
#' @export
#' @examples
#' # load very small example data
#' data(mini_counts, package = "chromVAR")
#' motifs <- getJasparMotifs()
#' library(motifmatchr)
#' motif_ix <- matchMotifs(motifs, mini_counts, 
#'   genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#'
#' # computing deviations
#' dev <- computeDeviations(object = mini_counts, 
#'                          annotations = motif_ix)
#'                          
#' # get covariability for just first three motifs                         
#' devcov <- deviationsCovariability(dev[1:3,])                         
deviationsCovariability <- function(object) {
  covs <- cov(t(assays(object)$z))
  vars <- row_sds(assays(object)$z)
  normed_covs <- covs/matrix(vars^2, nrow = nrow(covs), ncol = ncol(covs), 
                             byrow = FALSE)
  return(normed_covs)
}

#' @importFrom Biostrings nucleotideSubstitutionMatrix stringDist
get_kmer_dist <- function(kmers) {
  stopifnot(all_equal(nchar(kmers)))
  out <- as.matrix(stringDist(kmers, 
                              method = "substitutionMatrix", 
                              type = "overlap", 
                              gapOpening = Inf, 
                              substitutionMatrix = 
                                nucleotideSubstitutionMatrix(match = 1, 
                                                             mismatch = 0,
                                                             baseOnly = FALSE, 
                                                             type = "DNA")))
  out <- nchar(kmers[1]) - out
  diag(out) <- 0
  return(out)
}

get_mm_kmers <- function(kmer) {
  k <- nchar(kmer)
  nucs <- c("A", "C", "G", "T")
  kmer_nucs <- strsplit(kmer, "")[[1]]
  out <- c()
  for (i in seq_len(k)) {
    for (j in nucs[nucs != kmer_nucs[i]]) {
      kmod <- kmer_nucs
      kmod[i] <- j
      kmod <- paste(kmod, sep = "", collapse = "")
      out <- c(out, kmod)
    }
  }
  return(out)
}

expected_n_overlap_kmers <- function(max_extend){
  if (max_extend == 0){
    return(0)
  } else{
    return(4**(max_extend) + expected_n_overlap_kmers(max_extend - 1))
  }
}

get_overlap_kmers <- function(kmer, max_extend, dir = "both") {
  out <- c()
  if (max_extend <= 0) {
    return(out)
  }
  k <- nchar(kmer)
  nucs <- c("A", "C", "G", "T")
  if (max_extend == 1) {
    if (dir == "both") {
      out1 <- vapply(nucs, 
                     function(x) paste(x, substr(kmer, 1, k - 1), sep = "", 
                                             collapse = ""), 
                     "", 
                     USE.NAMES = FALSE)
      out2 <- vapply(nucs, 
                     function(x) paste(substr(kmer, 2, k), x, sep = "", 
                                             collapse = ""), 
                     "", 
                     USE.NAMES = FALSE)
      return(c(out1, out2))
    } else if (dir == "left") {
      out1 <- vapply(nucs, 
                     function(x) paste(x, substr(kmer, 1, k - 1), sep = "", 
                                             collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      return(out1)
    } else if (dir == "right") {
      out2 <- vapply(nucs, 
                     function(x) paste(substr(kmer, 2, k), x, sep = "", 
                                             collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      return(out2)
    }
  } else {
    if (dir == "both") {
      out1 <- vapply(nucs, 
                     function(x) paste(x, substr(kmer, 1, k - 1), sep = "", 
                                       collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      out2 <- vapply(nucs, 
                     function(x) paste(substr(kmer, 2, k), x, sep = "", 
                                       collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      return(c(out1, out2, 
               vapply(out1, 
                      get_overlap_kmers, 
                      rep("", expected_n_overlap_kmers(max_extend - 1)),
                      max_extend = max_extend - 1, 
                      dir = "left"), 
               vapply(out2, 
                      get_overlap_kmers, 
                      rep("", expected_n_overlap_kmers(max_extend - 1)),
                      max_extend = max_extend - 1, 
                      dir = "right")))
    } else if (dir == "left") {
      out1 <- vapply(nucs, 
                     function(x) 
                       paste(x, substr(kmer, 1, k - 1), sep = "", 
                                             collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      return(c(out1, vapply(out1, 
                            get_overlap_kmers, 
                            rep("", expected_n_overlap_kmers(max_extend - 1)),
                            max_extend = max_extend - 1, 
                            dir = "left")))
    } else if (dir == "right") {
      out2 <- vapply(nucs, 
                     function(x) 
                       paste(substr(kmer, 2, k), x, sep = "", collapse = ""), 
                     "",
                     USE.NAMES = FALSE)
      return(c(out2, 
               vapply(out2, 
                      get_overlap_kmers, 
                      rep("", expected_n_overlap_kmers(max_extend - 1)),
                      max_extend = max_extend - 1, 
                      dir = "right")))
    }
  }
}

kmers_to_names <- function(kmers, names) {
  ifelse(kmers %in% names, kmers, get_reverse_complement(kmers))
}

get_null_kmer_dist <- function(kmer, cov_mat, max_extend = 2) {
  kmers <- colnames(cov_mat)
  kmer <- kmers_to_names(kmer, kmers)
  o <- get_overlap_kmers(kmer, max_extend = max_extend)
  mm <- get_mm_kmers(kmer)
  
  omm <- 
    Reduce(union, 
           lapply(mm, 
                  function(k) 
                    unique(get_overlap_kmers(k,  max_extend = max_extend + 1))
           )
    )
  
  null_dist <- cov_mat[kmer, which(colnames(cov_mat) %ni% 
                                     kmers_to_names(unique(c(kmer, mm, o, omm)),
                                                    kmers))]
  return(null_dist)
}


get_covariable_kmers <- function(kmer, cov_mat, max_extend = 2) {
  
  # cov_mat <- deviationsCovariability(object)
  kmers <- colnames(cov_mat)
  kmer <- kmers_to_names(kmer, kmers)
  # Get Null Dist
  o <- get_overlap_kmers(kmer, max_extend = max_extend)
  mm <- get_mm_kmers(kmer)
  
  omm <- 
    Reduce(union, 
           lapply(mm, 
                  function(k) 
                    unique(get_overlap_kmers(k, max_extend = max_extend +  1))
           ))
  
  null_dist <- cov_mat[kmer, which(colnames(cov_mat) %ni%
                                     kmers_to_names(unique(c(kmer, 
                                                             mm, o, omm)), 
                                                    kmers))]
  
  candidates <- unique(c(mm, o))
  scores <- (cov_mat[kmer, unique(kmers_to_names(candidates, kmers))] - 
               mean(null_dist, 
                    na.rm = TRUE))/sd(null_dist, na.rm = TRUE)
  pvals <- pnorm(scores, lower.tail = FALSE)
  pvals.adj <- p.adjust(pvals)
  
  o_shifts <- if (max_extend >= 1) 
    do.call(c, lapply(seq_len(max_extend), 
                      function(x) c(rep(-x, 4^x), rep(x, 4^x)))) else NULL
  out <- data.frame(kmer = c(kmer, mm, o),
                    mismatch = c(NA, rep(seq_len(nchar(kmer)), 
                                         each = 3), 
                                 rep(NA, length(o))), 
                    shift = c(rep(0, length(mm) + 1), 
                              o_shifts), 
                    covariability = cov_mat[kmer, kmers_to_names(c(kmer, mm, o),
                                                                 kmers)], 
                    pval = pvals[kmers_to_names(c(kmer, 
                                                  mm, o), kmers)],
                    pval.adj = pvals.adj[kmers_to_names(c(kmer, mm, o), 
                                                        kmers)], 
                    stringsAsFactors = FALSE)
  
  return(out)
}


kmer_group_to_pwm <- function(kgroup, p = 0.01, threshold = 0.25) {
  min_shift <- min(kgroup$shift)
  max_shift <- max(kgroup$shift)
  k <- nchar(kgroup$kmer[1])
  nucs <- c("A", "C", "G", "T")
  out <- matrix(0, nrow = 4, ncol = k + abs(min_shift) + max_shift, 
                dimnames = list(nucs, 
                                NULL))
  if (min_shift < 0) {
    for (i in min_shift:-1) {
      ix <- which(kgroup$shift == i)
      nuc <- substr(kgroup$kmer[ix], 1, 1)
      covars <- vapply(nucs, 
                       function(x) 
                         max(c(0, kgroup$covariability[ix[which(nuc == x)]])),
                       0)
      covars[!is.finite(covars)] <- 0
      baseline <- rep(0.25, 4) #* (1 + covars)
      baseline <- baseline/sum(baseline)
      pvals <- vapply(nucs, 
                      function(x) min(kgroup$pval.adj[ix[which(nuc == x)]]),
                      0)
      covars[pvals > p] <- 0
      w <- max(covars)
      if (w > 1) 
        w <- 0 #1
      out[, i - min_shift + 1] <- 
        (1 - w) * baseline + 
        (w * covars^2/(sum(covars^2) +  all_true(covars == 0)))
    }
  }
  for (i in seq_len(k)) {
    ix <- c(1, which(kgroup$mismatch == i))
    nuc <- substr(kgroup$kmer[ix], i, i)
    covars <- vapply(nucs, 
                     function(x) 
                       max(c(0, kgroup$covariability[ix[which(nuc == x)]])),
                     0)
    covars[!is.finite(covars)] <- 0
    pvals <- vapply(nucs, 
                    function(x) min(kgroup$pval.adj[ix[which(nuc == x)]]),
                    0)
    covars[which(pvals > p & nucs != nuc[1])] <- 0
    out[, abs(min_shift) + i] <- covars^2/sum(covars^2)
  }
  if (max_shift > 0) {
    for (i in seq_len(max_shift)) {
      ix <- which(kgroup$shift == i)
      nuc <- substr(kgroup$kmer[ix], k, k)
      covars <- vapply(nucs, 
                       function(x) 
                         max(c(0, kgroup$covariability[ix[which(nuc == x)]])),
                       0)
      covars[!is.finite(covars)] <- 0
      baseline <- rep(0.25, 4) #* (1 + covars)
      baseline <- baseline/sum(baseline)
      pvals <- vapply(nucs, 
                      function(x) min(kgroup$pval.adj[ix[which(nuc ==  x)]]),
                      0)
      covars[which(pvals > p)] <- 0
      w <- max(covars)
      if (w > 1) 
        w <- 0 #1
      out[, abs(min_shift) + k + i] <- 
        (1 - w) * baseline + 
        (w * covars^2/(sum(covars^2) + all_true(covars == 0)))
    }
  }
  start <- 1
  bits <- function(x) {
    2 + sum(x * ifelse(x > 0, log2(x), 0))
  }
  if (min_shift < 0) {
    for (i in min_shift:-1) {
      if (bits(out[, i - min_shift + 1]) < threshold) {
        start <- start + 1
      } else {
        break
      }
    }
  }
  end <- ncol(out)
  if (max_shift > 0) {
    for (i in max_shift:1) {
      if (bits(out[, abs(min_shift) + k + i]) < threshold) {
        end <- end - 1
      } else {
        break
      }
    }
  }
  return(out[, start:end])
}




#' assembleKmers
#'
#' function to create de novo motifs from kmers based on deviations
#' @param object kmer chromVARDeviations object
#' @param threshold variability threshold
#' @param p p value threshold for inclusion of kmer
#' @param progress show progress bar?
#' @details function for assembling de novo kmers from kmer deviations
#' @return list with (1) motifs: de novo motif matrices, (2) seed: seed kmer
#'  for de novo motif 
#' @importFrom TFBSTools PWMatrix
#' @export
assembleKmers <- function(object, threshold = 1.5, p = 0.01, progress = TRUE) {
  devco <- deviationsCovariability(object)
  vars <- row_sds(assays(object)$z)
  cands <- 
    rownames(object)[order(vars, 
                           decreasing = TRUE)[seq_len(sum(vars > threshold, 
                                                          na.rm = TRUE))]]
  out <- list(motifs = list(), seed = list())
  nc <- length(cands)
  if (progress) pb <- txtProgressBar(min = 0, max = nc, style = 3)
  while (length(cands) > 1) {
    kgroup <- get_covariable_kmers(cands[1], devco)
    kmotif <- kmer_group_to_pwm(kgroup, p)
    if (max(kgroup$covariability, na.rm = TRUE) <= 1) {
      out$motifs <- c(out$motifs, list(kmotif))
      out$seed <- c(out$seed, cands[1])
    }
    nd <- get_null_kmer_dist(cands[1], devco)
    z <- (devco[cands[1], cands] - mean(nd, na.rm = TRUE)) / 
      sd(nd, na.rm = TRUE)
    pval <- pnorm(z, lower.tail = FALSE)
    d <- pwmDistance(kmotif, lapply(cands, seq_to_pwm))$d[1, ]
    exc <- intersect(which(d < 0.25), which(pval < p))
    cands <- cands[-exc]
    if (progress) setTxtProgressBar(pb, nc - length(cands))
  }
  if (progress) close(pb)
  denovo_motifs <- 
    do.call(PWMatrixList, 
            lapply(seq_along(out$motifs), 
                   function(x) PWMatrix(ID = paste0("denovo_", x), 
                                        name = paste0("denovo_",  x), 
                                        tags = list(seed = out$seed[[x]]), 
                                        profileMatrix = out$motifs[[x]])))
  names(denovo_motifs) <- name(denovo_motifs)
  return(denovo_motifs)
}

#' plotKmerMismatch
#'
#' @param kmer kmer, e.g. 'AAAAAAA'
#' @param cov_mat result from \code{\link{deviationsCovariability}}
#' @param pval p value threshold
#'
#' @return A plot
#' @export
plotKmerMismatch <- function(kmer, cov_mat, pval = 0.01) {
  
  kgroup <- get_covariable_kmers(kmer, cov_mat, max_extend = 0)
  ix <- which(!is.na(kgroup$mismatch))
  
  mm_df <- rbind(
    data.frame(val = ifelse(kgroup$covariability[ix] > 0, 
                            kgroup$covariability[ix]^2, 
                            0), pos = kgroup$mismatch[ix], 
               Nucleotide = substr(kgroup$kmer[ix], kgroup$mismatch[ix], 
                                   kgroup$mismatch[ix])), 
    data.frame(val = 1, pos = seq_len(nchar(kmer)), 
               Nucleotide = strsplit(kmer,  "")[[1]]))
  
  out <- ggplot(mm_df) + 
    geom_point(aes_string(x = "pos", y = "val", col = "Nucleotide"), 
               position = position_jitter(height = 0,  width = 0.1)) + 
    ylab("Shared variability with\nseed nucleotide") + xlab("Position") + 
    scale_x_continuous(breaks = c(seq_len(max(mm_df$pos)))) + 
    scale_y_continuous(breaks = c(0, 
                                  0.5, 1)) + 
    chromVAR_theme() + 
    scale_color_manual(name = "Nucleotide", breaks = c("A",   "C", "G", "T"), 
                       values = RColorBrewer::brewer.pal(4, "Dark2")) 
  return(out)
}



toIC <- function(mat) {
  mat * matrix(apply(mat, 2, function(x) 2 + sum(log(x) * x)), 
               nrow = nrow(mat), 
               ncol = ncol(mat), byrow = TRUE)
}
