
match_kmers <- function(kmers, seqs, tb.start = NA, tb.end = NA, max.mismatch = NA){
  stopifnot(all_false(duplicated(kmers)))
  if (is.character(kmers)) kmers = Biostrings::DNAStringSet(kmers)
  stopifnot(inherits(kmers,"DNAStringSet"))
  ##Make PDict
  if (is.na(tb.start)){
    pd <- Biostrings::PDict(kmers)
    indices  <- Biostrings::vwhichPDict(pd,seqs)
    indices_rc <- Biostrings::vwhichPDict(pd,Biostrings::reverseComplement(seqs))
  } else{
    pd <- Biostrings::PDict(kmers, tb.start = tb.start, tb.end = tb.end)
    indices  <- Biostrings::vwhichPDict(pd,
                                        seqs,
                                        fixed="subject",
                                        max.mismatch = max.mismatch)
    indices_rc <- Biostrings::vwhichPDict(pd,
                                          Biostrings::reverseComplement(seqs),
                                          fixed="subject",
                                          max.mismatch = max.mismatch)
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

# k-mer alignment --------------------------------------------------------------

kmer_dist <- function(kmers){
  out = as.matrix(Biostrings::stringDist(kmers))
  for (i in seq_along(kmers)){
    tmp = kmers
    tmp[[i]] <- Biostrings::reverseComplement(tmp[[i]])
    tmp_dist = as.matrix(Biostrings::stringDist(tmp))
    out[,i] = mapply(min, out[,i], tmp_dist[,i])
    out[i,] = mapply(min, out[i,], tmp_dist[i,])
  }
  as.dist(out)
}

kmer_overlap_dist <- function(kmers, kmer = NULL){
  if (is.character(kmers)) kmers = Biostrings::DNAStringSet(kmers)
  if (is.null(kmer)){
    out = as.matrix(Biostrings::stringDist(kmers, method = "substitutionMatrix",
                                           substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                         mismatch = 0,
                                                                                                         baseOnly = FALSE,
                                                                                                         type = "DNA"),
                                           type="overlap",gapOpening=-Inf))
    for (i in seq_along(kmers)){
      tmp = kmers
      tmp[[i]] <- Biostrings::reverseComplement(tmp[[i]])
      tmp_dist = as.matrix(Biostrings::stringDist(tmp, method = "substitutionMatrix",
                                                  substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                                mismatch = 0,
                                                                                                                baseOnly = FALSE,
                                                                                                                type = "DNA"),
                                                  type="overlap",gapOpening=-Inf))
      out[,i] = mapply(max, out[,i], tmp_dist[,i])
      out[i,] = mapply(max, out[i,], tmp_dist[i,])
      out[i,i] = Biostrings::nchar(kmers[[i]])
    }
  } else {
    if (is.character(kmer)) kmer = Biostrings::DNAStringSet(kmer)
    stopifnot(inherits(kmer,"DNAStringSet"))
    tmp1 = sapply(kmers, function(x)
      Biostrings::pairwiseAlignment(kmer,
                                    x,
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = 0,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf, scoreOnly = TRUE))
    tmp2 = sapply(kmers, function(x)
      Biostrings::pairwiseAlignment(Biostrings::reverseComplement(kmer),
                                    x,
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = 0,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf, scoreOnly = TRUE))
    out = mapply(max,tmp1,tmp2)
  }
  out
}





kmer_overlap_match <- function(kmers, kmer = NULL){
  if (is.character(kmers)) kmers = Biostrings::DNAStringSet(kmers)
  if (is.null(kmer)){
    out = as.matrix(Biostrings::stringDist(kmers, method = "substitutionMatrix",
                                         substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                       mismatch = -Inf,
                                                                                                       baseOnly = FALSE,
                                                                                                       type = "DNA"),
                                         type="overlap",gapOpening=-Inf))
    for (i in seq_along(kmers)){
      tmp = kmers
      tmp[[i]] <- Biostrings::reverseComplement(tmp[[i]])
      tmp_dist = as.matrix(Biostrings::stringDist(tmp, method = "substitutionMatrix",
                                                substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                              mismatch = -Inf,
                                                                                                              baseOnly = FALSE,
                                                                                                              type = "DNA"),
                                                type="overlap",gapOpening=-Inf))
      out[,i] = mapply(max, out[,i], tmp_dist[,i])
      out[i,] = mapply(max, out[i,], tmp_dist[i,])
      out[i,i] = Biostrings::nchar(kmers[[i]])
    }
  } else {
      if (is.character(kmer)) kmer = Biostrings::DNAStringSet(kmer)
      stopifnot(inherits(kmer,"DNAStringSet"))
      tmp1 = sapply(kmers, function(x)
                  Biostrings::pairwiseAlignment(kmer,
                                    x,
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = -Inf,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf, scoreOnly = TRUE))
      tmp2 = sapply(kmers, function(x)
      Biostrings::pairwiseAlignment(Biostrings::reverseComplement(kmer),
                                    x,
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = -Inf,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf, scoreOnly = TRUE))
      out = mapply(max,tmp1,tmp2)
  }
  out
}




extend_kmer <- function(kmer, sets, counts_mat, bg_peaks,
                        diff = 2,
                        niterations = 50,
                        BPPARAM = BiocParallel::bpparam()){

  ## Get kmers that overlap
  min_overlap = nchar(kmer) - diff
  kmer_dists = kmer_overlap_match(names(sets), kmer)
  close = which(kmer_dists >= min_overlap)
  if (length(close) >= 2){
    tmpsets = remove_nonoverlap(sets[close], kmer)
    tmpresults <- compute_variability(tmpsets, counts_mat, bg_peaks, niterations = niterations,
                                    BPPARAM = BPPARAM)
    tmpvar <- variability(tmpresults)

    setlen = sapply(tmpsets,length)
    nbg = 50
    bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(sets[[kmer]], x, replace=FALSE))), recursive = F)

    bgresults <- compute_variability(bgsets, counts_mat, bg_peaks, niterations = niterations,
                                   BPPARAM = BPPARAM)
    bgvar <- variability(bgresults)

    var_boost <- sapply(seq_along(tmpvar),
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
    var_boost <- pnorm(var_boost, lower.tail = FALSE)
    out_ix <- which(var_boost < 0.05)
    out = var_boost[out_ix]
    names(out) = names(tmpsets[out_ix])
    return(out)
  } else{
    return(NULL)
  }
}


add_left <- function(letter, kmer, gap = 0){
  return(BiocGenerics::paste(c(letter,rep("N",gap),kmer),collapse=""))
}

add_right <- function(letter, kmer, gap = 0){
  return(BiocGenerics::paste(c(kmer,rep("N",gap),letter),collapse=""))
}

#' @export
extend_kmer2 <- function(kmer, set, counts_mat, bg_peaks, seqs, 
                        max_extend = 3,
                        niterations = 50,
                        p.cutoff = 0.01,
                        BPPARAM = BiocParallel::bpparam()){

  nucs = c("A","C","G","T")
  ntest = max_extend * 2 * 4
  p.cutoff = p.cutoff / ntest
  i = 1
  kmer_head = rep("N",max_extend)
  kmer_tail = rep("N", max_extend)
  while (i <= max_extend){
    #left
    candidates = sapply(nucs, add_left, kmer, i-1)
    tmpsets = match_kmers(candidates, seqs, tb.start = i + 1, tb.end = i + nchar(kmer), max.mismatch = i-1)
    var_boost <- get_variability_boost(set,tmpsets,counts_mat,bg_peaks, niterations, BPPARAM)
    left_boost <- nucs[which(var_boost < p.cutoff)]
    if (length(left_boost)>0){
      kmer_head[max_extend - i + 1] = Biostrings::consensusString(left_boost,
                                                            ambiguityMap=Biostrings::IUPAC_CODE_MAP,
                                                            threshold = 0.25)
    }
    #right
    candidates = sapply(nucs, add_right, kmer, i-1)
    tmpsets = match_kmers(candidates, seqs, tb.start = 1, tb.end = nchar(kmer), max.mismatch = i-1)
    var_boost <- get_variability_boost(set,tmpsets,counts_mat,bg_peaks, niterations, BPPARAM)
    right_boost <- nucs[which(var_boost< p.cutoff)]
    if (length(right_boost)>0){
      kmer_tail[i] = Biostrings::consensusString(right_boost,
                                                            ambiguityMap=Biostrings::IUPAC_CODE_MAP,
                                                            threshold = 0.25)
    }
    i = i +1
  }
  if (!all_true(kmer_head == "N")){
    kmer_head = BiocGenerics::paste(kmer_head[(which(kmer_head != "N")[1]):(length(kmer_head))],
                                    collapse="")
  } else{
    kmer_head =""
  }
  if (!all_true(kmer_tail == "N")){
    kmer_tail = BiocGenerics::paste(kmer_tail[1:(rev(which(kmer_tail != "N"))[1])],
                                    collapse="")
  } else{
    kmer_tail =""
  }
  return(list(head = kmer_head, tail = kmer_tail))}


#' @export
extend_kmer3 <- function(kmer, set, counts_mat, bg_peaks, seqs, 
                         max_extend = 3,
                         niterations = 50,
                         p.cutoff = 0.01,
                         BPPARAM = BiocParallel::bpparam()){
  stopifnot(length(seqs) == length(set))
  nucs = c("A","C","G","T")
  ntest = max_extend * 2 * 4
  p.cutoff = p.cutoff / ntest
  i = 1
  kmer_head = rep("N",max_extend)
  kmer_tail = rep("N", max_extend)
  tmpsets = list()
  while (i <= max_extend){
    #left
    candidates = sapply(nucs, add_left, kmer, i-1)
    tmpsets = c(tmpsets,match_kmers(candidates, seqs, tb.start = i + 1, tb.end = i + nchar(kmer), max.mismatch = i-1))
    #right
    candidates = sapply(nucs, add_right, kmer, i-1)
    tmpsets = c(tmpsets,match_kmers(candidates, seqs, tb.start = 1, tb.end = nchar(kmer), max.mismatch = i-1))
    i = i + 1
  }
  tmpsets = lapply(tmpsets, function(x) set[x])
  var_boost <- get_variability_boost(set,tmpsets,counts_mat,bg_peaks, niterations, BPPARAM)
  i = 1
  while (i <= max_extend){
    left_boost <- nucs[which(var_boost[(8*(i-1)+1):(8*(i-1)+4)] < p.cutoff)]
    if (length(left_boost)>0){
      kmer_head[max_extend - i + 1] = Biostrings::consensusString(left_boost,
                                                                  ambiguityMap=Biostrings::IUPAC_CODE_MAP,
                                                                  threshold = 0.25)
    }
    right_boost <- nucs[which(var_boost[(8*(i-1)+5):(8*(i-1)+8)] < p.cutoff)]
    if (length(right_boost)>0){
      kmer_tail[i] = Biostrings::consensusString(right_boost,
                                                 ambiguityMap=Biostrings::IUPAC_CODE_MAP,
                                                 threshold = 0.25)
    }
    i = i +1
  }
  if (!all_true(kmer_head == "N")){
    kmer_head = BiocGenerics::paste(kmer_head[(which(kmer_head != "N")[1]):(length(kmer_head))],
                                    collapse="")
  } else{
    kmer_head =""
  }
  if (!all_true(kmer_tail == "N")){
    kmer_tail = BiocGenerics::paste(kmer_tail[1:(rev(which(kmer_tail != "N"))[1])],
                                    collapse="")
  } else{
    kmer_tail =""
  }
  return(list(head = kmer_head, tail = kmer_tail))}




#' @export
get_similar_kmers <- function(kmer, kmers, cutoff = 4){
  dists = kmer_overlap_dist(kmers,kmer)
  return(kmers[which(dists >= cutoff)])
}

#' @export
get_overlapping_kmers <- function(kmer, kmers, cutoff = 2){
  dists = kmer_overlap_match(kmers,kmer)
  return(kmers[which(dists >= cutoff)])
}

cor_helper <- function(a, b, ...){
  res = cor.test(a, b,...)
  return(c(r = res$estimate[[1]],p = res$p.value))
}

#' @export
get_correlated_results <- function(name, results, p.cutoff = 0.01, corr.cutoff = 0.5){
  mat = deviations(results)
  ix = which(names(results) == name)
  pcors = apply(mat, 1, function(x) cor_helper(x, mat[ix,],alternative = "greater"))
  sig = intersect(which(pcors["p",] < 0.01), which(pcors["r",] > corr.cutoff))
  return(names(results)[sig[sig != ix]])
}


#' @export
get_associated_kmers <- function(kmer,
                                 sets,
                                 results,
                                 counts_mat,
                                 bg_peaks,
                                 seqs,
                                 max_extend = 3,
                                 niterations = 50,
                                 p.cutoff = 0.01,
                                 k.cutoff = 4,
                                 BPPARAM = BiocParallel::bpparam()){

  #Step -1: get extended kmer
  extensions = extend_kmer3(kmer,
                            sets[[kmer]],
                            counts_mat,
                            bg_peaks,
                            seqs,
                            max_extend = max_extend,
                            niterations = niterations,
                            p.cutoff = p.cutoff,
                            BPPARAM = BPPARAM)
  ekmer = BiocGenerics::paste(extensions[["head"]],kmer,extensions[["tail"]], collapse="",sep="")

  if (length(sets) >1){
    #Step 0: Test sequence similarity
    sim_seqs = get_similar_kmers(ekmer, names(sets), cutoff = k.cutoff)
    #overlap_seqs = get_overlapping_kmers(kmer, names(sets), cutoff = 2)
    #candidates1 = sim_seqs[which(sim_seqs %ni% overlap_seqs)]
    #candidates2 = sim_seqs[which(sim_seqs %in% overlap_seqs)]

    #Step 1: Test Correlation for non-overlapping but similar kmers
    cor_res = get_correlated_results(kmer, results[sim_seqs], p.cutoff = p.cutoff)

    #Step 2: Test Correlation for overlapping kmers
    #tmpsets <- remove_nonoverlap(sets[candidates2], kmer)
    #var_boost <- get_variability_boost(sets[[kmer]],tmpsets, counts_mat, bg_peaks, niterations, BPPARAM)
    #var_boost <- pnorm(var_boost, lower.tail = FALSE)
    #boosters <- names(var_boost)[which(var_boost <= (p.cutoff / length(tmpsets)))]
    #cor_res2 = get_correlated_results(kmer, results[candidates2], p.cutoff = p.cutoff)
    
    #return...
    out <- list(kmer = kmer, ekmer = ekmer, similar = cor_res)
  } else {
    out <- list(kmer = kmer, ekmer = ekmer, similar = NULL)
  }

  return(out)

}

#' @export
group_kmers <- function(sets,
                                 results,
                                 counts_mat,
                                 bg_peaks,
                                 max_extend = 3,
                                 genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                                 niterations = 50,
                                 p.cutoff = 0.01,
                                 k.cutoff = 4,
                                 BPPARAM = BiocParallel::bpparam()){

  if (length(sets) != length(results) || !all_true(all.equal(sort(names(sets)),sort(names(results))))){
    warning("Names of sets and results not the same... using intersection")
    common_names <- intersect(names(sets),names(results))
    sets <- sets[common_names]
    results <- results[common_names]
  }
  
  p <- get_pvalues(results, adjust = TRUE)
  #order by pvalues...
  candidates <- names(p)[sort(p, index.return=T)$ix]
  out <- list()
  seqs <- Biostrings::getSeq(genome, counts_mat@peaks)
  while (length(candidates)>0){
    kmer_assoc = get_associated_kmers(candidates[1], sets[candidates], results[candidates], counts_mat, bg_peaks, seqs[sets[[candidates[1]]]], max_extend, niterations, p.cutoff, k.cutoff, BPPARAM)
    toremove <- c(kmer_assoc$kmer, kmer_assoc$similar, kmer_assoc$boost)
    candidates <- candidates[which(candidates %ni% toremove)]
    out = c(out, list(kmer_assoc))
    print(length(candidates))
    print(candidates)
  }
  names(out) <- sapply(out, function(obj) obj$kmer)
  return(out)
}




#' @export
make_motif_from_kmer_group <- function(kmer_group){
  shifts = c(get_kmer_alignment(kmer_group$kmer, kmer_group$ekmer)$shift)
  kmers = c(kmer_group$kmer)
  for (i in kmer_group$boost){
    al = get_kmer_alignment(i, kmer_group$ekmer)
    shifts = c(shifts, al$shift)
    kmers = c(kmers, al$kmer)
  }
  for (i in kmer_group$indep){
    al = get_kmer_alignment(i, kmer_group$ekmer)
    shifts = c(shifts, al$shift)
    kmers = c(kmers, al$kmer)
  }
  out = Biostrings::consensusMatrix(Biostrings::DNAStringSet(kmers),as.prob=FALSE,shift = shifts)[1:4,]
  out = out / max(colSums(out))
  out = out + matrix((1-colSums(out))/4, byrow = TRUE, ncol = ncol(out), nrow = nrow(out))
  return(out)}

#' @export
plot_kmer_group <- function(kmer_group, motif_list = NULL, top.motifs = 3){
  l = (length(kmer_group$similar) + length(kmer_group$boost) + 5)*2
  anno_df = data.frame(x = 0, y = l, label = "seed\nkmer")
  p  = ggplot() + ggmotif(kmer_group$kmer, y.pos = l)
  l = l - 3
  a = get_kmer_alignment(kmer_group$kmer, kmer_group$ekmer)$shift
  p  = p + ggmotif(kmer_group$ekmer, y.pos = l, x.pos = 1-a)
  anno_df = rbind(anno_df, data.frame(x = 0, y = l, label = "extended\nkmer"))
  l = l - 3
#   anno_df = rbind(anno_df, data.frame(x = 0, y = l, label = "overlapping kmers\nthat boost variability"))
#   for (i in kmer_group$boost){
#     al = get_kmer_alignment(i, kmer_group$ekmer)
#     p = p + ggmotif(al$kmer, y.pos = l, x.pos = al$shift - a)
#     l = l-1.5
#   }
#   l = l - 1.5
  anno_df = rbind(anno_df, data.frame(x = 0, y = l, label = "similar\nvariable kmers"))
  for (i in kmer_group$similar){
    al = get_kmer_alignment(i, kmer_group$ekmer)
    p = p + ggmotif(al$kmer, y.pos = l, x.pos = al$shift - a)
    l = l-1.5
  }
  if (!is.null(motif_list)){
    l = l - 1.5
    anno_df = rbind(anno_df, data.frame(x = 0, y = l, label = "similar motifs"))
    motif_matches <- TFBSTools::PFMSimilarity(motif_list, kmer_group$ekmer)
    scores <- sapply(motif_matches, function(x) x["relScore"])
    for (motif in names(motif_matches)[sort(scores, decreasing = TRUE, index.return=TRUE)$ix[1:top.motifs]]){
      m = as.matrix(motif_list[[motif]]@profileMatrix)
      m = m / matrix(colSums(m), byrow = FALSE, nrow = nrow(m), ncol=ncol(m))
      al = get_kmer_alignment(Biostrings::consensusString(m,ambiguityMap=Biostrings::IUPAC_CODE_MAP, threshold = 0.25), kmer_group$ekmer)
      p = p + ggmotif(m, y.pos = l, x.pos = al$shift - a)
      l = l - 1.5
    }
  }
  minxval  = min(sapply(ggplot_build(p)$data, function(obj) min(obj$x)))
  maxxval  = max(sapply(ggplot_build(p)$data, function(obj) max(obj$x)))
  xbreaks = ceiling(minxval):floor(maxxval)
  anno_df$x = minxval - 0.5
  #p = p + geom_text(data = anno_df, mapping = aes(x = x, y = y, label = label), hjust = 1, vjust =0)

  return(p + ggmotif_scale() + ggmotif_theme() +
           scale_x_continuous(breaks = xbreaks) +
           scale_y_continuous(breaks = anno_df$y, labels = anno_df$label) +
           xlab("position relative to start of seed kmer (bp)") +
           theme(axis.line.y = element_blank(),
                 axis.text.y = element_text(face="bold",size=12),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank()))
}


get_motif_list <- function(species = 'Homo sapiens', collection = "CORE", ...){
  opts = list()
  opts['species'] = species
  opts['collection'] = collection
  opts = c(opts, list(...))
  TFBSTools::getMatrixSet(JASPAR2014::JASPAR2014, opts)
}


#PFMatrixList = get_motif_list()



get_kmer_alignment<- function(kmer, reference){
  if (is.character(kmer)) kmer = Biostrings::DNAString(kmer)
        tmp1 = Biostrings::pairwiseAlignment(reference,
                                    kmer,
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = 0,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf)
        tmp2 = Biostrings::pairwiseAlignment(reference,Biostrings::reverseComplement(kmer),
                                    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                                  mismatch = 0,
                                                                                                  baseOnly = FALSE,
                                                                                                  type = "DNA"),
                                    type="overlap",gapOpening=-Inf)
        if (Biostrings::score(tmp1) >= Biostrings::score(tmp2)){
          al_start = Biostrings::start(Biostrings::aligned(tmp1,degap=TRUE))
          return(list(kmer = as.character(kmer), shift = al_start))
        } else if (Biostrings::score(tmp1) < Biostrings::score(tmp2)){
          al_start = Biostrings::start(Biostrings::aligned(tmp2,degap=TRUE))
          return(list(kmer = as.character(Biostrings::reverseComplement(kmer)),shift = al_start))
        }
}









kmer_cluster_number <-function(kmers,
                              max_k = 15,
                              out = c("plot","k"),
                              method = "globalSEmax"){

  out = match.arg(out)

  #Perform checks on arguments
  stopifnot(as.integer(max_k) == max_k && max_k > 1)


  tmpfun <- function(x, k){
    d = kmer_dist(x)
    h = hclust(d)
    cl = cutree(h, k)
    return(list("cluster" = cl))
  }

  g = cluster::clusGap(mat, tmpfun, K.max = max_k)
  rec = cluster::maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"], method=method)

  if (out == "plot"){

    tmp = data.frame(k = 1:max_k,
                     Gap = g$Tab[,"gap"],
                     ymax =  g$Tab[,"gap"] +  g$Tab[,"SE.sim"],
                     ymin = g$Tab[,"gap"] -  g$Tab[,"SE.sim"],
                     rec = (1:max_k == rec))

    tmp2 = data.frame(label = paste0("k = ",rec), k = tmp[tmp$rec,"k"],
                      Gap = tmp[tmp$rec,"ymax"] + 0.1 *(max(tmp$ymax) - min(tmp$ymin)),
                      ymax = tmp[tmp$rec,"ymax"], ymin = tmp[tmp$rec,"ymin"], rec = TRUE,
                      stringsAsFactors = FALSE)
    print(tmp2)
    p = ggplot2::ggplot(tmp,
                        ggplot2::aes_string(x = "k", y = "Gap", ymax = "ymax",
                                            ymin = "ymin", col = "rec")) +
      ggplot2::geom_point() + ggplot2::geom_errorbar() +
      geom_text(data = tmp2,
                ggplot2::aes_string(x = "k", y = "Gap", label = "label"),col="red") +
      scale_color_manual(values = c("black","red")) +
      theme(legend.position="none")

    print(p)
    invisible(list("plot" = p, "value" = rec))
  } else{
    return(rec)
  }
}





cluster_kmers <- function(kmers){

  d = sapply(kmers, function(x) sapply(kmers, kmer_dist, x))

}



align_pwms <- function(pwm1,pwm2){
  w1 = ncol(pwm1)
  w2 = ncol(pwm2)
  cc = sapply(1:(w1+w2-1), function(x) sum(pwm1[,max(w1-x+1,1):min(w1,w1+w2-x)]*pwm2[,max(1,x-w1+1):min(x,w2)]))
  pwm1_rev = pwm1[c(4,3,2,1),w1:1]
  cc_rev = sapply(1:(w1+w2-1), function(x) sum(pwm1_rev[,max(w1-x+1,1):min(w1,w1+w2-x)]*pwm2[,max(1,x-w1+1):min(x,w2)]))
  if (max(cc) >= max(cc_rev)){
    wm = which(cc == max(cc))
    pos = wm[length(wm)] - w1
    out = c(max(cc),pos,1)
  }else{
    wm = which(cc_rev == max(cc_rev))
    pos = wm[length(wm)] - w1
    out = c(max(cc_rev),pos,-1)
  }
  return(out)
}

align_cluster <- function(kmers){
  ##kmers should be sorted based on variability
  k = length(kmers[1])
  kmers = Biostrings::DNAStringSet(kmers)
  kmer_pwms = lapply(1:length(kmers), function(x) Biostrings::consensusMatrix(kmers[x])[1:4,])
  pwm = kmer_pwms[[1]]
  maxiter = 20
  while (i < maxiter){
    alignment = sapply(kmer_pwms, align_pwms, pwm)
    keep = which(alignment[1,] >= sum(apply(pwm,2,max))*(k-1)/ncol(pwm) )
    #if (length(keep) == 0) break
    pwm = consensusMatrix(c(kmers[keep][which(alignment[3,keep]==1)],Biostrings::reverseComplement(kmers[keep][which(alignment[3,keep]==-1)])),
                          shift = alignment[2,keep])[1:4,]
    kmers = kmers[-keep]

  }

}







align_kmers <- function(kmers){




}



# Gapped k-mers ----------------------------------------------------------------

#' @export
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

#Don't export
remove_rc <- function(seqs){
  temp <- cbind(as.character(seqs),
                as.character(Biostrings::reverseComplement(seqs)))
  temp[temp[,1] > temp[,2],1] <- temp[temp[,1] > temp[,2],2]
  Biostrings::DNAStringSet(unique(temp[,1]))
}

