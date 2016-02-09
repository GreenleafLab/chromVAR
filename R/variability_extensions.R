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
get_peak_dev_assoc <- function(counts_mat, dev, BPPARAM = BiocParallel::bpparam(), bg.reps = 100, progress.bar = TRUE){
  
  #break into groups of 500
  grpsize = 500
  grps <- lapply(1:(counts_mat@npeak %/% grpsize + ((counts_mat@npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,counts_mat@npeak)))
  
  helperfun <- function(grp){
    densemat = as.matrix(counts_mat@counts[grp,])
    cor_res = sapply(1:length(grp), 
                     function(x) pcaPP::cor.fk(densemat[x,],
                                                 dev))
    return(cor_res)
  }
  
  helperfun2 <- function(grp){
    densemat = as.matrix(counts_mat@counts[grp,])
    cor_res = sapply(1:length(grp), 
                     function(x) pcaPP::cor.fk(densemat[x,],
                                               dev[sample.int(length(dev))]))
    return(cor_res)
  }
  
  
  out <- do.call(c,BiocParallel::bplapply(grps, helperfun,  BPPARAM = BPPARAM))
   
  countab = tabulate(counts_mat@fragments_per_peak)
  sig = rep(1, length(out))
  if (progress.bar) pb = txtProgressBar()
  for (i in 1:length(countab)){
    if (countab[i] > 0){
      if (progress.bar && (i %% 10 ==0)) setTxtProgressBar(pb, i/(length(countab)/10))
      ix = which(counts_mat@fragments_per_peak == i)
      ix.sub = ix[sample.int(length(ix), size = bg.reps, replace = TRUE)]
      grps <- lapply(1:(bg.reps %/% grpsize + ((bg.reps %% grpsize)!=0)), function(x) ix.sub[((x-1)*grpsize +1):(min(x*grpsize,bg.reps))])
      bg <- do.call(c,BiocParallel::bplapply(grps, helperfun2,BPPARAM = BPPARAM))
      sig[ix] = sapply(out[ix], function(x) (x - mean(bg))/sd(bg))
    }
  }
  if (progress.bar) close(pb)
  sig = pnorm(sig, lower.tail = FALSE)
  sig[sig > 0.5] = 1 - sig[sig > 0.5]
  out <- data.frame(cor = out, pval = sig)
  
  return(out)
}



#' @export
get_variability_synergy <- function(sets, results, counts_mat, niterations = 50, BPPARAM = BiocParallel::bpparam(), nbg = 50, force = FALSE){
  stopifnot((length(sets) < 100) || (length(results) < 100) || force)
  #
  if (length(sets) != length(results) || !all_true(all.equal(sort(names(sets)),sort(names(results))))){
    warning("Names of sets and results not the same... using intersection")
    common_names <- intersect(names(sets),names(results))
    sets <- sets[common_names]
    results <- results[common_names]
  } 
  p <- get_pvalues(results, adjust = TRUE)
  #order by pvalues...
  setnames <- names(p)[sort(p, index.return=T)$ix]
  stopifnot(length(setnames)>1)
  l = length(setnames)
  outmat = matrix(nrow=l,ncol=l)
  for (i in 1:(l-1)){
    set <- sets[[setnames[i]]]
    tmpsets <- sets[setnames[(i+1):l]]
    outmat[i,(i+1):l] <- get_variability_boost(set, tmpsets, counts_mat, 
                                               niterations, BPPARAM, nbg)
  }
  outmat
}


#' @export
get_variability_boost <- function(set, sets, counts_mat, niterations= 50, BPPARAM = BiocParallel::bpparam(), nbg = 50){
  tmpsets <- remove_nonoverlap(sets, set)
  tmpresults <- compute_variability(tmpsets, counts_mat, niterations = niterations,
                                    BPPARAM = BPPARAM)
  tmpvar <- variability(tmpresults)
  setlen = sapply(tmpsets,length)
  
  bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(set, x, replace=FALSE))), recursive = F)
  
  bgresults <- compute_variability(bgsets, counts_mat, niterations = niterations,
                                   BPPARAM = BPPARAM)
  bgvar <- variability(bgresults)
  
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/ 
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
  #var_boost <- pnorm(var_boost, lower.tail = FALSE)
  return(var_boost)
}


#' @export
get_independent_variability <- function(to.remove, 
                                        sets, 
                                        counts_mat,
                                        nresult = length(sets),
                                        niterations = 50,
                                        BPPARAM = BiocParallel::bpparam()){
  tmpsets <- remove_overlap(sets, to.remove)
  tmpresults <- compute_variability(tmpsets, counts_mat, niterations = niterations,
                                    BPPARAM = BPPARAM)
  #tmpresults <- set_nresult(tmpresults, nresult)
  return(tmpresults)
}


#' @export
get_top_sets <- function(results, sets, counts_mat,
                         niterations = 50,
                         p_cutoff = 0.01,
                         max_iter = 25,
                         BPPARAM = BiocParallel::bpparam()){
  
  #Get only significant sets
  results <- subset_by_variability(results, cutoff = p_cutoff, adjusted = TRUE)
  tmpsets <- sets[names(results)]
  p <- get_pvalues(results, adjust = TRUE)
  
  #get max variable
  max_var <- which(variability(results) == max(variability(results)))
  max_var_name <- names(results)[max_var]
  min_p <- p[max_var]
  
  #initialize output
  out <- results[max_var]
  
  #Iterate...
  iter = 2
  candidates <- names(tmpsets)[names(tmpsets) != max_var_name]
  while( length(candidates) > 0 && iter < max_iter){
    cors <- get_peak_dev_assoc(counts_mat = fc[sets[[max_var_name]],], dev = deviations(results[[max_var_name]]), BPPARAM = BPPARAM, progress.bar = FALSE)
    to.remove <- sets[[max_var_name]][intersect(which(cors$cor > 0), which(cors$pval < p_cutoff))]
    tmpsets <- remove_overlap(tmpsets[candidates], to.remove)
    tmpresults <- compute_variability(tmpsets, counts_mat, niterations = niterations,
                                      BPPARAM = BPPARAM)
    tmpresults <- set_nresult(tmpresults, results@nresult)
    tmpresults <- subset_by_variability(tmpresults, cutoff = p_cutoff, adjusted = TRUE)
    if (length(tmpresults) ==0){
      break
    }
    # get most variable...
    max_var <- which(variability(tmpresults) == max(variability(tmpresults)))
    max_var_name <-  names(tmpresults)[max_var]
    min_p <- min(get_pvalues(tmpresults, adjust = TRUE))
    
    out <- c(out, tmpresults[max_var])
    candidates <- names(tmpresults)[-max_var]
    iter <- iter + 1
  }
  
  out <- set_nresult(out, results@nresult)
  
  return(out)
}



