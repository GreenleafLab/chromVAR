cor_helper <- function(a, b, ...){
  res = cor.test(a, b,...)
  return(c(r = res$estimate[[1]],p = res$p.value))
}

get_correlated_results <- function(name, results, p.cutoff = 0.01, corr.cutoff = 0.5){
  mat = deviations(results)
  ix = which(names(results) == name)
  pcors = apply(mat, 1, function(x) cor_helper(x, mat[ix,],alternative = "greater"))
  sig = intersect(which(pcors["p",] < 0.01), which(pcors["r",] > corr.cutoff))
  return(names(results)[sig[sig != ix]])
}


#' get_variable_peaks
#' 
#' function to find individual peaks with variability correlated with the variability 
#' of a set of peaks
#' @param counts_mat \code{\link{fragmentCounts}} object
#' @param deviation_result  \code{\link{deviationResult}} object
#' 
#' @export
get_variable_peaks <- function(counts_mat, deviation_result){
  
  dev = deviations(deviation_result)
  
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
                                               sample(dev)))
    return(cor_res)
  }
  
  out <- do.call(c,BiocParallel::bplapply(grps, helperfun))
   
  countab = tabulate(counts_mat@fragments_per_peak)
  sig = rep(1, length(out))
  for (i in 1:length(countab)){
    if (countab[i] > 0){
      print(i)
      ix = which(counts_mat@fragments_per_peak == i)
      ix.sub = ix[sample.int(length(ix), size = length(out), replace = TRUE)]
      grps <- lapply(1:(bg.reps %/% grpsize + ((bg.reps %% grpsize)!=0)), function(x) ix.sub[((x-1)*grpsize +1):(min(x*grpsize,bg.reps))])
      bg <- do.call(c,BiocParallel::bplapply(grps, helperfun2))
      sig[ix] = sapply(out[ix], function(x) ifelse(mean(bg) > x, mean(bg < x),mean(bg>x)))
    }
  }
  #sig = pnorm(sig, lower.tail = FALSE)
  #sig[sig > 0.5] = 1 - sig[sig > 0.5]
  out <- data.frame(cor = out, pval = sig)
  
  return(out)
}


#' @export
get_variability_synergy <- function(sets, results, counts_mat, niterations = 50, nbg = 50, force = FALSE){
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
                                               niterations, nbg)
  }
  outmat
}


#' @export
get_variability_boost <- function(set, sets, counts_mat, niterations= 50, nbg = 50){
  tmpsets <- remove_nonoverlap(sets, set)
  tmpresults <- compute_variability(c(list(set),tmpsets), counts_mat, niterations = niterations)
  tmpvar <- variability(tmpresults)[2:length(tmpresults)]
  basevar <- variability(tmpresults)[1]
  out <- match.arg(out)
  
  setlen = sapply(tmpsets,length)
  
  bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(set, x, replace=FALSE))), recursive = F)
  
  bgresults <- compute_variability(bgsets, counts_mat, niterations = niterations)
  bgvar <- variability(bgresults)
  
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/ 
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
  return(var_boost)
}


get_independent_variability <- function(to.remove, 
                                        sets, 
                                        counts_mat,
                                        nresult = length(sets),
                                        niterations = 50){
  tmpsets <- remove_overlap(sets, to.remove)
  tmpresults <- compute_variability(tmpsets, counts_mat, niterations = niterations)
  tmpresults <- set_nresult(tmpresults, nresult)
  return(tmpresults)
}


#' @export
get_top_sets <- function(results,
                          sets, 
                          counts_mat,
                         niterations = 50,
                         p_cutoff = 0.01,
                         max_iter = 25){
  
  #Get only significant sets
  results <- subset_by_variability(results, cutoff = p_cutoff, adjusted = TRUE)
  p <- get_pvalues(results, adjust = TRUE)
  
  #get max variable
  max_var <- which(variability(results) == max(variability(results)))
  max_var_name <- names(results)[max_var]
  
  #initialize output
  out <- c(max_var_name)
  
  #Iterate...
  iter = 2
  candidates <- names(tmpsets)[names(tmpsets) != max_var_name]
  while( length(candidates) > 0 && iter < max_iter){
    tmpsets <- remove_overlap(sets[candidates], sets[[max_var_name]])
    tmpresults <- compute_variability(tmpsets, counts_mat, niterations = niterations)
    tmpresults <- set_nresult(tmpresults, results@nresult)
    tmpresults <- subset_by_variability(tmpresults, cutoff = p_cutoff, adjusted = TRUE)
    if (length(tmpresults) ==0){
      break
    }
    # get most variable...
    max_var <- which(variability(tmpresults) == max(variability(tmpresults)))
    max_var_name <-  names(tmpresults)[max_var]
    
    out <- c(out, max_var_name)
    if (length(tmpresults > 1)){
      candidates <- names(tmpresults)[-max_var]
      iter <- iter + 1
    } else{
      break
    }
  }
  
  return(out)
}


