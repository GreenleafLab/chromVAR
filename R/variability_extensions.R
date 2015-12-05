
remove_overlap <- function(indices_list, index){
  if (is.character(index)){
    index = which(names(indices_list) == index)
  }
  toremove = indices_list[[index]]
  indices_list = indices_list[-index]
  indices_list = lapply(indices_list, function(x) x[which(x %ni% toremove)])
  return(indices_list)
} 

remove_nonoverlap <- function(indices_list, index){
  if (is.character(index)){
    index = which(names(indices_list) == index)
  }
  tokeep = indices_list[[index]]
  indices_list = indices_list[-index]
  indices_list = lapply(indices_list, function(x) x[which(x %in% tokeep)])
  return(indices_list)
} 

#' @export
get_variability_boost <- function(set, tmpsets, counts_mat, bg_peaks, niterations, BPPARAM, nbg = 50){
  tmpresults <- compute_variability(tmpsets, counts_mat, bg_peaks, niterations = niterations,
                                    BPPARAM = BPPARAM)
  tmpvar <- variability(tmpresults)
  setlen = sapply(tmpsets,length)
  
  bgsets <- unlist(lapply(setlen, function(x) lapply(1:nbg, function(y) sample(set, x, replace=FALSE))), recursive = F)
  
  bgresults <- compute_variability(bgsets, counts_mat, bg_peaks, niterations = niterations,
                                   BPPARAM = BPPARAM)
  bgvar <- variability(bgresults)
  
  var_boost <- sapply(seq_along(tmpvar), 
                      function(x) (tmpvar[x] - mean(bgvar[((x-1)*nbg+1):(x*nbg)]))/ 
                        sd(bgvar[((x-1)*nbg+1):(x*nbg)])  )
  var_boost <- pnorm(var_boost, lower.tail = FALSE)
  return(var_boost)
}


#' @export
get_independent_variability <- function(name, 
                                        sets, 
                                        counts_mat, 
                                        bg_peaks, 
                                        nresult = length(sets),
                                        niterations = 50,
                                        BPPARAM = BiocParallel::bpparam()){
  tmpsets <- remove_overlap(sets, name)
  tmpresults <- compute_variability(tmpsets, counts_mat, bg_peaks, niterations = niterations,
                                    BPPARAM = BPPARAM)
  tmpresults <- set_nresult(tmpresults, nresult)
  return(tmpresults)
}


#' @export
get_top_sets <- function(results, sets, counts_mat, bg_peaks, 
                         niterations = 50,
                         p_cutoff = 0.01,
                         max_iter = 25,
                         BPPARAM = BiocParallel::bpparam()){
  
  #Get only significant sets
  results <- subset_by_variability(results, cutoff = p_cutoff, adjusted = TRUE)
  sets <- sets[names(results)]
  p <- get_pvalues(results, adjust = TRUE)
  
  #get max variable
  max_var <- which(variability(results) == max(variability(results)))
  max_var_name <- names(results)[max_var]
  min_p <- p[max_var]
  
  #initialize output
  out <- results[max_var]
  
  #Iterate...
  iter = 2
  candidates = names(sets)[names(sets) != max_var_name]
  while( length(candidates) > 0 && iter < max_iter){
    tmpsets <- remove_overlap(c(sets[candidates],sets[max_var_name]), max_var_name)
    tmpresults <- compute_variability(tmpsets, counts_mat, bg_peaks, niterations = niterations,
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



