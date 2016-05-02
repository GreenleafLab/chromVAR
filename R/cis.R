
#'@importMethodsFrom GenomeInfoDb seqlevels seqnames
#'@export
get_cis_groups <- function(peaks, grpsize = 25, stepsize = 10){
  
  if (!isTRUE(all.equal(peaks, sort(peaks)))){
    stop("peaks must be sorted to be able to get cis groups")
  }
  
  chrs = seqlevels(peaks)
  out = list()
  out = do.call(c, BiocParallel::bplapply(seq_along(chrs), function(i) {
    chr_ix = which(seqnames(peaks) == chrs[i])
    if (length(chr_ix) > stepsize){
      tmp = lapply(1:(length(chr_ix) %/% stepsize), function(x) chr_ix[((x-1)*stepsize +1):(min((x-1)*stepsize + grpsize,length(chr_ix)))])
      names(tmp) = sapply(1:length(tmp), function(x) paste(chrs[i], x,sep="_",collapse=""))
      return(tmp)
    } else{
      return(list())
    }
  }))
  return(out)
}

deviations_to_cis_cor <- function(deviations, cis_groups, peaks, start, end, steps = 10000){
  
  #first correlation
  cis_cor1= cor(t(deviationss), use = "pairwise.complete.obs")
  
  #second correlation 
  cis_cor2 = cor(cis_cor2, use = "pairwise.complete.obs")
  
  #map to base pair
  cis_loc_means <- lapply(cis_groups, mean(start(peaks[cis_groups]) + width(peaks[cis_groups])/2))
  x <- rep(cis_loc_means, ncol(cis_cor2))
  y <- rep(cis_loc_means, each = ncol(cis_cor2))
  new_locs <- seq(start,end, steps)
  loc <- matrix(nrow = length(new_locs)**2, ncol = 2)
  loc[,1] = rep(new_locs, length(new_locs))
  loc[,2] = rep(new_locs, each = length(new_locs))
  z <- as.vector(cis_cor2)
  nx <- length(x)
  ny <- length(y)
  lx <- approx(x, 1:nx, loc[, 1])$y
  ly <- approx(y, 1:ny, loc[, 2])$y
  lx1 <- floor(lx)
  ly1 <- floor(ly)
  ex <- lx - lx1
  ey <- ly - ly1
  ex[lx1 == nx] <- 1
  ey[ly1 == ny] <- 1
  lx1[lx1 == nx] <- nx - 1
  ly1[ly1 == ny] <- ny - 1
  tmp <- z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + z[cbind(lx1 + 
                                                              1, ly1)] * ex * (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 - 
                                                                                                                     ex) * ey + z[cbind(lx1 + 1, ly1 + 1)] * ex * ey
  return(matrix(tmp, nrow = length(new_locs), ncol = length(new_locs)))
}

compute_cis_correlations <- function(counts, peaks){
  
  # Get bins
  
  # Compute deviations
  
  # First correlation
  
  # Second correlation
  
  # Map to base pair
  
}



