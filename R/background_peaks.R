#' get_gc
#' 
#' Computes GC content for peaks
#' @param peaks GenomicRanges
#' @param genome BSgenome object
#' @export
get_gc <- function(peaks, 
                   genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  seqs = Biostrings::getSeq(genome, peaks)
  nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
  gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
  return(gc)
}

#' get_background_peaks
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat counts matrix
#' @param bias vector with bias values for peaks, typically GC content
#' @param niterations number of background peaks to sample
#' @param window window size around peak in which to sample a background peak
#' @param with_replacement sample peaks with replacement? Default is TRUE
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by finding the window nearest neighbors to a peak in terms of 
#' GC content and # of fragments across samples using the Mahalanobis distance.  From those nearest
#' neighbors, niterations peaks are sampled from that background.
#' @export
get_background_peaks <- function(counts_mat, bias, niterations = 50, window = 500, with_replacement = TRUE, count = TRUE){

  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(countsum$fragments_per_peak + 1 ), bias)
  } else{
    norm_mat <- cbind(countsum$fragments_per_peak, bias)
  }
  reflected <- add_reflections(norm_mat, window = window)
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(reflected$data)))
  
  grpsize <- 10000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
  bghelper <- function(grp, reflected, tmp_vals, with_replacement, niterations){
    in1 = tmp_vals
    in2 = tmp_vals[grp,]
    tmp_nns <- FNN::get.knnx(in1, query = in2, k = window)$nn.index
    if (niterations == 1){
      return(matrix(apply(tmp_nns, 1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement))),
                    ncol = 1))
    } else{
      return(t(apply(tmp_nns,1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement)))))
    }
  }
  
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, reflected, tmp_vals, with_replacement, niterations))

  return(background_peaks)
}


# helper function, not exported
add_reflections <- function(x, window = 1000, as_ranks = FALSE){
  if (as_ranks){
    r = x
  } else{
    r <- apply(x,2, rank, ties.method="random")
  }
  half_window = window %/% 2
  n = nrow(x)
  mapping = c()
  out = x
  for (i in 1:ncol(x)){
    low = which(r[,i] <= half_window)
    high = which(r[,i] > n - half_window)
    low_add = x[low,]
    low_add[,i] = min(x[,i]) - (x[low,i] - min(x[,i]))
    high_add = x[high,]
    high_add[,i] = max(x[,i]) + (max(x[,i]) - x[high,i])
    out = rbind(out, low_add, high_add)
    mapping = c(mapping, low, high)
  }
  replace_reflected <- function(y){
    tmp <- which(y %in% (n+1):(n + half_window * ncol(x) *2))
    y <- replace(y, tmp, mapping[y[tmp]-n])
    return(y)
  }
  identify_same <- function(y){
    c(y, which(mapping == y) + n)
  }
  return(list("data" = out, "replace" = replace_reflected, "identify" = identify_same))
}

get_background_peaks_slow <- function(counts_mat, bias, niterations = 50, w = 0.1){
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")

  intensity = log10(countsum$fragments_per_peak)
  norm_mat = matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
  inv_cov_mat <- solve(cov(norm_mat))
  
  density <- do.call(c, BiocParallel::bplapply(1:nrow(counts_mat), maha_density, intensity, bias, inv_cov_mat, w))
  background_peaks <- do.call(rbind,BiocParallel::bplapply(1:nrow(counts_mat), bg_sample, intensity, bias, inv_cov_mat, density, niterations, w))
  
  return(background_peaks)
}


#' get_background_peaks_new
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat counts matrix
#' @param bias vector with bias values for peaks, typically GC content
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in  
#' GC content and # of fragments across samples using the Mahalanobis distance. 
#' The w paramter controls how similar background peaks should be.  The bs parameter
#' controls the precision with which the similarity is computed; increasing bs will
#' make the function run slower. Sensible default parameters are chosen for both.
#' @export
get_background_peaks_new <- function(counts_mat, bias, niterations = 50, w = 0.1, bs = 50){
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  
  intensity = log10(countsum$fragments_per_peak)
  norm_mat = matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
  
  chol_cov_mat <- chol(cov(norm_mat))
  trans_norm_mat <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))  
  
  #make bins
  bins1 = seq(min(trans_norm_mat[,1]), max(trans_norm_mat[,1]), length.out = bs)
  bins2 = seq(min(trans_norm_mat[,2]), max(trans_norm_mat[,2]), length.out = bs)
  
  bin_data = do.call(rbind,lapply(1:50, function(x) matrix(c(rep(bins1[x], bs),bins2), ncol = 2, byrow = FALSE)))
  
  bin_dist = euc_dist(bin_data)
  bin_p = dnorm(bin_dist, 0, w)
  
  bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx
  
  bin_density <- tabulate2(bin_membership, min_val = 1, max_val = bs**2)
  
  background_peaks <- bg_sample_helper(bin_membership-1, bin_p, bin_density, niterations)
  
  return(background_peaks)
}


get_background_peaks2 <- function(counts_mat, bias, niterations = 50, window = 500, with_replacement = TRUE, count = TRUE){
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(countsum$fragments_per_peak), bias)
  } else{
    norm_mat <- cbind(countsum$fragments_per_peak, bias)
  }
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
  bghelper <- function(grp, w, tmp_vals, with_replacement, niterations){
    in2 = tmp_vals[grp,]
    tmp_nns <- w$query(in2, k = window + 1, eps = 0)$nn.idx
    if (niterations == 1){
      return(matrix(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement)),
                    ncol = 1))
    } else{
      return(t(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement))))
    }
  }
  
  knn_kd_tree = nabor::WKNND(tmp_vals)
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, knn_kd_tree, tmp_vals, with_replacement, niterations))
  
  return(background_peaks)
}

get_background_peaks3 <- function(counts_mat, bias, niterations = 50, window = 1000, s = 0.1, with_replacement = TRUE, count = TRUE){
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(countsum$fragments_per_peak), bias)
  } else{
    norm_mat <- cbind(countsum$fragments_per_peak, bias)
  }
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
  densities <- get_2d_density(tmp_vals[,1], tmp_vals[,2], n = 50, s = 0.1)
  
  bghelper <- function(grp, w, tmp_vals, with_replacement, niterations, densities, s){
    in2 = tmp_vals[grp,]
    tmp_nns <- w$query(in2, k = window, eps = 0)
    if (niterations == 1){
      return(matrix(sapply(1:length(grp), function(x) sample(tmp_nns$nn.idx[x,], niterations, replace = with_replacement, 
                                                               prob = dnorm(tmp_nns$nn.dists[x,],0,s) / densities[tmp_nns$nn.idx[x,]] )),
                    ncol = 1))
    } else{
      return(t(sapply(1:length(grp), function(x) sample(tmp_nns$nn.idx[x,], niterations, replace = with_replacement, 
                                                          prob = dnorm(tmp_nns$nn.dists[x,],0,s) / densities[tmp_nns$nn.idx[x,]] ))))
    }
  }
  
  knn_kd_tree = nabor::WKNND(tmp_vals)
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, knn_kd_tree, tmp_vals, with_replacement, niterations, densities, s))
  
  return(background_peaks)
}
