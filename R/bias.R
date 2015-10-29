
compute_gc_bias <- function(seqs){
  nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
  gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
  return(gc)
}


add_bias_bins <- function(motif_mat, bias, counts, nbins = 20){
  #add bias bins
  bias_quantiles = quantile(bias, seq(0,1,1/nbins))
  bias_cut = cut(bias, breaks = bias_quantiles)
  bias_bins = split(1:length(bias), bias_cut)
  bias_bin_names = sapply(1:nbins, function(x) paste("bias_bin_",x,sep="",collapse=""))
  bias_mat = sparseMatrix(i = unlist(bias_bins), j = unlist(lapply(1:nbins, function(x) rep(x,length(bias_bins[[x]])))),
                          x = 1, dims = c(length(bias),nbins), dimnames = list(NULL,bias_bin_names))
  motif_mat = cBind(motif_mat, bias_mat)
  #make count bins
  pseudo_counts = counts + runif(length(counts),min = 0, max = 0.1)
  count_quantiles = quantile(pseudo_counts, seq(0,1,1/nbins))
  count_cut = cut(pseudo_counts, breaks = count_quantiles)
  count_bins = split(1:length(counts), count_cut)
  count_bin_names = sapply(1:nbins, function(x) paste("count_bin_",x,sep="",collapse=""))
  count_mat = sparseMatrix(i = unlist(count_bins), j = unlist(lapply(1:nbins, function(x) rep(x,length(count_bins[[x]])))),
                          x = 1, dims = c(length(bias),nbins), dimnames = list(NULL,count_bin_names))
  motif_mat = cBind(motif_mat, count_mat)
  return(motif_mat)
}
