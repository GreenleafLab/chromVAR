

mean_smooth <- function(X, window){
  if (is.na(as.integer(window)) || length(window) != 1 || window < 2 || window >= length(X)){
    stop("window must be an integer between 2 and length(X)")
  }
  pad_left = rev(X[1:(window %/% 2)])
  pad_right = rev(X[(length(X)-((window-1) %/% 2)):length(X)])
  cx <- c(0,cumsum(c(pad_left,X,pad_right)))
  return((cx[(window+1):(length(cx)-1)] - cx[1:(length(cx)-1-window)])/window)
}



calculateTFVariance <- function(counts_mat, tf_assign, bias, niterations = 30, window_size = 2500, cores = 10){
  
  niterations = niterations + 1
  
  n_tf = ncol(tf_assign)
  n_peaks = nrow(counts_mat)
  n_cell = ncol(counts_mat)
  
  counts = rowSums(counts_mat)
  reads_per_cell = colSums(counts_mat)
  total_reads = sum(counts_mat)
  
  #Get sorted counts
  
  counts_sort = sort(counts, index.return=T)$ix 
  counts_sort_rev = sort(counts_sort, index.return=T)$ix
  
  #Get sorted bias
  bias_sort = sort(bias, index.return=T)$ix
  bias_sort_rev = sort(bias_sort, index.return=T)$ix
  
  
  get_dev_tf <- function(tf_index){
    
    tf_vec = tf_assign[,tf_index]
    
    tf_count = sum(tf_vec)
    
    sampling_vec = mean_smooth(tf_vec[counts_sort],window_size)
    sampling_vec = sampling_vec[counts_sort_rev] / sum(sampling_vec)
    
    sampled_peaks = sparseMatrix(i= rep(1:niterations, each=tf_count),
                                 j = sample(1:n_peaks, size = tf_count * niterations, prob = sampling_vec, replace = TRUE),
                                 dims = c(niterations, n_peaks), 
                                 x = 1)
    
    observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
    expected =  sum(observed)*reads_per_cell/total_reads
    
    deviation = (observed - expected)**2
    
    bias_vec = mean_smooth(tf_vec[bias_sort],window_size)
    bias_vec = bias_vec[bias_sort_rev] * tf_count / sum(bias_vec)
    
    bias_peak_vec = mean_smooth(bias_vec[counts_sort],window_size)
    bias_peak_vec = bias_peak_vec[counts_sort_rev] * tf_count / sum(bias_peak_vec)
    
    bias_term = bias_vec %*% counts_mat
    
    corr_term = bias_peak_vec %*% counts_mat
    
    correction = (bias_term - corr_term) * sum(observed)/sum(bias_term)
    
    sampled_counts = t(sampled_peaks %*% counts_mat) + matrix(correction,byrow=F,nrow=n_cell,ncol=niterations)
    expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
    
    sampled_deviations = (sampled_counts - expected_sampled_counts)**2
    
    mean_sampled_deviations = rowMeans(sampled_deviations[,1:(niterations-1)])
    sd_sampled_deviataions = apply(sampled_deviations[,1:(niterations-1)],1, sd)
    
    extra_deviation = sampled_deviations[,niterations]
    
    normvar = sqrt(sum(deviation)/sum(mean_sampled_deviations))
    sd_var = sd(sqrt(sum(deviation)/colSums(sampled_deviations)))
    extravar = sqrt(sum(extra_deviation)/sum(mean_sampled_deviations))
    sd_extra = sd(sqrt(sum(extra_deviation)/colSums(sampled_deviations)))
    
    return(c(normvar,sd_var,extravar,sd_extra))
    
  }
  
  results = simplify2array(mclapply(1:n_tf, get_dev_tf, mc.cores = cores))
  
  results = cbind(colnames(tf_assign),as.data.frame(t(results)))
  colnames(results) = c("tfs","normvar","sd_var","extravar","sd_extra")
  
  return(results)
  
}


plotVar <- function(result){
  
  res_df = data.frame(var = c(result$normvar, result$extravar), 
                      min = c(result$normvar - result$sd_var, result$extravar - result$sd_extra),
                      max = c(result$normvar + result$sd_var, result$extravar + result$sd_extra),
                      tf = rep(result$tfs,2),
                      type = c(rep("Real",length(result$normvar)),rep("Permuted",length(result$extravar))),
                      ranks = c(rank(-1 * result$normvar),rank(-1 * result$extravar)))  
  
  ggplot(res_df, aes(x = ranks, y = var, color= type, ymin = min, ymax = max)) + geom_point() + geom_errorbar() + 
   xlab("Sorted TFs") + ylab("Variability") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))
  
  
}





