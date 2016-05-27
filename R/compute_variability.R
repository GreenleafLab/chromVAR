
quantile_helper <- function(values, quantiles, na.rm){
  if (na.rm){
    out <- quantile(values, quantiles, na.rm = TRUE)
  } else{
    if (!all_false(is.na(values))){
      out <- rep(NA,length(quantiles))
    } else{
      out <- quantile(values, quantiles)
    }
  }
  return(out)
}

#'compute_variability
#'
#'function to compute overall variability of motif sets across samples
#'@param deviations output from \code{\link{compute_deviations}}
#'@param bootstrap_error compute bootstrap confidence interval
#'@param bootstrap_samples number of bootstrap samples to take
#'@param bootstrap_quantiles quantiles for bootstrap 
#'@param na.rm remove NAs? default is true
#'@export
compute_variability <- function(deviations, 
                                bootstrap_error = TRUE, 
                                bootstrap_samples = 1000, 
                                bootstrap_quantiles = c(0.025,0.975),
                                na.rm = TRUE){
  

  if (is.matrix(deviations)){
    sd_deviations <- row_sds(deviations, na.rm)
    
    p_sd <- pchisq((ncol(deviations)-1) * (sd_deviations**2), 
                   df = (ncol(deviations)-1), 
                   lower.tail = FALSE)
    
    p_adj <- p.adjust(p = p_sd, method = "BH")  
    
    if (bootstrap_error){
      #check some arguments
      stopifnot(bootstrap_samples > 0 && all_whole(bootstrap_samples))
      stopifnot(all_true(bootstrap_quantiles > 0))
      stopifnot(all_true(bootstrap_quantiles < 1)) 
      stopifnot(bootstrap_quantiles[2] > bootstrap_quantiles[1])
      
      boot_sd <- do.call(rbind,BiocParallel::bplapply(1:bootstrap_samples, function(x) row_sds_perm(deviations, na.rm)))
      sd_error <- apply(boot_sd, 2, quantile_helper, quantiles = bootstrap_quantiles, na.rm = na.rm)
      
      out <- data.frame(variability = sd_deviations, 
                        bootstrap_lower_bound = sd_error[1,], 
                        bootstrap_upper_bound = sd_error[2,],
                        p_value = p_sd,
                        p_value_adj = p_adj,
                        row.names = rownames(deviations))    
    } else{
      out <- data.frame(variability = sd_deviations,
                        p_value = p_sd,
                        p_value_adj = p_adj,
                        row.names = rownames(deviations))    
    }    
  } else{
    sd_deviations <- sd(deviations, na.rm = na.rm)
    
    p_sd <- pchisq((length(deviations)-1) * (sd_deviations**2), 
                   df = (length(deviations)-1), 
                   lower.tail = FALSE)
    
    p_adj <- p_sd  
    
    if (bootstrap_error){
      #check some arguments
      stopifnot(bootstrap_samples > 0 && all_whole(bootstrap_samples))
      stopifnot(all_true(bootstrap_quantiles > 0))
      stopifnot(all_true(bootstrap_quantiles < 1)) 
      stopifnot(bootstrap_quantiles[2] > bootstrap_quantiles[1])
      
      boot_sd <- replicate(bootstrap_samples, sd(deviations[sample(1:length(deviations),length(deviations), replace = TRUE)], na.rm = na.rm))
      sd_error <- quantile_helper(boot_sd, quantiles = bootstrap_quantiles, na.rm = na.rm)
      
      out <- c(variability = sd_deviations, 
                        bootstrap_lower_bound = sd_error[1], 
                        bootstrap_upper_bound = sd_error[2],
                        p_value = p_sd,
                        p_value_adj = p_adj)    
    } else{
      out <- c(variability = sd_deviations,
                        p_value = p_sd,
                        p_value_adj = p_adj)    
    }    
  }
  return(out)
}



# Helper function --------------------------------------------------------------

compute_variability_single <- function(peak_set,
                                       counts_mat,
                                       background_peaks,
                                       expectation){

  tmp <- compute_deviations_single(peak_set, counts_mat, background_peaks, expectation)
  
  v <- sd(tmp$z, na.rm = TRUE)
  
  return(v)
}
