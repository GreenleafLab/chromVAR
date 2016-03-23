
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
  

  sd_deviations <- apply(deviations, 1, sd, na.rm = na.rm)
  
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
    
    bootstrap_indexes <- sample(ncol(deviations), 
                                ncol(deviations)*bootstrap_samples,
                                replace=TRUE)
    
    sd_error <- apply(deviations, 1, function(x){
      quantile_helper(sapply(1:bootstrap_samples, 
                      function(y) 
                        sd(x[bootstrap_indexes[(1 + (y-1)*ncol(deviations)):
                                                 (y*ncol(deviations))]], na.rm = na.rm)),
               bootstrap_quantiles, na.rm = na.rm)
    })
    
    
    
    
    out <- data.frame(variability = sd_deviations, 
                      bootstrap_lower_bound = sd_error[1,], 
                      bootstrap_upper_bound = sd_error[2,],
                      p_value = p_sd,
                      p_value_adj = p_adj)    
  } else{
    out <- data.frame(variability = sd_deviations,
                      p_value = p_sd,
                      p_value_adj = p_adj)    
  }
  return(out)
}



