# 
# 
# mean_smooth <- function(X, window){
#   if (is.na(as.integer(window)) || length(window) != 1 || window < 2 || window >= length(X)){
#     stop("window must be an integer between 2 and length(X)")
#   }
#   pad_left = rev(X[1:(window %/% 2)])
#   pad_right = rev(X[(length(X)-((window-1) %/% 2)):length(X)])
#   cx <- c(0,cumsum(c(pad_left,X,pad_right)))
#   return((cx[(window+1):(length(cx)-1)] - cx[1:(length(cx)-1-window)])/window)
# }
# 
# 
# 
# calculateTFVariance <- function(counts_mat, tf_assign, bias, niterations = 30, window_size = 2500, cores = 10){
# 
#   niterations = niterations + 1
# 
#   n_tf = ncol(tf_assign)
#   n_peaks = nrow(counts_mat)
#   n_cell = ncol(counts_mat)
# 
#   counts = rowSums(counts_mat)
#   reads_per_cell = colSums(counts_mat)
#   total_reads = sum(counts_mat)
# 
#   #Get sorted counts
# 
#   counts_sort = sort(counts, index.return=T)$ix
#   counts_sort_rev = sort(counts_sort, index.return=T)$ix
# 
#   #Get sorted bias
#   bias_sort = sort(bias, index.return=T)$ix
#   bias_sort_rev = sort(bias_sort, index.return=T)$ix
# 
# 
#   get_dev_tf <- function(tf_index){
# 
#     tf_vec = tf_assign[,tf_index]
# 
#     tf_count = sum(tf_vec)
# 
#     sampling_vec = mean_smooth(tf_vec[counts_sort],window_size)
#     sampling_vec = sampling_vec[counts_sort_rev] / sum(sampling_vec)
# 
#     sampled_peaks = sparseMatrix(i= rep(1:niterations, each=tf_count),
#                                  j = sample(1:n_peaks, size = tf_count * niterations, prob = sampling_vec, replace = TRUE),
#                                  dims = c(niterations, n_peaks),
#                                  x = 1)
# 
#     observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
#     expected =  sum(observed)*reads_per_cell/total_reads
# 
#     deviation = (observed - expected)**2
# 
#     bias_vec = mean_smooth(tf_vec[bias_sort],window_size)
#     bias_vec = bias_vec[bias_sort_rev] * tf_count / sum(bias_vec)
# 
#     bias_peak_vec = mean_smooth(bias_vec[counts_sort],window_size)
#     bias_peak_vec = bias_peak_vec[counts_sort_rev] * tf_count / sum(bias_peak_vec)
# 
#     bias_term = bias_vec %*% counts_mat
# 
#     corr_term = bias_peak_vec %*% counts_mat
# 
#     correction = (bias_term - corr_term) * sum(observed)/sum(bias_term)
# 
#     sampled_counts = t(sampled_peaks %*% counts_mat) + matrix(correction,byrow=F,nrow=n_cell,ncol=niterations)
#     expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
# 
#     sampled_deviations = (sampled_counts - expected_sampled_counts)**2
# 
#     mean_sampled_deviations = rowMeans(sampled_deviations[,1:(niterations-1)])
#     sd_sampled_deviataions = apply(sampled_deviations[,1:(niterations-1)],1, sd)
# 
#     extra_deviation = sampled_deviations[,niterations]
# 
#     normvar = sqrt(sum(deviation)/sum(mean_sampled_deviations))
#     sd_var = sd(sqrt(sum(deviation)/colSums(sampled_deviations)))
#     extravar = sqrt(sum(extra_deviation)/sum(mean_sampled_deviations))
#     sd_extra = sd(sqrt(sum(extra_deviation)/colSums(sampled_deviations)))
# 
#     return(c(normvar,sd_var,extravar,sd_extra))
# 
#   }
# 
#   results = simplify2array(mclapply(1:n_tf, get_dev_tf, mc.cores = cores))
# 
#   results = cbind(colnames(tf_assign),as.data.frame(t(results)))
#   colnames(results) = c("tfs","normvar","sd_var","extravar","sd_extra")
# 
#   return(results)
# 
# }
# 
# 
# 
# 
# plotVar <- function(result){
# 
#   res_df = data.frame(var = c(result$normvar, result$extravar),
#                       min = c(result$normvar - result$sd_var, result$extravar - result$sd_extra),
#                       max = c(result$normvar + result$sd_var, result$extravar + result$sd_extra),
#                       tf = rep(result$tfs,2),
#                       type = c(rep("Real",length(result$normvar)),rep("Permuted",length(result$extravar))),
#                       ranks = c(rank(-1 * result$normvar),rank(-1 * result$extravar)))
# 
#   ggplot(res_df, aes(x = ranks, y = var, color= type, ymin = min, ymax = max)) + geom_point() + geom_errorbar() +
#     xlab("Sorted TFs") + ylab("Variability") + scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$max)*1.05))
# 
# 
# }
# 
# 
# make_background_peak_set <- function(counts_mat, niterations = 30, window_size = 2500, cores = 10){
# 
#   niterations = niterations + 1
# 
#   n_peaks = nrow(counts_mat)
#   counts = rowSums(counts_mat)
# 
#   #Get sorted counts
# 
#   counts_sort = sort(counts, index.return=T)$ix
#   counts_sort_rev = sort(counts_sort, index.return=T)$ix
# 
#   get_matched_peak <- function(peak_index){
# 
#     #First define pool of peaks with similar peak intensiy
# 
#     peak_rank = counts_sort_rev[peak_index]
#     half_window = window_size %/% 2
# 
#     pool = counts_sort[max(1,peak_rank - half_window):min(n_peaks,peak_rank+half_window)]
# 
#     #Sample
# 
#     s = sample(pool, size = niterations, replace = T)
# 
#     return(s)
# 
#   }
# 
#   matched_peak_mat = simplify2array(mclapply(1:n_peaks, get_matched_peak, mc.cores = cores))
# 
#   return(matched_peak_mat)
# }
# 
# make_bias_correction_table <- function(counts_mat, bias, niterations = 30, window_size = 2500, cores = 10){
# 
#   niterations = niterations + 1
# 
#   n_peaks = nrow(counts_mat)
# 
#   counts = rowSums(counts_mat)
# 
#   #Get sorted counts
# 
#   counts_sort = sort(counts, index.return=T)$ix
#   counts_sort_rev = sort(counts_sort, index.return=T)$ix
# 
#   #Get sorted bias
#   bias_sort = sort(bias, index.return=T)$ix
#   bias_sort_rev = sort(bias_sort, index.return=T)$ix
# 
#   get_bias_correction <- function(peak_index){
# 
#     #get expected reads based on similar Tn5 bias
#     bias_rank = bias_sort_rev[peak_index]
#     half_window = window_size %/% 2
# 
#     bias_pool = bias_sort[max(1,bias_rank - half_window):min(n_peaks,bias_rank+half_window)]
#     bias_weights = rep(0, n_peaks)
#     bias_weights[bias_pool] = length(bias_pool)/(half_window * 2 + 1)
# 
#     bias_reads = bias_weights %*% counts_mat
# 
#     #get expected reads based on counts given Tn5 bias
# 
#     weights = mean_smooth(bias_weights[counts_sort], window_size)[counts_sort_rev]
# 
#     bias_peak_reads = weights %*% counts_mat
# 
#     #get correction
# 
#     correction = (bias_reads - bias_peak_reads) * counts[peak_index]/sum(bias_reads)
# 
#     return(as.vector(correction))
#   }
# 
#   correction_mat = simplify2array(mclapply(1:n_peaks, get_bias_correction, mc.cores = cores))
# 
#   return(correction_mat)
# 
# }
# 
# 
# 
# calculateTFVariance2 <- function(counts_mat, tf_indices, bias, niterations = 30, window_size = 2500, cores = 10){
# 
#   matched_peak_mat <- make_background_peak_set(counts_mat, niterations = niterations, window_size = window_size, cores = cores)
#   correction_mat <- make_bias_correction_table(counts_mat, bias, niterations = niterations, window_size = window_size, cores = cores)
# 
#   niterations = niterations + 1
# 
#   n_tf = ncol(tf_assign)
#   n_peaks = nrow(counts_mat)
#   n_cell = ncol(counts_mat)
# 
#   counts = rowSums(counts_mat)
#   reads_per_cell = colSums(counts_mat)
#   total_reads = sum(counts_mat)
# 
#   expected_per_peak = outer(counts, reads_per_cell/total_reads)
# 
#   get_tf_deviations <- function(tf_index){
# 
#     counts_correction = rowSums(correction_mat[,tf_indices[[tf_index]]])
# 
#     background_mat = sparseMatrix(i = rep(1:niterations,length(tf_indices[[tf_index]])),j = as.vector(matched_peak_mat[,tf_indices[[tf_index]]]),  dims = c(niterations, n_peaks), x = 1)
# 
#     sampled_counts = t(background_mat %*% counts_mat + matrix(counts_correction, byrow = T, ncol = n_cell, nrow = niterations))
# 
#     expected_sampled_counts = t(outer(colSums(sampled_counts), reads_per_cell/total_reads))
# 
#     #sampled_deviation = sampled_counts - expected_sampled_counts
# 
#     tf_vec = rep(0, n_peaks)
#     tf_vec[tf_indices[[tf_index]]] = 1
# 
#     observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
#     expected =  sum(observed)*reads_per_cell/total_reads
# 
#     deviation = (observed - expected)
# 
# #     get_z_loo <- function(ind){
# #       obs_loo = (observed - counts_mat[ind,]) - (expected - expected_per_peak[ind,])
# #       s_loo =  (sampled_counts - (counts_mat[matched_peak_mat[,ind],] + matrix(correction_mat[,ind], byrow=T, nrow = niterations, ncol = n_cell))) - (expected_sampled_counts  - (expected_per_peak[matched_peak_mat[,ind],]))
# #       z_loo = (obs_loo - colMeans(s_loo))/apply(s_loo, 2, sd)
# #       return(z_loo)
# #     }
# #
# #     observed_loo = (matrix(observed, nrow = sum(tf_vec), ncol = length(observed), byrow=T) - counts_mat[tf_vec,]) - (matrix(expected, nrow = sum(tf_vec), ncol = length(observed), byrow=T) - expected_per_peak[tf_vec,])
# #
# #     sampled_loo = sapply(1:niterations, function(x) as.matrix((matrix(sampled_counts[,x], nrow = sum(tf_vec), ncol = length(observed), byrow=T) - (counts_mat[matched_peak_mat[x,tf_indices[[tf_index]]],] + t(correction_mat[,tf_indices[[tf_index]]]))) -
# #                            (matrix(expected_sampled_counts[,x], nrow = sum(tf_vec), ncol = length(observed), byrow=T) - expected_per_peak[matched_peak_mat[x,tf_indices[[tf_index]]],])))
# #
# #     z_loo = sapply(1:sum(tf_vec), function(x) (observed_loo[x,]- sapply(1:niterations, ))
# #
#     result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]),
#                              extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]),
#                              observed_deviations = as.vector(deviation),
#                              tf = names(tf_indices)[tf_index])
# 
#     result = compute_z_score(result)
#     result = compute_variability(result)
# 
# 
#     return(result)
# 
#   }
# 
#   sampled_deviations = simplify2array(mclapply(1:n_tf,get_tf_background, mc.cores = cores))
# 
#   observed_per_tf = sparseMatrix(i = unlist(lapply(1:n_tf,function(x) rep(x,length(tf_indices[[x]])))), j = unlist(tf_indices), x= 1, dims=c(n_tf,n_peaks)) %*% counts_mat
#   expected_per_tf = outer(rowSums(observed_per_tf),reads_per_cell/total_reads)
#   deviation = observed_per_tf - expected_per_tf
# 
#   #Compute Z-score Based Metrics
# 
#   z_scores = (deviation - t(apply(sampled_deviations[1:(niterations-1),,],2:3,mean))) / t(apply(sampled_deviations[1:(niterations-1),,],2:3,sd))
#   sd_z = apply(z_scores,1, sd)
#   z_sig = sapply((n_cell-1)*sd_z**2,pchisq,n_cell-1,lower.tail=F)
#   z_sig_adj = p.adjust(z_sig,method="BH")
# 
#   #Compute same metrics for extra iteration
# 
#   sim_z_scores = (t(sampled_deviations[niterations,,]) - t(apply(sampled_deviations[1:(niterations-1),,],2:3,mean))) / t(apply(sampled_deviations[1:(niterations-1),,],2:3,sd))
#   sim_sd_z = apply(sim_z_scores,1, sd)
#   sim_z_sig = sapply((n_cell-1)*sim_sd_z**2,pchisq,n_cell-1,lower.tail=F)
#   sim_z_sig_adj = p.adjust(sim_z_sig,method="BH")
# 
#   #Compute Metrics from original paper
# 
#   norm_deviation = deviation/sqrt(t(apply(sampled_deviations[1:(niterations-1),,]**2,2:3,mean)))
# 
#   norm_variability = sqrt(apply(deviation**2,1,sum) / apply(t(apply(sampled_deviations[1:(niterations-1),,]**2,2:3,mean)),1,sum))
#   error_norm_variability = apply(sapply(1:(niterations-1), function(x) sqrt(apply(deviation**2,1,sum) / apply(sampled_deviations[x,,]**2,2,sum))),1,sd)
# 
# 
# }
# 
# 
# 
# calculate_deviations <- function(counts_mat, motif_mat, bias, cutoff = 0, cutoff.type = 'lower',
#                                  niterations = 100, window_size = 2500, BPPARAM=bpparam()){
# 
#   motif_mat[motif_mat > cutoff] =1
# 
#   niterations = niterations + 1
# 
#   n_tf = ncol(motif_mat)
#   n_peaks = nrow(counts_mat)
#   n_cell = ncol(counts_mat)
# 
#   counts = rowSums(counts_mat)
#   reads_per_cell = colSums(counts_mat)
#   total_reads = sum(counts_mat)
# 
#   #Get sorted counts
# 
#   counts_sort = sort(counts, index.return=T)$ix
#   counts_sort_rev = sort(counts_sort, index.return=T)$ix
# 
#   #Get sorted bias
#   bias_sort = sort(bias, index.return=T)$ix
#   bias_sort_rev = sort(bias_sort, index.return=T)$ix
# 
# 
#   get_dev_tf <- function(tf_index){
# 
#     tf_vec = motif_mat[,tf_index]
#     tf_count = sum(tf_vec)
# 
#     sampling_vec = mean_smooth(tf_vec[counts_sort],window_size)
#     sampling_vec = sampling_vec[counts_sort_rev] / sum(sampling_vec)
# 
#     sampled_peaks = sparseMatrix(i= rep(1:niterations, each=tf_count),
#                                  j = sample(1:n_peaks, size = tf_count * niterations, prob = sampling_vec, replace = TRUE),
#                                  dims = c(niterations, n_peaks),
#                                  x = 1)
# 
#     observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
#     expected =  sum(observed)*reads_per_cell/total_reads
# 
#     deviation = (observed - expected)
# 
#     bias_vec = mean_smooth(tf_vec[bias_sort],window_size)
#     bias_vec = bias_vec[bias_sort_rev] * tf_count / sum(bias_vec)
# 
#     bias_peak_vec = mean_smooth(bias_vec[counts_sort],window_size)
#     bias_peak_vec = bias_peak_vec[counts_sort_rev] * tf_count / sum(bias_peak_vec)
# 
#     bias_term = bias_vec %*% counts_mat
# 
#     corr_term = bias_peak_vec %*% counts_mat
# 
#     correction = (bias_term - corr_term) * sum(observed)/sum(bias_term)
# 
#     sampled_counts = t(sampled_peaks %*% counts_mat) + matrix(correction,byrow=F,nrow=n_cell,ncol=niterations)
#     expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
# 
# 
#     result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]),
#                              extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]),
#                              observed_deviations = as.vector(deviation),
#                              tf = names(tf_indices)[tf_index])
# 
#     result = compute_z_score(result)
#     result = compute_variability(result)
# 
#     return(result)
#   }
# 
# 
#   results = deviationResultSet(results = bplapply(1:n_tf, get_dev_tf, BPPARAM = BPPARAM), tfs = names(tf_indices))
#   results = adjust_p_values(results)
# 
#   return(results)
# 
# }
# 
# 
# toFragmentEnds <- function(x, atac.offsets = FALSE) {
#   ##Get Left and Right
#   x_left <- left(x)
#   strand(x_left) <- "+"
#   x_right <- right(x)
#   strand(x_right) <- "-"
#   ##Convert to GRangesList
#   grl <- GenomicAlignments::grglist(c(x_left, x_right)[S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))],
#                                     order.as.in.query=TRUE,
#                                     drop.D.ranges=FALSE)
#   out <- GenomicAlignments:::shrinkByHalf(grl)
#   names(out) <- names(x)
#   ##Change sizes
#   out <- resize(out, width = 1)
#   ##ATAC-seq correction
#   if (atac.offsets){
#     out <- resize(promoters(out,upstream = 0, downstream = 5), width=1, fix="end")
#   }
#   return(out)
# }
# 
# 
# getRG <- function(bam, n = 100000){
# 
#   bamfile = Rsamtools::BamFile(bam, yieldSize = n)
# 
#   open(bamfile)
#   read_groups <- unique(scanBam(bamfile, param = Rsamtools::ScanBamParam(tag="RG"))[[1]]$tag$RG)
#   close(bamfile)
# 
#   return(read_groups)
# }
# 
# 
# getFragmentCountsByRG <- function(bam, bed, BPPARAM = bpparam()){
# 
#   bamfile = Rsamtools::BamFile(bam, asMates = TRUE)
# 
#   #split by RG tags
#   RG_tags = unique(scanBam(bamfile, param = Rsamtools::ScanBamParam(tag="RG"))[[1]]$tag$RG)
# 
#   paired = GenomicAlignments::readGAlignmentPairs(bamfile, param = Rsamtools::ScanBamParam(tag="RG"))
#   fe = toFragmentEnds(paired, atac = F)
# 
#   o = as.data.frame(GenomicRanges::findOverlaps(bed, fe, type= "any", ignore.strand = T))
# 
#   o = cbind(o, data.frame("RG" = mcols(left(paired))$RG[o$subjectHits], stringsAsFactors = F))
# 
#   get_rg_peaks <- function(RG){
# 
#     o$queryHits[which(o$RG == RG)]
# 
#   }
# 
#   rg_peaks <- lapply(RG_tags, get_rg_peaks)
# 
#   n = length(RG_tags)
#   out <- sparseMatrix( j = unlist(lapply(1:n, function(x) rep(x, length(rg_peaks[[x]])))),
#                        i = unlist(lapply(1:n, function(x) rg_peaks[[x]])),
#                        x = 1,
#                        dims = c(length(bed),n),
#                        dimnames = list(NULL, RG_tags))
#   return(out)
# }
# 
# 
# 
# getBackgroundPeakSets <- function(peak_counts, bias, niterations = 51, window = 2500, BPPARAM = BiocParallel::bpparam()){
#   #Standardise peak counts and bias
#   peak_counts_norm = (log(peak_counts) - mean(log(peak_counts)))/sd(log(peak_counts))
#   bias_norm = (bias - mean(bias)) / sd(bias)
#   norm_mat = cbind(bias_norm, peak_counts_norm)
#   sample_peaks <- function(indices){
#     nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
#     s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
#     return(s)
#   }
#   n = length(peak_counts)
#   chunks = split(1:n, cut(1:n, breaks = seq(0, n+1, 1000) ))
#   sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
#   return(sampled_peaks)
# }
# 
# 
# getBackgroundPeakSets2 <- function(peak_counts, bias, niterations = 51, window = 2500, BPPARAM = BiocParallel::bpparam()){
#   #Standardise peak counts and bias
#   peak_counts_norm = rank(peak_counts, ties.method="average")
#   bias_norm = rank(bias, ties.method="average")
#   norm_mat = cbind(bias_norm, peak_counts_norm)
#   sample_peaks <- function(indices){
#     nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
#     s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
#     return(s)
#   }
#   n = length(peak_counts)
#   chunks = split(1:n, cut(1:n, breaks = seq(0, n+1, 1000) ))
#   sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
#   return(sampled_peaks)
# }
# 
# getBackgroundPeakSets3 <- function(peak_counts, bias, niterations = 51, window = 2500, BPPARAM = BiocParallel::bpparam()){
#   #Standardise peak counts and bias
#   peak_counts_norm = rank(peak_counts, ties.method="random")
#   bias_norm = rank(bias, ties.method="random")
#   norm_mat = cbind(bias_norm, peak_counts_norm)
#   ##Reflect across boundaries
#   half_window = window %/% 2
#   lowcounts = which(peak_counts_norm > length(peak_counts) - half_window)
#   lowmat = norm_mat[lowcounts,]
#   lowmat[,2] = length(peak_counts)*2 -lowmat[,2]
#   highcounts = which(peak_counts_norm < half_window)
#   highmat = norm_mat[highcounts,]
#   highmat[,2] = -highmat[,2]
#   lowbias = which(bias_norm > length(peak_counts) - half_window)
#   lowbiasmat = norm_mat[lowbias,]
#   lowbiasmat[,1] = length(peak_counts)*2 -lowbiasmat[,1]
#   highbias = which(bias_norm < half_window)
#   highbiasmat = norm_mat[highbias,]
#   highbiasmat[,1] = -highbiasmat[,1]
#   norm_mat = rbind(norm_mat, lowmat, highmat, lowbiasmat, highbiasmat)
#   sample_peaks <- function(indices){
#     nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
#     s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
#     return(s)
#   }
#   n = length(peak_counts)
#   chunks = split(1:n, cut(1:n, breaks = seq(0, n+1, 1000) ))
#   sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
#   #replace "reflected" paks
#   to_be_replaced = BiocParallel::bplapply((length(peak_counts)+1):nrow(norm_mat), function(x) which(sampled_peaks == x), BPPARAM = BPPARAM)
#   mapping = c(lowcounts, highcounts, lowbias, highbias)
#   replace_values = BiocParallel::bplapply(1:length(to_be_replaced), function(x) rep(mapping[x],length(to_be_replaced[[x]]) ), BPPARAM = BPPARAM)
#   sampled_peaks = replace(sampled_peaks, unlist(to_be_replaced), unlist(replace_values))
#   return(sampled_peaks)
# }
# 
# 
# setMethod("getBackgroundPeakSets", "fragmentCounts",
#           function(counts_mat, niterations = 50, window = 2500, BPPARAM = BiocParallel::bpparam()){
#             #if bias not available...
#             if ("bias" %ni% colnames(S4Vectors::mcols(counts_mat@peaks))){
#               stop("Peaks must have metadata column named 'bias'. Compute using compute_bias.")
#             }
#             #Standardise peak counts and bias
#             peak_counts_norm = rank(counts_mat@fragments_per_peak, ties.method="random")
#             bias_norm = rank(counts_mat@peaks$bias, ties.method="random")
#             norm_mat = cbind(bias_norm, peak_counts_norm)
#             n = nrow(norm_mat)
#             ##Reflect across boundaries
#             half_window = window %/% 2
#             lowcounts = which(peak_counts_norm > n - half_window)
#             lowmat = norm_mat[lowcounts,]
#             lowmat[,2] = n*2 -lowmat[,2]
#             highcounts = which(peak_counts_norm < half_window)
#             highmat = norm_mat[highcounts,]
#             highmat[,2] = -highmat[,2]
#             lowbias = which(bias_norm > n - half_window)
#             lowbiasmat = norm_mat[lowbias,]
#             lowbiasmat[,1] = n*2 -lowbiasmat[,1]
#             highbias = which(bias_norm < half_window)
#             highbiasmat = norm_mat[highbias,]
#             highbiasmat[,1] = -highbiasmat[,1]
#             norm_mat = rbind(norm_mat, lowmat, highmat, lowbiasmat, highbiasmat)
#             sample_peaks <- function(indices){
#               nn = FNN::get.knnx(norm_mat, norm_mat[indices,], algorithm = "kd_tree", k = window)$nn.index
#               s = t(sapply(1:nrow(nn), function(x) sample(nn[x,], size = niterations, replace = TRUE)))
#               return(s)
#             }
#             chunks = split(1:n, cut(1:n, breaks = seq(0, n + (n %% 1000), 1000) ))
#             if (niterations > 1){
#               sampled_peaks = do.call(rbind,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM))
#             }
#             else if (niterations == 1){
#               sampled_peaks = matrix(do.call(c,BiocParallel::bplapply(chunks, sample_peaks, BPPARAM = BPPARAM)), nrow = n, ncol =1)
#             }
#             #replace "reflected" paks
#             to_be_replaced = BiocParallel::bplapply((n+1):nrow(norm_mat), function(x) which(sampled_peaks == x), BPPARAM = BPPARAM)
#             mapping = c(lowcounts, highcounts, lowbias, highbias)
#             replace_values = BiocParallel::bplapply(1:length(to_be_replaced), function(x) rep(mapping[x],length(to_be_replaced[[x]]) ), BPPARAM = BPPARAM)
#             sampled_peaks = replace(sampled_peaks, unlist(to_be_replaced), unlist(replace_values))
#             background_peaks = backgroundPeaks(background_peaks = sampled_peaks, peaks = counts_mat@peaks)
#             return(background_peaks)
#           })
# 



