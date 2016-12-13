# Collection of internal utility functions

# To be deleted? ---------------------------------------------------------------

filterSeqs <- function(seqs){
  freqs <- Biostrings::alphabetFrequency(seqs)
  seqs <- seqs[which(rowSums(freqs[,-(1:4)])==0)]
  return(seqs)
}

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# Not-in operator --------------------------------------------------------------

'%ni%' = Negate('%in%')

# last element--------------------------------------------------------------------
last_element <- function(x) x[length(x)]

# check for installed package --------------------------------------------------

is.installed <- function(pkg) is.element(pkg, installed.packages()[,1])

# Mathematical shortcuts -------------------------------------------------------

mean_smooth <- function(X, window){
  if (is.na(as.integer(window)) || length(window) != 1 || window < 2 ||
        window >= length(X)){
    stop("window must be an integer between 2 and length(X)")
  }
  pad_left = rev(X[1:(window %/% 2)])
  pad_right = rev(X[(length(X)-((window-1) %/% 2)):length(X)])
  cx <- c(0,cumsum(c(pad_left,X,pad_right)))
  return((cx[(window+1):(length(cx)-1)] - cx[1:(length(cx)-1-window)])/window)
}

rms <- function(x) sqrt(mean(x**2))

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}

sd.along <- function(a, i) {
  n <- dim(a)[i]
  tmp.var <- means.along(a*a,i) - (means.along(a,i)**2)
  return(sqrt(tmp.var * n/(n-1)))
}




# Functions to test whether a vector is all TRUE or all FALSE ------------------
all_true <- function(x){
  stopifnot(inherits(x,"logical"))
  ifelse(sum(x)==length(x), TRUE, FALSE)
}

all_false <- function(x){
  stopifnot(inherits(x,"logical"))
  ifelse(sum(x)==0, TRUE, FALSE)
}

all_equal <- function(x){
  stopifnot(is.vector(x))
  all_true(x == x[1])
}



# Function to test whether vector is all whole number


all_whole <- function(x, tol = .Machine$double.eps^0.5){
  all_true(abs(x - round(x)) < tol)
}

# Function to merge lists, either by name or order -----------------------------

merge_lists <- function(..., by = c("order","name")){
  by <- match.arg(by)
  if (by == "order"){
    lx = sapply(list(...),length)
    if (!all_true(lx == lx[1])){
      max_lx = max(lx)
      fixed_lx = lapply(list(...), function(x){
        tmp_lx = length(x)
        if (tmp_lx < max_lx){
          return(c(x, sapply(1:(max_lx - tmp_lx), function(y) list(NULL))))
        } else{
          return(x)
        }
      })
      return(do.call(merge_lists,fixed_lx))
    }
    return(Map(c,...))
  } else{
    input <- list(...)
    if (length(input) ==1){
      output = input[[1]]
      if (all_false(duplicated(names(output)))){
        return(output)
      } else{
        tmp = input[[1]]
        unames = unique(names(output))
        output = lapply(unames, function(x) unlist(tmp[names(tmp) == x],
                                                   recursive = FALSE,
                                                   use.names=FALSE))
        names(output) = unames
        return(output)
      }
    } else{
      output = merge_lists(input[[1]], by = "name")
      for (j in 2:length(input)){
        tmp = merge_lists(input[[j]], by = "name")
        common_names = names(output)[which(names(output) %in% names(tmp))]
        unique_old = names(output)[which(names(output) %ni% common_names)]
        unique_new = names(tmp)[which(names(tmp) %ni% common_names)]
        output = c(Map(c,output[common_names],tmp[common_names]),
                   output[unique_old],tmp[unique_new])
      }
      return(output)
    }
  }
}

## remove overlapping or non-overlapping items
remove_overlap <- function(x, to.remove){
  #removes items in vector to.remove from x
  if (inherits(x, "list")){
    return(lapply(x, function(y) y[which(y %ni% to.remove)]))
  } else if (inherits(x, "vector")){
    return(x[which(x %ni% to.remove)])
  } else{
    stop("x must be list or vector")
  }
}

remove_nonoverlap <- function(x, to.keep){
  #keep only items in x that are in to.keep
  if (inherits(x, "list")){
    return(lapply(x, function(y) y[which(y %in% to.keep)]))
  } else if (inherits(x, "vector")){
    return(x[which(x %in% to.keep)])
  } else{
    stop("x must be list or vector")
  }
}

# Modified tabulate funcion ----------------------------------------------------

tabulate2<-function(x,min_val,max_val){
  if (max_val <= min_val){
    stop("max_val must be greater than min_val")
  }
  if (min_val<0 && max_val >0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    p = tabulate(x)[1:max_val]
    z = length(which(x == 0))
    out = c(n,z,p)
    out[which(is.na(out))] = 0
    names(out)=min_val:max_val
    return(out)}
  else if (min_val==0 && max_val >0){
    p = tabulate(x)[1:max_val]
    z = length(which(x == 0))
    out = c(z,p)
    out[which(is.na(out))] = 0
    names(out)=min_val:max_val
    return(out)}
  else if (min_val > 0 && max_val >0){
    out = tabulate(x)[min_val:max_val]
    out[which(is.na(out))] = 0
    names(out)=min_val:max_val
    return(out)
  }
  else if (min_val <0 && max_val == 0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    z = length(which(x == 0))
    out = c(n,z)
    out[which(is.na(out))] = 0
    names(out)=min_val:max_val
    return(out)}
  else if (min_val <0 && max_val < 0){
    n = rev(tabulate(-1*(x))[1:(-min_val)])
    out = n
    out[which(is.na(out))] = 0
    names(out)=min_val:max_val
    return(out)}
  else{
    stop("something may be amiss with min_val or max_val")
  }
}


# Functions with ranges --------------------------------------------------------

getstrand<- function(ranges){
  return(as.character(BiocGenerics::strand(ranges)))
}


# sorting alphanumeric ---------------------------------------------------------

split_alpha_numeric <- function(x){
  is.nonnumeric <- function(x) {
    is.na(suppressWarnings(as.numeric(x)))
  }
  y <- sapply(strsplit(x,"")[[1]], is.nonnumeric)
  diffs <- diff(c(TRUE,y))
  splitpoints = which(diffs != 0)
  if (length(splitpoints) == 0){
    return(x)
  } else if (length(splitpoints)==1){
    return(c(substr(x,0,splitpoints-1),substr(x,splitpoints,nchar(x))))
  }
  out <- substr(x,0,splitpoints[1]-1)
  for (i in 2:length(splitpoints)){
    out <- c(out, substr(x,splitpoints[i-1],splitpoints[i]-1))
  }
  out <- c(out, substr(x,splitpoints[length(splitpoints)],nchar(x)))
  return(out)
}

mxsort <- function(x){
  if (!is.character(x))
    return(sort(x))
  which.nas <- which(is.na(x))
  which.blanks <- which(x == "")
  split_x <- lapply(x, split_alpha_numeric)
  lx = sapply(split_x,length)
  if (!all_true(lx == lx[1])){
    max_lx = max(lx)
    split_x = lapply(split_x, function(y){
      tmp_lx = length(y)
      if (tmp_lx < max_lx){
        return(c(y, rep(NA, max_lx - tmp_lx)))
      } else{
        return(y)
      }
    })
  }
  vecs <- do.call(merge_lists,split_x)
  make.numeric <- function(x){
     ifelse(!(is.na(suppressWarnings(as.numeric(x)))), suppressWarnings(as.numeric(x)), x)
  }
  vecs <- lapply(vecs, make.numeric)
  ord <- do.call(order, c(vecs,list(na.last=FALSE)))
  return(x[ord])
}

# making permutations ----------------------------------------------------------

#' @export
make_bias_bins <- function(object, nbins = 25, frac = 0.3){
  bias = rowData(object)$bias
  if (is.null(bias)){
    stop("bias column in rowData of object must be set!")
  }
  npeaks = nrow(object)
  fragments_per_peak <- get_fragments_per_peak(object)
  #make bias bins
  bias_quantiles = quantile(bias, seq(0,1,1/nbins))
  bias_cut = cut(bias, breaks = bias_quantiles)
  bias_bins = split(1:npeaks, bias_cut)
  names(bias_bins) = sapply(1:nbins, function(x) paste("bias_bin_",x,sep="",collapse=""))
  #make count bins
  pseudo_counts = fragments_per_peak + runif(npeaks,min = 0, max = 0.1)
  count_quantiles = quantile(pseudo_counts, seq(0,1,1/nbins))
  count_cut = cut(pseudo_counts, breaks = count_quantiles)
  count_bins = split(1:npeaks, count_cut)
  names(count_bins) = sapply(1:nbins, function(x) paste("count_bin_",x,sep="",collapse=""))
  #make bias / count bins
  nbins = round(sqrt(nbins))
  bias_quantiles = quantile(bias, seq(0,1,1/nbins))
  bias_cut = cut(bias, breaks = bias_quantiles)
  tmp_bias_bins = split(1:npeaks, bias_cut)
  count_quantiles = quantile(pseudo_counts, seq(0,1,1/nbins))
  count_cut = cut(pseudo_counts, breaks = count_quantiles)
  tmp_count_bins = split(1:npeaks, count_cut)
  bias_count_bins = sapply(1:nbins, function(x) sapply(1:nbins, function(y) intersect(tmp_bias_bins[[y]], tmp_count_bins[[x]])))
  names(bias_count_bins) = sapply(1:nbins, function(x) sapply(1:nbins, function(y) paste("bias_count_bin_",x,"_",y,sep="",collapse="")))
  tmp = c(bias_bins, count_bins, bias_count_bins)
  sets = lapply(tmp, function(x) sample(x, size = length(x)* frac))
  mean_bias <- sapply(sets, function(x) mean(bias[x]))
  mean_counts <- sapply(sets, function(x) mean(fragments_per_peak[x]))
  rd <- DataFrame(bias = mean_bias, counts = mean_counts)
  out <- SummarizedExperiment(assays = list(matches = convert_from_ix_list(sets, nrow(object))),
                              colData = rd)
  return(out)
}

#' @export
make_permuted_sets <- function(object, annotations, window = 10){
  if (inherits(annotations, "SummarizedExperiment")){
    peak_indices <- convert_to_ix_list(assays(annotations)$match)
  } else {
    stop("annotations must be given as SummarizedExperiment with a 'matches' slot")
  }
  bg <- get_background_peaks_alternative(object, niterations = 1, window = window)
  sets <- lapply(seq_along(peak_indices), function(x) bg[peak_indices[[x]],1])
  names(sets) <- sapply(names(peak_indices), function(x) paste("permuted.",x,collapse="",sep=""))
  cd <- colData(annotations)
  if ("name" %in% colnames(cd)){
    cd$name <- sapply(cd$name, function(x) paste("permuted.",x,collapse="",sep=""))
  } else{
    cd$name <- names(sets)
  }
  rownames(cd) <- names(sets)
  out <- SummarizedExperiment(assays = list(matches = convert_from_ix_list(sets, nrow(object))),
                              colData = cd)
  return(out)
}

get_background_peaks_alternative <- function(object, niterations = 50, window = 500){
  
  bias = rowData(object)$bias
  
  if (is.null(bias)){
    stop("bias column in rowData must be set!")
  }
  
  fragments_per_peak <- get_fragments_per_peak(object)
  if (min(fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  npeak <- nrow(object)
  
  norm_mat <- cbind(fragments_per_peak, bias)
  
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(1:(npeak %/% grpsize + ((npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,npeak)))
  
  bghelper <- function(grp, tmp_vals, niterations){
    tmp_nns <- nabor::knn(tmp_vals, tmp_vals[grp,], window + 1, eps = 0)$nn.idx
    if (niterations == 1){
      return(matrix(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = TRUE)),
                    ncol = 1))
    } else{
      return(t(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = TRUE))))
    }
  }
  
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, tmp_vals, niterations))
  
  return(background_peaks)
}



# Convert matrix of indices to list of lists -----------------------------------

convert_to_ix_list <- function(ix){
  stopifnot(inherits(ix,"Matrix") || inherits(ix,"matrix"))
  if (inherits(ix,"sparseMatrix")){
    tmp <- summary(ix)
    tmp$j <- factor(tmp$j, levels = 1:ncol(ix), ordered = TRUE)
    out <- split(tmp$i, tmp$j)    
  } else{
    out <- lapply(1:ncol(ix), function(x) which(ix[,x] != 0))
  }
  names(out) = colnames(ix)
  return(out)  
}

# Convert list of indices to matrix  -------------------------------------------
convert_from_ix_list <- function(ix, npeaks){
  stopifnot(inherits(ix,"list"))
  sparseMatrix(j = unlist(lapply(seq_along(ix), function(x) rep(x, length(ix[[x]]))),use.names = FALSE),
                      i = unlist(ix, use.names = FALSE),
                      dims = c(npeaks, length(ix)),
                      x = 1,
               dimnames = list(NULL, names(ix)))
  
}


# Jaccard ----------------------------------------------------------------------

get_jaccard <- function(x){
  u <- crossprod(x)
  t <- colSums(x)
  i <- outer(t,t,FUN = "+")
  return(as.matrix(u/(i - u)))
}







