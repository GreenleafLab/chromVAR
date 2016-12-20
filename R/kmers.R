get_reverse_complement <- function(x){
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

ambiguity_mapping <- lapply(Biostrings::IUPAC_CODE_MAP, function(x){
  letters = strsplit(x,"")[[1]]
  out = rep(0, length(letters))
  out[which(letters == "A")] = 1
  out[which(letters == "C")] = 2
  out[which(letters == "G")] = 3
  out[which(letters == "T")] = 4
  return(out)  
})


seq_to_pwm <- function(in_seq, mismatch = 0){
  mat <- matrix(mismatch, ncol = nchar(in_seq), nrow = 4)
  rownames(mat) <- c("A","C","G","T")
  in_seq <- strsplit(in_seq,"")[[1]]
  for (x in seq_along(in_seq)){
    mat[ambiguity_mapping[[in_seq[x]]],x] = 1
  }
  return(mat) 
}



#' deviations_covariability
#'
#' @param object deviations result
#'
#' @return "covariability" matrix
#' @details Returns the "covariability" between motifs/kmers/peaksets.  "Covariability"
#' is defined as covariance between Z-scores divided by variance of Z-scores for
#' one motif/kmer/peakset (the row).
#' @export
deviations_covariability <- function(object){
  covs <- cov(t(assays(object)$z))
  vars <- row_sds(assays(object)$z)
  normed_covs <- covs / matrix(vars**2, nrow = nrow(covs), ncol = ncol(covs), byrow = FALSE)
  return(normed_covs)
}

get_kmer_dist <- function(kmers){
  stopifnot(all_equal(nchar(kmers)))
  out <- as.matrix(Biostrings::stringDist(kmers, method = "substitutionMatrix", type="overlap",
                         gapOpening = Inf,
                         substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                                                       mismatch = 0,
                                                                                       baseOnly = FALSE,
                                                                                 type = "DNA")))
  out <- nchar(kmers[1]) - out
  diag(out)<- 0
  return(out)
}

get_mm_kmers <- function(kmer){
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  kmer_nucs = strsplit(kmer,"")[[1]]
  out <- c()
  for (i in 1:k){
    for (j in nucs[nucs != kmer_nucs[i]]){
      kmod = kmer_nucs
      kmod[i] = j
      kmod = paste(kmod, sep="",collapse="")
      out <- c(out,kmod)
    }
  }
  return(out)
}


get_overlap_kmers <- function(kmer, max_extend, dir = "both"){
  out = c()
  if (max_extend <= 0){
    return(out)
  }
  k <- nchar(kmer)
  nucs <- c("A","C","G","T")
  if (max_extend == 1){
    if (dir == "both"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)
      return(c(out1,out2))
    } else if (dir == "left"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      return(out1)
    } else if (dir == "right"){
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)
      return(out2)
    }
  } else{
    if (dir == "both"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)
      return(c(out1,out2, sapply(out1, get_overlap_kmers, max_extend = max_extend - 1, dir = "left"),
               sapply(out2, get_overlap_kmers, max_extend = max_extend - 1, dir = "right")))
    } else if (dir == "left"){
      out1 <- sapply(nucs, function(x) paste(x,substr(kmer,1,k-1),sep="",collapse=""), USE.NAMES = FALSE)
      return(c(out1,sapply(out1, get_overlap_kmers, max_extend = max_extend - 1, dir = "left")))
    } else if (dir == "right"){
      out2 <- sapply(nucs, function(x) paste(substr(kmer,2,k),x,sep="",collapse=""), USE.NAMES = FALSE)
      return(c(out2,
               sapply(out2, get_overlap_kmers, max_extend = max_extend - 1, dir = "right")))
    }
  }
}

kmers_to_names <- function(kmers, names){
  ifelse(kmers %in% names, kmers, get_reverse_complement(kmers))
}

get_null_kmer_dist <- function(kmer, cov_mat, max_extend = 2){
  kmers <- colnames(cov_mat)
  kmer <- kmers_to_names(kmer, kmers)
  o <- get_overlap_kmers(kmer, max_extend = max_extend)
  mm <- get_mm_kmers(kmer)

  omm <- Reduce(union, lapply(mm, function(k) unique(get_overlap_kmers(k, max_extend = max_extend + 1))))

  null_dist <- cov_mat[kmer,which(colnames(cov_mat) %ni% kmers_to_names(unique(c(kmer,mm,o,omm)),kmers))]
  return(null_dist)
}


get_covariable_kmers <- function(kmer, cov_mat,  max_extend = 2){

  #cov_mat <- deviations_covariability(object)
  kmers <- colnames(cov_mat)
  kmer <- kmers_to_names(kmer, kmers)
  #Get Null Dist
  o <- get_overlap_kmers(kmer, max_extend = max_extend)
  mm <- get_mm_kmers(kmer)

  omm <- Reduce(union, lapply(mm, function(k) unique(get_overlap_kmers(k, max_extend = max_extend + 1))))

  null_dist <- cov_mat[kmer,which(colnames(cov_mat) %ni% kmers_to_names(unique(c(kmer,mm,o,omm)),kmers))]

  candidates <- unique(c(mm,o))
  scores <- (cov_mat[kmer,unique(kmers_to_names(candidates,kmers))] - mean(null_dist,na.rm=TRUE)) / sd(null_dist,na.rm=TRUE)
  pvals <- pnorm(scores, lower.tail = FALSE)
  pvals.adj <- p.adjust(pvals)

  o_shifts <- if (max_extend >=1) do.call(c, lapply(1:max_extend, function(x) c(rep(-x,4^x),rep(x,4^x)))) else NULL
  out <- data.frame(kmer = c(kmer, mm, o),
                    mismatch = c(NA, rep(1:nchar(kmer),each = 3), rep(NA,length(o))),
                    shift = c(rep(0,length(mm)+1),o_shifts),
                    covariability = cov_mat[kmer, kmers_to_names(c(kmer, mm, o), kmers)],
                    pval = pvals[kmers_to_names(c(kmer, mm, o), kmers)],
                    pval.adj = pvals.adj[kmers_to_names(c(kmer, mm, o), kmers)],
                    stringsAsFactors = FALSE)

  return(out)
}


kmer_group_to_pwm <- function(kgroup, p = 0.01, threshold = 0.25){
  min_shift = min(kgroup$shift)
  max_shift = max(kgroup$shift)
  k = nchar(kgroup$kmer[1])
  nucs = c("A","C","G","T")
  out <- matrix(0, nrow = 4, ncol = k + abs(min_shift) + max_shift,
                dimnames = list(nucs,NULL))
  if (min_shift < 0 ){
    for (i in min_shift:-1){
      ix <- which(kgroup$shift == i)
      nuc <- substr(kgroup$kmer[ix], 1,1)
      covars <- sapply(nucs, function(x) max(c(0,kgroup$covariability[ix[which(nuc == x)]])))
      covars[!is.finite(covars)] = 0
      baseline <- rep(0.25,4) * (1 + covars)
      baseline <- baseline / sum(baseline)
      pvals <- sapply(nucs, function(x) min(kgroup$pval.adj[ix[which(nuc == x)]]))
      covars[pvals > p] = 0
      w <- max(covars)
      if (w > 1) w = 1
      out[,i - min_shift  + 1] = (1 - w) * baseline + (w * covars**2 / (sum(covars**2) + all_true(covars == 0)))
    }
  }
  for (i in 1:k){
    ix <- c(1,which(kgroup$mismatch == i))
    nuc <- substr(kgroup$kmer[ix], i,i)
    covars <- sapply(nucs, function(x) max(c(0,kgroup$covariability[ix[which(nuc == x)]])))
    covars[!is.finite(covars)] = 0
    pvals <- sapply(nucs, function(x) min(kgroup$pval.adj[ix[which(nuc == x)]]))
    covars[pvals > p] = 0
    out[,abs(min_shift) + i] = covars**2 / sum(covars**2)
  }
  if (max_shift > 0 ){
    for (i in 1:max_shift){
      ix <- which(kgroup$shift == i)
      nuc <- substr(kgroup$kmer[ix], k,k)
      covars <- sapply(nucs, function(x) max(c(0,kgroup$covariability[ix[which(nuc == x)]])))
      covars[!is.finite(covars)] = 0
      baseline <- rep(0.25,4) * (1 + covars)
      baseline <- baseline / sum(baseline)
      pvals <- sapply(nucs, function(x) min(kgroup$pval.adj[ix[which(nuc == x)]]))
      covars[pvals > p] = 0
      w <- max(covars)
      if (w > 1) w = 1
      out[,abs(min_shift) + k + i] = (1 - w) * baseline + (w * covars**2 / (sum(covars**2) + all_true(covars == 0)))
    }
  }
  start = 1
  bits <- function(x){
    2 + sum(x*ifelse(x > 0, log2(x),0))
  }
  if (min_shift < 0 ){
    for (i in min_shift:-1){
      if (bits(out[,i - min_shift  + 1]) < threshold){
        start = start + 1
      } else{
        break
      }
    }
  }
  end = ncol(out)
  if (max_shift > 0 ){
    for (i in max_shift:1){
      if (bits(out[,abs(min_shift) + k + i]) < threshold){
        end = end - 1
      } else{
        break
      }
   }
  }
  return(out[,start:end])
}


#' plot_kmer_mismatch
#'
#' @param kmer kmer, e.g. "AAAAAAA"
#' @param cov_mat result from \code{\link{deviations_covaraibility}}
#' @param pval p value threshold
#'
#' @return A plot
#' @export
plot_kmer_mismatch <- function(kmer, cov_mat, pval = 0.01){

  kgroup <- get_covariable_kmers(kmer, cov_mat, max_extend = 0)
  ix <- which(!is.na(kgroup$mismatch))

  mm_df <- rbind(data.frame(val = ifelse(kgroup$pval.adj[ix] < 0.01,kgroup$covariability[ix]**2,0), pos = kgroup$mismatch[ix],
                       Nucleotide = substr(kgroup$kmer[ix], kgroup$mismatch[ix],kgroup$mismatch[ix])),
                 data.frame(val = 1, pos = 1:nchar(kmer), Nucleotide = strsplit(kmer,"")[[1]]))

  ggplot(mm_df) +
    geom_point( aes(x = pos, y = val, col = Nucleotide), position = position_jitter(height = 0, width = 0.2))  +
     ylab("Relative Nucleotide\nVariability Score") +
    xlab("Position") + scale_x_continuous(breaks = c(1:max(mm_df$pos)))+ scale_y_continuous(breaks = c(0,0.5,1))+
    chromVAR_theme() + scale_color_manual(name="Nucleotide",
                                          breaks = c("A","C","G","T"),
                                          values = RColorBrewer::brewer.pal(4,"Dark2")) +
    ggmotif(kmer, x.pos = 0, y.pos = max(c(1.1, max(mm_df$val) + 0.1)), ht = 0.15, wt = 1) + ggmotif_scale()


}


#' assemble_kmers
#'
#' @param object kmer deviations object
#' @param threshold variability threshold
#' @param p p value threshold for inclusion of kmer
#' @details function for assembling de novo kmers from kmer deviations
#' @return list with (1) motifs: de novo motif matrices, (2) seed: seed kmer for de novo motif 
#' @export
assemble_kmers <- function(object, threshold = 1.5, p = 0.01){
  devco = deviations_covariability(object)
  vars = row_sds(assays(object)$z)
  cands = rownames(object)[order(vars, decreasing = TRUE)[1:sum(vars > threshold,na.rm=TRUE)]]
  out = list(motifs = list(), seed = list())
  nc <- length(cands)
  pb <- txtProgressBar(min = 0, max = nc, style = 3)
  while (length(cands) >1){
    kgroup = get_covariable_kmers(cands[1], devco)
    kmotif = kmer_group_to_pwm(kgroup, p)
    if (max(kgroup$covariability, na.rm = TRUE) <= 1){
      out$motifs = c(out$motifs, list(kmotif))
      out$seed = c(out$seed, cands[1])
    }
    nd = get_null_kmer_dist(cands[1], devco)
    z <- (devco[cands[1],cands] - mean(nd, na.rm = TRUE)) / sd(nd, na.rm = TRUE)
    pval <- pnorm(z, lower.tail = FALSE)
    d <- pwm_distance(kmotif, lapply(cands,seq_to_pwm))$d[1,]
    exc <- intersect(which(d < 0.25), which(pval < p))
    cands <- cands[-exc]
    setTxtProgressBar(pb, nc -length(cands))
  }
  close(pb)
  denovo_motifs <- do.call(TFBSTools::PWMatrixList,lapply(seq_along(out$motifs), function(x) TFBSTools::PWMatrix(ID = paste0("denovo_",x),
                                                                                                                 name = paste0("denovo_",x),
                                                                                                                 tags = list(seed = out$seed[[x]]),
                                                                                                                 profileMatrix = out$motifs[[x]])))
  names(denovo_motifs) <- name(denovo_motifs)
  return(denovo_motifs)
}



#' @import plotly shiny
motif_explorer <- function(motifs,
                          kmers,
                          tsne,
                          annotation_column){


  #if ("tsne" %in% names(tsne)) tsne <- tsne$tsne

  type <- c(rep("motif",nrow(motifs)),rep("kmer",nrow(kmers)))

  vars <- row_sds(rbind(assays(motifs)$z,assays(kmers)$z), FALSE)
  mat <- rbind(assays(motifs)$deviations,assays(kmers)$deviations)[tsne$ix,]

  tmp1 <- motifs[tsne$ix[which(type[tsne$ix] == "motif")],]
  rowData(tmp1) <- NULL
  tmp2 <- kmers[tsne$ix[which(type[tsne$ix] == "kmer")] - nrow(motifs),]
  rowData(tmp2) <- NULL
  covars <- deviations_covariability(rbind(tmp1,tmp2))
  pwm_dist <- pwm_distance(c(rowData(motifs)$pwm,
                             do.call(TFBSTools::PWMatrixList, lapply(rownames(kmers), function(x) TFBSTools::PWMatrix(ID = x, profileMatrix = seq_to_pwm(x)))))[tsne$ix])

  if (!is.null(motifs)){
    cd <- colData(motifs)
  } else{
    cd <- colData(kmers)
  }

  if (is.null(annotation_column)){
    annotation_column <- "none"
  } else {
    stopifnot(annotation_column %in% colnames(cd))
  }


  ui <- fluidPage(
    fluidRow(column(3,plotlyOutput("plot1")),
             column(3,
                    fluidRow(plotOutput("plot2")),
                    fluidRow(plotlyOutput("plot3"))),
             column(3,
                    fluidRow(plotOutput("plot4")),
                    fluidRow(plotOutput("plot5"))))
  )


  server <- function(input, output, session) {


    # Render the plot
    output$plot1 <- renderPlotly({
      if (annotation_column == "none"){
        p1 = ggplot(data.frame(x = tsne$tsne$Y[,1], y = tsne$tsne$Y[,2],
                               text = c(rownames(motifs),rownames(kmers))[tsne$ix],
                               Variability = vars[tsne$ix],
                               type = type[tsne$ix]),
                    aes(x = x, y = y,  text = text, size = Variability, shape = type))  +
          geom_point(size = 2) + chromVAR_theme(12) + scale_size(range=c(0.5,2.5)) + scale_shape_manual(values = c(1,16)) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2")
      } else {
        mean_dev_anno <- sapply(unique(cd[,annotation_column]),
                                function(x) rowMeans(mat[,which(cd[,annotation_column] == x)]))
        motif_anno <- apply(mean_dev_anno,1,function(x) colnames(mean_dev_anno)[which.max(x)])

        p1 = ggplot(data.frame(x = tsne$tsne$Y[,1], y = tsne$tsne$Y[,2], color = motif_anno,
                               text = c(rownames(motifs),rownames(kmers))[tsne$ix],
                               Variability = vars[tsne$ix],
                               type = type[tsne$ix]),
                    aes(x = x, y = y, col = color, text = text, size = Variability, shape = type))  +
          geom_point(size = 2) + chromVAR_theme(12) + scale_size(range=c(0.5,2.5)) + scale_shape_manual(values = c(1,16)) +
          labs(col = "Max Deviations:\n") +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid::unit(0.5,"lines"))
      }
      ggplotly(p1)
    })

    s <- reactive({event_data(event = "plotly_click")})

    output$plot2 <- renderPlot({
      sel = s()
      if (length(sel) == 0) return(NULL)
      print(sel$pointNumber)
      print(sel$curveNumber)
      print(type[tsne$ix[sel$pointNumber+1]])
      if (type[tsne$ix[sel$pointNumber+1]] == "motif"){
        pwm = pwm_to_prob(rowData(motifs)[tsne$ix[sel$pointNumber+1],"pwm"][[1]])[[1]]
      } else{
        pwm = seq_to_pwm(rownames(kmers)[tsne$ix[sel$pointNumber+1] - nrow(motifs)])
      }
      ggmotif::ggmotif_plot(pwm)
    })

    output$plot3 <- renderPlotly({
      sel = s()
      if (length(sel) == 0) return(NULL)
      p = ggplot(data.frame(Similarity = pwm_dist$dist[sel$pointNumber,],
                        Covariability = covars[sel$pointNumber,],
                        type = type[tsne$ix],
                        text = c(rownames(motifs),rownames(kmers))[tsne$ix])) +
        geom_point(aes(x=Similarity,
                       y = Covariability,
                       shape = type,
                       text = text)) +
        scale_shape_manual(values = c(1,16)) + xlab("Distance between PWM")
      ggplotly(p)
    })

    output$plot4 <- renderPlot({
      NULL
    })

    output$plot5 <- renderPlot({
      NULL
    })
  }

  shinyApp(ui, server, options = list(width = 750))
}


toIC <- function(mat){
  mat * matrix(apply(mat,2, function(x) 2 + sum(log(x)*x)),nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
}
