# Plotting ---------------------------------------------------------------------            

#' plot_variability
#' 
#' plot variability of motifs/etc
#' @param variability output from \code{\link{compute_variability}}
#' @param xlab label for x-axis (default is "Sorted TFs")
#' @param n  number of toppoints to label?
#' @param labels names of sets. if not given, uses rownames of variability
#' @import ggplot2
#' @export
plot_variability <- function(variability, xlab = "Sorted TFs", 
                   n = 3, labels = NULL){
  
  if (is.null(labels)) labels = rownames(variability)
  res_df = cbind(variability,ranks = rank(-1 * variability$variability,
                                          ties.method="random"),
                 tf = labels)
  
  ylab = "Variability"
  
  if ("bootstrap_lower_bound" %ni% colnames(variability)){
    
    out = ggplot2::ggplot(res_df, ggplot2::aes_string(x = "ranks", 
                                                      y = "variability")) + 
      geom_point() +
      ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + 
      ggplot2::scale_y_continuous(expand=c(0,0),
                                  limits=c(0,max(res_df$variability, na.rm =TRUE)*1.05)) +
      chromVAR_theme() 
  } else{
    
    out = ggplot2::ggplot(res_df, ggplot2::aes_string(x = "ranks", 
                                                    y = "variability",
                                                    min = "bootstrap_lower_bound", 
                                                    max = "bootstrap_upper_bound")) + 
    ggplot2::geom_point()+ ggplot2::geom_errorbar() +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + 
    ggplot2::scale_y_continuous(expand=c(0,0),
                                limits=c(0,max(res_df$bootstrap_upper_bound, na.rm =TRUE)*1.05)) +
    chromVAR_theme()    
  }
  
  if (n >=1){
    top_df = res_df[res_df$ranks <= n,]
    out = out + ggplot2::geom_text(data = top_df, 
                                   ggplot2::aes_string(x = "ranks", 
                                                       y = "variability", 
                                                       label = "tf"),
                                   size = 3, hjust=-0.45,col = "Black")
  } 
  
  return(out)
  
}


## helper functions (not exported) ----------------------------------------------

chromVAR_theme <-function(base_size = 12, base_family="Helvetica"){
  theme(line = element_line(colour = "black", size = 0.5, linetype = 1, 
                            lineend = "butt"), 
        rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1), 
        text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, 
                            vjust = 0.5, angle = 0, lineheight = 0.9), 
        strip.text = element_text(size = rel(0.8)), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = rel(0.8), colour = "black"),
        axis.text.x = element_text(vjust = 1), 
        axis.text.y = element_text(hjust = 1), 
        axis.ticks = element_line(colour = "black"), 
        axis.title.x = element_text(), 
        axis.title.y = element_text(angle = 90), 
        axis.ticks.length = grid::unit(0.15, "cm"), 
        axis.ticks.margin = grid::unit(0.1, "cm"), 
        legend.background = element_rect(colour = NA), 
        legend.margin = grid::unit(0.2, "cm"), 
        legend.key = element_blank(), 
        legend.key.size = grid::unit(1.2, "lines"), 
        legend.key.height = NULL, 
        legend.key.width = NULL, 
        legend.text = element_text(size = rel(0.8)), 
        legend.text.align = NULL, 
        legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0), 
        legend.title.align = NULL, 
        legend.position = "right", 
        legend.direction = NULL, 
        legend.justification = "center", 
        legend.box = NULL, 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.margin = grid::unit(0.25, "lines"), 
        panel.margin.x = NULL, 
        panel.margin.y = NULL, 
        strip.background = element_blank(), 
        strip.text.x = element_text(), 
        strip.text.y = element_text(angle = -90), 
        plot.background = element_blank(), 
        plot.title = element_text(size = rel(1.2)), 
        plot.margin = grid::unit(c(0.5, 0.5, 0.25, 0.25), "lines"), complete = TRUE)
}



pretty_scientific <- function(l) {
  # format as scientific
  l <- format(l, nsmall = 0, scientific = TRUE)
  # remove + sign
  l <- gsub("+", "", l, fixed=T)
  # break into prefix and suffix
  pre <- sapply(l, function(x) substr(x,1,gregexpr("e",x)[[1]][1]-1))
  post <- format(as.numeric(sapply(l, function(x) substr(x,gregexpr("e",x)[[1]][1]+1,nchar(x)))))
  # combine prefix and suffix with plotmath
  out <- sapply(1:length(l), function(x) paste(pre[x],"%*%10^",post[x],sep="",collapse=""))
  out[which(pre=="")]=NA
  # return as expression
  return(parse(text=out))
}


order_of_magnitude <- function(x){
  if (x==0){
    return(0)
  }
  else if (x< 0){
    x = -1 * x
  }
  return(floor(log10(x)))
}


pretty_scale_format <- function(l){
  digits = order_of_magnitude(max(l)) - order_of_magnitude(min(diff(l))) + 2
  l = signif(l, digits = digits)
  if (max(l)>1000){
    return(pretty_scientific(l))
  }
  else if (max(l)<0.001){
    return(pretty_scientific(l))
  }
  else{return(format(l, nsmall = 0))}
}

cor_dist <- function(x){
  as.dist(1 - cor(t(x), use = "pairwise.complete.obs"))
}


#' plot_deviations
#' 
#' plot heatmap of deviations
#' @param X deviations matrix (either z-score or fold change)
#' @param xlabel label for x axis
#' @param ylabel label for y axis
#' @param name title for plot
#' @param cluster_row cluster rows?
#' @param cluster_col cluster columns?
#' @param cluster_row_dist distance function to use for clustering rows
#' @param cluster_col_dist distance function to use for clustering columns
#' @param sample_annotation annotation of sample
#' @return ggplot object
#' @export
plot_deviations <- function(X, 
                            xlabel = NA, 
                            ylabel = NA,  
                            name = NA, 
                            cluster_row = TRUE, 
                            cluster_col = TRUE, 
                            cluster_row_dist = cor_dist, 
                            cluster_col_dist = cor_dist,
                            sample_annotation = NULL,
                            set_names = NULL){
  
  if (cluster_row){
    rowclust = hclust(cluster_row_dist(X))
    X = X[rowclust$order,]
    if (!is.null(set_names)){
      set_names = set_names[rowclust$order]
    }
  }
  if (cluster_col){
    colclust = hclust(cluster_col_dist(t(X)))
    X = X[,colclust$order]
    if (!is.null(sample_annotation)){
      sample_annotation = sample_annotation[colclust$order]
    }
  }
  
  df = cbind(data.frame("y" = factor(rownames(X), levels = rownames(X), ordered=T)),X)
  mdf = reshape2::melt(df, id = "y")
  p = ggplot(mdf, aes(x=variable, y=y)) + 
    geom_raster(aes(fill=value)) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 6),
      plot.margin = grid::unit(c(0.1,0.1,0.1,0.1),"cm")
    )
  if (cluster_row){
    dendro_row = ggdendro::dendro_data(rowclust)
    dendro_segments_row = dendro_row$segments
    dendro_segments_row[,c("y","yend")] = (dendro_segments_row[,c("y","yend")] * sqrt(nrow(X)*ncol(X)) / max(dendro_segments_row[,c("y","yend")]) * 0.4) + ncol(X) + 1
    p = p + geom_segment(data = dendro_segments_row,  mapping = aes(x=y,xend=yend, y = x, yend = xend, col=NULL))
  }
  if (cluster_col){
    dendro_col = ggdendro::dendro_data(colclust)
    dendro_segments_col = dendro_col$segments
    dendro_segments_col[,c("y","yend")] = (dendro_segments_col[,c("y","yend")]* sqrt(nrow(X)*ncol(X)) / max(dendro_segments_col[,c("y","yend")])  * 0.4) + nrow(X) + 1
    p = p + geom_segment(data = dendro_segments_col,  mapping = aes(x=x,xend=xend, y = y, yend = yend, col=NULL))
  }
  
  colors = c(scales::muted("blue"),"white",scales::muted("red"))
  limits = c(min(mdf$value,na.rm=T), max(mdf$value,na.rm=T))
  guidebreaks = c(limits[1],0,limits[2])

  p = p + scale_fill_gradient2(name = "Deviation\nScore",limits = limits,low = colors[1], mid = colors[2], high = colors[3], 
                               breaks = guidebreaks, label = pretty_scale_format, expand=c(0,0), midpoint = 0)
  if (!is.null(set_names)){
    p = p + scale_y_discrete(breaks = levels(mdf$y), labels = set_names)  
  }
  if (!is.na(xlabel)){
    p = p + xlab(xlabel)
  } else{
    p = p + theme(axis.title.x = element_blank())
  }
  if (!is.na(ylabel)){
    p = p + ylab(ylabel)
  } else{
    p = p + theme(axis.title.y = element_blank())
  }  
  if (!is.na(name)){
    p = p + ggtitle(name)
  }
  if (!is.null(sample_annotation)){  
    cols = RColorBrewer::brewer.pal(n  = length(levels(sample_annotation)), "Dark2")
    anno_col = cols[as.integer(sample_annotation)]
    tmp_df = data.frame(x = 1:ncol(X),  y = -2.75,  z = anno)
    p = p +  annotate("rect",xmin = 1:ncol(X) -0.5, xmax = 1:ncol(X) + 0.5, ymax = -1.5, ymin = min(-2 - nrow(X) *0.025,-3), col = anno_col,
                      fill = anno_col) + geom_point(data = tmp_df, 
                        mapping= aes(x = x, y =y, colour = z, fill = NULL), size = 0) + 
      scale_color_brewer(name = "Sample\nAnnotation", palette = "Dark2") + theme(legend.box = "horizontal", legend.background=element_blank(), axis.text.x = element_blank())
  }
  return(p)

}

# Plotting kmer group ----------------------------------------------------------

#'@export
plot_kmer_group <- function(a, similar_motifs = NULL, plot.consensus = TRUE){
  p = ggplot()
  for (i in 1:nrow(a)){
    p = p + ggmotif(a$kmer[i], y.pos = (nrow(a)-i)*1.25, x.pos = a$shift[i])
  }  
  anno_df = data.frame(y = (nrow(a)-1)*1.25+0.5, label = "K-mers")
  if (plot.consensus){  
    consensus = Biostrings::consensusMatrix(Biostrings::DNAStringSet(a$kmer),
                                            shift = a$shift - min(a$shift))[1:4,]
    consensus_shift = align_pwms(consensus, seq_to_pwm(a$kmer[1]), both_strands = FALSE)$pos
    tmp_y = min(sapply(ggplot_build(p)$data, function(obj) min(obj$y))) - 2
    p = p + ggmotif(consensus / matrix(apply(consensus,2,sum),nrow = 4, ncol=ncol(consensus),byrow=TRUE), 
                    y.pos = tmp_y, 
                    x.pos = consensus_shift)
    anno_df = rbind(anno_df, data.frame(y = tmp_y + 0.5, label = "Consensus"))
  }
  if (!is.null(similar_motifs)){
    if (!plot.consensus){
      consensus = Biostrings::consensusMatrix(Biostrings::DNAStringSet(a$kmer),
                                              shift = a$shift - min(a$shift))[1:4,]
      consensus_shift = align_pwms(consensus, seq_to_pwm(a$kmer[1]), both_strands = FALSE)$pos
    }
    tmp_y = min(sapply(ggplot_build(p)$data, function(obj) min(obj$y))) - 2
    if (inherits(similar_motifs,"PFMatrixList")){
      similar_motifs = TFBSTools::PWMatrixList(lapply(similar_motifs, TFBSTools::toPWM))
    }    
    for (i in 1:length(similar_motifs)){
      m = as.matrix(similar_motifs[[i]])
      tmp_a = align_pwms(m, consensus, both_strands = TRUE)
      tmp_shift = tmp_a$pos[1] + consensus_shift
      if (tmp_a$strand[1] == -1) m = Biostrings::reverseComplement(m)
      p = p + ggmotif(m, 
                      y.pos = tmp_y - max(m), 
                      x.pos = tmp_shift)
      anno_df = rbind(anno_df, data.frame(y = tmp_y - (i-1) * 1.25 + 0.5, label = name(similar_motifs[[i]])))
      tmp_y = tmp_y - max(apply(m, 2, function(x) sum(abs(x))))
    }
  }
  
  out = p + ggmotif_scale() + ggmotif:::ggmotif_theme() +
    scale_x_continuous(breaks = 0:max(sapply(1:nrow(a), function(x) nchar(a$kmer[x]) + a$shift[x]))) +
    scale_y_continuous(breaks = anno_df$y, labels = anno_df$label) +
    xlab("position relative to start of seed kmer (bp)") +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face="bold",size=12),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  return(out)
}





