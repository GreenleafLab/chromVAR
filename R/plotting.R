# Plotting ---------------------------------------------------------------------

#' plot_variability
#'
#' plot variability of motifs/etc
#' @param variability output from \code{\link{compute_variability}}
#' @param xlab label for x-axis (default is "Sorted TFs")
#' @param n  number of toppoints to label?
#' @param labels names of sets. if not given, uses rownames of variability
#' @import ggplot2 plotly
#' @export
plot_variability <- function(variability, xlab = "Sorted TFs",
                   n = 3, labels = variability$name, interactive = TRUE){

  res_df = cbind(variability, rank = rank(-1 * variability$variability,
                ties.method="random"), annotation = labels)

  ylab = "Variability"

  if ("bootstrap_lower_bound" %ni% colnames(variability)){

    out = ggplot(res_df,aes_string(x = "rank", y = "variability")) +
      geom_point() + xlab(xlab) + ylab(ylab) +
      scale_y_continuous(expand=c(0,0),  limits=c(0,max(res_df$variability, na.rm =TRUE)*1.05)) +
      chromVAR_theme()
    
  } else {

    out = ggplot(res_df, aes_string(x = "rank", y = "variability",
            min = "bootstrap_lower_bound",max = "bootstrap_upper_bound",
            label = "annotation")) + geom_point()+ geom_errorbar() +
            xlab(xlab) + ylab(ylab) + chromVAR_theme() +
            scale_y_continuous(expand=c(0,0),limits=c(0,max(res_df$bootstrap_upper_bound, na.rm =TRUE)*1.05))
  }

  if (interactive){
    return(ggplotly(out))
  } else{
    if (n >=1){
      top_df = res_df[res_df$rank <= n,]
      out = out + geom_text(data = top_df,
                            size = 3, hjust=-0.45,col = "Black")
    }
    return(out)
  }

}

#' chromVAR_theme
#'
#' @param base_size base font size
#' @param base_family base font family
#'
#' @return ggplot2 theme
#' @export 
#'

chromVAR_theme <- function(base_size = 12, base_family="Helvetica"){
  theme_bw(base_size = base_size, base_family = base_family)  %+replace%
    theme(panel.border = element_blank(),
          axis.line.x = element_line(colour = "black", size = 0.5, linetype = 1,
                                   lineend = "butt"),
          axis.line.y = element_line(colour = "black", size = 0.5, linetype = 1,
                                     lineend = "butt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          plot.background = element_blank())
}

## helper functions (not exported) ----------------------------------------------


cor_dist <- function(x){
  as.dist(1 - cor(t(x), use = "pairwise.complete.obs"))
}


#' plot_deviations
#'
#' plot heatmap of deviations
#' @param object chromVAR backend object
#' @param what "z" or "deviations"
#' @param cluster_rows cluster rows?
#' @param cluster_cols cluster columns?
#' @param sample_annotation annotation of sample
#' @param show_rownames show row names?
#' @param show_colnames show column names?
#' @return pheatmap output
#' @import pheatmap
#' @export
plot_deviations_heatmap <- function(object, top = 100, what = c("z","deviations"),
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            sample_annotation = colData(object),
                            show_rownames = FALSE, show_colnames = FALSE, ...){

  what <- match.arg(what)

  x <- if (what =="z") assays(object)$z else assays(object)$deviations

  if (top < nrow(object)){
    vars <- compute_variability(object, bootstrap_error = FALSE)$variability
    x <- x[which(vars > quantile(vars, 1 - (top/length(vars)))),]
  }

  if (isTRUE(cluster_cols)){
    colpc = prcomp(x)
    colpci = summary(colpc)$importance[2,]
    cluster_cols = hclust(dist(colpc$rotation[,which(colpci > 0.01)]))
  }

  pheatmap(x, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
           annotation_col = as.data.frame(sample_annotation), show_rownames = show_rownames,
           show_colnames = show_colnames, ...)

}




variability_table <- function(var_df){
  DT::datatable(var_df, options = list(order = list(2,'desc'), pageLength = 5),
                selection = list(mode = "single", selected = which.max(var_df$variability))) %>% DT::formatSignif(colnames(var_df), 3)
}


#' plot_deviations_tsne
#'
#' @param object deviations result object
#' @param tsne result from \code{\link{deviations_tsne}}
#' @param var_df variability result
#' @param annotation_column column name for colData(object) of annotation to be used
#' for coloring samples
#'
#' @return shiny widget
#' @export
plot_deviations_tsne <- function(object,
                                 tsne,
                                 var_df = NULL,
                                 annotation_column = NULL,
                                 motif = NULL,
                                 shiny = interactive()){

  if (shiny) return(plot_deviations_tsne_shiny(object, tsne, var_df, annotation_column))

  stopifnot(annotation_column %in% colnames(colData(object)))
  anno = colData(object)[,annotation_column]

  if ("tsne" %in% names(tsne)) tsne <- tsne$tsne

  out <- list()
  if (!is.null(annotation_column)){
    for (i in annotation_column){
      anno = colData(object)[,i]
      out[[i]] =ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = anno, text = colnames(object)),
             aes(x = x, y = y, col = color, text = text))  +
        geom_point(size = 2) + chromVAR_theme(12) +
        scale_color_brewer(palette = "Dark2", name = i) +
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
    }
  }
  if (!is.null(motif)){
    for (i in motif){
      if (i %in% rownames(object)){
        ix <- motif
      } else if (i %in% rowData(object)$name){
        ix <- which(rowData(object)$name == i)
        if (length(ix) > 1) ix <- ix[which.max(row_sds(assays(object[ix,])$z))]
      } else if (is.numeric(i)){
        ix <- i
      } else{
        stop("motif name invalid")
      }
      out[[i]] <- ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = assays(object)[["deviations"]][ix,], text = colnames(object)),
                         aes(x = x, y = y, col = color, text = text)) + geom_point(size= 2) +
        scale_color_gradient2(name = rowData(object)[ix,"name"],mid = "lightgray", low = "blue", high = "red") + chromVAR_theme(12) +
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
    }
  }
  return(out)
}

#' @import plotly shiny
plot_deviations_tsne_shiny <- function(object,
                                 tsne,
                                 var_df,
                                 annotation_column){

  if (is.null(var_df)) var_df <- compute_variability(object)

  if ("tsne" %in% names(tsne)) tsne <- tsne$tsne
  if (is.null(annotation_column)){
    annotation_column <- "none"
  } else {
    stopifnot(annotation_column %in% colnames(colData(object)))
  }


  ui <- fluidPage(
    fluidRow(column(3,
                    h4("Select options for coloring plots:"),
                    br(),
                    selectInput("color", "Color first plot by:",
                                choices = c("none",colnames(colData(object))),
                                selected = annotation_column),
                    br(),
                    p(strong("Color second plot by selecting from table:"))),
            column(9,DT::dataTableOutput('tbl', width = 350))),
    fluidRow(column(6,plotlyOutput("plot1")),
             column(6,plotlyOutput("plot2")))
  )

  var_tab <- variability_table(var_df)

  server <- function(input, output, session) {


    output$tbl = DT::renderDataTable(
      var_tab
    )

    # Render the plot
    output$plot1 <- renderPlotly({
      if (input$color == "none"){
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], text = colnames(object)),
                    aes(x = x, y = y, text = text))  +
          geom_point(size = 2) + chromVAR_theme(12) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      } else{
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = colData(object)[,input$color], text = colnames(object)),
                  aes(x = x, y = y, col = color, text = text))  +
          geom_point(size = 2) + chromVAR_theme(12) +
          labs(col = input$color)  +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      }
      ggplotly(p1)
    })

    output$plot2 <- renderPlotly({
      s = input$tbl_rows_selected
      if (length(s) == 0) return(NULL)
      p2 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = assays(object)[["deviations"]][s,], text = colnames(object)),
                  aes(x = x, y = y, col = color, text = text)) + geom_point(size= 2) +
        scale_color_gradient2(name = rowData(object)[s,"name"],mid = "lightgray", low = "blue", high = "red") + chromVAR_theme(12) +
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      ggplotly(p2)
    })

  }

  shinyApp(ui, server, options = list(width = 750))
}

